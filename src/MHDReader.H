#include "Proto.H"
#include "hdf5.h"
#include "MHD_Input_Parsing.H"
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <sstream>
extern Parsefrominputs inputs;

using namespace Proto;

class MHDReader
{
    public:
    inline MHDReader(){}

    inline static std::string resolvedPath(const std::string& a_filename)
    {
        std::error_code error;
        auto absolutePath = std::filesystem::absolute(a_filename, error);
        return error ? a_filename : absolutePath.lexically_normal().string();
    }

    inline static void abortFileError(
            const std::string& a_filename,
            const std::string& a_reason)
    {
        int rank = 0;
#ifdef PR_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        std::ostringstream message;
        message << "\nHelioCubed ERROR [MPI rank " << rank << "]: "
                << a_reason << "\n"
                << "  Requested HDF5 file: " << a_filename << "\n"
                << "  Resolved HDF5 path:  " << resolvedPath(a_filename)
                << "\n";
        std::fputs(message.str().c_str(), stderr);
        std::fflush(stderr);
#ifdef PR_MPI
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
        std::abort();
    }

    inline static void requireFileExists(const std::string& a_filename)
    {
        std::error_code error;
        bool exists = std::filesystem::exists(a_filename, error);
        if (error || !exists)
        {
            std::string reason = "HDF5 input file does not exist or is inaccessible.";
            if (error)
            {
                reason += " Filesystem error: " + error.message();
            }
            abortFileError(a_filename, reason);
        }
    }

    inline static void requireFileOpened(hid_t a_file, const std::string& a_filename)
    {
        if (a_file < 0)
        {
            abortFileError(a_filename, "HDF5 was unable to open the input file.");
        }
    }
   
     
    template<typename... Args>
    inline void readFile(hid_t* a_file, std::string a_filename, Args... a_params)
    {
        char fname[100];

        // The following pragmas suppress some GCC warnings that are overly conservative
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
        sprintf(fname, a_filename.c_str(), a_params...);
#pragma GCC diagnostic pop

        requireFileExists(fname);

#ifdef PR_MPI
        auto mpi_comm = MPI_COMM_WORLD;
        auto mpi_info = MPI_INFO_NULL;

        auto p_access = H5Pcreate(H5P_FILE_ACCESS);
        assert(H5Pset_fapl_mpio(p_access, mpi_comm, mpi_info) >= 0);
        *a_file = H5Fopen(fname, H5F_ACC_RDONLY, p_access);
        assert(H5Pclose(p_access) >= 0);
#else
        *a_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
#endif
        requireFileOpened(*a_file, fname);
    }
   
    template<unsigned int C>
    inline Box readDomain(hid_t* a_file)
    {
        // READ NUM COMPONENTS
        int numComps;
        auto att_comps = H5Aopen(*a_file, "num_components", H5P_DEFAULT);
        assert(H5Aread(att_comps, H5T_STD_I32LE, &numComps) >= 0); 
        assert(H5Aclose(att_comps) >= 0);

        PROTO_ASSERT(numComps == C,
                "MHDReader::readLevel | Error: Incorrect number of components. %i != %i",
                C, numComps);

        // READ DOMAIN SIZES
        int domainSize[2];
        auto att_domain = H5Aopen(*a_file, "domain", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
        
        return Box(Point(1, domainSize[0], domainSize[1]));
    }
    
    template<typename T, unsigned int C, MemType MEM, typename... Args>
    inline void readDataset(hid_t* a_file,
            BoxData<T, C, MEM>& a_data, 
            int a_index)
    {
        hid_t H5T_T;
        HDF5Handler::getH5DataType<T>(&H5T_T);
        auto ds_data = H5Dopen2(*a_file, ("/data" + std::to_string(a_index)).c_str(), H5P_DEFAULT);
        if (MEM == HOST)
        {
            assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data.data()) >= 0);
        } else {
            int dataSize = a_data.size();
            T* rawData = (T*)proto_malloc<HOST>(dataSize);
            assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData) >= 0);
            proto_memcpy<HOST, MEM>(rawData, a_data.data(), dataSize);
            proto_free<HOST>(rawData);
        }
        assert(H5Dclose(ds_data) >= 0);
    }

    template<typename T, unsigned int C, MemType MEM, typename... Args>
    inline void readData(
            std::vector<BoxData<T, C, MEM>>& a_data,
            std::string a_filename,
            Args... a_params)
    {
        PROTO_ASSERT(DIM == 3,
                "MHDReader::readLevel | Error: Invalid for DIM != 3");

        PR_TIME("MHDReader::readLevel");
        // PARSE THE FILE NAME
        char fname[100];

        // The following pragmas suppress some GCC warnings that are overly conservative
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
        sprintf(fname, a_filename.c_str(), a_params...);
#pragma GCC diagnostic pop

        requireFileExists(fname);

#ifdef PR_MPI
        auto mpi_comm = MPI_COMM_WORLD;
        auto mpi_info = MPI_INFO_NULL;

        auto p_access = H5Pcreate(H5P_FILE_ACCESS);
        assert(H5Pset_fapl_mpio(p_access, mpi_comm, mpi_info) >= 0);
        auto file = H5Fopen(fname, H5F_ACC_RDONLY, p_access);
        assert(H5Pclose(p_access) >= 0);
#else
        auto file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
#endif
        requireFileOpened(file, fname);
        hid_t H5T_T;
        HDF5Handler::getH5DataType<T>(&H5T_T);

        // READ DOMAIN SIZES
        int domainSize[2];
        auto att_domain = H5Aopen(file, "domain", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
        
        // READ NUM COMPONENTS
        int numComps;
        auto att_comps = H5Aopen(file, "num_components", H5P_DEFAULT);
        assert(H5Aread(att_comps, H5T_STD_I32LE, &numComps) >= 0); 
        assert(H5Aclose(att_comps) >= 0);

        PROTO_ASSERT(numComps == C,
                "MHDReader::readLevel | Error: Incorrect number of components. %i != %i",
                C, numComps);

        PROTO_ASSERT(inputs.domainSizez > 0 && inputs.domainSizey > 0,
                "MHDReader::readData | Error: theta and phi domain sizes must be positive.");
        PROTO_ASSERT(domainSize[0] % inputs.domainSizez == 0,
                "MHDReader::readData | Error: HDF5 phi size (%i) is not divisible by domainSizez (%i).",
                domainSize[0], inputs.domainSizez);
        PROTO_ASSERT(domainSize[1] % inputs.domainSizey == 0,
                "MHDReader::readData | Error: HDF5 theta size (%i) is not divisible by domainSizey (%i).",
                domainSize[1], inputs.domainSizey);

        int coarse_ratio_phi = domainSize[0]/inputs.domainSizez;
        int coarse_ratio_theta = domainSize[1]/inputs.domainSizey;
        int coarse_phi_size = domainSize[0]/coarse_ratio_phi;
        int coarse_theta_size = domainSize[1]/coarse_ratio_theta;
        int coarse_frame_size = coarse_phi_size*coarse_theta_size;

        int numDatasets;
        att_domain = H5Aopen(file, "num_datasets", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, &numDatasets) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
        
        // READ DATA
        int dataSize = domainSize[0]*domainSize[1]*numComps;
        
        auto ds_data = H5Dopen2(file, "/data0", H5P_DEFAULT);

        a_data.clear();
        a_data.resize(numDatasets);
        for (int ii = 0; ii < numDatasets; ii++)
        {
            ds_data = H5Dopen2(file, ("/data" + std::to_string(ii)).c_str(), H5P_DEFAULT);
            // DEFINE BOX DATA
            Box domainBox(Point::Zeros(),Point(0, coarse_theta_size-1, 2*coarse_phi_size-1));
            a_data[ii].define(domainBox);
            
            Box domainBox2(Point::Zeros(),Point(0, domainSize[0]-1, domainSize[1]-1));
            BoxData<double,C> a_data2(domainBox2);
            
            
            if (MEM == HOST)
            {
            //     assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data.data()) >= 0);
                assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data2.data()) >= 0);
                // This is needed to make transpose. h5 file is (phi, theta). We need (theta, phi).
                BoxData<double,C> a_data3 = Stencil<double>::AvgDown(
                        Point(1, coarse_ratio_phi, coarse_ratio_theta))(a_data2);
                for (int c = 0; c < numComps; c++){  
                    for (int i = 0; i < coarse_phi_size; i++ )
                    {
                        for (int j = 0; j < coarse_theta_size; j++ )
                        {
                            // Index in the final matrix.
                            int index1 = 2*c*coarse_frame_size + i*coarse_theta_size+j;
                            // Index in the original matrix.
                            int index2 = c*coarse_frame_size + (j*coarse_phi_size+i);
                            a_data[ii].data()[index1] = a_data3.data()[index2];
                            a_data[ii].data()[index1+coarse_frame_size] = a_data[ii].data()[index1];
                        }
                    }
                }

            } else {
                T* rawData = (T*)proto_malloc<HOST>(dataSize);
                assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData) >= 0);
                proto_memcpy<HOST, MEM>(rawData, a_data[ii].data(), dataSize);
                proto_free<HOST>(rawData);
            }
        }
        assert(H5Dclose(ds_data) >= 0);
        barrier();
        // assert(H5Fclose(file) >= 0);
    }


    
    template<typename T, typename... Args>
    inline void readGeom(
            std::vector<T>& a_data,
            std::string a_filename,
            Args... a_params)
    {
        PROTO_ASSERT(DIM == 3,
                "MHDReader::readGeom | Error: Invalid for DIM != 3");

        PR_TIME("MHDReader::readGeom");
        hid_t file;
        readFile(&file, a_filename, a_params...);
        
        hid_t H5T_T;
        HDF5Handler::getH5DataType<T>(&H5T_T);

        // READ DOMAIN SIZES
        int domainSize[2];
        auto att_domain = H5Aopen(file, "domain", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
       
        // READ DATA
        int dataSize = domainSize[1];
        auto ds_data = H5Dopen2(file, "/geometry/dtheta", H5P_DEFAULT);

        // DEFINE VECTOR
        a_data.resize(dataSize);
        
        // READ DATA
        assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_data[0]) >= 0);
        
        // CLEAN UP
        assert(H5Dclose(ds_data) >= 0);
        barrier();
        assert(H5Fclose(file) >= 0);
    }
};
