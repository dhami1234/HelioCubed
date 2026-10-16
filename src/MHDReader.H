#pragma once

#include "Proto.H"
#include "hdf5.h"
#include "MHD_Error.H"
#include "MHD_Input_Parsing.H"
#include <cstdio>
#include <filesystem>
#include <type_traits>
#include <vector>

extern Parsefrominputs inputs;
using namespace Proto;

class MHDReader
{
public:
    inline static std::string resolvedPath(const std::string& filename)
    {
        std::error_code error;
        auto path = std::filesystem::absolute(filename, error);
        return error ? filename : path.lexically_normal().string();
    }

    inline static void abortFileError(const std::string& filename,
                                     const std::string& reason)
    {
        MHD::require(false, "%s\n  Requested HDF5 file: %s\n  Resolved HDF5 path: %s",
                     reason.c_str(), filename.c_str(), resolvedPath(filename).c_str());
    }

    inline static void requireFileExists(const std::string& filename)
    {
        std::error_code error;
        bool exists = std::filesystem::exists(filename, error);
        if (error || !exists)
        {
            std::string reason = "HDF5 input file does not exist or is inaccessible.";
            if (error) reason += " Filesystem error: " + error.message();
            abortFileError(filename, reason);
        }
    }

    inline static void requireFileOpened(hid_t file, const std::string& filename)
    {
        if (file < 0) abortFileError(filename, "HDF5 was unable to open the input file.");
    }

    template<typename... Args>
    inline void readFile(hid_t* file, std::string filename, Args... params)
    {
        filename = formatFilename(filename, params...);
        requireFileExists(filename);
#ifdef PR_MPI
        auto access = checked(H5Pcreate(H5P_FILE_ACCESS), filename, "create file access properties");
        checked(H5Pset_fapl_mpio(access, MPI_COMM_WORLD, MPI_INFO_NULL), filename,
                "set MPI file access");
        *file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, access);
        checked(H5Pclose(access), filename, "close file access properties");
#else
        *file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
#endif
        requireFileOpened(*file, filename);
    }

    template<unsigned int C>
    inline Box readDomain(hid_t* file)
    {
        int sizes[2];
        readDomainSizes<C>(*file, sizes);
        return Box(Point(1, sizes[0], sizes[1]));
    }

    template<typename T, unsigned int C, MemType MEM>
    inline void readDataset(hid_t* file, BoxData<T, C, MEM>& data, int index)
    {
        const auto name = "/data" + std::to_string(index);
        if constexpr (MEM == HOST)
        {
            readValues(*file, name, data.data(), data.size());
        }
        else
        {
            std::vector<T> host(data.size());
            readValues(*file, name, host.data(), host.size());
            proto_memcpy<HOST, MEM>(host.data(), data.data(), host.size()*sizeof(T));
        }
    }

    template<typename T, unsigned int C, MemType MEM, typename... Args>
    inline void readData(std::vector<BoxData<T, C, MEM>>& data,
                         std::string filename, Args... params)
    {
        static_assert(DIM == 3, "MHDReader requires DIM=3");
        PR_TIME("MHDReader::readData");
        hid_t file;
        readFile(&file, filename, params...);
        filename = fileName(file);
        int sizes[2];
        readDomainSizes<C>(file, sizes);
        require(inputs.domainSizez > 0 && inputs.domainSizey > 0, filename,
                "Theta and phi domain sizes must be positive.");
        require(sizes[0] % inputs.domainSizez == 0 && sizes[1] % inputs.domainSizey == 0,
                filename, "HDF5 phi/theta sizes must be divisible by domainSizez/domainSizey.");
        const int phi = inputs.domainSizez, theta = inputs.domainSizey;
        const Point ratio(1, sizes[0]/phi, sizes[1]/theta);
        int frames;
        readIntAttribute(file, "num_datasets", &frames, 1);
        require(frames > 0, filename, "num_datasets must be positive.");

        data.clear();
        data.resize(frames);
        const Box rawBox(Point::Zeros(), Point(0, sizes[0]-1, sizes[1]-1));
        const Box outBox(Point::Zeros(), Point(0, theta-1, 2*phi-1));
        for (int frame = 0; frame < frames; ++frame)
        {
            BoxData<T, C, HOST> raw(rawBox);
            readDataset(&file, raw, frame);
            BoxData<T, C, HOST> coarse = Stencil<T>::AvgDown(ratio)(raw);
            BoxData<T, C, HOST> transposed(outBox);
            // Files store (phi, theta). Boundary interpolation expects
            // (theta, phi), with a repeated phi period for solar rotation.
            for (unsigned int c = 0; c < C; ++c)
                for (int i = 0; i < phi; ++i)
                    for (int j = 0; j < theta; ++j)
                    {
                        const T value = coarse(Point(0, i, j), c);
                        transposed(Point(0, j, i), c) = value;
                        transposed(Point(0, j, i+phi), c) = value;
                    }
            data[frame].define(outBox);
            proto_memcpy<HOST, MEM>(transposed.data(), data[frame].data(),
                                    transposed.size()*sizeof(T));
        }
        checked(H5Fclose(file), filename, "close boundary file");
    }

    template<typename T, typename... Args>
    inline void readGeom(std::vector<T>& data, std::string filename, Args... params)
    {
        static_assert(DIM == 3, "MHDReader requires DIM=3");
        PR_TIME("MHDReader::readGeom");
        hid_t file;
        readFile(&file, filename, params...);
        filename = fileName(file);
        int sizes[2];
        readIntAttribute(file, "domain", sizes, 2);
        require(sizes[0] > 0 && sizes[1] > 0, filename, "HDF5 domain sizes must be positive.");
        data.resize(sizes[1]);
        readValues(file, "/geometry/dtheta", data.data(), data.size());
        checked(H5Fclose(file), filename, "close geometry file");
    }

private:
    inline static void require(bool condition, const std::string& filename,
                               const std::string& reason)
    {
        if (!condition) abortFileError(filename, reason);
    }

    template<typename Status>
    inline static Status checked(Status status, const std::string& filename,
                                  const std::string& operation)
    {
        require(status >= 0, filename, "HDF5 failed to " + operation + ".");
        return status;
    }

    template<typename... Args>
    inline static std::string formatFilename(const std::string& pattern, Args... params)
    {
        if constexpr (sizeof...(Args) == 0) return pattern;
        const int size = std::snprintf(nullptr, 0, pattern.c_str(), params...);
        require(size >= 0, pattern, "Unable to format HDF5 filename.");
        std::vector<char> buffer(static_cast<size_t>(size)+1);
        std::snprintf(buffer.data(), buffer.size(), pattern.c_str(), params...);
        return std::string(buffer.data());
    }

    inline static std::string fileName(hid_t file)
    {
        const auto size = H5Fget_name(file, nullptr, 0);
        require(size >= 0, "<open HDF5 file>", "Unable to get HDF5 filename.");
        std::vector<char> buffer(static_cast<size_t>(size)+1);
        checked(H5Fget_name(file, buffer.data(), buffer.size()), "<open HDF5 file>", "get filename");
        return std::string(buffer.data());
    }

    inline static void readIntAttribute(hid_t file, const std::string& name,
                                        int* values, size_t count)
    {
        const auto filename = fileName(file);
        auto attr = checked(H5Aopen(file, name.c_str(), H5P_DEFAULT), filename, "open attribute " + name);
        auto space = checked(H5Aget_space(attr), filename, "get attribute space " + name);
        require(H5Sget_simple_extent_npoints(space) == static_cast<hssize_t>(count), filename,
                "Incorrect number of values in attribute " + name + ".");
        checked(H5Aread(attr, H5T_NATIVE_INT, values), filename, "read attribute " + name);
        checked(H5Sclose(space), filename, "close attribute space " + name);
        checked(H5Aclose(attr), filename, "close attribute " + name);
    }

    template<unsigned int C>
    inline static void readDomainSizes(hid_t file, int* sizes)
    {
        int components;
        readIntAttribute(file, "num_components", &components, 1);
        const auto filename = fileName(file);
        require(components == C, filename, "Incorrect number of HDF5 components: expected " +
                std::to_string(C) + ", got " + std::to_string(components) + ".");
        readIntAttribute(file, "domain", sizes, 2);
        require(sizes[0] > 0 && sizes[1] > 0, filename, "HDF5 domain sizes must be positive.");
    }

    template<typename T>
    inline static void readValues(hid_t file, const std::string& name, T* values, size_t count)
    {
        const auto filename = fileName(file);
        auto dataset = checked(H5Dopen2(file, name.c_str(), H5P_DEFAULT), filename, "open dataset " + name);
        auto space = checked(H5Dget_space(dataset), filename, "get dataset space " + name);
        require(H5Sget_simple_extent_npoints(space) == static_cast<hssize_t>(count), filename,
                "Incorrect number of values in dataset " + name + ".");
        hid_t type = -1;
        // Proto's generic type helper currently has no float specialization.
        if constexpr (std::is_same<T, float>::value) type = H5T_NATIVE_FLOAT;
        else HDF5Handler::getH5DataType<T>(&type);
        checked(type, filename, "select dataset memory type");
        checked(H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values), filename, "read dataset " + name);
        checked(H5Sclose(space), filename, "close dataset space " + name);
        checked(H5Dclose(dataset), filename, "close dataset " + name);
    }
};
