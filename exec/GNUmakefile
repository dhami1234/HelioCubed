VARS_OLD := $(.VARIABLES)
#DEPXX=clang++
#CXX=clang++
DEPXX=mpicxx
CXX=mpicxx
DIM=3
TURB=0
USE_HDF5=TRUE
HDF5_LIB=/opt/homebrew/Cellar/hdf5-mpi/1.14.2/lib
HDF5_INC=/opt/homebrew/Cellar/hdf5-mpi/1.14.2/include
# HDF5_LIB=/opt/homebrew/Cellar/hdf5-mpi/1.12.2_1/lib
# HDF5_INC=/opt/homebrew/Cellar/hdf5-mpi/1.12.2_1/include
# HDF5_LIB=/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/intel-19.1.1.217/hdf5-1.10.6-ow5p5eh25vi4yxdpeenmtinbjed3xt7e/lib
# HDF5_INC=/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/intel-19.1.1.217/hdf5-1.10.6-ow5p5eh25vi4yxdpeenmtinbjed3xt7e/include
# HDF5_LIB=$(TACC_HDF5_DIR)/lib
# HDF5_INC=$(TACC_HDF5_DIR)/include

# BLIS_INC=/opt/homebrew/Cellar/blis/0.9.0/include/blis
PROTO= ../../proto

# LAPACK_LIB=-lblas -llapack
CXXOPTFLAGS=-O3 
CXXOPTFLAGS+=-DmemoryFootPrint
CXXOPTFLAGS+= -DPR_HDF5=TRUE
CXXOPTFLAGS+= -DPR_MPI=TRUE
CXXOPTFLAGS+= -DPR_STACK_ALLOC=4294967296
CXXOPTFLAGS+= -DPR_OPS=TRUE
TARGET:=HelioCubed.exe

SRC=$(patsubst ../src/%,%,$(wildcard ../src/*.cpp)) $(wildcard *.cpp) $(wildcard $(PROTO)/examples/_common/*.cpp) 
NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))

NODIRDEP = $(subst .cpp,.d,$(NODIRSRC))
DEPENDS=$(patsubst %.d,d/%.d, $(NODIRDEP))

VPATH=o. ../src  $(PROTO)/examples/_common



CPPFLAGS:=-DDIM=$(DIM) -DTURB=$(TURB) -I$(PROTO)/include -I$(PROTO)/include/base -I$(PROTO)/include/amr -I$(PROTO)/EBProto/include -I. -I../src -I$(PROTO)/examples/_common 
CXXFLAGS:= $(CXXOPTFLAGS) -std=c++14 -w

LIBS += $(LAPACK_LIB)

ifdef HDF5_LIB
LIBS+= -L$(HDF5_LIB) -lhdf5 -lz
CPPFLAGS+= -I$(HDF5_INC)
endif

ifdef BLIS_INC
CPPFLAGS+= -I$(BLIS_INC)
endif

ifdef FFTW3_DIR
LIBS+= -L$(FFTW3_DIR)/lib -lfftw3
CPPFLAGS+= -I$(FFTW3_DIR)/include
endif

DEP=$(DEPXX) -MM $(CXXFLAGS)
ifeq ($(DEPXX),nvcc)
DEP=$(DEPXX) -M $(CXXFLAGS)
LIBS+= -lnvToolsExt
endif

ifeq ($(DEPXX),hipcc)
CPPFLAGS+=-I$(ROCM_PATH)/rocthrust/include/ -I$(ROCM_PATH)/rocprim/include/ -I$(ROCM_PATH)/include/roctracer
CXXFLAGS+=--amdgpu-target=gfx906 -Wno-non-c-typedef-for-linkage -Wno-unused-command-line-argument
LIBS+= -L/$(ROCM_PATH)/lib/
LIBS+= -lroctx64
endif


LINK=$(CXX)
OBJDIR = o
DEPDIR = d


$(TARGET): $(OBJ)
	$(LINK) -o $(TARGET) $(OBJ) $(LIBS)

$(OBJDIR):
	@if [ ! -d "./$(OBJDIR)" ];then     \
		mkdir $(OBJDIR);           \
	fi

$(DEPDIR):
	@if [ ! -d "./$(DEPDIR)" ];then     \
		mkdir $(DEPDIR);           \
	fi

#compile non-main sources
$(OBJDIR)/%.o : %.cpp   $(OBJDIR) $(DEPDIR)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(CPPFLAGS)
	$(DEP) $<   $(CPPFLAGS) > $(DEPDIR)/$*.d



.PHONY: clean dataclean realclean print dirs


print-%: ; @echo $* = $($*)

-include $(DEPENDS)

vars:
	$(foreach v,$(filter-out $(VARS_OLD) VARS_OLD,$(.VARIABLES)), $(info $(v) = $($(v))))

clean:
	rm -rf d o *.o *.exe *.d $(PROTO)/examples/_common/*.o ../src/*.o $(PROTO)/examples/_common/*.d ../src/*.d

dataclean:
	rm -rf *.hdf5 *.vtk *.table pout*



