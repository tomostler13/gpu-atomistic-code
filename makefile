# Compilers
SHELL:=/bin/bash
GCC=g++  -DCOMP='"GNU C++ Compiler $(shell g++ --version | head -n 1 | cut -b 5-)"' -DHOSTNAME='"$(shell hostname)"'
NVCC=nvcc  -DCOMP='"NVIDIA C++ Compiler $(shell nvcc --version | tail -n 2 | head -n 1)"' -DHOSTNAME='"$(shell hostname)"'
GITINFO=-DGIT_SHA1='"$(shell git rev-parse HEAD)"' -DGITDIRTY='"$(shell git status -s | wc -l)"'
export LANG=C
export LC_ALL=C
# LIBS
DEFS=-DNDEBUG
CUDEFS=-DCUDA
LIBS= -lfftw3_omp -lfftw3 -lfftw3f -lm -lconfig++ -lstdc++ -llapack -lblas
CPULIBS= -fopenmp -lpthread
CUDALIBS= -L/usr/local/cuda/lib64/ -lcurand -lcudart -lcufft
OPT_LEVEL=-O3
GCC_FLAGS= $(OPT_LEVEL)
#NVCC_FLAGS= -g $(OPT_LEVEL) -I/usr/local/cuda/include -m64 -ccbin /usr/bin/g++-4.4 --ptxas-options=-v -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20 
NVCC_FLAGS= $(OPT_LEVEL) -I/usr/local/cuda/include -m64 -ccbin /usr/bin/g++-4.4 --ptxas-options=-v -gencode=arch=compute_13,code=sm_13 -gencode=arch=compute_13,code=compute_13 


# Objects
OBJECTS= \
obj/config.o \
obj/error.o \
obj/random.o \
obj/geom.o \
obj/util.o \
obj/mat.o \
obj/intmat.o \
obj/spins.o \
obj/fields.o \
obj/exch.o \
obj/anis.o \
obj/llgCPU.o \
obj/maths.o \
obj/sim.o

SWITCHOBJ= \
obj/main.o \
obj/llg.o

NVCCOBJ= \
obj/cuda.o \
obj/cufields.o \
obj/cuint.o

CUDA_OBJECTS=$(OBJECTS:.o=_cuda.o)
NVCC_OBJECTS=$(NVCCOBJ:.o=_cuda.o)
SWITCH_OBJECTS=$(SWITCHOBJ:.o=_cuda.o)

EXECUTABLE=ASD

all: $(OBJECTS) gcc

# Serial Targets
gcc: $(OBJECTS) $(SWITCHOBJ)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTS) $(SWITCHOBJ) -o $(EXECUTABLE) $(LIBS) $(CPULIBS)

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<
$(SWITCHOBJ): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<


# cuda targets
gcc-cuda: $(SWITCH_OBJECTS) $(NVCC_OBJECTS) $(CUDA_OBJECTS)
	$(NVCC) $(CUDA_OBJECTS) $(SWITCH_OBJECTS) $(NVCC_OBJECTS) $(CUDALIBS) $(LIBS) -o $(EXECUTABLE) $(GITINFO) $(DEFS)

$(CUDA_OBJECTS): obj/%_cuda.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_FLAGS) $(DEFS) $(GITINFO) $<

$(NVCC_OBJECTS) : obj/%_cuda.o: src/%.cu
	$(NVCC) $(NVCC_FLAGS) $(GITINFO) $(CUDEFS) $(DEFS) -c $< -o $@

$(SWITCH_OBJECTS) : obj/%_cuda.o: src/%.cpp
	$(NVCC) $(NVCC_FLAGS) $(GITINFO) $(CUDEFS) $(DEFS) -c $< -o $@

clean:
	@rm -f obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE)

tidy:	
	@rm -f *~
	@rm -f src/*~
