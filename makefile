# Compilers
SHELL:=/bin/bash
HOSTNAME=$(shell hostname)
GITINFO=-DGIT_SHA1='"$(shell git rev-parse HEAD)"' -DGITDIRTY='"$(shell git status -s | grep -v ? | wc -l)"'
GCC=g++  -DCOMP='"GNU C++ Compiler $(shell g++ --version | head -n 1 | cut -b 5-)"' -DHOSTNAME='"$(shell hostname)"' ${GITINFO}
NVCC=nvcc  -DCOMP='"NVIDIA C++ Compiler $(shell nvcc --version | tail -n 2 | head -n 1)"' -DHOSTNAME='"$(shell hostname)"' ${GITINFO}
export LANG=C
export LC_ALL=C
# LIBS
DEFS=-DNDEBUG
CUDEFS=-DCUDA
#This part is hostname dependent (library paths etc)
include host.args
include files.in

CUDA_OBJECTS=$(OBJECTS:.o=_cuda.o)
NVCC_OBJECTS=$(NVCCOBJ:.o=_cuda.o)
SWITCH_OBJECTS=$(SWITCHOBJ:.o=_cuda.o)

EXECUTABLE=ASD

all: $(OBJECTS) gcc

# Serial Targets
gcc: $(OBJECTS) $(SWITCHOBJ)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTS) $(SWITCHOBJ) $(STATIC_LINK) -o $(EXECUTABLE) $(LIBS)

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $<
$(SWITCHOBJ): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $<


# cuda targets
gcc-cuda: $(SWITCH_OBJECTS) $(NVCC_OBJECTS) $(CUDA_OBJECTS)
	$(NVCC) $(CUDA_OBJECTS) $(SWITCH_OBJECTS) $(NVCC_OBJECTS) $(STATIC_LINK) $(CUDALIBS) $(LIBS) -o $(EXECUTABLE) $(DEFS)

$(CUDA_OBJECTS): obj/%_cuda.o: src/%.cpp
	$(GCC) $(GCC_FLAGS) $(CUDEFS) $(DEFS) -c $< $(CPULIBS) $(LIBS) -o $@

$(NVCC_OBJECTS) : obj/%_cuda.o: src/%.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDEFS) $(DEFS) -c $< $(CUDALIBS) $(LIBS) -o $@

$(SWITCH_OBJECTS) : obj/%_cuda.o: src/%.cpp
	$(NVCC) $(NVCC_FLAGS) $(CUDEFS) $(DEFS) -c $< $(CUDALIBS) $(LIBS) -o $@

clean:
	@rm -f obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE)

tidy:	
	@rm -f *~
	@rm -f src/*~
