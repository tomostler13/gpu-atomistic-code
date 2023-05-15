# Compilers
SHELL:=/bin/bash
LIBS=-L/home/731392/lib/libconfig/1.7.3/lib -lconfig++ -I/home/731392/lib/fftw/3.3.10/include -L/home/731392/lib/fftw/3.3.10/include
INC=-I/home/731392/lib/libconfig/1.7.3/include -I/home/731392/lib/fftw/3.3.10/include -L/home/731392/lib/fftw/3.3.10/include
HOSTNAME=$(shell hostname)
OPT_LEVEL=-O3
GCC=g++ -g $(OPTLEVEL)
export LANG=C
export LC_ALL=C
# LIBS
DEFS=-DNDEBUG
#This part is hostname dependent (library paths etc)
include files.in


EXECUTABLE=gen_dw_grid

all: $(OBJECTS) gcc

# Serial Targets
gcc: $(OBJECTS)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTS) -o $(EXECUTABLE) $(LIBS) $(INC)

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(INC) $<


clean:
	@rm -f obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE)

tidy:	
	@rm -f *~
	@rm -f src/*~
