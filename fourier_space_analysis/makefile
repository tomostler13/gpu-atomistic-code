# Compilers
SHELL:=/bin/bash
GCC=g++
export LANG=C
export LC_ALL=C
# LIBS
DEFS=-DNDEBUG
LIBS= -lfftw3 -lfftw3f -lm  -lstdc++
STATIC_LINK=
CPULIBS=
STATIC_LINK=
OPT_LEVEL=-O3
GCC_FLAGS= $(OPT_LEVEL)

# Objects
OBJECTS= \
obj/timeseries.o

EXECUTABLE=timeseries

all: $(OBJECTS) gcc

# Serial Targets
gcc: $(OBJECTS)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTS) $(SWITCHOBJ) $(STATIC_LINK) -o $(EXECUTABLE) $(LIBS)


$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<

clean:
	@rm -f obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE)

tidy:	
	@rm -f *~
	@rm -f src/*~
