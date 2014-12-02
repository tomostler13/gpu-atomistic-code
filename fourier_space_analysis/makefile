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
OBJECTSTS= \
obj/timeseries.o
OBJECTSSSF= \
obj/analyse_structure_factor.o

EXECUTABLETS=timeseries
EXECUTABLESSF=ssf

all: $(OBJECTS) ts ssf

# Time series targets
ts: $(OBJECTSTS)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTSTS) $(STATIC_LINK) -o $(EXECUTABLETS) $(LIBS)


$(OBJECTSTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<

# Time series targets
ssf: $(OBJECTSSSF)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTSSSF) $(STATIC_LINK) -o $(EXECUTABLESSF) $(LIBS)


$(OBJECTSSSF): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<

clean:
	@rm -f obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE)

tidy:	
	@rm -f *~
	@rm -f src/*~
