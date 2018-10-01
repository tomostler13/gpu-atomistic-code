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
OBJECTSTSALL= \
obj/timeseries_all.o
OBJECTSTSTM= \
obj/timeseries_tm.o
OBJECTSTS= \
obj/timeseries.o
OBJECTSSSF= \
obj/analyse_structure_factor.o
OBJECTS2SLSSF= \
obj/twosublatssf.o

EXECUTABLETS=timeseries
EXECUTABLETSALL=timeseries_all
EXECUTABLETSTM=timeseries_tm
EXECUTABLESSF=ssf
EXECUTABLE2SLSSF=2slssf

all: $(OBJECTS) tsa ts ssf tm

# Time series all targets
tsa: $(OBJECTSTSALL)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTSTSALL) $(STATIC_LINK) -o $(EXECUTABLETS) $(LIBS)


$(OBJECTSTSALL): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<

# Time series two magnon
tm: $(OBJECTSTSTM)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTSTSTM) $(STATIC_LINK) -o $(EXECUTABLETM) $(LIBS)


$(OBJECTSTSTM): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<

# Time series targets
ts: $(OBJECTSTS)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTSTS) $(STATIC_LINK) -o $(EXECUTABLETSALL) $(LIBS)


$(OBJECTSTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<

# Time series targets
ssf: $(OBJECTSSSF)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTSSSF) $(STATIC_LINK) -o $(EXECUTABLESSF) $(LIBS)


$(OBJECTSSSF): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<

# Time series targets
ssf: $(OBJECTS2SLSSF)
	$(GCC) $(DEFS) $(GCC_FLAGS) $(OBJECTS2SLSSF) $(STATIC_LINK) -o $(EXECUTABLE2SLSSF) $(LIBS)


$(OBJECTS2SLSSF): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(DEFS) $(GCC_FLAGS) $(GITINFO) $<

clean:
	@rm -f obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE)

tidy:	
	@rm -f *~
	@rm -f src/*~
