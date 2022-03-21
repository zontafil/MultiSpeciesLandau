# CC=g++
CC=mpiCC
CFLAGS=-lm
DEPS = coulombStructurePreserving.h coulomb_kernel.h
OBJ = coulomb_kernel.o
OBJTEST = test.o

ifdef SILLY
	DEFINES+=-DVERBOSE_LEVEL=3
else ifdef VERBOSE
	DEFINES+=-DVERBOSE_LEVEL=2
else
	DEFINES+=-DVERBOSE_LEVEL=1
endif

ifdef DEBUG
	CFLAGS+=-O0 -g
else
	CFLAGS+=-O3
endif

ifdef CUDA
	CC=nvcc
	DEFINES+=-DCUDA=1
else
	CC=g++
endif

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) $(DEFINES) -c -o $@ $<
%.o: %.cu $(DEPS)
	$(CC) $(CFLAGS) $(DEFINES) -c -o $@ $<

coulombStructurePreserving: coulombStructurePreserving.o $(OBJ)
	$(info Using compiler ${CC})
	$(CC) $(CFLAGS) -o $@ $^
test: test.o
	$(info Using compiler ${CC})
	$(CC) $(CFLAGS) -o test $(OBJTEST)

clean:
	rm -f *.o coulombStructurePreserving test