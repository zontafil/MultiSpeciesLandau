# CC=g++
CC=mpiCC
CFLAGS=-lm
DEPS = coulombStructurePreserving.h coulomb_kernel.h
OBJ = coulomb_kernel.o coulombStructurePreserving.o

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

ifeq ($(CUDA), 1)
	CC=nvcc
	CFLAGS+=--expt-relaxed-constexpr
	DEFINES+=-DCUDA=1
	DEFINES+=-DEIGEN_MAX_STATIC_ALIGN_BYTES=16
else
	CC=g++
endif

ifeq ($(CUDA), 1)
%.o: %.cu $(DEPS)
	$(CC) $(CFLAGS) $(DEFINES) -c -o $@ $<
endif
%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) $(DEFINES) -c -o $@ $<

coulomb: main.o $(OBJ)
	$(info Using compiler ${CC})
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f *.o coulomb