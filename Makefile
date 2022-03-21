# CC=g++
CC=mpiCC
CFLAGS=-fopenmp -lm
DEPS = coulombStructurePreserving.h
OBJ = coulombStructurePreserving.o
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

CC=g++

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) $(DEFINES) -c -o $@ $<

coulombStructurePreserving: coulombStructurePreserving.o
	$(info Using compiler ${CC})
	$(CC) $(CFLAGS) -o coulombStructurePreserving $(OBJ)
test: test.o
	$(info Using compiler ${CC})
	$(CC) $(CFLAGS) -o test $(OBJTEST)

clean:
	rm -f *.o coulombStructurePreserving test