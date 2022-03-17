# CC=g++
CC=mpiCC
CFLAGS=-fopenmp -lm
DEPS = piccoulomb.h
OBJ = piccoulomb.o
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

piccoulomb: piccoulomb.o
	$(info Using compiler ${CC})
	$(CC) $(CFLAGS) -o piccoulomb $(OBJ)
test: test.o
	$(info Using compiler ${CC})
	$(CC) $(CFLAGS) -o test $(OBJTEST)

clean:
	rm -f *.o piccoulomb test