# CC=g++
CC=mpiCC
CFLAGS=-fopenmp -O3 -lm
DEPS = piccoulomb.h
OBJ = piccoulomb.o

ifeq ($(MPI),1)
	DEFINES+=-DENABLEMPI
	CC=mpiCC
else
	CC=gcc
endif

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) $(DEFINES) -c -o $@ $<

bmc: piccoulomb.o
	$(info Using compiler ${CC})
	$(CC) $(CFLAGS) -o piccoulomb $(OBJ)

clean:
	rm -f *.o piccoulomb