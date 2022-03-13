# CC=g++
CC=mpiCC
CFLAGS=-fopenmp -lm -O0 -g
# CFLAGS=-fopenmp -O3 -lm
DEPS = piccoulomb.h
OBJ = piccoulomb.o
OBJTEST = test.o

ifdef VERBOSE
	DEFINES+=-DVERBOSE_LEVEL=2
else
	DEFINES+=-DVERBOSE_LEVEL=1
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