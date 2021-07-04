CC = gcc

CCFLAGS = -g

LIB_PATH = -L/mnt/scratch/papi-3.6.2/lib
INCLUDE = -I/mnt/scratch/papi-3.6.2/include
LIBS = -lpapi -lcblas -latlas -lpthread

EXECUTABLES = matvector
MATMULT_OBJS = matvector.o
NUMBERGEN_OBJS = numbergen.o

all: $(EXECUTABLES)

.c.o:
	$(CC) $(CCFLAGS) $(INCLUDE) -c $<

matvector: $(MATMULT_OBJS)
	$(CC) -o  matvector $(CCFLAGS) $(LIB_PATH) $(LIBS) $(INCLUDE) $(MATMULT_OBJS)

clean:
	rm -f *.o $(EXECUTABLES)

