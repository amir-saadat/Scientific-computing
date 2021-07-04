
CC=mpicc
#CFLAGS=-Wall -g

all: cholMPI1
	echo "Done"

INCLUDE := -I/usr/include

LIB := -lcblas -llapack -lm

cholMPI1 : cholMPI1.c
	${CC} ${CFLAGS} ${INCLUDE} ${LIB} -o cholMPI1 cholMPI1.c


clean:
	${RM} *.out core cholMPI cholMPI1

