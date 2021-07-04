
CC=g++
#CFLAGS=-Wall -g

all: chol
	echo "Done"

LIB := -lpapi -lcblas -llapack

chol : chol.cpp
	${CC} ${CFLAGS} ${LIB} -o chol.out chol.cpp


clean:
	${RM} *.out core

