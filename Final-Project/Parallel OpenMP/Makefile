################################################################################
CC = gcc -fopenmp

all: example
	echo "Done"

# Add source files here
#EXECUTABLE	:= cholOMP9

# C/C++ source files (compiled with gcc / c++)
CCFILES		:= cholOMP9.c
# Additional libraries needed by the project

LIB := -lm -lpapi -L/share/apps/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

LIB_HOME=/share/apps/mkl/lib/intel64

export LD_LIBRARY_PATH=/share/apps/mkl/lib/intel64 

# include files
INCLUDE := -I/share/apps/mkl/include

example : 
	${CC} ${LIB} ${INCLUDE} ${CCFILES} -o cholOMP9

clean:
	${RM} cholOMP8 cholOMP9 *.out *.a

################################################################################

