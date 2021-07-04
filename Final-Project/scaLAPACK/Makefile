################################################################################

export LD_LIBRARY_PATH=/share/apps/mkl/lib/intel64 
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib 

CC = mpicc

all: example
	echo "Done"

# Add source files here
#EXECUTABLE	:= LUtest

# C/C++ source files (compiled with gcc / c++)
CCFILES		:= LUtest.c
# Additional libraries needed by the project

#LIB := -lm -lpapi -L/share/apps/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
LIB := -lcblas -llapack

#LIB_HOME :=/share/apps/mkl/lib/intel64

#LIBBLACS := -L/share/apps/mkl/lib/intel64 -lmkl_blacs_ilp64 -lmkl_blacs_intelmpi_ilp64 -lmkl_blacs_intelmpi_lp64 -lmkl_blacs_lp64 -lmkl_blacs_openmpi_ilp64 -lmkl_blacs_openmpi_lp64 -lmkl_blacs_sgimpt_ilp64 -lmkl_blacs_sgimpt_lp64 -lmkl_scalapack_ilp64 -lmkl_scalapack_lp64

LIBBLACS :=-L/home/mnafarse/lib/lib -lscalapack

#LIBMPI := -L/usr/lib64/openmpi/lib -lmpi 

#LIBBLACS := -L/share/apps/mkl/lib/intel64 -lmkl_blacs_ilp64 -lmkl_blacs_intelmpi_ilp64 -lmkl_blacs_intelmpi_lp64 -lmkl_blacs_lp64  -lmkl_scalapack_ilp64 -lmkl_scalapack_lp64


# include files
#INCLUDE := -I/share/apps/mkl/include

example : 
	${CC} ${CCFILES} ${INCLUDE} ${LIBMPI} ${LIBBLACS} ${LIB} -o LUtest

clean:
	${RM} LUtest *.out *.a

################################################################################

