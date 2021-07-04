################################################################################


export INCLUDE_PATH=/share/apps/mkl/include 
export INCLUDE_PATH=/home/asaadat/project/plasma-installer_2.4.5/install/include

export LD_LIBRARY_PATH=/share/apps/mkl/lib/intel64 
export LD_LIBRARY_PATH=/home/asaadat/project/plasma-installer_2.4.5/install/lib 


CC = gcc -fopenmp

CCFLAGS = -g

all: example
	echo "Done"

# Add source files here
#EXECUTABLE	:= example_dpotrf

# C/C++ source files (compiled with gcc / c++)
CCFILES		:= example_dpotrf.c
# Additional libraries needed by the project


LIB := -lm -lpapi -L/share/apps/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

PLASMA_LIB := -L/home/asaadat/project/plasma-installer_2.4.5/install/lib -lplasma -lcoreblas -lplasma -lquark -llapacke 
		 
LIB_HOME=/share/apps/mkl/lib/intel64


# include files
INCLUDE := -I/share/apps/mkl/include

PLASMA_INCLUDE := -I/home/asaadat/project/plasma-installer_2.4.5/install/include

.c.o:
	$(CC) $(CCFLAGS) $(INCLUDE) $(PLASMA_INCLUDE) -c $<

example : 
	${CC} ${CCFILES} ${CCFLAGS} ${INCLUDE} ${PLASMA_INCLUDE} ${LIB} ${PLASMA_LIB}  -o dpotrf

clean:
	${RM} dpotrf *.out *.a

################################################################################

