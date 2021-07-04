################################################################################
CC = /usr/local/cuda/bin/nvcc

all: example
	echo "Done"

# Add source files here
#EXECUTABLE	:= chol_qr_it
# Cuda source files (compiled with cudacc)
CUFILES		:= chol_qr_main.cu chol_qr_it.cu
# C/C++ source files (compiled with gcc / c++)
CCFILES		:= 
# Additional libraries needed by the project
ifeq ($(emu),1)
        LIB := -lcublasemu
else
        LIB := -lcublas -lf77blas -latlas -lcblas -llapack -lblas -L/usr/local/cuda/sdk/C/lib -lcutil_x86_64 -lm

endif

# include files
INCLUDE := -I/usr/local/cuda/sdk/C/common/inc

example : 
	${CC} ${LIB} ${INCLUDE} ${CUFILES} -o chol_qr_it
################################################################################

