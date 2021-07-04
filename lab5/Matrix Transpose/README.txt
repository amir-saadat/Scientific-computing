myfile.dat is the input file as matrix C.

data.dat contains Matrix dimension, Number of processors in row, Number of processors in column, 
number of local row and number of local column, respectively. (n, Pr, Pc, nLocRow, nLocCol)


(Note that n is not necessarily Pr*nLocRow or Pc*nLocCol. Because the code is generally written 
even for prime number dimensions. So, it is possible that the rows or columns of most right and 
bottom processors are not equal to nLocRow or nLocCol. nLocRow, nLocCol are the original order 
and may be disobeyed by right and bottom borders)

please do the following:
1.clean make
2.make
3.(Let's say you have 2 row processes and 3 col processes)
mpirun -np (2*3) ./a.out
or
mpirun.openmpi -np (2*3) ./a.out

 
