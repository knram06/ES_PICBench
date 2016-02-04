Building the code:
The code consists of a single main.c file which includes all other relevant headers. Can be compiled

with gcc:
gcc main.c -Wall -fopenmp -lm

with icc:
icc main.c -Wall -openmp -lm

to get a default a.out executable.

Running the program:
The executable takes 3 arguments:
./a.out <coarse grid points on one side> <number of levels> <gauss seidel iterations>

So a typical run might be:
./a.out 3 6 3

which corresponds to having 3 nodal points on the coarset grid and the finest level would have (2*2^5)+1=65 nodal points.

Please obtain timing runs by running the script run_timings.sh which carries
out timing runs for
./a.out 3 7 3 and
./a.out 3 8 3

up to 64 threads

