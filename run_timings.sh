
# remove the previous version
rm a.out

# build the program
icc main.c -Wall -openmp -lm

# run number of times
# run 1
for i in 1 2 4 8 16 32 64; do
    export OMP_NUM_THREADS=$i
    ./a.out 3 7 3
done

# run 2
for i in 1 2 4 8 16 32 64; do
    export OMP_NUM_THREADS=$i
    ./a.out 3 8 3
done

