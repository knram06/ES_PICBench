icc main.c -O3 -march=native -Wall -L$PARDISO_LIC_PATH -lpardiso500-INTEL1301-X86-64 -L/usr/global/intel/composer_xe/2013.0.079/composer_xe_2013.0.079/compiler/lib/intel64  -L/usr/global/intel/mkl/11.0.0.079/lib/intel64 -lifcore -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lm -lgfortran -fopenmp
