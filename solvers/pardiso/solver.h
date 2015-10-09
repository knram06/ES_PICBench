#ifndef PARDISO_SOLVER_H
#define PARDISO_SOLVER_H

/* PARDISO prototype. */
// taken from pardiso_unsym.c
void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);

int     n;
int     *ia;
int     *ja;
double  *a;
int      nnz;
int      mtype; //= 11;        /* Real unsymmetric matrix */

/* RHS and solution vectors. */
double   *b, *x;
int      nrhs;          /* Number of right hand sides. */

/* Internal solver memory pointer pt,                  */
/* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
/* or void *pt[64] should be OK on both architectures  */ 
void    *pt[64];

/* Pardiso control parameters. */
int      iparm[64];
double   dparm[64];
int      solver;
int      maxfct, mnum, phase, error, msglvl;

/* Number of processors. */
int      num_procs;

/* Auxiliary variables. */
char    *var;
int      i, k;

double   ddum;              /* Double dummy */
int      idum;              /* Integer dummy. */

int SolverInitialize(int *argc, char ***argv)
{
    mtype = 11;         // real unsymmetric matrix
    solver = 0;         // use the recursive iterative solver
    error = 0;
    nrhs = 1;           // number of right hand sides

    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         return 1;
    }
    else
        printf("[PARDISO]: License check was successful ... \n");

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;

    iparm[10] = 0; /* no scaling  */
    iparm[12] = 0; /* no matching */
    
    maxfct = 1;         /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */

    return 0;
}

int buildSolverMatCSRAndVec(int *rowOffsets, int *colIndices, double *vals, double *rhs, const int numRows)
{
    /* -------------------------------------------------------------------- */    
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */ 
    nnz = rowOffsets[numRows];

    for(i = 0; i < numRows+1; i++)
        rowOffsets[i] += 1;
    for(i = 0; i < nnz; i++)
        colIndices[i] += 1;

    // let the global vars point to the data
    n = numRows;
    ia = rowOffsets;
    ja = colIndices;
    a = vals;
    b = rhs;

    return 0;
}

void checkMatAndVec()
{
    // only for debugging
    pardiso_chkmatrix(&mtype, &n, a, ia, ja, &error);
    if (error != 0) {
        printf("\nERROR in consistency of matrix: %d", error);
        exit(1);
    }

    pardiso_chkvec (&n, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR  in right hand side: %d", error);
        exit(1);
    }
    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR right hand side: %d", error);
        exit(1);
    }
}

int initSolverParameters()
{
    phase = 13;
    iparm[3] = 61;
    return 0;
}

int setSolutionVector(double *sol)
{
    x = sol;
    return 0;
}

int SolverLinSolve()
{
    do
    {
        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs,
                iparm, &msglvl, b, x, &error, dparm);
    } while(dparm[33] > 1e-6);
   
    if (error != 0)
    {
        printf("\nERROR during solve?: %d\n", error);
        exit(1);
    }

    return 0;
}


#endif
