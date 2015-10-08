
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
double   b[8], x[8], diag[8];
int nrhs;

int fun()
{
    /* Matrix data. */
    int    n = 8;
    int    ia[ 9] = { 0, 4, 7, 9, 11, 12, 15, 17, 20 };
    int    ja[20] = { 0,    2,       5, 6, 
                         1, 2,    4,
                            2,             7,
                               3,       6,
                         1,
                            2,       5,    7,
                         1,             6,
                            2,          6, 7 };
    double  a[20] = { 7.0,      1.0,           2.0, 7.0, 
                          -4.0, 8.0,      2.0,
                                1.0,                     5.0,
                                     7.0,           9.0,
                          -4.0,
                                7.0,           3.0,      8.0,
                           1.0,                    11.0,
                               -3.0,                2.0, 5.0 };

    int      nnz = ia[n];
    int      mtype = 11;        /* Real unsymmetric matrix */

    /* RHS and solution vectors. */
    nrhs = 1;          /* Number of right hand sides. */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 

/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters and initialize the solvers      */
/*     internal adress pointers. This is only necessary for the FIRST   */
/*     call of the PARDISO solver.                                      */
/* ---------------------------------------------------------------------*/
      
    error = 0;
    solver = 1;
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

    iparm[11] = 0; /* no scaling  */
    iparm[13] = 0; /* no matching */
    
    maxfct = 1;         /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */

/* -------------------------------------------------------------------- */    
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */ 
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }

    /* Set right hand side to i. */
    for (i = 0; i < n; i++) {
        b[i] = i;
        x[i] = 0;
    }

/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
    
    pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
    if (error != 0) {
        printf("\nERROR in consistency of matrix: %d", error);
        exit(1);
    }

/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

    pardiso_chkvec (&n, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR  in right hand side: %d", error);
        exit(1);
    }

/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */

    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR right hand side: %d", error);
        exit(1);
    }
 
/* -------------------------------------------------------------------- */    
/* ..  Prepare the preconditioner                                       */
/* -------------------------------------------------------------------- */

    //iparm[3] = 0;
    phase = 12;
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);
    //}
   
    if (error != 0) {
        printf("\nERROR during preconditioner: %d", error);
        exit(3);
    }


/* -------------------------------------------------------------------- */    
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */    
    iparm[3] = 81;
    phase = 33;

    //iparm[7] = 150;       /* Max numbers of iterative refinement steps. */

    //for(i = 0; i < 100; i++)
    //{
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);
    //}
   
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    printf("%lf\n", dparm[34]);

    printf("\nSolve completed ... ");
    printf("\nThe solution of the system is: ");
    for (i = 0; i < n; i++) {
        printf("\n x [%d] = % f", i, x[i] );
    }
    printf ("\n");

/* -------------------------------------------------------------------- */    
/* ... Inverse factorization.                                           */                                       
/* -------------------------------------------------------------------- */  
   
/* -------------------------------------------------------------------- */    
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */ 
    for (i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */ 
    phase = -1;                 /* Release internal memory. */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    return 0;
} 
