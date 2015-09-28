#ifndef SOLVER_H
#define SOLVER_H

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

// store relevant data as global for better access
Mat A;
Vec b, x;
KSP ksp;

// buffer pointers which need to freed AFTER destroying arrays
PetscInt *rows, *cols;
PetscScalar *values, *rhsVals;

#undef __FUNCT__
#define __FUNCT__ "SolverInitialize"
PetscErrorCode SolverInitialize(int *argc, char ***argv)
{
    PetscFunctionBegin;
    PetscInitialize(argc, argv, (char*)"petsc_opts.cfg", NULL);
    MatCreate(PETSC_COMM_WORLD, &A);
    VecCreate(PETSC_COMM_WORLD, &b);
    VecCreate(PETSC_COMM_WORLD, &x);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "buildSolverMatCSRAndVec"
PetscErrorCode buildSolverMatCSRAndVec(const int *rowOffsets, const int *colIndices, const double *vals, const double *rhs, const int numRows)
{
    PetscFunctionBegin;

    // the last entry in the rowOffsets will actually be the total nonzeros
    // so just use that
    const int totalNonZeros = rowOffsets[numRows]; // takes care of one off indexing
    const int numRowsOffset = numRows+1;

    VecSetFromOptions(b);

    // declare solver specific stuff here
    rows = malloc(sizeof(PetscInt) * numRowsOffset);
    cols = malloc(sizeof(PetscInt) * totalNonZeros);
    values = malloc(sizeof(PetscScalar) * totalNonZeros);
    rhsVals = malloc(sizeof(PetscScalar) * numRows);

    // now loop through and copy over stuff - costly?
    int i;
    for(i = 0; i < numRowsOffset; i++)
        rows[i] = (PetscInt)(rowOffsets[i]);

    for(i = 0; i < totalNonZeros; i++)
    {
        cols[i]   = (PetscInt)(colIndices[i]);
        values[i] = (PetscScalar)(vals[i]);
    }

    // build up the matrix
    MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, (PetscInt)numRows, (PetscInt)numRows, rows, cols, values, &A);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // fill in the b vector also
    for(i = 0; i < numRows; i++)
        rhsVals[i] = (PetscScalar)rhs[i];

    VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, numRows, rhsVals, &b);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    // allocate the same storage for x
    VecDuplicate(b, &x);

    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    PetscFunctionReturn(0);
} // end of buildSolverMatCSRAndVec

#undef __FUNCT__
#define __FUNCT__ "initSolverParameters"
PetscErrorCode initSolverParameters()
{
    PetscFunctionBegin;

    KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);
    KSPSetFromOptions(ksp);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SolverLinSolve"
PetscErrorCode SolverLinSolve()
{
    PetscFunctionBegin;
    KSPSolve(ksp, b, x);
    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "getSolution"
PetscErrorCode getSolution(double *sol)
{
    PetscFunctionBegin;

    PetscInt sizeX;
    VecGetLocalSize(x, &sizeX);

    PetscScalar *s;
    VecGetArray(x, &s);

    int i;
    for(i = 0; i < sizeX; i++)
        sol[i] = (double)s[i];

    VecRestoreArray(x, &s);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SolverFinalize"
PetscErrorCode SolverFinalize()
{
    PetscFunctionBegin;
    KSPDestroy(&ksp);
    //VecDestroy(&x);
    VecDestroy(&b);
    free(rhsVals);

    MatDestroy(&A);
    free(values); free(cols); free(rows);
    PetscFinalize();
    PetscFunctionReturn(0);
}

#endif
