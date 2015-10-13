#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

// store relevant data as global for easier access
Mat A;
MatNullSpace nullspace;
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

    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullspace);

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

    MatSetNearNullSpace(A, nullspace);

    // fill in the b vector also
    for(i = 0; i < numRows; i++)
        rhsVals[i] = (PetscScalar)rhs[i];

    VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, numRows, rhsVals, &b);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    // allocate the same storage for x
    VecDuplicate(b, &x);
    VecCopy(b, x);

    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    PetscFunctionReturn(0);
} // end of buildSolverMatCSRAndVec

#undef __FUNCT__
#define __FUNCT__ "initSolverParameters"
PetscErrorCode initSolverParameters()
{
    PetscFunctionBegin;

    KSPSetOperators(ksp, A, A);
    KSPSetFromOptions(ksp);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SolverLinSolve"
PetscInt SolverLinSolve()
{
    PetscFunctionBegin;
    PetscInt it;
    KSPSolve(ksp, b, x);

    KSPGetIterationNumber(ksp, &it);
    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    return it;
}


// TODO: optimize these calls?
// seems inefficient
#undef __FUNCT__
#define __FUNCT__ "updateRHS"
PetscErrorCode updateRHS(const double *rhs, const int* indices, const int valCount)
{
    PetscFunctionBegin;

    PetscInt i;
    PetscInt ni = (PetscInt)(valCount);
    PetscScalar *y;

    VecGetArray(b, &y);

    // modify array
    for(i = 0; i < valCount; i++)
        y[ indices[i] ] = (PetscScalar)rhs[i];

    VecRestoreArray(b, &y);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "getSolution"
PetscErrorCode getSolution(double *sol)
{
    PetscFunctionBegin;

    PetscInt i,sizeX;
    VecGetLocalSize(x, &sizeX);

    PetscScalar *s;
    VecGetArray(x, &s);

    for(i = 0; i < sizeX; i++)
        sol[i] = (double)s[i];

    VecRestoreArray(x, &s);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "getRHS"
PetscErrorCode getRHS(double *rhs)
{
    PetscFunctionBegin;

    PetscInt i,sizeX;
    VecGetLocalSize(b, &sizeX);

    PetscScalar *s;
    VecGetArray(b, &s);

    for(i = 0; i < sizeX; i++)
        rhs[i] = (double)s[i];

    VecRestoreArray(b, &s);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SolverFinalize"
PetscErrorCode SolverFinalize()
{
    PetscFunctionBegin;

    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);

    MatNullSpaceDestroy(&nullspace);

    free(rhsVals);
    free(values); free(cols); free(rows);
    PetscFinalize();
    PetscFunctionReturn(0);
}

#endif
