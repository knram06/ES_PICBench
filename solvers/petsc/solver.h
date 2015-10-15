#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscmat.h>
#include <petscvec.h>

#include <math.h>

// store relevant data as global for easier access
DM da;
Vec x;
KSP ksp;
PetscInt numRows;

// buffer pointers which need to freed AFTER destroying arrays
PetscInt *rows, *cols;
PetscScalar *values, *rhsVals;

#undef __FUNCT__
#define __FUNCT__ "SolverInitialize"
PetscErrorCode SolverInitialize(int *argc, char ***argv)
{
    PetscFunctionBegin;
    PetscInitialize(argc, argv, (char*)"petsc_opts.cfg", NULL);
    //MatCreate(PETSC_COMM_WORLD, &A);
    //VecCreate(PETSC_COMM_WORLD, &b);
    VecCreate(PETSC_COMM_WORLD, &x);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(KSP ksp, Mat J, Mat A, void* ctx)
{
    PetscFunctionBeginUser;

    // TODO: refine this approch later
    // Just use createseq to copy into local
    Mat temp;
    // AVOID THIS! - repositions the array to a new location, not the previous location
    MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, numRows, numRows, rows, cols, values, &temp);

    // now copy this over into the relevant location
    //MatCopy(temp, A, DIFFERENT_NONZERO_PATTERN);

    // TODO: Use MatSetValue and set using the CSR format
    PetscInt i, j;
    for(i = 0; i < numRows; i++)
    {
        for(j = rows[i]; j < rows[i+1]; j++)
            MatSetValue(A, i, cols[j], values[j], INSERT_VALUES);
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    PetscFunctionReturn(0);

} // end of ComputeMatrix

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(KSP ksp, Vec b, void* ctx)
{
    PetscFunctionBeginUser;
    PetscScalar *array;
    PetscInt i;

    DMDAVecGetArray(da, b, &array);
    for(i = 0; i < numRows; i++)
        array[i] = rhsVals[i];

    DMDAVecRestoreArray(da, b, &array);

    //VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, numRows, rhsVals, &b);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    //VecView(b, PETSC_VIEWER_STDOUT_SELF);

    // allocate the same storage for x
    VecDuplicate(b, &x);
    VecCopy(b, x);

    PetscFunctionReturn(0);
} // end of ComputeRHS

#undef __FUNCT__
#define __FUNCT__ "buildSolverMatCSRAndVec"
PetscErrorCode buildSolverMatCSRAndVec(const int *rowOffsets, const int *colIndices, const double *vals, const double *rhs, const int numberOfRows)
{
    PetscFunctionBegin;

    // update the global variable
    numRows = numberOfRows;

    // take the cube root to get the one side length
    //PetscInt nr = (PetscInt)round(pow(numRows, 1./3));

    // the last entry in the rowOffsets will actually be the total nonzeros
    // so just use that
    const int totalNonZeros = rowOffsets[numRows]; // takes care of one off indexing
    const int numRowsOffset = numRows+1;

    // create the DM object
    DMDACreate1d(
            PETSC_COMM_WORLD,
            DM_BOUNDARY_NONE,
            numRows,
            1, // dof
            1, // stencil width
            PETSC_NULL,
            &da
            );
    //DMDASetInterpolationType(da, DMDA_Q1);
    KSPSetDM(ksp, da);

    KSPSetComputeRHS(ksp, ComputeRHS, NULL);
    KSPSetComputeOperators(ksp, ComputeMatrix, NULL);

    //VecSetFromOptions(b);

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

    for(i = 0; i < numRows; i++)
        rhsVals[i] = (PetscScalar)rhs[i];

    PetscFunctionReturn(0);
} // end of buildSolverMatCSRAndVec

#undef __FUNCT__
#define __FUNCT__ "initSolverParameters"
PetscErrorCode initSolverParameters()
{
    PetscFunctionBegin;
    KSPSetFromOptions(ksp);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SolverLinSolve"
PetscInt SolverLinSolve()
{
    PetscFunctionBegin;

    Vec t;
    Mat m;

    PetscInt it;
    KSPSolve(ksp, NULL, NULL);

    KSPGetRhs(ksp, &t);
    //VecView(t, PETSC_VIEWER_STDOUT_SELF);

    KSPGetOperators(ksp, &m, NULL);
    //MatView(m, PETSC_VIEWER_STDOUT_SELF);

    KSPGetIterationNumber(ksp, &it);
    KSPGetSolution(ksp, &x);
    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    return it;
}


//// TODO: optimize these calls?
//// seems inefficient
//#undef __FUNCT__
//#define __FUNCT__ "updateRHS"
//PetscErrorCode updateRHS(const double *rhs, const int* indices, const int valCount)
//{
//    PetscFunctionBegin;
//
//    PetscInt i;
//    PetscInt ni = (PetscInt)(valCount);
//    PetscScalar *y;
//
//    VecGetArray(b, &y);
//
//    // modify array
//    for(i = 0; i < valCount; i++)
//        y[ indices[i] ] = (PetscScalar)rhs[i];
//
//    VecRestoreArray(b, &y);
//    PetscFunctionReturn(0);
//}

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

//#undef __FUNCT__
//#define __FUNCT__ "getRHS"
//PetscErrorCode getRHS(double *rhs)
//{
//    PetscFunctionBegin;
//
//    PetscInt i,sizeX;
//    VecGetLocalSize(b, &sizeX);
//
//    PetscScalar *s;
//    VecGetArray(b, &s);
//
//    for(i = 0; i < sizeX; i++)
//        rhs[i] = (double)s[i];
//
//    VecRestoreArray(b, &s);
//    PetscFunctionReturn(0);
//}

#undef __FUNCT__
#define __FUNCT__ "SolverFinalize"
PetscErrorCode SolverFinalize()
{
    PetscFunctionBegin;

    KSPDestroy(&ksp);
    DMDestroy(&da);
    VecDestroy(&x);
    //VecDestroy(&b);
    //MatDestroy(&A);

    free(rhsVals);
    free(values); free(cols); free(rows);
    PetscFinalize();
    PetscFunctionReturn(0);
}

#endif
