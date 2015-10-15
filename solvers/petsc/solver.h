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
int numRows;

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

    Mat temp;
    // build up the matrix
    MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, numRows, numRows, rows, cols, values, &temp);

    MatCopy(temp, A, DIFFERENT_NONZERO_PATTERN); 

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    PetscFunctionReturn(0);

} // end of ComputeMatrix

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(KSP ksp, Vec b, void* ctx)
{
    PetscFunctionBeginUser;
    PetscScalar ***array;
    PetscInt i,j,k;
    PetscInt xs, ys, zs, xm, ym, zm;

    DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
    DMDAVecGetArray(da, b, &array);
    for(k = zs; k < zs+zm; k++)
    {
        for(j = ys; j < ys + ym; j++)
        {
            for(i = xs; i < xs+xm; i++)
            {
                array[k][j][i] = GRID_1D(rhsVals, i, j, k);
            }
        } // end of j loop
    } // end of k loop
    DMDAVecRestoreArray(da, b, &array);

    //VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, numRows, rhsVals, &b);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

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
    numRows = numberOfRows;

    // take the cube root to get the one side length
    PetscInt nr = (PetscInt)round(pow(numRows, 1./3));

    // create the DM object
    DMDACreate3d(
            PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
            DMDA_STENCIL_STAR,
            nr, nr, nr,
            PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
            1, // dof
            1, // stencil width
            PETSC_NULL, PETSC_NULL, PETSC_NULL,
            &da
            );
    //DMDASetInterpolationType(da, DMDA_Q0);
    KSPSetDM(ksp, da);

    KSPSetComputeRHS(ksp, ComputeRHS, (void*)&numRows);
    KSPSetComputeOperators(ksp, ComputeMatrix, (void*)&numRows);

    // the last entry in the rowOffsets will actually be the total nonzeros
    // so just use that
    const int totalNonZeros = rowOffsets[numRows]; // takes care of one off indexing
    const int numRowsOffset = numRows+1;

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
    PetscInt it;
    KSPSolve(ksp, NULL, NULL);

    Mat temp;
    KSPGetOperators(ksp, &temp, NULL);
    MatView(temp, PETSC_VIEWER_STDOUT_SELF);

    KSPGetIterationNumber(ksp, &it);
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
