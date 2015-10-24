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
PetscErrorCode ComputeMatrix(KSP ksp, Mat J,Mat jac, void *ctx)
{
    PetscErrorCode ierr;
    PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,num, numi, numj, numk;
    PetscScalar    v[7],Hx,Hy,Hz,HyHzdHx,HxHzdHy,HxHydHz;
    MatStencil     row, col[7];
    DM             da;
    MatNullSpace   nullspace;

    PetscFunctionBeginUser;
    ierr    = KSPGetDM(ksp,&da);CHKERRQ(ierr);
    ierr    = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
    Hx      = GRID_LENGTH / (PetscReal)(mx-1);
    Hy      = GRID_LENGTH / (PetscReal)(my-1);
    Hz      = GRID_LENGTH / (PetscReal)(mz-1);
    HyHzdHx = Hy*Hz/Hx;
    HxHzdHy = Hx*Hz/Hy;
    HxHydHz = Hx*Hy/Hz;
    ierr    = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

    const double center[2]       = {GRID_LENGTH / 2., GRID_LENGTH / 2.};
    const double capillaryRadius = CAPILLARY_RADIUS;

    // extractor dimensions
    const double extractorInner  = EXTRACTOR_INNER_RADIUS;
    const double extractorOuter  = EXTRACTOR_OUTER_RADIUS;
    for (k=zs; k<zs+zm; k++)
    {
        PetscScalar z = k*Hz;
        for (j=ys; j<ys+ym; j++)
        {
            PetscScalar y = j*Hy;
            for (i=xs; i<xs+xm; i++)
            {
                row.i = i; row.j = j; row.k = k;
                num = 0; numi=0; numj=0; numk=0;

                PetscScalar ty = y - center[0];
                PetscScalar tz = z - center[1];

                PetscScalar rr = ty*ty + tz*tz;

                // if on boundary points
                if(   i == 0 || i == mx-1
                   || j == 0 || j == my-1
                   || k == 0 || k == mz-1)
                {

                    // if on X-Faces
                    if(i==0 || i == mx-1)
                    {
                        if(i == 0)
                        {
                            if( rr > capillaryRadius*capillaryRadius)
                            {
                                v[num]   = -1;
                                col[num].i = i+1;
                                col[num].j = j;
                                col[num].k = k;
                                //ierr = MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                                num++; numi++;
                            }
                        }
                        else if (i == mx-1)
                        {
                            if( (rr < extractorInner*extractorInner)
                                    ||
                                    (rr > extractorOuter*extractorOuter))
                            {
                                v[num]   = -1;
                                col[num].i = i-1;
                                col[num].j = j;
                                col[num].k = k;
                                //ierr = MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                                num++; numi++;
                            }
                        }

                        v[num] = 1;
                        col[num].i = i;
                        col[num].j = j;
                        col[num].k = k;
                        num++;
                        MatSetValuesStencil(jac, 1, &row, num, col, v, INSERT_VALUES);
                    } // end of if on X-Faces check

                    // if on other boundary faces
                    else
                    {
                        if (k == mz-1) {
                            v[num]     = -1;
                            col[num].i = i;
                            col[num].j = j;
                            col[num].k = k-1;
                            num++; numk++;
                        }
                        else if (k == 0) {
                            v[num]     = -1;
                            col[num].i = i;
                            col[num].j = j;
                            col[num].k = k+1;
                            num++; numk++;
                        }

                        if (j == my-1) {
                            v[num]     = -1;
                            col[num].i = i;
                            col[num].j = j-1;
                            col[num].k = k;
                            num++; numj++;
                        }
                        else if (j == 0) {
                            v[num]     = -1;
                            col[num].i = i;
                            col[num].j = j+1;
                            col[num].k = k;
                            num++; numj++;
                        }
                        //v[num]     = (PetscReal)(numk)*HxHydHz + (PetscReal)(numj)*HxHzdHy + (PetscReal)(numi)*HyHzdHx;
                        v[num]     = (PetscReal)(numk) + (PetscReal)(numj) + (PetscReal)(numi);
                        col[num].i = i;   col[num].j = j;   col[num].k = k;
                        num++;
                        ierr = MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    }
                }
                else
                {
                    //v[0] = -HxHydHz;                          col[0].i = i;   col[0].j = j;   col[0].k = k-1;
                    //v[1] = -HxHzdHy;                          col[1].i = i;   col[1].j = j-1; col[1].k = k;
                    //v[2] = -HyHzdHx;                          col[2].i = i-1; col[2].j = j;   col[2].k = k;
                    //v[3] = 2.0*(HxHydHz + HxHzdHy + HyHzdHx); col[3].i = i;   col[3].j = j;   col[3].k = k;
                    //v[4] = -HyHzdHx;                          col[4].i = i+1; col[4].j = j;   col[4].k = k;
                    //v[5] = -HxHzdHy;                          col[5].i = i;   col[5].j = j+1; col[5].k = k;
                    //v[6] = -HxHydHz;                          col[6].i = i;   col[6].j = j;   col[6].k = k+1;

                    v[0] = -1;                          col[0].i = i;   col[0].j = j;   col[0].k = k-1;
                    v[1] = -1;                          col[1].i = i;   col[1].j = j-1; col[1].k = k;
                    v[2] = -1;                          col[2].i = i-1; col[2].j = j;   col[2].k = k;
                    v[3] = 2.0*(3);                     col[3].i = i;   col[3].j = j;   col[3].k = k;
                    v[4] = -1;                          col[4].i = i+1; col[4].j = j;   col[4].k = k;
                    v[5] = -1;                          col[5].i = i;   col[5].j = j+1; col[5].k = k;
                    v[6] = -1;                          col[6].i = i;   col[6].j = j;   col[6].k = k+1;
                    ierr = MatSetValuesStencil(jac,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);
                }
            }
        }
    }
    ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    //ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
    //ierr = MatSetNullSpace(J,nullspace);CHKERRQ(ierr);
    //ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/*#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(KSP ksp, Mat J,Mat jac, void *ctx)
{
    PetscErrorCode ierr;
    PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,num, numi, numj, numk;
    PetscScalar    v[7],Hx,Hy,Hz,HyHzdHx,HxHzdHy,HxHydHz;
    MatStencil     row, col[7];
    DM             da;
    const int nr = (int)round(pow(numRows, 1./3));
    const int nrSq = nr*nr;

    PetscFunctionBeginUser;
    ierr    = KSPGetDM(ksp,&da);CHKERRQ(ierr);
    ierr    = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
    Hx      = GRID_LENGTH / (PetscReal)(mx-1);
    Hy      = GRID_LENGTH / (PetscReal)(my-1);
    Hz      = GRID_LENGTH / (PetscReal)(mz-1);
    HyHzdHx = Hy*Hz/Hx;
    HxHzdHy = Hx*Hz/Hy;
    HxHydHz = Hx*Hy/Hz;
    ierr    = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
    for (k=zs; k<zs+zm; k++)
    {
        for (j=ys; j<ys+ym; j++)
        {
            for (i=xs; i<xs+xm; i++)
            {
                row.i = i; row.j = j; row.k = k;
                int p = INDEX_1D(nr, i, j, k);
                int num = 0;

                int q;
                for(q = rows[p]; q < rows[p+1]; q++)
                {
                    int tempJ = cols[q];
                    int rk = tempJ / nrSq; tempJ = tempJ % nrSq;
                    int rj = tempJ / nr;
                    int ri = tempJ % nr;

                    col[num].i = ri;
                    col[num].j = rj;
                    col[num].k = rk;

                    v[num] = (PetscScalar)values[q];
                    num++;
                }

                MatSetValuesStencil(jac, 1, &row, num, col, v, INSERT_VALUES);
            } // end of i loop
        } // end of j loop
    } // end of k loop

    ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    //ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
    //ierr = MatSetNullSpace(J,nullspace);CHKERRQ(ierr);
    //ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
*/

//#undef __FUNCT__
//#define __FUNCT__ "ComputeMatrix"
//PetscErrorCode ComputeMatrix(KSP ksp, Mat jac, Mat A, void* ctx)
//{
//    PetscFunctionBeginUser;
//    PetscInt i,j;
//    PetscScalar v[7];
//    MatStencil row, col[7];
//
//    //Mat temp;
//    // build up the matrix
//    //MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, numRows, numRows, rows, cols, values, &temp);
//    //MatCopy(temp, A, DIFFERENT_NONZERO_PATTERN); 
//
//    const int nr = (int)round(pow(numRows, 1./3));
//    const int nrSq = nr*nr;
//
//    for(i = 0; i < numRows; i++)
//    {
//        PetscInt tempI = i;
//
//        PetscInt rk = tempI / nrSq; tempI = tempI % nrSq;
//        PetscInt rj = tempI / nr;
//        PetscInt ri = tempI % nr;
//
//        row.i = ri; row.j = rj; row.k = rk;
//
//        int count = 0;
//        for(j = rows[i]; j < rows[i+1]; j++)
//        {
//            PetscInt tempJ = cols[j];
//
//            rk = tempJ / nrSq; tempJ = tempJ % nrSq;
//            rj = tempJ / nr;
//            ri = tempJ % nr;
//
//            col[count].i = ri;
//            col[count].j = rj;
//            col[count].k = rk;
//
//            v[count] = (PetscScalar)values[j];
//            count++;
//        } // end of J loop where rowwise non zeros are accumulated into v
//
//        MatSetValuesStencil(A, 1, &row, count, col, v, INSERT_VALUES);
//    }
//
//    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
//    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
//
//    PetscFunctionReturn(0);
//
//} // end of ComputeMatrix

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(KSP ksp, Vec b, void* ctx)
{
    PetscFunctionBeginUser;
    PetscScalar ***array;
    DM dm;
    PetscInt i,j,k, mx, my, mz;
    PetscInt xs, ys, zs, xm, ym, zm;
    PetscScalar Hx, Hy, Hz;

    KSPGetDM(ksp, &dm);
    DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);
    Hx      = GRID_LENGTH / (PetscReal)(mx-1);
    Hy      = GRID_LENGTH / (PetscReal)(my-1);
    Hz      = GRID_LENGTH / (PetscReal)(mz-1);

    const double center[2]       = {GRID_LENGTH / 2., GRID_LENGTH / 2.};
    const double capillaryRadius = CAPILLARY_RADIUS;

    // extractor dimensions
    const double extractorInner  = EXTRACTOR_INNER_RADIUS;
    const double extractorOuter  = EXTRACTOR_OUTER_RADIUS;

    DMDAGetCorners(dm, &xs, &ys, &zs, &xm, &ym, &zm);
    DMDAVecGetArray(dm, b, &array);

    for (k=zs; k<zs+zm; k++)
    {
        PetscScalar z = k*Hz;
        for (j=ys; j<ys+ym; j++)
        {
            PetscScalar y = j*Hy;
            for (i=xs; i<xs+xm; i++)
            {
                // supply a default
                array[k][j][i] = 0;

                // if on boundary points
                if( i == 0 || i == mx-1 )
                {
                    PetscScalar ty = y - center[0];
                    PetscScalar tz = z - center[1];

                    PetscScalar rr = ty*ty + tz*tz;

                    if(i==0)
                    {
                        if( rr <= capillaryRadius*capillaryRadius)
                            array[k][j][i] = CAPILLARY_VOLTAGE;
                    }
                    else if (i == mx-1)
                    {
                        if( (rr >= extractorInner*extractorInner)
                                &&
                            (rr <= extractorOuter*extractorOuter) )
                            array[k][j][i] = EXTRACTOR_VOLTAGE;
                    }

                } // end of if on X-Faces check
            } // end of i loop
        } // end of j loop
    } // end of k loop

    DMDAVecRestoreArray(dm, b, &array);

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
    numRows = numberOfRows;

    // take the cube root to get the one side length
    PetscInt nr = (PetscInt)round(pow(numRows, 1./3));

    // create the DM object
    DMDACreate3d(
            PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
            DMDA_STENCIL_STAR,
            -nr, -nr, -nr,
            PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
            1, // dof
            1, // stencil width
            PETSC_NULL, PETSC_NULL, PETSC_NULL,
            &da
            );
    DMDASetUniformCoordinates(da, 0, GRID_LENGTH, 0, GRID_LENGTH, 0, GRID_LENGTH);
    //DMView(da, PETSC_VIEWER_STDOUT_SELF);
    //DMDASetInterpolationType(da, DMDA_Q1);
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

    //printf("Mat: \n");
    //Mat temp;
    //KSPGetOperators(ksp, &temp, NULL);
    //MatView(temp, PETSC_VIEWER_STDOUT_SELF);

    KSPGetIterationNumber(ksp, &it);
    //printf("solution: \n");
    KSPGetSolution(ksp, &x);
    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);

    //Vec b;
    //printf("rhs: \n");
    //KSPGetRhs(ksp, &b);
    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    return it;
}


//// TODO: optimize these calls? SEEMS wrong too!!
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
    VecDestroy(&x);
    //VecDestroy(&b);
    //MatDestroy(&A);
    DMDestroy(&da);

    free(rhsVals);
    free(values); free(cols); free(rows);
    PetscFinalize();
    PetscFunctionReturn(0);
}

#endif
