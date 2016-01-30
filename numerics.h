#ifndef NUMERICS_H
#define NUMERICS_H

#include "solvers/multigrid/mg_3d.h"
//double single_step_solve(double* grid, const int numNodes, const double sorOmega)
//{
//    int i, j, k;
//    double temp, norm = 0.;
//
//    // iterate across all points
//    for (i = 1; i < numNodes - 1; i++)
//    {
//        for(j = 1; j < numNodes - 1; j++)
//        {
//            for(k = 1; k < numNodes - 1; k++)
//            {
//                temp = GRID_1D(grid, i, j, k);
//
//                // obtain the value got by averaging
//                double val = (1./6) * (GRID_1D(grid,i+1,j,k) + GRID_1D(grid,i-1,j,k) + GRID_1D(grid,i,j-1,k) + GRID_1D(grid,i,j+1,k) + GRID_1D(grid,i,j,k-1) + GRID_1D(grid,i,j,k+1));
// 
//                // obtain the value with SOR
//                GRID_1D(grid, i,j,k) = val*sorOmega + (1-sorOmega)*temp;
//                //std::cout << grid[i][j][k].potential << " " << i << " " << j << " " << k << std::endl;
// 
//                // store the difference in values at same grid points in temp itself
//                // for norm calculation
//                temp = GRID_1D(grid,i,j,k) - temp;
//                norm += temp * temp;
//            }
//        }
//    }
//
//    return norm;
//}

void resetRHSInteriorPoints(double *rhs, GridInfo *gInfo)
{
    int i, j, k;
    const int N = gInfo->numNodes;
    const int NN = N*N;

    #pragma omp for schedule(static)
    for(i = 1; i < N-1; i++)
    {
        const int nni = NN*i;
        for(j = 1; j < N-1; j++)
        {
            const int nj = N*j;
            for(k = 1; k < N-1; k++)
            {
                int pos = nni + nj + k;
                rhs[pos] = 0;
            }
        }
    } // end of outermost i loop
}

void Solve(double toler, int maxIter, double *threadNorm)
{
    int iterCount = 0;
    int maxThreads = omp_get_max_threads();
    int tid = omp_get_thread_num();

    int i; double norm = 0;
    //threadNorm[tid] = SolverGetResidual();
    //#pragma omp barrier
    //#pragma omp single
    //{
    // let all threads do this so that their 'norm' is the same
    for(i = 0; i < maxThreads; i++)
    {
        // square and sum it to get the l2-norm
        // at the end
        norm += threadNorm[i]*threadNorm[i];
    }
    norm = sqrt(norm);
    //}
    double cmpNorm = toler*norm;

    #pragma omp single
    {
        printf("%10s %20s\n", "Iter_Count", "Norm");
        printf("%10d %20.8e\n", iterCount, norm);
    }

    //printf("Max threads: %d\n", maxThreads);

    //#pragma omp for schedule(static)
    //{
    for(iterCount = 1; norm > cmpNorm && (iterCount < maxIter); iterCount++)
    {
        threadNorm[tid] = SolverLinSolve();

        #pragma omp barrier // VERY IMPORTANT!!
        // let one thread calculate the actual norm
        #pragma omp single
        {
            int i;
            norm = 0;
            for(i = 0; i < maxThreads; i++)
            {
                // square and sum it to get the l2-norm
                // at the end
                norm += threadNorm[i]*threadNorm[i];
            }
            norm = sqrt(norm);

            if(!(iterCount % ITER_HEADER_INTERVAL))
                printf("%10s %20s\n", "Iter_Count", "Norm");

            if(!(iterCount % ITER_INTERVAL) )
                printf("%10d %20.8e\n", iterCount, norm);
        }
    } // end of iteration loop
    #pragma omp single
    {
    printf("%10d %20.8e\n", iterCount, norm);

    // if broke out due to MAX_ITER, warn so
    if(iterCount == maxIter)
        fprintf(stderr, "Stopped iterations due to MAX_ITER limit of %d\n", maxIter);
    }
    //} // end of OMP PRAGMA loop
}

// calculate the ElectricField once Node values are known
void calcElectricField(EField* ElectricField, const double* grid, GridInfo* gInfo)
{
    const int    N = gInfo->numNodes;
    const int   NN = N*N;
    const double invSpacing = gInfo->invSpacing;
    int i, j, k;

    // electric field is the NEGATIVE gradient of the potential field
    // so caching the multFactor we use just before filling the value
    const double multFactor = (-invSpacing * 0.5);

    // boundaries need to be handled separately
    #pragma omp for schedule(static)
    for(i = 0; i < N; i++)
    {
        bool isI_0 = (i == 0);
        bool isI_LEN = (i == (N - 1));
        bool isI_0_OR_LEN = isI_0 || isI_LEN;

        const int nni = NN*i;
        for(j = 0; j < N; j++)
        {
            bool isJ_0 = (j == 0);
            bool isJ_LEN = (j == (N - 1));
            bool isJ_0_OR_LEN = isJ_0 || isJ_LEN;

            const int nj = N*j;
            for(k = 0; k < N; k++)
            {
                bool isK_0 = (k == 0);
                bool isK_LEN = (k == (N - 1));
                bool isK_0_OR_LEN = isK_0 || isK_LEN;

                EField* elecField = &ElectricField[nni + nj + k];

                /**************************************/
                /****** X = 0 || X = GRID_LENGTH ******/
                double valX = 0.;
                if(isI_0_OR_LEN)
                {
                    if(isI_0)
                        valX = ( -grid[NN*(i+2) + nj + k]
                             + 4*grid[NN*(i+1) + nj + k]
                             - 3*grid[   nni   + nj + k] );

                    else if (isI_LEN)
                        valX = -( -grid[NN*(i-2) + nj + k]
                              + 4*grid[NN*(i-1) + nj + k]
                              - 3*grid[  nni    + nj + k] );
                }
                else
                    valX = (grid[NN*(i+1) + nj + k] - grid[NN*(i-1) +nj + k] );

                // by now ElectricField[i][j][k] needs to be filled
                elecField->components[0] = multFactor * valX;


                /**************************************/
                /****** Y = 0 || Y = GRID_LENGTH ******/
                double valY = 0;
                if(isJ_0_OR_LEN)
                {
                    if (isJ_0)
                        valY = ( -grid[nni + N*(j+2) + k]
                             + 4*grid[nni + N*(j+1) + k]
                             - 3*grid[nni +   nj    + k]);
                    else if (isJ_LEN)
                        valY = -(-grid[nni + N*(j-2) + k]
                             + 4*grid[nni + N*(j-1) + k]
                             - 3*grid[nni +   nj    + k] );
                }
                else
                    valY = (grid[nni + N*(j+1) + k] - grid[nni + N*(j-1) + k]);

                // by now ElectricField[i][j][k] needs to be filled
                elecField->components[1] = multFactor * valY;


                /**************************************/
                /****** Z = 0 || Z = GRID_LENGTH ******/
                double valZ = 0;
                if(isK_0_OR_LEN)
                {
                    if (isK_0)
                        valZ = ( -grid[nni + nj + k+2]
                             + 4*grid[nni + nj + k+1]
                             - 3*grid[nni + nj + k  ] );
                    else if (isK_LEN)
                        valZ = -( -grid[nni + nj + k-2]
                              + 4*grid[nni + nj + k-1]
                              - 3*grid[nni + nj + k  ] );
                }
                else
                    valZ = (grid[nni + nj + k+1] - grid[nni + nj + k-1]);

                // by now ElectricField[i][j][k] needs to be filled
                elecField->components[2] = multFactor * valZ;
            } // end of k loop
        } // end of j loop
    } // end of i loop
}

// setup the boundary conditions
//int setupBoundaryConditions(double* grid, GridInfo* gInfo, BoundaryNode *bNodes)
//{
//    const int numNodes           = gInfo->numNodes;
//    const double spacing         = gInfo->spacing;
//
//    const double center[2]       = {GRID_LENGTH / 2., GRID_LENGTH / 2.};
//    const double capillaryRadius = CAPILLARY_RADIUS;
//
//    // extractor dimensions
//    const double extractorInner  = EXTRACTOR_INNER_RADIUS;
//    const double extractorOuter  = EXTRACTOR_OUTER_RADIUS;
//
//    // store the extents in terms of indices for better comparisons
//    //const int extentIndices[2]   = { (int)( fabs(center[0] - capillaryRadius) * gInfo->invSpacing ), (int)( fabs( center[1] + capillaryRadius) * gInfo->invSpacing ) };
//
//    // store the voltage parameters
//    //const double capillaryVoltage = CAPILLARY_VOLTAGE;
//    //const double extractorVoltage = EXTRACTOR_VOLTAGE;
//
//    int i, j, k;
//    int nodeCount = 0;
//
//    // set boundary conditions
//    /*********************************************************/
//    // loop across all X-Z faces on the near and farther side
//    // j = 0 and j = SIZE_Y - 1 facesa
//    double x = 0., y = 0., z = 0.;
//
//    for(i = 0; i < numNodes; i++)
//    {
//        x = spacing * i;
//        for(k = 0; k < numNodes; k++)
//        {
//            //z = spacing * k;
//            // checking with x2 - 2y2 + z2
//            //y = 0.;
//            //GRID_1D(grid, i, 0, k) =TEST_FUNCTION;
//
//            // enforce Neumann boundary nodes
//            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, 0, k);
//            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, 1, k);
//            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, 2, k);
//            nodeCount++;
//
//            //y = GRID_LENGTH;
//            //GRID_1D(grid, i, numNodes-1, k).potential=TEST_FUNCTION;
//
//            // Y = GRID_LENGTH
//            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, numNodes-1, k);
//            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, numNodes-2, k);
//            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, numNodes-3, k);
//            nodeCount++;
//        }
//
//        // on X-Y faces
//        for(j = 0; j < numNodes; j++)
//        {
//            //y = spacing*j;
//            //z = 0.;
//            //GRID_1D(grid, i, j, 0).potential=TEST_FUNCTION;
//
//            // enforce Neumann boundary nodes
//            // Z = 0
//            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, j, 0);
//            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, j, 1);
//            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, j, 2);
//            nodeCount++;
//
//            //z = GRID_LENGTH;
//            //GRID_1D(grid, i, j, numNodes-1).potential=TEST_FUNCTION;
//
//            // Z = GRID_LENGTH
//            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, j, numNodes-1);
//            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, j, numNodes-2);
//            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, j, numNodes-3);
//            nodeCount++;
//        }
//    }
//
//    // on Y-Z faces
//    for(j = 0; j < numNodes; j++)
//    {
//        double ty = (spacing*j - center[0]);
//
//        for(k = 0; k < numNodes; k++)
//        {
//            double tz = (spacing*k - center[1]);
//            double sumSqs = ty*ty + tz*tz;
//
//            // CAPILLARY side
//            // we need points INSIDE the capillary for Neumann BC
//            if( sumSqs < capillaryRadius*capillaryRadius )
//            {
//                GRID_1D(grid, 0, j, k) = CAPILLARY_VOLTAGE;
//            }
//            else
//            {
//                bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, 0, j, k);
//                bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, 1, j, k);
//                bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, 2, j, k);
//
//                nodeCount++;
//            }
//
//            // EXTRACTOR side
//            if( (sumSqs > extractorInner*extractorInner) && (sumSqs < extractorOuter*extractorOuter ) )
//            {
//                GRID_1D(grid, numNodes-1, j, k) = EXTRACTOR_VOLTAGE;
//            }
//            else
//            {
//                bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, numNodes-1, j, k);
//                bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, numNodes-2, j, k);
//                bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, numNodes-3, j, k);
//
//                nodeCount++;
//            }
//
//        }
//    }
//
//    return nodeCount;
//    //std::cout << grid[1][1][0].potential << " " << grid[2][3][4].potential << std::endl;
//}

void enforceNeumannBC(BoundaryNode* bNodes, const int nodeCount)
{
    int i;

    for(i = 0; i < nodeCount; i++)
        *(bNodes[i].bndryNodes[0]) = (1./3) * (4 * *(bNodes[i].bndryNodes[1]) - *(bNodes[i].bndryNodes[2]));

}
#endif
