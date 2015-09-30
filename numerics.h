#ifndef NUMERICS_H
#define NUMERICS_H

#include <stdlib.h>
#include <assert.h>

typedef struct __csr_t
{
    double *mat;
    int *rowOffsets;
    int *colIndices;
    int numRows;
} MatCSR;

// simple typedef to hold col index
// and corresponding value
typedef struct __col_index_t
{
    int index;
    double val;
} ColVal;

// compare function
int cv_compare(const void *a, const void *b)
{
    ColVal *cv_a = (ColVal*)a;
    ColVal *cv_b = (ColVal*)b;
    if((*cv_a).index < (*cv_b).index)
        return -1;
    else if((*cv_a).index == (*cv_b).index)
        return 0;
    else
        return 1;
}

void allocCSRForm(MatCSR* matcsr, const int numRows, const int maxNZPerRow)
{
    const int maxNZ = numRows * maxNZPerRow;

    // have matrix prefilled with zeros
    matcsr->mat = calloc(maxNZ, sizeof(double));

    matcsr->rowOffsets = malloc(sizeof(int) * (numRows+1));

    // first starting point has to be zero
    matcsr->rowOffsets[0] = 0; 
    matcsr->colIndices = malloc(sizeof(int) * maxNZ);
    matcsr->numRows = numRows;
}

void deallocCSRForm(MatCSR *matcsr)
{
    free(matcsr->mat);
    free(matcsr->rowOffsets);
    free(matcsr->colIndices);
}

void buildSparseMatAndRHSVec(int *rowOffsets, int *colIndices, double *vals,
                             double *rhs, GridInfo* gridInfo)
{
    int i, j, k;
    const int numNodes = gridInfo->numNodes;

    // for convenient indexing
    int *rowPtr = &rowOffsets[1];

    const int selfCoeff = -6;
    const int nonSelfCoeff = 1;

    int runningIndex = 0;
    for(i = 0; i < numNodes; i++)
    {
        bool isI_0 =   (i == 0);
        bool isI_LEN = (i == (numNodes-1));
        bool isI_0_OR_LEN = isI_0 || isI_LEN;

        for(j = 0; j < numNodes; j++)
        {
            bool isJ_0 =   (j == 0);
            bool isJ_LEN = (j == (numNodes-1));
            bool isJ_0_OR_LEN = isJ_0 || isJ_LEN;

            for(k = 0; k < numNodes; k++)
            {
                bool isK_0 =   (k == 0);
                bool isK_LEN = (k == (numNodes-1));
                bool isK_0_OR_LEN = isK_0 || isK_LEN;

                const int I = INDEX_1D(numNodes, i, j, k);

                // fill in with some sentinel values
                // assuming repeated values 9 entries could be made here
                // i.e. corner points, self coefficients will add up
                //int cols[9]      = {-1, -1, -1,
                //                    -1, -1, -1,
                //                    -1, -1, -1};
                //double values[9] = {0, 0, 0,
                //                    0, 0, 0,
                //                    0, 0, 0};
                ColVal cvals[9];
                int elementsInserted = 0;

                // supply a default
                //rhs[I] = 0;

                // if on boundary
                if(isI_0_OR_LEN || isJ_0_OR_LEN || isK_0_OR_LEN)
                {
                    // if on Dirichlet boundary
                    if(isI_0_OR_LEN)
                    {
                        // for now keep all dirichlet
                        cvals[elementsInserted].index = I;
                        cvals[elementsInserted].val = 1;
                        elementsInserted++;

                        if(isI_0)
                            rhs[I] = 0; // (say)
                        else
                            rhs[I] = 10;
                    }
                    // on Neumann boundaries
                    // ORDERING IS VERY IMPORTANT - here
                    // other option is to sort cols and vals
                    // and then INSERT
                    else
                    {
                        // if on J = 0 or J = LEN boundaries
                        if(isJ_0)
                        {
                            // at point i,j,k
                            cvals[elementsInserted].index = I;
                            cvals[elementsInserted].val = 3;
                            elementsInserted++;

                            // at point i,j+1,k
                            cvals[elementsInserted].index = INDEX_1D(numNodes, i, j+1, k);
                            cvals[elementsInserted].val = -4;
                            elementsInserted++;

                            // at point i,j+2,k
                            cvals[elementsInserted].index = INDEX_1D(numNodes, i, j+2, k);
                            cvals[elementsInserted].val = 1;
                            elementsInserted++;
                        }
                        else if(isJ_LEN)
                        {
                            // at point i,j-2,k
                            cvals[elementsInserted].index = INDEX_1D(numNodes, i, j-2, k);
                            cvals[elementsInserted].val = 1;
                            elementsInserted++;

                            // at point i,j-1,k
                            cvals[elementsInserted].index = INDEX_1D(numNodes, i, j-1, k);
                            cvals[elementsInserted].val = -4;
                            elementsInserted++;

                            // at point i,j,k
                            cvals[elementsInserted].index = I;
                            cvals[elementsInserted].val = 3;
                            elementsInserted++;
                        }
                        if(isK_0)
                        {
                            // at point i,j,k
                            cvals[elementsInserted].index = I;
                            cvals[elementsInserted].val = 3;
                            elementsInserted++;

                            // at point i,j,k+1
                            cvals[elementsInserted].index = INDEX_1D(numNodes, i, j, k+1);
                            cvals[elementsInserted].val = -4;
                            elementsInserted++;

                            // at point i,j,k+2
                            cvals[elementsInserted].index = INDEX_1D(numNodes, i, j, k+2);
                            cvals[elementsInserted].val = 1;
                            elementsInserted++;
                        }
                        else if(isK_LEN)
                        {
                            // at point i,j,k-2
                            cvals[elementsInserted].index = INDEX_1D(numNodes, i, j, k-2);
                            cvals[elementsInserted].val = 1;
                            elementsInserted++;

                            // at point i,j,k-1
                            cvals[elementsInserted].index = INDEX_1D(numNodes, i, j, k-1);
                            cvals[elementsInserted].val = -4;
                            elementsInserted++;

                            // at point i,j,k
                            cvals[elementsInserted].index = I;
                            cvals[elementsInserted].val = 3;
                            elementsInserted++;
                        }
                        // if on K = 0 or K = LEN boundaries
                    } // end of else check, i.e. if node is on Neumann boundaries
                } // end of if check, i.e. if node is on any boundary

                // else
                // fill in coefficients for interior point behaviour
                else
                {
                    // all these points are guaranteed to be in the
                    // interior, so we can fill in the coefficients as needed
                    int J = 0;

                    // explicitly writing this - to make it look like
                    // LOOP UNROLLING?
                    // We expect this convention to be in increasing order
                    // of I, since we define I as k+n*j+n*n*i)
                    // CONVENTION: i-1, j-1, k-1, k, k+1, j+1, i+1
                    J = INDEX_1D(numNodes, i-1, j  , k  );
                    cvals[elementsInserted].index = J;
                    cvals[elementsInserted].val = nonSelfCoeff;
                    elementsInserted++;

                    J = INDEX_1D(numNodes, i  , j-1, k  );
                    cvals[elementsInserted].index = J;
                    cvals[elementsInserted].val = nonSelfCoeff;
                    elementsInserted++;

                    J = INDEX_1D(numNodes, i  , j  , k-1);
                    cvals[elementsInserted].index = J;
                    cvals[elementsInserted].val = nonSelfCoeff;
                    elementsInserted++;

                    /********** SELF COEFF ***************/
                    /*************************************/
                    //J = INDEX_1D(numNodes, i  , j  , k  );
                    cvals[elementsInserted].index = I;
                    cvals[elementsInserted].val = selfCoeff;
                    elementsInserted++;
                    /*************************************/
                    /*************************************/

                    J = INDEX_1D(numNodes, i  , j  , k+1);
                    cvals[elementsInserted].index = J;
                    cvals[elementsInserted].val = nonSelfCoeff;
                    elementsInserted++;

                    J = INDEX_1D(numNodes, i  , j+1, k  );
                    cvals[elementsInserted].index = J;
                    cvals[elementsInserted].val = nonSelfCoeff;
                    elementsInserted++;

                    J = INDEX_1D(numNodes, i+1, j  , k  );
                    cvals[elementsInserted].index = J;
                    cvals[elementsInserted].val = nonSelfCoeff;
                    elementsInserted++;

                } // end of else condition - i.e. at interior points

                // now sort based on column index - most external solvers
                // might rely on the column indices being sorted - at
                // least PETSc does http://lists.mcs.anl.gov/pipermail/petsc-users/2012-August/014856.html
                qsort(cvals, elementsInserted, sizeof(ColVal), cv_compare);

                // now all the repeated col indices (which should be I) will occur
                // next to each other
                int p;
                int marker_I = -1; // marker to note where to write to later
                double sumCoeff_I = 0.;
                for(p = 0; p < elementsInserted; p++)
                {
                    if (cvals[p].index != I)
                    {
                        colIndices[runningIndex] = cvals[p].index;
                        vals[runningIndex] = cvals[p].val;
                        runningIndex++;
                    }
                    else
                    {
                        if(marker_I == -1)
                        {
                            marker_I = runningIndex;
                            runningIndex++;
                        }
                        sumCoeff_I += cvals[p].val;
                    }
                } // end of for loop for p
                
                // now write out the accumulated value in the stored index
                assert(marker_I != -1);
                colIndices[marker_I] = I;
                vals[marker_I] = sumCoeff_I;

                // update the row offsets vector
                // using the one off adjusted pointer
                rowPtr[I] = runningIndex;
            } // end of k loop
        } // end of j loop
    } // end of i loop

}

double single_step_solve(Node* grid, const int numNodes, const double sorOmega)
{
    int i, j, k;
    double temp, norm = 0.;

    // iterate across all points
    for (i = 1; i < numNodes - 1; i++)
    {
        for(j = 1; j < numNodes - 1; j++)
        {
            for(k = 1; k < numNodes - 1; k++)
            {
                temp = GRID_1D(grid, i, j, k).potential;

                // obtain the value got by averaging
                double val = (1./6) * (GRID_1D(grid,i+1,j,k).potential + GRID_1D(grid,i-1,j,k).potential + GRID_1D(grid,i,j-1,k).potential + GRID_1D(grid,i,j+1,k).potential + GRID_1D(grid,i,j,k-1).potential + GRID_1D(grid,i,j,k+1).potential);
 
                // obtain the value with SOR
                GRID_1D(grid, i,j,k).potential = val*sorOmega + (1-sorOmega)*temp;
                //std::cout << grid[i][j][k].potential << " " << i << " " << j << " " << k << std::endl;
 
                // store the difference in values at same grid points in temp itself
                // for norm calculation
                temp = GRID_1D(grid,i,j,k).potential - temp;
                norm += temp * temp;
            }
        }
    }

    return norm;
}

// calculate the ElectricField once Node values are known
void calcElectricField(EField* ElectricField, double* grid, GridInfo* gInfo)
{
    const int   numNodes = gInfo->numNodes;
    const double invSpacing = gInfo->invSpacing;
    int i, j, k;

    // electric field is the NEGATIVE gradient of the potential field
    // so caching the multFactor we use just before filling the value
    const double multFactor = (-invSpacing / 2);

    // boundaries need to be handled separately
    for(i = 0; i < numNodes; i++)
    {
        bool isI_0 = (i == 0);
        bool isI_LEN = (i == (numNodes - 1));
        bool isI_0_OR_LEN = isI_0 || isI_LEN;

        for(j = 0; j < numNodes; j++)
        {
            bool isJ_0 = (j == 0);
            bool isJ_LEN = (j == (numNodes - 1));
            bool isJ_0_OR_LEN = isJ_0 || isJ_LEN;

            for(k = 0; k < numNodes; k++)
            {
                bool isK_0 = (k == 0);
                bool isK_LEN = (k == (numNodes - 1));
                bool isK_0_OR_LEN = isK_0 || isK_LEN;

                EField* elecField = &GRID_1D(ElectricField, i, j, k);

                /**************************************/
                /****** X = 0 || X = GRID_LENGTH ******/
                double val = 0.;
                if(isI_0_OR_LEN)
                {
                    if(isI_0)
                        val = ( -GRID_1D(grid,i+2,j,k)
                             + 4*GRID_1D(grid,i+1,j,k)
                             - 3*GRID_1D(grid,i  ,j,k) );

                    else if (isI_LEN)
                        val = -( -GRID_1D(grid,i-2,j,k) 
                              + 4*GRID_1D(grid,i-1,j,k)
                              - 3*GRID_1D(grid,i  ,j,k) );
                }
                else
                    val = (GRID_1D(grid,i+1,j,k) - GRID_1D(grid,i-1,j,k));

                // by now ElectricField[i][j][k] needs to be filled
                elecField->components[0] = multFactor * val;


                /**************************************/
                /****** Y = 0 || Y = GRID_LENGTH ******/
                if(isJ_0_OR_LEN)
                {
                    if (isJ_0)
                        val = ( -GRID_1D(grid,i,j+2,k)
                             + 4*GRID_1D(grid,i,j+1,k)
                             - 3*GRID_1D(grid,i,j  ,k) );
                    else if (isJ_LEN)
                        val = -(-GRID_1D(grid,i,j-2,k) 
                             + 4*GRID_1D(grid,i,j-1,k)
                             - 3*GRID_1D(grid,i,j  ,k) );
                }
                else
                    val = (GRID_1D(grid,i,j+1,k) - GRID_1D(grid,i,j-1,k));

                // by now ElectricField[i][j][k] needs to be filled
                elecField->components[1] = multFactor * val;


                /**************************************/
                /****** Z = 0 || Z = GRID_LENGTH ******/
                if(isK_0_OR_LEN)
                {
                    if (isK_0)
                        val = ( -GRID_1D(grid,i,j,k+2)
                             + 4*GRID_1D(grid,i,j,k+1)
                             - 3*GRID_1D(grid,i,j,k  ) );
                    else if (isK_LEN)
                        val = -( -GRID_1D(grid,i,j,k-2) 
                              + 4*GRID_1D(grid,i,j,k-1)
                              - 3*GRID_1D(grid,i,j,k  ) );
                }
                else
                    val = (GRID_1D(grid,i,j,k+1) - GRID_1D(grid,i,j,k-1));

                // by now ElectricField[i][j][k] needs to be filled
                elecField->components[2] = multFactor * val;

            } // end of k loop
        } // end of j loop
    } // end of i loop
}

// setup the boundary conditions
int setupBoundaryConditions(Node* grid, GridInfo* gInfo, BoundaryNode *bNodes)
{
    const int numNodes           = gInfo->numNodes;
    const double spacing         = gInfo->spacing;

    const double center[2]       = {GRID_LENGTH / 2., GRID_LENGTH / 2.};
    const double capillaryRadius = CAPILLARY_RADIUS;

    // extractor dimensions
    const double extractorInner  = EXTRACTOR_INNER_RADIUS;
    const double extractorOuter  = EXTRACTOR_OUTER_RADIUS;

    // store the extents in terms of indices for better comparisons
    //const int extentIndices[2]   = { (int)( fabs(center[0] - capillaryRadius) * gInfo->invSpacing ), (int)( fabs( center[1] + capillaryRadius) * gInfo->invSpacing ) };

    // store the voltage parameters
    //const double capillaryVoltage = CAPILLARY_VOLTAGE;
    //const double extractorVoltage = EXTRACTOR_VOLTAGE;

    int i, j, k;
    int nodeCount = 0;

    // set boundary conditions
    /*********************************************************/
    // loop across all X-Z faces on the near and farther side
    // j = 0 and j = SIZE_Y - 1 facesa
    double x = 0., y = 0., z = 0.;

    for(i = 0; i < numNodes; i++)
    {
        x = spacing * i;
        for(k = 0; k < numNodes; k++)
        {
            z = spacing * k;
            // checking with x2 - 2y2 + z2
            //grid[i][0][k] = x*x + z*z;          // y is zero here
            y = 0.;
            GRID_1D(grid, i, 0, k).potential=TEST_FUNCTION;

            // enforce Neumann boundary nodes
            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, 0, k);
            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, 1, k);
            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, 2, k);
            nodeCount++;

            y = GRID_LENGTH;
            GRID_1D(grid, i, numNodes-1, k).potential=TEST_FUNCTION;

            // Y = GRID_LENGTH
            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, numNodes-1, k);
            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, numNodes-2, k);
            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, numNodes-3, k);
            nodeCount++;
        }

        // on X-Y faces
        for(j = 0; j < numNodes; j++)
        {
            y = spacing*j;
            z = 0.;
            GRID_1D(grid, i, j, 0).potential=TEST_FUNCTION;

            // enforce Neumann boundary nodes
            // Z = 0
            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, j, 0);
            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, j, 1);
            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, j, 2);
            nodeCount++;

            z = GRID_LENGTH;
            GRID_1D(grid, i, j, numNodes-1).potential=TEST_FUNCTION;

            // Z = GRID_LENGTH
            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, j, numNodes-1);
            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, j, numNodes-2);
            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, j, numNodes-3);
            nodeCount++;
        }
    }

    // on Y-Z faces
    for(j = 0; j < numNodes; j++)
    {
        double ty = (spacing*j - center[0]);

        for(k = 0; k < numNodes; k++)
        {
            double tz = (spacing*k - center[1]);
            double sumSqs = ty*ty + tz*tz;

            // CAPILLARY side
            // we need points INSIDE the capillary for Neumann BC
            if( sumSqs < capillaryRadius*capillaryRadius )
            {
                GRID_1D(grid, 0, j, k).potential = CAPILLARY_VOLTAGE;
            }
            else
            {
                bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, 0, j, k);
                bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, 1, j, k);
                bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, 2, j, k);

                nodeCount++;
            }

            // EXTRACTOR side
            if( (sumSqs > extractorInner*extractorInner) && (sumSqs < extractorOuter*extractorOuter ) )
            {
                GRID_1D(grid, numNodes-1, j, k).potential = EXTRACTOR_VOLTAGE;
            }
            else
            {
                bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, numNodes-1, j, k);
                bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, numNodes-2, j, k);
                bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, numNodes-3, j, k);

                nodeCount++;
            }

        }
    }

    return nodeCount;
    //std::cout << grid[1][1][0].potential << " " << grid[2][3][4].potential << std::endl;
}

void enforceNeumannBC(BoundaryNode* bNodes, const int nodeCount)
{
    int i;

    for(i = 0; i < nodeCount; i++)
        (bNodes[i].bndryNodes[0])->potential = (1./3) * (4 * (bNodes[i].bndryNodes[1])->potential - (bNodes[i].bndryNodes[2])->potential);

}
#endif
