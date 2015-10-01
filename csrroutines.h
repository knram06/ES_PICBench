#ifndef CSRROUTINES_H
#define CSRROUTINES_H

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
    matcsr->mat = malloc(sizeof(double) * maxNZ);

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
    const double center[2]       = {GRID_LENGTH / 2., GRID_LENGTH / 2.};
    const double capillaryRadius = CAPILLARY_RADIUS;

    // extractor dimensions
    const double extractorInner  = EXTRACTOR_INNER_RADIUS;
    const double extractorOuter  = EXTRACTOR_OUTER_RADIUS;

    int i, j, k;
    const int numNodes = gridInfo->numNodes;
    const double spacing = gridInfo->spacing;

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
            const double y = j*spacing;

            for(k = 0; k < numNodes; k++)
            {
                bool isK_0 =   (k == 0);
                bool isK_LEN = (k == (numNodes-1));
                bool isK_0_OR_LEN = isK_0 || isK_LEN;
                const double z = k*spacing;

                const int I = INDEX_1D(numNodes, i, j, k);

                // assuming repeated values 9 entries could be made here
                // i.e. corner points, self coefficients will add up
                ColVal cvals[9];
                int elementsInserted = 0;

                // supply a default
                rhs[I] = 0;

                // TODO: improve this using ghost cells?
                // if on boundary faces
                if(isI_0_OR_LEN || isJ_0_OR_LEN || isK_0_OR_LEN)
                {
                    // if on X-Faces
                    if(isI_0_OR_LEN)
                    {
                        double ty = y - center[0];
                        double tz = z - center[1];

                        double rr = ty*ty + tz*tz;
                        if(isI_0)
                        {
                            // check if within capillary
                            if (rr <= capillaryRadius*capillaryRadius)
                            {
                                // if so impose Dirichlet
                                cvals[elementsInserted].index = I;
                                cvals[elementsInserted].val = 1;
                                elementsInserted++;

                                rhs[I] = CAPILLARY_VOLTAGE;
                            }
                            // else impose Neumann
                            {
                                // at point i,j,k
                                cvals[elementsInserted].index = I;
                                cvals[elementsInserted].val = 3;
                                elementsInserted++;

                                // at point i+1,j,k
                                cvals[elementsInserted].index = INDEX_1D(numNodes, i+1, j, k);
                                cvals[elementsInserted].val = -4;
                                elementsInserted++;

                                // at point i+2,j,k
                                cvals[elementsInserted].index = INDEX_1D(numNodes, i+2, j, k);
                                cvals[elementsInserted].val = 1;
                                elementsInserted++;
                            }
                        } // end of check for if X == 0 face
                        // else on Extractor face then
                        else
                        {
                            // if within Extractor
                            if((rr >= EXTRACTOR_INNER_RADIUS*EXTRACTOR_INNER_RADIUS)
                                                &&
                               (rr <= EXTRACTOR_OUTER_RADIUS*EXTRACTOR_OUTER_RADIUS))
                            {
                                // impose Dirichlet condition
                                cvals[elementsInserted].index = I;
                                cvals[elementsInserted].val = 1;
                                elementsInserted++;

                                rhs[I] = EXTRACTOR_VOLTAGE;
                            }
                            // if outside the extractor
                            else
                            {
                                // at point i-2,j,k
                                cvals[elementsInserted].index = INDEX_1D(numNodes, i-2, j, k);
                                cvals[elementsInserted].val = 1;
                                elementsInserted++;

                                // at point i-1,j,k
                                cvals[elementsInserted].index = INDEX_1D(numNodes, i-1, j, k);
                                cvals[elementsInserted].val = -4;
                                elementsInserted++;

                                // at point i,j,k
                                cvals[elementsInserted].index = I;
                                cvals[elementsInserted].val = 3;
                                elementsInserted++;
                            }
                        } // end of else check - to check if on Extractor face
                    } // end of check, if on X-Faces

                    // on Neumann boundaries
                    // we  sort cols and vals
                    // and then INSERT
                    // if on Y and Z boundary faces
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



#endif
