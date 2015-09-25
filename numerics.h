#ifndef NUMERICS_H
#define NUMERICS_H

void build_A_matrix(double* A, const int sizeA)
{
    const int mat_side_size = sizeA*sizeA*sizeA;
    int i, j, k;
    int I, J;

    // define coeff in this manner
    // to make this positive definite
    const int selfCoeff = 6;
    const int nonSelfCoeff = -1;

    const int steps[2] = {-1, 1};
    int c;

    /*******************************************/
    /************All - Nodes********************/
    for(i = 0; i < sizeA; i++)
    {
        for(j = 0; j < sizeA; j++)
        {
            for(k = 0; k < sizeA; k++)
            {
                I = INDEX_1D(sizeA, i, j, k);
                const int NI = mat_side_size * I;
                A[NI+I] = selfCoeff;

                for(c = 0; c < 2; c++)
                {
                    int step = steps[c];
                    int iStep = i + step;
                    int jStep = j + step;
                    int kStep = k + step;

                    if((iStep >= 0) && (iStep < sizeA))
                    {
                        J = INDEX_1D(sizeA, iStep, j, k);
                        A[NI+J] = nonSelfCoeff;
                    }
                    if((jStep >= 0) && (jStep < sizeA))
                    {
                        J = INDEX_1D(sizeA, i, jStep, k);
                        A[NI+J] = nonSelfCoeff;
                    }
                    if((kStep >= 0) && (kStep < sizeA))
                    {
                        J = INDEX_1D(sizeA, i, j, kStep);
                        A[NI+J] = nonSelfCoeff;
                    }
                }
            } // end of k loop
        } // end of j loop
    }
}

void build_b_vec(double* b, GridInfo* gridInfo, const int sideSize)
{
    int i, j, k;
    int I, J;
    const double spacing = gridInfo->spacing;
    const int steps[2] = {-1, 1};
    int c;

    /*******************************************/
    /************All - Nodes********************/
    for(i = 0; i < sideSize; i++)
    {
        const double x = spacing*i;
        for(j = 0; j < sideSize; j++)
        {
            const double y = spacing*j;
            for(k = 0; k < sideSize; k++)
            {
                const double z = spacing*k;
                I = INDEX_1D(sideSize, i, j, k);
                b[I] = 0.;

                for(c = 0; c < 2; c++)
                {
                    int step = steps[c];
                    int iStep = i + step;
                    int jStep = j + step;
                    int kStep = k + step;

                    // X-Face boundary
                    if(iStep < 0)
                        b[I] += TEST_FUNCTION(0, y, z);
                    else if (iStep == sideSize)
                        b[I] += TEST_FUNCTION(GRID_LENGTH, y, z);

                    // Y-Face boundary
                    if(jStep < 0)
                        b[I] += TEST_FUNCTION(x, 0,           z);
                    else if (jStep == sideSize)
                        b[I] += TEST_FUNCTION(x, GRID_LENGTH, z);

                    // Z-Face boundary
                    if(kStep < 0)
                        b[I] += TEST_FUNCTION(x, y,           0);
                    else if (kStep == sideSize)
                        b[I] += TEST_FUNCTION(x, y, GRID_LENGTH);
                }
            }
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
void calcElectricField(EField* ElectricField, Node* grid, GridInfo* gInfo)
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
                        val = ( -GRID_1D(grid,i+2,j,k).potential
                                + 4*GRID_1D(grid,i+1,j,k).potential
                                - 3*GRID_1D(grid,i  ,j,k).potential );

                    else if (isI_LEN)
                        val = -( -GRID_1D(grid,i-2,j,k).potential 
                                + 4*GRID_1D(grid,i-1,j,k).potential
                                - 3*GRID_1D(grid,i  ,j,k).potential );
                }
                else
                    val = (GRID_1D(grid,i+1,j,k).potential - GRID_1D(grid,i-1,j,k).potential);

                // by now ElectricField[i][j][k] needs to be filled
                elecField->components[0] = multFactor * val;


                /**************************************/
                /****** Y = 0 || Y = GRID_LENGTH ******/
                if(isJ_0_OR_LEN)
                {
                    if (isJ_0)
                        val = ( -GRID_1D(grid,i,j+2,k).potential
                                + 4*GRID_1D(grid,i,j+1,k).potential
                                - 3*GRID_1D(grid,i,j  ,k).potential );
                    else if (isJ_LEN)
                        val = -(-GRID_1D(grid,i,j-2,k).potential 
                                + 4*GRID_1D(grid,i,j-1,k).potential
                                - 3*GRID_1D(grid,i,j  ,k).potential );
                }
                else
                    val = (GRID_1D(grid,i,j+1,k).potential - GRID_1D(grid,i,j-1,k).potential);

                // by now ElectricField[i][j][k] needs to be filled
                elecField->components[1] = multFactor * val;


                /**************************************/
                /****** Z = 0 || Z = GRID_LENGTH ******/
                if(isK_0_OR_LEN)
                {
                    if (isK_0)
                        val = ( -GRID_1D(grid,i,j,k+2).potential
                                + 4*GRID_1D(grid,i,j,k+1).potential
                                - 3*GRID_1D(grid,i,j,k  ).potential );
                    else if (isK_LEN)
                        val = -( -GRID_1D(grid,i,j,k-2).potential 
                                + 4*GRID_1D(grid,i,j,k-1).potential
                                - 3*GRID_1D(grid,i,j,k  ).potential );
                }
                else
                    val = (GRID_1D(grid,i,j,k+1).potential - GRID_1D(grid,i,j,k-1).potential);

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
        //x = spacing * i;
        for(k = 0; k < numNodes; k++)
        {
            //z = spacing * k;
            //y = 0.;
            //GRID_1D(grid, i, 0, k).potential=TEST_FUNCTION;

            // enforce Neumann boundary nodes
            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, 0, k);
            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, 1, k);
            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, 2, k);
            nodeCount++;

            //y = GRID_LENGTH;
            //GRID_1D(grid, i, numNodes-1, k).potential=TEST_FUNCTION;

            // Y = GRID_LENGTH
            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, numNodes-1, k);
            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, numNodes-2, k);
            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, numNodes-3, k);
            nodeCount++;
        }

        // on X-Y faces
        for(j = 0; j < numNodes; j++)
        {
            //y = spacing*j;
            //z = 0.;
            //GRID_1D(grid, i, j, 0).potential=TEST_FUNCTION;

            // enforce Neumann boundary nodes
            // Z = 0
            bNodes[nodeCount].bndryNodes[0] = &GRID_1D(grid, i, j, 0);
            bNodes[nodeCount].bndryNodes[1] = &GRID_1D(grid, i, j, 1);
            bNodes[nodeCount].bndryNodes[2] = &GRID_1D(grid, i, j, 2);
            nodeCount++;

            //z = GRID_LENGTH;
            //GRID_1D(grid, i, j, numNodes-1).potential=TEST_FUNCTION;

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


/********** CG-SOLVER Routines ***********/

// TODO: add and check restrict usage here
void matvec(const double* const mat, double *vec, double *result, const int n)
{
    int i, j;
    for(i = 0; i < n; i++)
    {
        double sum = 0.;
        for(j = 0; j < n; j++)
            sum += GRID_1D(mat, i, j) * vec[j];

        result[i] = sum;
    }
}

void vec_sub(double *a, double *b, double *res, const int n)
{
    int i;
    for(i = 0; i < n; i++)
        res[i] = a[i] - b[i];
}

void vec_add(double* a, double *b, double *res, const int n)
{
    int i;
    for(i = 0; i < n; i++)
        res[i] = a[i] + b[i];
}

void vec_mult(double val, double *v, double *res, const int n)
{
    int i;
    for(i = 0; i < n; i++)
        res[i] = val*v[i];
}

void daxpy(const double alpha, const double *x, double *y, const int n)
{
    int i;
    for(i = 0; i < n; i++)
        y[i] += alpha*x[i];
}


void vec_set_equal(double *dst, const double* const src, const int n)
{
    int i;
    for(i = 0; i < n; i++)
        dst[i] = src[i];
}

double vec_transp_mult(const double* const v, const double* const r, const int n)
{
    double ret = 0.;
    int i;

    for(i = 0; i < n; i++)
        ret += v[i] * r[i];

    return ret;
}

int cg_solve(const double* A, const double* b, double tolerance,
                   double* x, const int n)
{
    double *r = malloc(n * sizeof(double));
    double *d = malloc(n * sizeof(double));

    // some ops to setup initial conditions and values
    matvec(A, x, r, n);
    vec_sub(b, r, r, n);
    vec_set_equal(d, r, n);

    // store some variables used in the while looping
    double residual = vec_transp_mult(r, r, n);
    double *Ad = malloc(n*sizeof(double));
    double *temp = malloc(n*sizeof(double));

    // residual is the square, so tolerance is also squared
    int iterCount = 0;
    tolerance = tolerance*tolerance;
    while(residual > tolerance)
    {
        iterCount++;
        matvec(A, d, Ad, n);

        double alpha = residual / vec_transp_mult(d, Ad, n);
        daxpy(alpha, d, x, n);

        // multiply alpha and d - store in temp
        //vec_mult(alpha, d, temp, n);

        // add x and alpha*d - store in x
        //vec_add(x, temp, x, n);

        // now do d_transp * Ad and store in temp
        daxpy(-alpha, Ad, r, n);
        //vec_mult(alpha, Ad, temp, n);
        //vec_sub(r, temp, r, n);

        double residualOld = residual;
        residual = vec_transp_mult(r, r, n);

        double beta = residual / residualOld;

        // multiply Beta and d
        vec_mult(beta, d, temp, n);
        vec_add(r, temp, d, n);
    }

    // IMPORTANT!! - free all malloc'ed memory
    free(temp); free(Ad);
    free(d); free(r);

    return iterCount;
}

#endif
