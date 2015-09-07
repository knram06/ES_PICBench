#ifndef NUMERICS_H
#define NUMERICS_H

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

void enforceNeumannBC(BoundaryNode* bNodes, const int nodeCount)
{
    int i;

    for(i = 0; i < nodeCount; i++)
        (bNodes[i].bndryNodes[0])->potential = (1./3) * (4 * (bNodes[i].bndryNodes[1])->potential - (bNodes[i].bndryNodes[2])->potential);

}
#endif
