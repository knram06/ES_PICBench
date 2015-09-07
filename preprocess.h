#ifndef PREPROCESS_H
#define PREPROCESS_H

void parseMDFileToParticles(Particle particleData[], FILE* fp)
{
    // buffer to read in line data
    char line[256];
    
    // keep track of particleCount
    int particleCount = 0;

    // fgets stops reading when it encounters the newline character
    while(fgets(line, sizeof(line), fp))
    {
        char* token;
        int fieldIndex = 1;

        //printf("line: %s\n", line);
        token = strtok(line, " ");
        particleData[particleCount].x = atof(token);

        while(fieldIndex < 8)
        {
            token = strtok(NULL, " ");

            switch(fieldIndex)
            {
                case 1:
                    particleData[particleCount].y = atof(token);//+ (GRID_LENGTH / 2.);
                    break;

                case 2:
                    particleData[particleCount].z = atof(token);//+ (GRID_LENGTH / 2.);
                    break;

                case 3:
                    particleData[particleCount].Vx = atof(token);
                    break;

                case 4:
                    particleData[particleCount].Vy = atof(token);
                    break;

                case 5:
                    particleData[particleCount].Vz = atof(token);
                    break;

                case 6:
                    particleData[particleCount].mass = atof(token);
                    break;

                case 7:
                    particleData[particleCount].charge= atof(token);
                    break;
            }

            /*if(fieldIndex == 1)
                particleData[particleCount].y = atof(token) + (GRID_LENGTH / 2.);
            else if(fieldIndex == 2)
                particleData[particleCount].z = atof(token) + (GRID_LENGTH / 2.);
            else if(fieldIndex == 3)
                particleData[particleCount].Vx = atof(token);
            else if(fieldIndex == 4)
                particleData[particleCount].Vy = atof(token);
            else if(fieldIndex == 5)
                particleData[particleCount].Vz = atof(token);
            else if(fieldIndex == 6)
                particleData[particleCount].mass = atof(token);
            else if(fieldIndex == 7)
                particleData[particleCount].charge= atof(token);
            */

            // increment fieldIndex
            fieldIndex++;
        }

        // increment particleCounter
        particleCount++;
    }
}

void allocateGrid(Node** grid, GridInfo* gInfo)
{
    const int totalNodes = gInfo->totalNodes;
    int i, j, k;

    (*grid) = malloc(totalNodes * sizeof(Node));
    assert((*grid) != NULL);
}

void deallocGrid(Node** grid)
{
    free((*grid));
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

int accumulateNeumannBCNodes(Node*** grid, GridInfo* gInfo, BoundaryNode* bNodes)
{
    const int numNodes           = gInfo->numNodes;
    const double spacing         = gInfo->spacing;

    const double center[2]       = {GRID_LENGTH / 2., GRID_LENGTH / 2.};
    const double capillaryRadius = CAPILLARY_RADIUS;

    // extractor dimensions
    const double extractorInner  = EXTRACTOR_INNER_RADIUS;
    const double extractorOuter  = EXTRACTOR_OUTER_RADIUS;

    int i, j, k;
    int nodeCount = 0;

    /**** Y-Z faces ****/
    // Y-Z face with capillary
    // i.e. X = 0 position with center at corner
    for(j = 0; j < numNodes; j++)
    {
        //double y = spacing * j;
        double ty = (spacing*j - center[0]);

        for(k = 0; k < numNodes; k++)
        {
            //double z = spacing * k;
            double tz = (spacing*k - center[1]);
            double sumSqs = ty*ty + tz*tz;

            // CAPILLARY side
            // we need points OUTSIDE the capillary for Neumann BC
            if( sumSqs >= capillaryRadius*capillaryRadius )
            {
                bNodes[nodeCount].bndryNodes[0] = &grid[0][j][k];
                bNodes[nodeCount].bndryNodes[1] = &grid[1][j][k];
                bNodes[nodeCount].bndryNodes[2] = &grid[2][j][k];

                nodeCount++;
            }

            // EXTRACTOR side
            if( (sumSqs <= extractorInner*extractorInner) || (sumSqs >= extractorOuter*extractorOuter ) )
            {
                bNodes[nodeCount].bndryNodes[0] = &grid[numNodes-1][j][k];
                bNodes[nodeCount].bndryNodes[1] = &grid[numNodes-2][j][k];
                bNodes[nodeCount].bndryNodes[2] = &grid[numNodes-3][j][k];

                nodeCount++;

            }

        }
    }

    /**** X-Z faces ****/
    for(i = 0; i < numNodes; i++)
    {
        for(k = 0; k < numNodes; k++)
        {
            // Y = 0
            bNodes[nodeCount].bndryNodes[0] = &grid[i][0][k];
            bNodes[nodeCount].bndryNodes[1] = &grid[i][1][k];
            bNodes[nodeCount].bndryNodes[2] = &grid[i][2][k];

            nodeCount++;

            // Y = GRID_LENGTH
            bNodes[nodeCount].bndryNodes[0] = &grid[i][numNodes-1][k];
            bNodes[nodeCount].bndryNodes[1] = &grid[i][numNodes-2][k];
            bNodes[nodeCount].bndryNodes[2] = &grid[i][numNodes-3][k];

            nodeCount++;
        }
    }

    /**** X-Y faces ****/
    for(i = 0; i < numNodes; i++)
    {
        for(j = 0; j < numNodes; j++)
        {
            // Z = 0
            bNodes[nodeCount].bndryNodes[0] = &grid[i][j][0];
            bNodes[nodeCount].bndryNodes[1] = &grid[i][j][1];
            bNodes[nodeCount].bndryNodes[2] = &grid[i][j][2];
            nodeCount++;

            // Z = GRID_LENGTH
            bNodes[nodeCount].bndryNodes[0] = &grid[i][j][numNodes-1];
            bNodes[nodeCount].bndryNodes[1] = &grid[i][j][numNodes-2];
            bNodes[nodeCount].bndryNodes[2] = &grid[i][j][numNodes-3];
            nodeCount++;
        }
    }

    return nodeCount;
}

#endif
