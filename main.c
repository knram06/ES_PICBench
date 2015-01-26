#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
// strtok issues resolved with including this header
#include <string.h>   

#define GRID_LENGTH (3e-4)
#define NUM_NODES 101


/* geometry dimension */
// capillary centered on YZ face
// origin at Corner of domain
// --> Particle data MUST be corrected for new origin
#define CAPILLARY_RADIUS (1.3e-5)
#define EXTRACTOR_INNER_RADIUS (1e-4)
#define EXTRACTOR_OUTER_RADIUS (1.4e-4)
#define CAPILLARY_VOLTAGE 0.
#define EXTRACTOR_VOLTAGE (-1350.)

// define all problem parameters in terms of macros
#define MD_FILE "./input_data/MD_data/10_input_Pos_Q488_20130318.inp"

#define PARTICLE_SIZE ((int)5e4)
#define TIMESTEPS ((int)8e3)

//#define TEST_FUNCTION (x*x - 2*y*y + z*z)
#define TEST_FUNCTION 0.

// declare the Particle data type here
typedef struct
{
    double x, y, z;
    double Vx, Vy, Vz;
    double mass, charge;
} Particle;

// TODO: just use a double array instead?
typedef struct
{
    //double pos[3];
    double potential;
} Node;

typedef struct
{
    Node* bndryNodes[3];
} BoundaryNode;

typedef struct
{
    double components[3];
} EField;

typedef struct
{
    int    numNodes;
    double spacing, invSpacing;
} GridInfo;

// utility function
unsigned int countLinesInFile(FILE* fp);

// pre-process
void parseMDFileToParticles(Particle particleData[], FILE* fp);
void allocateGrid(Node**** grid, GridInfo* gInfo);
void deallocGrid(Node**** grid, GridInfo* gInfo);
void setupBoundaryConditions(Node*** grid, GridInfo* gInfo);
int accumulateNeumannBCNodes(Node*** grid, GridInfo* gInfo, BoundaryNode* bNodes);

// numerics related
void solve(Node*** grid, GridInfo* gInfo, const double tolerance, const double sorOmega);
double single_step_solve(Node*** grid, const int numNodes, const double sorOmega);
void enforceNeumannBC(BoundaryNode* bNodes, const int nodeCount);

// post process
void writeOutputData(const char* fileName, Node*** grid, EField*** ElectricField, GridInfo* gInfo);


// calculate the ElectricField once Node values are known
void calcElectricField(EField*** ElectricField, Node*** grid, GridInfo* gInfo)
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

                /**************************************/
                /****** X = 0 || X = GRID_LENGTH ******/
                double val = 0.;
                if(isI_0_OR_LEN)
                {
                    if(isI_0)
                        val = ( -grid[i+2][j][k].potential
                             + 4*grid[i+1][j][k].potential
                             - 3*grid[i  ][j][k].potential );

                    else if (isI_LEN)
                        val = -( -grid[i-2][j][k].potential 
                              + 4*grid[i-1][j][k].potential
                              - 3*grid[i  ][j][k].potential );
                }
                else
                    val = (grid[i+1][j][k].potential - grid[i-1][j][k].potential);

                // by now ElectricField[i][j][k] needs to be filled
                ElectricField[i][j][k].components[0] = multFactor * val;


                /**************************************/
                /****** Y = 0 || Y = GRID_LENGTH ******/
                if(isJ_0_OR_LEN)
                {
                    if (isJ_0)
                        val = ( -grid[i][j+2][k].potential
                             + 4*grid[i][j+1][k].potential
                             - 3*grid[i][j  ][k].potential );
                    else if (isJ_LEN)
                        val = -(-grid[i][j-2][k].potential 
                             + 4*grid[i][j-1][k].potential
                             - 3*grid[i][j  ][k].potential );
                }
                else
                    val = (grid[i][j+1][k].potential - grid[i][j-1][k].potential);

                // by now ElectricField[i][j][k] needs to be filled
                ElectricField[i][j][k].components[1] = multFactor * val;


                /**************************************/
                /****** Z = 0 || Z = GRID_LENGTH ******/
                if(isK_0_OR_LEN)
                {
                    if (isK_0)
                        val = ( -grid[i][j][k+2].potential
                             + 4*grid[i][j][k+1].potential
                             - 3*grid[i][j][k  ].potential );
                    else if (isK_LEN)
                        val = -( -grid[i][j][k-2].potential 
                              + 4*grid[i][j][k-1].potential
                              - 3*grid[i][j][k  ].potential );
                }
                else
                    val = (grid[i][j][k+1].potential - grid[i][j][k-1].potential);

                // by now ElectricField[i][j][k] needs to be filled
                ElectricField[i][j][k].components[2] = multFactor * val;

            } // end of k loop
        } // end of j loop
    } // end of i loop

}

void allocateEField(EField**** grid, GridInfo* gInfo)
{
    const int numNodes = gInfo->numNodes;
    int i, j, k;

    (*grid) = malloc(numNodes * sizeof(EField**));
    assert((*grid) != NULL);

    for(i = 0; i < numNodes; i++)
    {
        (*grid)[i] = malloc(numNodes * sizeof(EField*));
        for(j = 0; j < numNodes; j++)
        {
            (*grid)[i][j] = malloc(numNodes * sizeof(EField));
            for(k = 0; k < numNodes; k++)
            {
                // initialize the struct values
                (*grid)[i][j][k].components[0] = 0.;
                (*grid)[i][j][k].components[1] = 0.;
                (*grid)[i][j][k].components[2] = 0.;
                //val->pos[0] = 0.;
                //val->pos[1] = 0.;
                //val->pos[2] = 0.;
            }
        }
    }
}

void deallocEField(EField**** grid, GridInfo* gInfo)
{
    const int numNodes = gInfo->numNodes;
    int i, j;

    for(i = 0; i < numNodes; i++)
    {
        for(j = 0; j < numNodes; j++)
            free((*grid)[i][j]);

        free((*grid)[i]);
    }
    free((*grid));
}

int main()
{
    // read in MD data
    FILE *fp = fopen(MD_FILE, "r");
    if(fp == NULL)
    {
        fprintf(stderr, "Error in opening file %s \n", MD_FILE);
        return EXIT_FAILURE;
    }

    // count the number of lines in the input file
    // so that we can preallocate later
    unsigned int particleCount = countLinesInFile(fp);
    printf("Number of lines read was %d.\n", particleCount);

    // fill in GridInfo data
    GridInfo gridInfo;
    gridInfo.numNodes = NUM_NODES;
    gridInfo.spacing = GRID_LENGTH / (NUM_NODES - 1);
    gridInfo.invSpacing = 1./(gridInfo.spacing);        // just caching this to avoid divisions?

    // now preallocate the particles data array
    Particle* MD_data = malloc(particleCount * sizeof(Particle));

    // loop through the file and tokenize entries
    rewind(fp);         // rewind file to the beginning
    parseMDFileToParticles(MD_data, fp);
    fclose(fp);

    // allocate the grid
    Node*** grid = NULL;
    allocateGrid(&grid, &gridInfo);

    EField*** ElectricField = NULL;
    allocateEField(&ElectricField, &gridInfo);

    // setup boundary conditions
    setupBoundaryConditions(grid, &gridInfo);

    // preallocate NeumannBC nodes in terms of MAX possible
    // i.e. num of sides of cube * num nodes per side
    BoundaryNode* bNodes = malloc( 6 * (NUM_NODES)*(NUM_NODES) * sizeof(BoundaryNode) );

    // calculate the actual BoundaryNodes count and resize the array to that instead
    // to avoid wastage
    int nodeCount = accumulateNeumannBCNodes(grid, &gridInfo, bNodes);
    bNodes = realloc(bNodes, nodeCount * sizeof(BoundaryNode));

    // enforce boundary conditions
    // impose Neumann BCs
    // solve and at each step, impose Neumann BCs?
    double tolerance = 1e-9, sorOmega = 1.9;
    double norm = 100.;

    while(norm >= tolerance)
    {
        norm = sqrt(single_step_solve(grid, gridInfo.numNodes, sorOmega));

        // TODO: norm should be updated with this calculation no?
        // as it changes the values in the grid?
        enforceNeumannBC(bNodes, nodeCount);

        printf("norm: %e\n", norm);
    }

    calcElectricField(ElectricField, grid, &gridInfo);

    // allocate for particles
    Particle* domainParticles = malloc(PARTICLE_SIZE * sizeof(Particle));

    // for required number of timesteps,
    int i;
    for(i = 0; i < TIMESTEPS; i++)
    {

    }

    // introduce the particles
    // then move them

    // write out data for post processing
    writeOutputData("out.vtk", grid, ElectricField, &gridInfo);



    free(domainParticles);
    free(bNodes);
    deallocEField(&ElectricField, &gridInfo);
    deallocGrid(&grid, &gridInfo);
    free(MD_data);

    return 0;
}

unsigned int countLinesInFile(FILE* fp)
{
    unsigned int lineCount = 0;
    char ch;

    while(!feof(fp))
    {
        // get a character
        ch = fgetc(fp);

        // if it is a newline, then increment linecount
        if(ch == '\n')
            lineCount++;
    }

    return lineCount;
}

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
                    particleData[particleCount].y = atof(token) + (GRID_LENGTH / 2.);
                    break;

                case 2:
                    particleData[particleCount].z = atof(token) + (GRID_LENGTH / 2.);
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

void allocateGrid(Node**** grid, GridInfo* gInfo)
{
    const int numNodes = gInfo->numNodes;
    int i, j, k;

    (*grid) = malloc(numNodes * sizeof(Node**));
    assert((*grid) != NULL);

    for(i = 0; i < numNodes; i++)
    {
        (*grid)[i] = malloc(numNodes * sizeof(Node*));
        for(j = 0; j < numNodes; j++)
        {
            (*grid)[i][j] = malloc(numNodes * sizeof(Node));
            for(k = 0; k < numNodes; k++)
            {
                // initialize the struct values
                Node* val = &(*grid)[i][j][k];
                //val->pos[0] = 0.;
                //val->pos[1] = 0.;
                //val->pos[2] = 0.;
                val->potential = 0.;
            }
        }
    }
}

void deallocGrid(Node**** grid, GridInfo* gInfo)
{
    const int numNodes = gInfo->numNodes;
    int i, j;

    for(i = 0; i < numNodes; i++)
    {
        for(j = 0; j < numNodes; j++)
            free((*grid)[i][j]);

        free((*grid)[i]);
    }
    free((*grid));
}

// setup the boundary conditions
void setupBoundaryConditions(Node*** grid, GridInfo* gInfo)
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
            //grid[i][0][k] = TEST_FUNCTION;
            grid[i][0][k].potential=TEST_FUNCTION;
            y = GRID_LENGTH;
            //grid[i][SIZE_Y-1][k] = TEST_FUNCTION;
            grid[i][numNodes - 1][k].potential=TEST_FUNCTION;
        }

        // on X-Y faces
        for(j = 0; j < numNodes; j++)
        {
            y = spacing*j;
            z = 0.;
            grid[i][j][0].potential=TEST_FUNCTION;
            z = GRID_LENGTH;
            grid[i][j][numNodes - 1].potential=TEST_FUNCTION;
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
                grid[0][j][k].potential = CAPILLARY_VOLTAGE;
            }

            // EXTRACTOR side
            if( (sumSqs > extractorInner*extractorInner) && (sumSqs < extractorOuter*extractorOuter ) )
            {
                grid[numNodes - 1][j][k].potential = EXTRACTOR_VOLTAGE;
            }

        }
    }
    //std::cout << grid[1][1][0].potential << " " << grid[2][3][4].potential << std::endl;
}

/*!
 * Solve takes in a tolerance
 * 
 * Solves the Laplace Equation using Finite Difference and 
 * using SOR (Successive Over Relaxation) for the resulting
 * Linear system.
 *
 * \todo Improve this to print iteration count and headers
 * Headers in a periodic manner can be printed.
 *
 */
void solve(Node*** grid, GridInfo* gInfo, const double tolerance, const double sorOmega)
{
    const int numNodes = gInfo->numNodes;

    double norm = 100.;
    int count = 0;

    // solve using criteria
    // solution in the "inner" cube
    while (norm >= tolerance)
    {
        count++;
        norm = sqrt(single_step_solve(grid, numNodes, sorOmega));

        printf("norm: %e\n", norm);
    }
}

double single_step_solve(Node*** grid, const int numNodes, const double sorOmega)
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
                temp = grid[i][j][k].potential;

                // obtain the value got by averaging
                double val = (1./6) * (grid[i+1][j][k].potential + grid[i-1][j][k].potential + grid[i][j-1][k].potential + grid[i][j+1][k].potential + grid[i][j][k-1].potential + grid[i][j][k+1].potential);
 
                // obtain the value with SOR
                grid[i][j][k].potential = val*sorOmega + (1-sorOmega)*temp;
                //std::cout << grid[i][j][k].potential << " " << i << " " << j << " " << k << std::endl;
 
                // store the difference in values at same grid points in temp itself
                // for norm calculation
                temp = grid[i][j][k].potential - temp;
                norm += temp * temp;
            }
        }
    }

    return norm;
}

// function for writing out values
void writeOutputData(const char* fileName, Node*** grid, EField*** ElectricField, GridInfo* gInfo)
{
    FILE* fileValues = fopen(fileName, "w");
    int i, j, k;

    const int numNodes = gInfo->numNodes;
    const double spacing = gInfo->spacing;

    const int totalNodes = numNodes*numNodes*numNodes;

    // write the VTK header
    fprintf(fileValues, "# vtk DataFile Version 2.0\n"
                        "Potential data\n"
                        "ASCII\n"
                        "DATASET STRUCTURED_GRID\n"
                        "DIMENSIONS %d %d %d\n"
                        "POINTS %d float\n", numNodes, numNodes, numNodes, totalNodes
                        );

    // cache some of the data which will be written later
    double* potentialValues = malloc(totalNodes * sizeof(double));
    double* electricField   = malloc(3 * totalNodes * sizeof(double));

    int count = 0;
    for(i = 0; i < numNodes; i++)
    {
        double x = spacing * i;
        for(j = 0; j < numNodes; j++)
        {
            double y = spacing * j;
            for(k = 0; k < numNodes; k++)
            {
                double z = spacing * k;

                fprintf(fileValues, "%10.8e %10.8e %10.8e\n", x, y, z);

                // update the potential data array
                potentialValues[count] = grid[i][j][k].potential;

                // update the electric field array
                const int pos = 3*count;
                electricField[pos]     = ElectricField[i][j][k].components[0];
                electricField[pos + 1] = ElectricField[i][j][k].components[1];
                electricField[pos + 2] = ElectricField[i][j][k].components[2];

                count++;
            }
        }
    }

    // now write out the potential values
    fprintf(fileValues, "\n"
                        "POINT_DATA %d\n"
                        "SCALARS potential float 1\n"
                        "LOOKUP_TABLE default\n", totalNodes
            );
    for(count = 0; count < totalNodes; count++)
        fprintf(fileValues, "%10.8e\n", potentialValues[count]);

    // write out the ElectricField values
    fprintf(fileValues, "\n"
                        "VECTORS ElectricField float\n");
    for(count = 0; count < totalNodes; count++)
    {
        const int pos = 3*count;
        fprintf(fileValues, "%10.8e %10.8e %10.8e\n", electricField[pos], electricField[pos+1], electricField[pos+2] );
    }


    free(electricField);
    free(potentialValues);
    fclose(fileValues);
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

void enforceNeumannBC(BoundaryNode* bNodes, const int nodeCount)
{
    int i;

    for(i = 0; i < nodeCount; i++)
        (bNodes[i].bndryNodes[0])->potential = (1./3) * (4 * (bNodes[i].bndryNodes[1])->potential - (bNodes[i].bndryNodes[2])->potential);

}
