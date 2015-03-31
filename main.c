#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
// strtok issues resolved with including this header
#include <string.h>
#include <time.h>

#define GRID_LENGTH (3e-4)
#define NUM_NODES 65

/* geometry dimension */
// capillary centered on YZ face
// origin at Corner of domain
// --> Particle data MUST be corrected for new origin
#define CAPILLARY_RADIUS (1.40625e-5)
#define EXTRACTOR_INNER_RADIUS (9.375e-5)
#define EXTRACTOR_OUTER_RADIUS (1.40625e-4)
#define CAPILLARY_VOLTAGE 10.
#define EXTRACTOR_VOLTAGE (-1350.)

// timesteps info
#define T_MD (230e-12)
#define T_PIC (30e-12)

// Physical constants
#define ELECTRONIC_CHARGE (1.60217657e-19)

// array bounds and margins
#define LOST_PARTICLES_MARGIN 50

// define all problem parameters in terms of macros
#define MD_FILE "./input_data/MD_data/10_input_Pos_Q488_20130318.inp"

#define PARTICLE_SIZE ((int)5e4)
#define TIMESTEPS ((int)0)
#define ITER_INTERVAL (200)
#define ITER_HEADER_INTERVAL (5000)
#define POST_WRITE_FILES (false)
#define POST_INTERVAL (1000)
#define POST_WRITE_PATH ("output/")

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

// release of particles functions
//void calculateReleaseRate(int* Nrel, double* Nfrac, int particleCount);
void releaseParticles(const int numParticlesToRelease,
                      const Particle* inputData, const int inputCount,
                      Particle* domainParticles, int* domainParticleBound,
                      int* lostParticlesArray, int* lostParticleBound);
Particle randomizeParticleAttribs(Particle inputParticle);
bool isParticleInDomain(const Particle p)
{
    return (p.x <= GRID_LENGTH)
        && (p.y <= GRID_LENGTH)
        && (p.z <= GRID_LENGTH);
}
void swapGapsWithEndParticles(Particle* domainParticles, int* domainParticleBound,
                              int* lostParticlesArray, int* lostParticleBound);

// warning when using const EField***, omitting for now
int moveParticlesInField(Particle* domainParticles, int domainParticleBound,
                         int* lostParticlesArray,
                         //int* lostParticlesBound,
                         EField*** ElectricField, GridInfo* gInfo);

// numerics related
void solve(Node*** grid, GridInfo* gInfo, const double tolerance, const double sorOmega);
double single_step_solve(Node*** grid, const int numNodes, const double sorOmega);
double enforceNeumannBC(BoundaryNode* bNodes, const int nodeCount, double sorOmega);

// post process
void writeOutputData(const char* fileName, Node*** grid, EField*** ElectricField, GridInfo* gInfo);
void writeOutputDataXML(const char* fileName, Node*** grid, EField*** ElectricField, GridInfo* gInfo);
void writeParticleData(const char* fileName, const Particle* particleData, int particleCount);


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
    printf("Consolidating Neumann BC nodes into a different structure....");
    int nodeCount = accumulateNeumannBCNodes(grid, &gridInfo, bNodes);
    bNodes = realloc(bNodes, nodeCount * sizeof(BoundaryNode));
    printf("done\n");

    writeOutputData("initial.vtk", grid, ElectricField, &gridInfo);
    // enforce boundary conditions
    // impose Neumann BCs
    // solve and at each step, impose Neumann BCs?
    double tolerance = 1e+3, sorOmega = 1.9;
    double norm = 1e5;

    int iterCount = 1;
    //while(norm >= tolerance)
    clock_t start = clock(), diff;
    for(iterCount = 1; norm >= tolerance; iterCount++)
    {
        norm = sqrt( single_step_solve(grid, gridInfo.numNodes, sorOmega) + enforceNeumannBC(bNodes, nodeCount, sorOmega) );

        // TODO: norm should be updated with this calculation no?
        // as it changes the values in the grid?
        //enforceNeumannBC(bNodes, nodeCount);

        if(!(iterCount % ITER_HEADER_INTERVAL))
            printf("%10s %20s\n", "Iter_Count", "Norm");

        if(!(iterCount % ITER_INTERVAL) )
            printf("%10d %20.8e\n", iterCount, norm);

    }
    diff = clock() - start;
    double solveTime = diff /CLOCKS_PER_SEC;

    printf("\nCalculating Electric Field.....");
    calcElectricField(ElectricField, grid, &gridInfo);
    printf("done\n");

    // write out data for post processing
    writeOutputData("out.vtk", grid, ElectricField, &gridInfo);


    // allocate for particles
    Particle* domainParticles = malloc(PARTICLE_SIZE * sizeof(Particle));
    int totalParticlesCount = 0;

    // calculate the release rate
    int Nrel;
    double Nfrac, runningNfrac, temp;
    double particleReleaseRate = particleCount * (T_PIC/ T_MD);

    // set the integer and fractional parts
    Nfrac = modf(particleReleaseRate, &temp);
    Nrel = (int)(temp);
    runningNfrac = Nfrac;

    //calculateReleaseRate(&Nrel, &Nfrac, particleCount);

    // array to keep track of lost particles
    int* lostParticles = malloc( (Nrel + LOST_PARTICLES_MARGIN) * sizeof(int) );
    int lostParticleBound = -1;

    start = clock();
    // for required number of timesteps
    int i, lostParticleCount = 0;
    for(i = 1; i <= TIMESTEPS; i++)
    {
        printf("\nTimestep %d:\n", i);

        // calculate the current timestep's release rate
        runningNfrac = modf(runningNfrac, &temp);
        const int numParticlesToRelease = Nrel + (int)(temp);

        //totalParticlesBound += numParticlesToRelease - lostParticleCount;
        lostParticleBound = (lostParticleCount - 1);        // adjust for one off issue
        // introduce the particles
        //releaseParticles(1,
        releaseParticles(numParticlesToRelease,
                         MD_data, particleCount,
                         domainParticles, &totalParticlesCount,
                         lostParticles,
                         &lostParticleBound);          // adjust for one off issue

        //totalParticlesBound = (totalParticlesCount - 1);
        swapGapsWithEndParticles(domainParticles, &totalParticlesCount,
                                 lostParticles, &lostParticleBound);

        printf("Total Number of Particles: %d\n", totalParticlesCount); // need +1 for the one-off offset

        // then move them
        lostParticleCount = moveParticlesInField(domainParticles, totalParticlesCount,
                                                 lostParticles,//&lostParticleBound,
                                                 ElectricField, &gridInfo);
        printf("%d Particles left the domain\n", lostParticleCount);
        //printf("%d empty slots in the domain particles\n", lostParticleBound+1);

        // IMPORTANT: add Nfrac to runningNfrac to adjust correctly
        // for the fractional part
        runningNfrac += Nfrac;

        char outputPath[50];
        if(POST_WRITE_FILES && !(i % POST_INTERVAL) )
        {
            sprintf(outputPath, "%s/particleOutput_%d.txt", POST_WRITE_PATH, i);
            writeParticleData(outputPath, domainParticles, totalParticlesCount);
        }
    }

    diff = clock() - start;
    double timeStepsTime = diff /CLOCKS_PER_SEC;

    printf("\nTiming Info\n%10s %10.8e\n%10s %10.8e\n", "Solve", solveTime, "TimeSteps", timeStepsTime);


    free(lostParticles);
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
        double ty = fabs(spacing*j - center[0]);

        for(k = 0; k < numNodes; k++)
        {
            double tz = fabs(spacing*k - center[1]);
            //double sumSqs = ty*ty + tz*tz;

            // CAPILLARY side
            // we need points INSIDE the capillary for Neumann BC
            if( (ty <= capillaryRadius) && (tz <= capillaryRadius) )
            {
                grid[0][j][k].potential = CAPILLARY_VOLTAGE;
            }

            // EXTRACTOR side
            if( ((ty > extractorInner) || (tz > extractorInner))
                                &&
                ((ty <= extractorOuter) && (tz <= extractorOuter))
              )
            {
                grid[numNodes - 1][j][k].potential = EXTRACTOR_VOLTAGE;
            }

        }
    }
    //std::cout << grid[1][1][0].potential << " " << grid[2][3][4].potential << std::endl;
}

//void calculateReleaseRate(int* Nrel, double* Nfrac, int particleCount)
//{
//    const double tMD = T_MD;
//    const double tPIC = T_PIC;
//
//    const int nMD = particleCount;
//
//    double particleReleaseRate = nMD * tPIC / tMD;
//    double intPart;
//
//    (*Nfrac) = modf(particleReleaseRate, &intPart);
//    (*Nrel) = (int)(intPart);
//}

// TODO: document the args
void releaseParticles(const int numParticlesToRelease,
                      const Particle* inputData, const int inputCount,
                      Particle* domainParticles, int* domainParticleCount,
                      int* lostParticlesArray, int* lostParticleBound)
{
    int i;

    // TODO: move this one level above?
    // increment the domainParticle count in preparation for insertion
    //(*domainParticleBound) += (numParticlesToRelease - (*lostParticlesBound));
    //int particleBoundCopy = (*domainParticleBound);

    // assert just in case
    assert(*domainParticleCount <= PARTICLE_SIZE);

    // loop through MD_data and randomize particle attributes
    for(i = 0; i < numParticlesToRelease; i++)
    {
        Particle releasedParticle = randomizeParticleAttribs(inputData[rand() % inputCount] );

        // if bound is not 0
        if(*lostParticleBound >= 0)
        {
            // insert particles from the end
            domainParticles[ lostParticlesArray[(*lostParticleBound)] ] = releasedParticle;
            (*lostParticleBound)--;
            //(*domainParticleBound)--;
        }
        // if it is 0
        // then insert remaining elements at the end of the domainParticle array
        else
        {
            // store a copy of domainParticleBound
            domainParticles[ *domainParticleCount] = releasedParticle;

            // increment counter since we are inserting at the end of the
            // domainParticles array
            (*domainParticleCount)++;
        }
    }

}

void swapGapsWithEndParticles(Particle* domainParticles, int* domainParticleCount,
                              int* lostParticlesArray, int* lostParticleBound)
{
    /* if at the end of this, lostParticlesBound is still not zero,
     * that means we lost a large enough number of particles, that the
     * releasedParticles did not fill up the gaps.
     * So swap all gaps with the last few particles and decrement bounds
     * appropriately */
    while(*lostParticleBound >= 0)
    {
        // swap lostParticlesBound with domainParticleBound value
        (*domainParticleCount)--;       // important to decrement this first?
        domainParticles[ lostParticlesArray[*lostParticleBound] ] = domainParticles[ (*domainParticleCount) ];
        (*lostParticleBound)--;
    }

    // at the end of this, lostParticleBound HAS to be -1 again no?
    //assert(*lostParticleBound == -1);
}

int moveParticlesInField(Particle* domainParticles, int domainParticleBound,
                         int* lostParticlesArray,//int* lostParticlesBound,
                         EField*** ElectricField, GridInfo* gInfo)
{
    const double invSpacing = gInfo->invSpacing;
    const double halfTimeStep = 0.5 * T_PIC;

    // FOR NOW
    // Approximate each particle to nearest node
    // TODO: Weight it based on distances to other nodes?

    // loop through all particles
    int i, lostParticleCount = 0;
    //int lostParticlesBound = -1;
    for(i = 0; i < domainParticleBound; i++)
    {
        Particle* p = &domainParticles[i];

        // round off and get the nearest indices
        int iPos = (int)(round( (p->x)*invSpacing) );
        int jPos = (int)(round( (p->y)*invSpacing) );
        int kPos = (int)(round( (p->z)*invSpacing) );

        // now store the ElectricField value
        EField* elecField = &ElectricField[iPos][jPos][kPos];

        // TODO: generalize as a loop through the number of components?
        const double multFactor = (p->charge * ELECTRONIC_CHARGE * T_PIC/ p->mass);

        // update the velocity
        // save the old velocities
        double old_Vx = p->Vx,
               old_Vy = p->Vy,
               old_Vz = p->Vz;

        p->Vx += elecField->components[0] * multFactor;
        p->Vy += elecField->components[1] * multFactor;
        p->Vz += elecField->components[2] * multFactor;

        // update the position
        p->x += halfTimeStep * (old_Vx + p->Vx);
        p->y += halfTimeStep * (old_Vy + p->Vy);
        p->z += halfTimeStep * (old_Vz + p->Vz);

        // check if the particle is out of bounds
        if(!isParticleInDomain(*p))
        {
            // if so, increment lost particle counter
            lostParticlesArray[lostParticleCount] = i;
            lostParticleCount++;
        }
    }

    return lostParticleCount;
}

Particle randomizeParticleAttribs(Particle inputParticle)
{
    //Particle ret;                                                                       

    double randNum = (rand() / (double)(RAND_MAX)) * 2. * M_PI;
    double rYZ = sqrt(inputParticle.y*inputParticle.y + inputParticle.z*inputParticle.z);

    // adjust for the fact that our origin is at corner
    // whereas MD data origin is at capillary center
    // TODO: enforce this in a better place?
    inputParticle.y = rYZ * cos(randNum) + (GRID_LENGTH/2.);
    inputParticle.z = rYZ * sin(randNum) + (GRID_LENGTH/2.);

    // choose a different random number for velocity - just for more randomization
    randNum = (rand() / (double)(RAND_MAX)) * 2. * M_PI;
    double velYZ = sqrt(inputParticle.Vy*inputParticle.Vy + inputParticle.Vz*inputParticle.Vz);
    inputParticle.Vy = velYZ * cos(randNum);
    inputParticle.Vz = velYZ * sin(randNum);

    return inputParticle;
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

    // temp adjustment to match grid positions
    const double halfGridLength = GRID_LENGTH/2.;
    int count = 0;
    for(i = 0; i < numNodes; i++)
    {
        double x = spacing * i;
        for(j = 0; j < numNodes; j++)
        {
            double y = (spacing * j) - halfGridLength;
            for(k = 0; k < numNodes; k++)
            {
                double z = (spacing * k) - halfGridLength;

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

void writeOutputDataXML(const char* fileName, Node*** grid, EField*** ElectricField, GridInfo* gInfo)
{
    FILE* fileValues = fopen(fileName, "w");
    int i, j, k;

    const int numNodes = gInfo->numNodes;
    const double spacing = gInfo->spacing;

    const int totalNodes = numNodes*numNodes*numNodes;

    // write the VTK header
    fprintf(fileValues,
            "<VTKFile type = \"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            "  <StructuredGrid WholeExtent=\"%f %f %f %f %f %f\">\n",
                                            0., GRID_LENGTH, 0., GRID_LENGTH, 0., GRID_LENGTH);

    fprintf(fileValues,
            "    <Piece Extent=\"%f %f %f %f %f %f\">\n",
                               0., GRID_LENGTH, 0., GRID_LENGTH, 0., GRID_LENGTH);

    fprintf(fileValues,
            "    </Piece>\n"
            "  </StructuredGrid>\n"
            "</VTKFile>\n"
            );
    fclose(fileValues);
}

void writeParticleData(const char* fileName, const Particle* particleData, int particleCount)
{
    FILE* fileValues = fopen(fileName, "w");
    int i;

    // write the particle data header
    fprintf(fileValues, "%20s %20s %20s %20s %20s %20s %20s %20s\n", "x", "y", "z", "Vx", "Vy", "Vz", "mass", "charge");

    for(i = 0; i < particleCount; i++)
    {
        const Particle* p = &particleData[i];

        fprintf(fileValues, "%20.8e %20.8e %20.8e %20.8e %20.8e %20.8e %20.8e %20.8e\n",
                              p->x,  p->y,  p->z, p->Vx, p->Vy,  p->Vz, p->mass, p->charge);
    }

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
        double ty = fabs(spacing*j - center[0]);

        for(k = 0; k < numNodes; k++)
        {
            //double z = spacing * k;
            double tz = fabs(spacing*k - center[1]);
            //double sumSqs = ty*ty + tz*tz;

            // CAPILLARY side
            // we need points OUTSIDE the capillary for Neumann BC
            if( (ty > capillaryRadius) || (tz > capillaryRadius) )
            {
                bNodes[nodeCount].bndryNodes[0] = &grid[0][j][k];
                bNodes[nodeCount].bndryNodes[1] = &grid[1][j][k];
                bNodes[nodeCount].bndryNodes[2] = &grid[2][j][k];

                nodeCount++;
            }

            // EXTRACTOR side
            if( ((ty <= extractorInner) && (tz <= extractorInner))
                               || 
                ((ty > extractorOuter) || (tz > extractorOuter))
              )
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

double enforceNeumannBC(BoundaryNode* bNodes, const int nodeCount, double sorOmega)
{
    int i;
    double norm = 0.;

    for(i = 0; i < nodeCount; i++)
    {
        double temp = (bNodes[i].bndryNodes[0])->potential;
        double val = (1./3) * (4 * (bNodes[i].bndryNodes[1])->potential - (bNodes[i].bndryNodes[2])->potential);

        (bNodes[i].bndryNodes[0])->potential = val*sorOmega + (1-sorOmega)*temp;

        temp = (bNodes[i].bndryNodes[0])->potential - temp;
        norm += temp*temp;
    }

    return norm;

}
