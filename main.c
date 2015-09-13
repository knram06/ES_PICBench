#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
// strtok issues resolved with including this header
#include <string.h>
#include <time.h>

#define GRID_LENGTH (3e-4)
#define NUM_NODES 101

/*Macro for 3D to 1D indexing */
#define GRID_1D(grid, i, j, k) ( grid[(k) + NUM_NODES*(j) + NUM_NODES*NUM_NODES*(i)] )

/* geometry dimension */
// capillary centered on YZ face
// origin at Corner of domain
// --> Particle data MUST be corrected for new origin
#define CAPILLARY_RADIUS (1.326e-5)
#define EXTRACTOR_INNER_RADIUS (1e-4)
#define EXTRACTOR_OUTER_RADIUS (1.4e-4)
#define CAPILLARY_VOLTAGE 0.
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
#define TIMESTEPS ((int)8e3)
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
    int    numNodes, totalNodes;
    double spacing, invSpacing;
} GridInfo;

#include "utilities.h"
#include "preprocess.h"
#include "particleFns.h"
#include "numerics.h"
#include "postprocess.h"

void allocateEField(EField** grid, GridInfo* gInfo)
{
    const int totalNodes = gInfo->totalNodes;
    int i, j, k;

    (*grid) = malloc(totalNodes * sizeof(EField));
    assert((*grid) != NULL);

}

void deallocEField(EField** grid)
{
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
    gridInfo.totalNodes = gridInfo.numNodes * gridInfo.numNodes * gridInfo.numNodes;
    gridInfo.spacing = GRID_LENGTH / (NUM_NODES - 1);
    gridInfo.invSpacing = 1./(gridInfo.spacing);        // just caching this to avoid divisions?

    // now preallocate the particles data array
    Particle* MD_data = malloc(particleCount * sizeof(Particle));

    // loop through the file and tokenize entries
    rewind(fp);         // rewind file to the beginning
    parseMDFileToParticles(MD_data, fp);
    fclose(fp);

    // allocate the grid
    Node* grid = NULL;
    allocateGrid(&grid, &gridInfo);

    EField* ElectricField = NULL;
    allocateEField(&ElectricField, &gridInfo);

    // preallocate NeumannBC nodes in terms of MAX possible
    // i.e. num of sides of cube * num nodes per side
    BoundaryNode *bNodes = malloc( 6 * (NUM_NODES)*(NUM_NODES) * sizeof(BoundaryNode) );

    // calculate the actual BoundaryNodes count and resize the array to that instead
    // to avoid wastage
    printf("Consolidating Neumann BC nodes into a different structure....");
    // setup boundary conditions
    int nodeCount = setupBoundaryConditions(grid, &gridInfo, bNodes);

    bNodes = realloc(bNodes, nodeCount * sizeof(BoundaryNode));
    printf("done\n");

    // enforce boundary conditions
    // impose Neumann BCs
    // solve and at each step, impose Neumann BCs?
    double tolerance = 1e-9, sorOmega = 1.9;
    double norm = 100.;

    int iterCount = 1;
    //while(norm >= tolerance)
    clock_t start = clock(), diff;
    for(iterCount = 1; norm >= tolerance; iterCount++)
    {
        norm = sqrt(single_step_solve(grid, gridInfo.numNodes, sorOmega));

        // TODO: norm should be updated with this calculation no?
        // as it changes the values in the grid?
        enforceNeumannBC(bNodes, nodeCount);

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

    /*! Array to keep track of lost particles
     * The idea behind is that, particles may leave the domain and this is an index
     * of such particles.
     * Using this, new particles can be inserted into those locations appropriately.
     */
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


    // write out data for post processing
    writeOutputData("out.vtk", grid, ElectricField, &gridInfo);

    free(lostParticles);
    free(domainParticles);
    free(bNodes);
    deallocEField(&ElectricField);
    deallocGrid(&grid);
    free(MD_data);

    return 0;
}


