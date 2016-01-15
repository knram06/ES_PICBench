#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
// strtok issues resolved with including this header
#include <string.h>
#include <limits.h>
#include <time.h>

#define GRID_LENGTH (3e-4)
//#define NUM_NODES 101

/*Macro for 3D to 1D indexing */
//#define GRID_1D(grid, i, j, k) ( grid[(k) + NUM_NODES*(j) + NUM_NODES*NUM_NODES*(i) ] )
#define INDEX_1D(n, i, j, k) ( (k) + (n)*(j) + (n)*(n)*(i) )
#define GRID_1D(grid, i, j, k) (grid[ INDEX_1D(NUM_NODES, i, j, k) ] )

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
#define FREE_SPACE_PERMITTIVITY (8.8542e-12)

// array bounds and margins
#define LOST_PARTICLES_MARGIN 50

// define all problem parameters in terms of macros
#define MD_FILE "./input_data/MD_data/10_input_Pos_Q488_20130318.inp"

//#define PARTICLE_SIZE ((int)5e4)
#define PARTICLE_SORT_INTERVAL (20)

#define MAX_ITER (200)
#define TIMESTEPS ((int)0)
#define ITER_INTERVAL (200)
#define ITER_HEADER_INTERVAL (1500)
#define POST_WRITE_FILES (true)
#define POST_INTERVAL (200)
#define POST_WRITE_PATH ("output/")

#define POISSON_TIMESTEPS ((int)2000)
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
    double* bndryNodes[3];
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
#include "csrroutines.h"
#include "particleFns.h"
#include "numerics.h"
#include "postprocess.h"

#include "solvers/multigrid/mg_3d.h"

void allocateEField(EField** grid, GridInfo* gInfo)
{
    const int totalNodes = gInfo->totalNodes;

    (*grid) = malloc(totalNodes * sizeof(EField));
    assert((*grid) != NULL);

}

void deallocEField(EField** grid)
{
    free((*grid));
}

int main(int argc, char **argv)
{
    // read in MD data
    FILE *fp = fopen(MD_FILE, "r");
    if(fp == NULL)
    {
        fprintf(stderr, "Error in opening file %s \n", MD_FILE);
        return EXIT_FAILURE;
    }

    SolverInitialize(argc, argv);

    // count the number of lines in the input file
    // so that we can preallocate later
    unsigned int particleCount = countLinesInFile(fp);
    printf("Number of lines read was %d.\n", particleCount);

    // allocate the grid
    double *grid = NULL, *rhs = NULL;
    double h;
    int finestGridNum = SolverGetDetails(&grid, &rhs, &h);

    // fill in GridInfo data
    GridInfo gridInfo;
    gridInfo.numNodes = finestGridNum;
    gridInfo.totalNodes = finestGridNum*finestGridNum*finestGridNum;
    gridInfo.spacing = h;
    gridInfo.invSpacing = 1./h;

    // now preallocate the particles data array
    Particle* MD_data = malloc(particleCount * sizeof(Particle));

    // loop through the file and tokenize entries
    rewind(fp);         // rewind file to the beginning
    parseMDFileToParticles(MD_data, fp);
    fclose(fp);

    // TODO: Uncomment and fix this
    EField* ElectricField = NULL;
    allocateEField(&ElectricField, &gridInfo);

    //// preallocate NeumannBC nodes in terms of MAX possible
    //// i.e. num of sides of cube * num nodes per side
    //BoundaryNode *bNodes = malloc( 6 * (NUM_NODES)*(NUM_NODES) * sizeof(BoundaryNode) );

    //// calculate the actual BoundaryNodes count and resize the array to that instead
    //// to avoid wastage
    //printf("Consolidating Neumann BC nodes into a different structure....");
    //// setup boundary conditions
    //int nodeCount = setupBoundaryConditions(grid, &gridInfo, bNodes);

    //bNodes = realloc(bNodes, nodeCount * sizeof(BoundaryNode));
    //printf("done\n");

    SolverSetupBoundaryConditions();
    // solve and at each step, impose Neumann BCs?
    // FOR TESTING ONLY
    double norm = SolverGetInitialResidual(), tolerance = 1e-6;

    // FMG Initialization
    //printf("Carrying out FMG Initialization.... ");
    //SolverFMGInitialize();
    //printf("done\n");

    SolverResetTimingInfo();
    Solve(norm, tolerance, MAX_ITER);
    SolverPrintTimingInfo();

    printf("\nCalculating Electric Field.....");
    calcElectricField(ElectricField, grid, &gridInfo);
    printf("done\n");
    //writeOutputData("laplace.vtk", grid, ElectricField, &gridInfo);
    //return 0;

    // calculate the release rate
    int Nrel;
    double Nfrac, runningNfrac, temp;
    double particleReleaseRate = particleCount * (T_PIC/ T_MD);

    // set the integer and fractional parts
    Nfrac = modf(particleReleaseRate, &temp);
    Nrel = (int)(temp);
    runningNfrac = Nfrac;

    // Array to keep track of lost particles
    //The idea behind is that, particles may leave the domain and this is an index
    //of such particles.
    //Using this, new particles can be inserted into those locations appropriately.
    int* lostParticles = malloc( (Nrel + LOST_PARTICLES_MARGIN) * sizeof(int) );
    int lostParticleBound = -1;

    // SINGLE PARTICLE release check for preallocating
    int iterCountTillExit = simulateSingleParticleTillDomainExit(MD_data, particleCount, ElectricField, &gridInfo);
    int particlePreallocCount = 2*Nrel * iterCountTillExit;
    printf("%d Iterations were done for Particle to leave domain with Laplace solution.\nReleasing %d particles per timesteps, approx total particles to preallocate is %d (doubled - covering twice the domain length) \n", iterCountTillExit, Nrel, particlePreallocCount);


    // allocate for particles
    Particle* domainParticles = malloc(particlePreallocCount * sizeof(Particle));
    int totalParticlesCount = 0;

    clock_t start = clock(), diff;
    // for required number of timesteps
    int i;
    int lostParticleCount = 0;
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
            sprintf(outputPath, "%s/particleOutput/particleOutput_%d.txt", POST_WRITE_PATH, i);
            writeParticleData(outputPath, domainParticles, totalParticlesCount);
        }
    } // end of Laplace loop
    //return 0;

    diff = clock() - start;
    double timeStepsTime = diff /CLOCKS_PER_SEC;

    // TEMP - remove later
    resortParticles(domainParticles, totalParticlesCount);

    /*********************************************/
    /***********POISSON SOLVER********************/
    /*********************************************/
    //double *rhs = malloc(sizeof(double) * 8 * totalParticlesCount);
    //int *rhsIndices = malloc(sizeof(int) * 8 * totalParticlesCount);
    //ColVal *indVals = malloc(sizeof(ColVal) * 8 * totalParticlesCount);

    //FILE *iterData = fopen("iter_data.txt", "w");
    start = clock();
    for(i = 1; i <= POISSON_TIMESTEPS; i++)
    {
        clock_t tstart = clock();
        printf("\nPoisson Timestep %d:\n", i);

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

        // Re-sort Particles based on their distance from the capillary center
        // improves locality and some searches
        // BUT don't do this every iteration, since we will use qsort
        // and QUICKSORT WORST CASE performance is when array is ALMOST SORTED - O(n^2)
        if( !(i % PARTICLE_SORT_INTERVAL) )
            resortParticles(domainParticles, totalParticlesCount);

        // reset the rhs vector
        resetRHSInteriorPoints(rhs, &gridInfo);       // resets only interior points
        // based on the new particle positions, update charge fractions at nodes
        updateChargeFractions(domainParticles, totalParticlesCount, rhs, &gridInfo);

        // SOLVE here
        norm = SolverGetResidual();
        Solve(norm, tolerance, MAX_ITER);

        // update the electric field
        calcElectricField(ElectricField, grid, &gridInfo);

        // then move the particles
        lostParticleCount = moveParticlesInField(domainParticles, totalParticlesCount,
                                                 lostParticles,//&lostParticleBound,
                                                 ElectricField, &gridInfo);
        printf("%d Particles left the domain\n", lostParticleCount);
        //printf("%d empty slots in the domain particles\n", lostParticleBound+1);

        // IMPORTANT: add Nfrac to runningNfrac to adjust correctly
        // for the fractional part
        runningNfrac += Nfrac;

        // avoid timing the IO portion
        diff = clock() - tstart;
        double timeTaken = (double)diff/CLOCKS_PER_SEC;
        //fprintf(iterData, "%10d %10lf\n", iterNum, timeTaken);

        char outputPath[50];
        if(POST_WRITE_FILES && !(i % POST_INTERVAL) )
        {
            sprintf(outputPath, "%s/particleOutput/particleOutput_%d.txt", POST_WRITE_PATH, i);
            writeParticleData(outputPath, domainParticles, totalParticlesCount);

            sprintf(outputPath, "%s/solutions/out_%d.vtk", POST_WRITE_PATH, i);
            writeOutputData(outputPath, grid, ElectricField, &gridInfo);
        }
    } // end of Poisson timestep loop
    //fclose(iterData);

    writeOutputData("poisson.vtk", grid, ElectricField, &gridInfo);
    //getRHS(rhs);
    //writeVectorToFile("poisson_v.txt", rhs, gridInfo.totalNodes);

    diff = clock() - start;
    double poissonStepsTime = diff /CLOCKS_PER_SEC;
    //printf("\nTiming Info\n%10s %10.8e\n%10s %10.8e\n%10s %10.8e\n", "Solve", solveTime, "TimeSteps", timeStepsTime, "Poisson TimeSteps", poissonStepsTime);


    //free(rhsIndices);
    //free(rhs);
    //deallocCSRForm(&mcsr);
    SolverFinalize();

    free(lostParticles);
    free(domainParticles);
    //free(bNodes);
    deallocEField(&ElectricField);
    //deallocGrid(&grid);
    free(MD_data);

    return 0;
}

