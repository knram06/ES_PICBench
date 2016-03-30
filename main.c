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
#define TIMESTEPS ((int)2000)
#define ITER_INTERVAL (200)
#define ITER_HEADER_INTERVAL (1500)
#define POST_WRITE_FILES (false)
#define POST_INTERVAL (200)
#define POST_WRITE_PATH ("output/")

#define POISSON_TIMESTEPS ((int)300)
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
//#include "csrroutines.h"
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

    SolverSetupBoundaryConditions();
    // solve and at each step, impose Neumann BCs?
    double tolerance = 1e-6;

    // FMG Initialization
    //printf("Carrying out FMG Initialization.... ");
    //SolverFMGInitialize();
    //printf("done\n");

    const int maxThreads = omp_get_max_threads(); // test with 8 in gcc -g mode in gdb
    printf("Max threads: %d\n", maxThreads);

    double *threadNorm = calloc(maxThreads, sizeof(double));
    // HACKY WAY - fill in the first element of threadNorm to seed the initial norm
    // calculation
    threadNorm[0] = SolverGetInitialResidual();
    SolverResetTimingInfo();
    #pragma omp parallel
    { Solve(tolerance, MAX_ITER, threadNorm); }
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
    int* lostParticles = malloc( (2*Nrel) * sizeof(int) );
    int lostParticleBound = -1;

    // SINGLE PARTICLE release check for preallocating
    int iterCountTillExit = simulateSingleParticleTillDomainExit(MD_data, particleCount, ElectricField, &gridInfo);
    int particlePreallocCount = 2*Nrel * iterCountTillExit;
    printf("%d Iterations were done for Particle to leave domain with Laplace solution.\nReleasing %d particles per timesteps, approx total particles to preallocate is %d (doubled - covering twice the domain length) \n", iterCountTillExit, Nrel, particlePreallocCount);

    // allocate for particles
    Particle* domainParticles = malloc(particlePreallocCount * sizeof(Particle));
    int totalParticlesCount = 0;

    // for required number of timesteps
    int i;
    int lostParticleCount = 0;
    int numParticlesToRelease;

    // store for one extra space, i.e. with zero index
    int *localLostParticlesCount = calloc(maxThreads+1, sizeof(int));
    int *threadOffsetLostParticles = &(localLostParticlesCount[1]);

    double timingTemp;
    TimingInfo *tInfo = NULL;
    const char* stageNames[5] = {"ReleaseParticles", "SwapGaps", "UpdateChargeFrns", "MoveParticles", "ThreadUpdates"};
    allocTimingInfo(&tInfo, stageNames, 4);

    double start = omp_get_wtime();
    #pragma omp parallel private(i) //num_threads(8)
    {
        int t;                                  // per thread loop counter
        int tid = omp_get_thread_num();         // thread id

        // per thread space for storing lost particles
        int *localLostParticles = calloc(2*Nrel, sizeof(int));

        // SEEDS for thread-safe random number generation
        // using malloc so that we get a random seed - based on whatever
        // the value is at startup
        unsigned int *randSeeds = malloc(maxThreads*sizeof(unsigned int));

    for(i = 1; i <= TIMESTEPS; i++)
    {
        // Implicit barrier is NECESSARY! So that releaseParticles sees the calculation update!!
        #pragma omp single
        {
        printf("\nTimestep %d:\n", i);

        // calculate the current timestep's release rate
        runningNfrac = modf(runningNfrac, &temp);
        numParticlesToRelease = Nrel + (int)(temp);
        lostParticleBound = (lostParticleCount - 1);        // adjust for one off issue
        }

        #pragma omp master
        timingTemp = omp_get_wtime();

        // introduce the particles
        releaseParticles(numParticlesToRelease,
                         MD_data, particleCount,
                         domainParticles, &totalParticlesCount,
                         lostParticles,
                         &lostParticleBound, // adjust for one off issue
                         randSeeds);
        #pragma omp master
        {
        tInfo->timeTaken[0] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[0]++;

        timingTemp = omp_get_wtime();
        }

        swapGapsWithEndParticles(domainParticles, &totalParticlesCount,
                                 lostParticles, &lostParticleBound);
        #pragma omp master
        {
        tInfo->timeTaken[1] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[1]++;

        //# pragma omp single
        printf("Total Number of Particles: %d\n", totalParticlesCount); // need +1 for the one-off offset
        timingTemp = omp_get_wtime();
        }

        resetRHSInteriorPoints(rhs, &gridInfo);       // resets only interior points

        #pragma omp master
        timingTemp = omp_get_wtime();

        updateChargeFractions(domainParticles, totalParticlesCount, rhs, &gridInfo);

        #pragma omp master
        {
        tInfo->timeTaken[2] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[2]++;

        timingTemp = omp_get_wtime();
        }

        // then move them
        threadOffsetLostParticles[tid] = moveParticlesInField(
                                domainParticles, totalParticlesCount,
                                localLostParticles,//&lostParticleBound,
                                ElectricField, &gridInfo);
        #pragma omp master
        {
        tInfo->timeTaken[3] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[3]++;
        }

        // form a cumulative sum for the localLostParticlesCount
        // may be negligible for small number of threads - not timing this
        #pragma omp barrier     // IMPORTANT!!
        #pragma omp single
        {
            for(t = 1; t < maxThreads; t++)
                threadOffsetLostParticles[t] += threadOffsetLostParticles[t-1];
        }

        #pragma omp master
        timingTemp = omp_get_wtime();

        // now using the cumulative sum array, copy over from shared local lost particles
        // arrays to the global one
        {
            int localCount = localLostParticlesCount[tid];
            for(t = localLostParticlesCount[tid]; t < localLostParticlesCount[tid+1]; t++)
                lostParticles[t] = localLostParticles[t - localCount];

            #pragma omp barrier // ensure all threads finish the update, done for the timing
            #pragma omp master
            {
            tInfo->timeTaken[4] += (omp_get_wtime() - timingTemp);
            tInfo->numCalls[4]++;
            }
        }

        #pragma omp single // just so we can use the implicit barrier
        {
            // update the lostParticleCount - which is the last entry in the
            // cumulative array
            lostParticleCount = threadOffsetLostParticles[maxThreads-1];

            printf("%d Particles left the domain\n", lostParticleCount);
            // IMPORTANT: add Nfrac to runningNfrac to adjust correctly
            // for the fractional part
            runningNfrac += Nfrac;

            char outputPath[50];
            if(POST_WRITE_FILES && !(i % POST_INTERVAL) )
            {
                sprintf(outputPath, "%s/particleOutput/particleOutput_%d.txt", POST_WRITE_PATH, i);
                writeParticleData(outputPath, domainParticles, totalParticlesCount);
            }
        }
    } // end of Laplace loop
    free(randSeeds);
    free(localLostParticles);
    }
    free(localLostParticlesCount);
    printTimingInfo(tInfo);

    double diff = omp_get_wtime() - start;
    printf("Laplace steps time: %10.8lf\n", diff);


    /*********************************************/
    /***********POISSON SOLVER********************/
    /*********************************************/
    /*
    const char *stageNames[8] = {"Release Particles", "SwapGaps", "ReSort", "ResetRHS", "UpdateChargeFrns", "Solve","calcElecField", "moveParticlesInField"};
    allocTimingInfo(&tInfo, stageNames, 8);

    // reset the Solver Timing Info
    SolverResetTimingInfo();

    double start = omp_get_wtime();
    #pragma omp parallel private(i)
    {
        int t;                                  // per thread loop counter
        int tid = omp_get_thread_num();         // thread id

        // per thread space for storing lost particles
        int *localLostParticles = calloc(Nrel + LOST_PARTICLES_MARGIN, sizeof(int));

        // SEEDS for thread-safe random number generation
        unsigned int *randSeeds = calloc(maxThreads, sizeof(unsigned int));

    for(i = 1; i <= POISSON_TIMESTEPS; i++)
    {
        //clock_t tstart = clock();
        #pragma omp single
        {
        printf("\nPoisson Timestep %d:\n", i);

        // calculate the current timestep's release rate
        runningNfrac = modf(runningNfrac, &temp);
        numParticlesToRelease = Nrel + (int)(temp);

        //totalParticlesBound += numParticlesToRelease - lostParticleCount;
        lostParticleBound = (lostParticleCount - 1);        // adjust for one off issue
        }

        #pragma omp master
        timingTemp = omp_get_wtime();

        // introduce the particles
        releaseParticles(numParticlesToRelease,
                         MD_data, particleCount,
                         domainParticles, &totalParticlesCount,
                         lostParticles,
                         &lostParticleBound,          // adjust for one off issue
                         randSeeds);

        #pragma omp master
        {
        tInfo->timeTaken[0] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[0]++;

        //totalParticlesBound = (totalParticlesCount - 1);
        timingTemp = omp_get_wtime();
        }
        swapGapsWithEndParticles(domainParticles, &totalParticlesCount,
                                 lostParticles, &lostParticleBound);
        #pragma omp master
        {
        tInfo->timeTaken[1] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[1]++;

        printf("Total Number of Particles: %d\n", totalParticlesCount); // need +1 for the one-off offset

        // Re-sort Particles based on their distance from the capillary center
        // improves locality and some searches
        // BUT don't do this every iteration, since we will use qsort
        // and QUICKSORT WORST CASE performance is when array is ALMOST SORTED - O(n^2)
        if( !(i % PARTICLE_SORT_INTERVAL) )
        {
            timingTemp = omp_get_wtime();
            resortParticles(domainParticles, totalParticlesCount);
            tInfo->timeTaken[2] += (omp_get_wtime() - timingTemp);
            tInfo->numCalls[2]++;
        }

        // reset the rhs vector
        timingTemp = omp_get_wtime();
        } // end of OMP MASTER

        resetRHSInteriorPoints(rhs, &gridInfo);       // resets only interior points

        #pragma omp master
        {
        tInfo->timeTaken[3] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[3]++;

        // based on the new particle positions, update charge fractions at nodes
        timingTemp = omp_get_wtime();
        }

        //#pragma omp single
        updateChargeFractions(domainParticles, totalParticlesCount, rhs, &gridInfo);

        #pragma omp master
        {
        tInfo->timeTaken[4] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[4]++;
        }

        threadNorm[tid] = SolverGetResidual();
        #pragma omp barrier     // VERY IMPORTANT! to ensure threadNorm array is correctly filled by all threads before moving on to calculating the norm
        #pragma omp master
        timingTemp = omp_get_wtime();

        Solve(tolerance, MAX_ITER, threadNorm);

        #pragma omp master
        {
        tInfo->timeTaken[5] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[5]++;

        // update the electric field
        timingTemp = omp_get_wtime();
        }
        calcElectricField(ElectricField, grid, &gridInfo);

        #pragma omp master
        {
        tInfo->timeTaken[6] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[6]++;

        // then move the particles
        timingTemp = omp_get_wtime();
        }

        threadOffsetLostParticles[tid] = moveParticlesInField(
                                domainParticles, totalParticlesCount,
                                localLostParticles,//&lostParticleBound,
                                ElectricField, &gridInfo);

        #pragma omp master
        {
        tInfo->timeTaken[7] += (omp_get_wtime() - timingTemp);
        tInfo->numCalls[7]++;
        }

        // form a cumulative sum for the localLostParticlesCount
        #pragma omp barrier     // IMPORTANT!!
        #pragma omp single
        {
            for(t = 1; t < maxThreads; t++)
                threadOffsetLostParticles[t] += threadOffsetLostParticles[t-1];
        }

        // now using the cumulative sum array, copy over from shared local lost particles
        // arrays to the global one
        {
            int localCount = localLostParticlesCount[tid];
            for(t = localLostParticlesCount[tid]; t < localLostParticlesCount[tid+1]; t++)
                lostParticles[t] = localLostParticles[t - localCount];
        }

        #pragma omp single // just so we can use the implicit barrier
        {
            // update the lostParticleCount - which is the last entry in the
            // cumulative array
            lostParticleCount = threadOffsetLostParticles[maxThreads-1];

            printf("%d Particles left the domain\n", lostParticleCount);
            // IMPORTANT: add Nfrac to runningNfrac to adjust correctly
            // for the fractional part
            runningNfrac += Nfrac;

            // avoid timing the IO portion
            //diff = clock() - tstart;
            //double timeTaken = (double)diff/CLOCKS_PER_SEC;
            //fprintf(iterData, "%10d %10lf\n", iterNum, timeTaken);

            char outputPath[50];
            if(POST_WRITE_FILES && !(i % POST_INTERVAL) )
            {
                sprintf(outputPath, "%s/particleOutput/particleOutput_%d.txt", POST_WRITE_PATH, i);
                writeParticleData(outputPath, domainParticles, totalParticlesCount);

                sprintf(outputPath, "%s/solutions/out_%d.vtk", POST_WRITE_PATH, i);
                writeOutputData(outputPath, grid, ElectricField, &gridInfo);
            }
        }
    } // end of Poisson timestep loop
    free(randSeeds);
    free(localLostParticles);
    }
    free(localLostParticlesCount);
    free(threadNorm);
    //fclose(iterData);
    double diff = omp_get_wtime() - start;
    //double poissonStepsTime = diff /CLOCKS_PER_SEC;
    printf("Poisson steps time: %10.8lf\n", diff);

    SolverPrintTimingInfo();
    printTimingInfo(tInfo);
    */
    //writeOutputData("poisson.vtk", grid, ElectricField, &gridInfo);
    //writeVectorToFile("poisson_v.txt", rhs, gridInfo.totalNodes);

    //printf("\nTiming Info\n%10s %10.8e\n%10s %10.8e\n%10s %10.8e\n", "Solve", solveTime, "TimeSteps", timeStepsTime, "Poisson TimeSteps", poissonStepsTime);

    //free(rhsIndices);
    //free(rhs);
    //deallocCSRForm(&mcsr);
    deAllocTimingInfo(&tInfo);
    SolverFinalize();

    free(lostParticles);
    free(domainParticles);
    //free(bNodes);
    deallocEField(&ElectricField);
    //deallocGrid(&grid);
    free(MD_data);

    return 0;
}

