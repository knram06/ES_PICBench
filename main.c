#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
// strtok issues resolved with including this header
#include <string.h>   

#define GRID_LENGTH 1.
#define NUM_NODES 101
//#define SIZE_X 101
//#define SIZE_Y 101
//#define SIZE_Z 101


/* geometry dimension */
// capillary centered on YZ face
// origin at Corner of domain
// --> Particle data MUST be corrected for new origin
#define CAPILLARY_RADIUS 0.04
#define EXTRACTOR_INNER_RADIUS 0.2
#define EXTRACTOR_OUTER_RADIUS (1.2 * EXTRACTOR_INNER_RADIUS)

// define all problem parameters in terms of macros
#define MD_FILE "./input_data/MD_data/10_input_Pos_Q488_20130318.inp"

#define TEST_FUNCTION (x*x - 2*y*y + z*z)
//#define TEST_FUNCTION 0.

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
    int    numNodes;
    double spacing, invSpacing;
} GridInfo;

unsigned int countLinesInFile(FILE* fp);
void parseMDFileToParticles(Particle particleData[], FILE* fp);
void allocateGrid(Node**** grid, GridInfo* gInfo);
void deAllocGrid(Node**** grid, GridInfo* gInfo);
void setupBoundaryConditions(Node*** grid, GridInfo* gInfo);
void solve(Node*** grid, GridInfo* gInfo, double tolerance, double sorOmega);

void enforceNeumannBC(Node*** grid, GridInfo* gInfo)
{
    const int numNodes               = gInfo->numNodes;
    const double spacing             = gInfo->spacing;

    const double capillaryRadius     = CAPILLARY_RADIUS;

    // extractor dimensions
    const double extractorInner      = EXTRACTOR_INNER_RADIUS;
    const double extractorOuter      = EXTRACTOR_OUTER_RADIUS;

    int i, j, k;
    // Y-Z face with capillary
    // i.e. X = 0 position with center at corner
    for(j = 0; j < numNodes; j++)
    {
        double y = spacing * j;

        for(k = 0; k < numNodes; k++)
        {
            double z = spacing * k;
        }
    }

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
    unsigned int lineCount = countLinesInFile(fp);
    printf("Number of lines read was %d.\n", lineCount);

    // fill in GridInfo data
    GridInfo gridInfo;
    gridInfo.numNodes = NUM_NODES;
    gridInfo.spacing = GRID_LENGTH / (NUM_NODES - 1);
    gridInfo.invSpacing = 1./(gridInfo.spacing);        // just caching this to avoid divisions?

    // now preallocate the particles data array
    Particle* MD_data = malloc(lineCount * sizeof(Particle));

    // loop through the file and tokenize entries
    parseMDFileToParticles(MD_data, fp);
    fclose(fp);

    // allocate the grid
    Node*** grid = NULL;
    allocateGrid(&grid, &gridInfo);

    // setup boundary conditions
    setupBoundaryConditions(grid, &gridInfo);

    // enforce boundary conditions
    // impose Neumann BCs
    // solve and at each step, impose Neumann BCs?
    solve(grid, &gridInfo, 1e-12, 1.9);


    deAllocGrid(&grid, &gridInfo);
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
    int count = 0;
    // buffer to read in line data
    char line[256];
    
    // rewind to the beginning
    rewind(fp);

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
            if(fieldIndex == 1)
                particleData[particleCount].y = atof(token);
            else if(fieldIndex == 2)
                particleData[particleCount].z = atof(token);
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

            // increment fieldCounter
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

void deAllocGrid(Node**** grid, GridInfo* gInfo)
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
    const int numNodes = gInfo->numNodes;
    const double spacing = gInfo->spacing;
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
        y = spacing*j;
        for(k = 0; k < numNodes; k++)
        {
            z = spacing * k;
            x = 0.;
            grid[0][j][k].potential=TEST_FUNCTION;
            x = GRID_LENGTH;
            grid[numNodes - 1][j][k].potential=TEST_FUNCTION;
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
void solve(Node*** grid, GridInfo* gInfo, double tolerance, double sorOmega)
{
    const int numNodes = gInfo->numNodes;

    double norm = 100.; double temp = 0., val = 0.;
    int i, j, k;
    int count = 0;

    const double inv6 = 1./6;

    // solve using criteria
    // solution in the "inner" cube
    while (norm >= tolerance)
    {
        count++;
        norm = 0.;

        // iterate across all points
        for (i = 1; i < numNodes - 1; i++)
        {
            for(j = 1; j < numNodes - 1; j++)
            {
                for(k = 1; k < numNodes - 1; k++)
                {
                    temp = grid[i][j][k].potential;

                    // obtain the value got by averaging
                    val = inv6 * (grid[i+1][j][k].potential + grid[i-1][j][k].potential + grid[i][j-1][k].potential + grid[i][j+1][k].potential + grid[i][j][k-1].potential + grid[i][j][k+1].potential);

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
        norm = sqrt(norm);
        printf("norm: %e\n", norm);
    }
    //imposeNeumannBCs();
}

