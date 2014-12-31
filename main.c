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

// define all problem parameters in terms of macros
#define MD_FILE "./input_data/MD_data/10_input_Pos_Q488_20130318.inp"


// declare the Particle data type here
typedef struct
{
    double x, y, z;
    double Vx, Vy, Vz;
    double mass, charge;
} Particle;

typedef struct
{
    double pos[3];
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

    // enforce boundary conditions
    // impose Neumann BCs
    // solve and at each step, impose Neumann BCs?


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
    int i, j;

    (*grid) = malloc(numNodes * sizeof(Node**));
    for(i = 0; i < numNodes; i++)
    {
        (*grid)[i] = malloc(numNodes * sizeof(Node*));
        for(j = 0; j < numNodes; j++)
            (*grid)[i][j] = malloc(numNodes * sizeof(Node));
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

