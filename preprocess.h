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

void allocateGrid(double** grid, GridInfo* gInfo)
{
    const int totalNodes = gInfo->totalNodes;

    (*grid) = malloc(totalNodes * sizeof(double));
    assert((*grid) != NULL);
}

void deallocGrid(double** grid)
{
    free((*grid));
}

#endif
