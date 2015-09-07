#ifndef POSTPROCESS_H
#define POSTPROCESS_H

// function for writing out values
void writeOutputData(const char* fileName, Node* grid, EField* ElectricField, GridInfo* gInfo)
{
    FILE* fileValues = fopen(fileName, "w");
    int i, j, k;

    const int numNodes = gInfo->numNodes;
    const double spacing = gInfo->spacing;

    const int totalNodes = gInfo->totalNodes;

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
                potentialValues[count] = GRID_1D(grid, i,j,k).potential;

                EField* elecField = &GRID_1D(ElectricField, i, j, k);
                // update the electric field array
                const int pos = 3*count;
                electricField[pos]     = elecField->components[0];
                electricField[pos + 1] = elecField->components[1];
                electricField[pos + 2] = elecField->components[2];

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
#endif
