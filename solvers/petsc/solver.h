#ifndef SOLVER_H
#define SOLVER_H

#include <petscksp.h>
#include <petscmat.h>

// store relevant data as global for better access
Mat A;
Vec b, x;

void SolverInitialize(int *argc, char ***argv)
{
    PetscInitialize(argc, argv, NULL, NULL);
    MatCreate(PETSC_COMM_WORLD, &A);
}

void buildSolverMatrixFromCSR(const int *rowOffsets, const int *colIndices, const double *vals, const int numRows)
{
    // the last entry in the rowOffsets will actually be the total nonzeros
    // so just use that
    const int totalNonZeros = rowOffsets[numRows]; // takes care of one off indexing
    const int numRowsOffset = numRows+1;

    // declare solver specific stuff here
    PetscInt *rows = malloc(sizeof(PetscInt) * numRowsOffset);
    PetscInt *cols = malloc(sizeof(PetscInt) * totalNonZeros);
    PetscScalar *values = malloc(sizeof(PetscScalar) * totalNonZeros);

    // now loop through and copy over stuff - costly?
    int i;
    for(i = 0; i < numRowsOffset; i++)
        rows[i] = (PetscInt)(rowOffsets[i]);

    for(i = 0; i < totalNonZeros; i++)
    {
        cols[i]   = (PetscInt)(colIndices[i]);
        values[i] = (PetscScalar)(vals[i]);
    }

    // build up the matrix
    MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, (PetscInt)numRows, (PetscInt)numRows, rows, cols, values, &A);

    MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    free(values);
    free(cols); free(rows);
}

void SolverFinalize()
{
    MatDestroy(&A);
    PetscFinalize();
}

#endif
