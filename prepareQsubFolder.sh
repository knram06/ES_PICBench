#!/bin/sh

SCRATCH=~/scratch
QSUB_SCRATCH=$SCRATCH/qsub_runs/PICC_Code
CP_PATH=$QSUB_SCRATCH/$1

MAIN_DIR_FILES="sub.pbs main.c *.h"

# create the source and include folders
mkdir $CP_PATH/

# copy over the relevant files
cp -r solvers/ $CP_PATH/
cp -r output/ $CP_PATH/
cp -r input_data/ $CP_PATH/
cp $MAIN_DIR_FILES $CP_PATH/

