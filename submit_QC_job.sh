#!/bin/bash
# Written by Nina Gilshteyn, M.S.
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=4:00:00,h_data=40G,h_vmem=infinity

## Modify the parallel environment
## and the number of cores as needed:

#$ -pe shared 1
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea


# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/4.9.3

module load anaconda3
. $CONDA_DIR/etc/profile.d/conda.sh
conda activate project

## $VAR is the raw gene expression CSV filename passed in from bash_loop_submitter.sh

python3 NEW_get_QC_plots.py $VAR


# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
