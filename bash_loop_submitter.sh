#!/bin/bash

# Written by Nina Gilshteyn, M.S.
# Loop through all raw gene expression CSV files and submit a QC job for each
for i in *csv; do
  echo "Submitting job for $i"
  qsub -v VAR="$i" submit_QC_job.sh  #this tells the computer to use the $i within the submit_QC_job.sh since I am not using it directly in shell
done
