#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=6:00:00
#SBATCH --job-name ClusterBarcodes

# Script to find each of the individual read count files, and process
# to collapse barcodes with a Hamming distance of two or less.
#
# 11 June 2019
#

# cd into working directory
cd $SLURM_SUBMIT_DIR

# Load required modules
module load CCEnv
module load StdEnv/2020
module load gcc/9.3.0
module load r/4.0.2
module load openmpi/4.0.3

# Run R script
ls *.txt | parallel -j 4 Rscript --no-restore collapseCloseBarcodesParallel3.R {} {.}_collapsed.txt

