#!/bin/bash
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user ethayer3@illinois.edu
#SBATCH -o pt2aggr_logging.out.log
#SBATCH -e pt2aggr_error.err.log
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G

module load cellranger/7.1.0

cellranger aggr --id=aggr_et02b04_092424 --csv=aggr_timepoints.csv
