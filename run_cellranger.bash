#!/bin/bash
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user ethayer3@illinois.edu
#SBATCH -o run_logging.out.log
#SBATCH -e run_error.err.log
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G

module purge
module load cellranger/7.1.0

cellranger multi --id pool_1 --csv input_csv/pool1_input.csv --localcores 32 --localmem 96  --localvmem 96 --disable-ui 
cellranger multi --id pool_2 --csv input_csv/pool2_input.csv --localcores 32 --localmem 96  --localvmem 96 --disable-ui 
cellranger multi --id pool_3 --csv input_csv/pool3_input.csv --localcores 32 --localmem 96  --localvmem 96 --disable-ui 
