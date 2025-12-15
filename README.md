# scRNAseq_pICtimecourse
_scRNA-seq experiment of A549 cells either unstimulated with pIC or stimulated with pIC for 4, 8, 12, or 16 hours._
Scripts associated with the ATAC-seq data in Single-cell heterogeneity in interferon induction potential is heritable and governed by variation in cell state, Thayer et al. (2025).

## Overview
This repository contains the scripts used to process and analyze the scRNA-seq data collected from A549s wither unstimulated or stimulated for 4, 8, 12, or 16 hours with pIC. The goal of this work is similar to https://github.com/BROOKELAB/Temporal-scRNA-seq-H3N2, but instead of using virus as a stimulus, we use pIC. Processing of the data is slightly different than previous sequencing experiments done by our group because we used the Chromium Next GEM Single Cell FixedRNA Sample Preparation Kit, and this required pooling the samples after fixation and hybridization. 

## Preliminary Processing
Raw reads were demultiplexed using bcl-convert v4.1.7 Conversion Software (Illumina). The reads could then be processed using the following:
1. Cell Ranger **multi**
  - Input: .fastq files with locations in pool1_input.csv, pool2_input.csv, pool3_input.csv
  - Output: .h5 files for each pool
  - Script: run_cellranger.bash

2. Cell Ranger **aggr**
  - Input: .h5 files with locations in aggr_timepoints.csv
  - Output: One .h5 file with all of the experimental conditions data combined
  - Script: aggr_timepoints_092424.bash

## scRNA-seq Analysis Requiring a Supercomputer
### et02b04_initialanalysis_12.12.25.R
This script performs:\
:white_check_mark: Initial upload of the .h5 file and processing including normalization and transformation\
:white_check_mark: Saving the .rds file after processing and export of most data to a .csv to bypass the time-consuming steps\
:white_check_mark: Upload the supplemantary file, barcodes_umap_probabilities_info.csv, and join the two csv files for further analysis\
:white_check_mark: More data reconfiguration to then perform the Mann-Whitney analysis to find top geneswhose expression influences the probability of transitioning to Terminal 1\
:white_check_mark: Cutoff determination and data visualization to make sense of the distributions
### Requirements
**R (v4.5.2)**
- Seurat (v5.3.1)
- dittoSeq (v1.22.0)
- tidyverse (v2.0.0)
- sctransform (v0.4.2)
- glmGamPoi (v1.22.0)
- DropletUtils (v1.30.0)
- simpleSingleCell (v1.34.0)
- scater (v1.38.0)
- scran (v1.38.0)
- BiocSingular (v1.26.0)
- magrittr (v2.0.4)
- dplyr (v1.1.4)
- ggplot2 (v4.0.1)
- ggrepel (v0.9.6)
- data.table (v1.17.8)

## scRNA-seq Analysis Post-Export as a Smaller Datafile (_no need for supercomputer_)













