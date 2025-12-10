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
### Script ...
