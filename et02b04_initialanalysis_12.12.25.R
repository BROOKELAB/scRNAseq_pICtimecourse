## ------------------------------------------------------------------------------------------------------------
library(Seurat)
library(dittoSeq)
library(tidyverse)
library(sctransform)
library(glmGamPoi)
library(DropletUtils)
library(simpleSingleCell)
library(scater)
library(scran)
library(BiocSingular)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)


## ------------------------------------------------------------------------------------------------------------
boxplot_theme <- theme(legend.title = element_text(color = "black", size = 8, 
                                                  face = "bold"),
                       panel.background = element_rect(color ="black", linewidth = 0.5, fill = "white"),
                       legend.text = element_text(color = "black", size = 8),
                       legend.key = element_blank(),
                       axis.title.x = element_text(color = "black", size = 8, 
                                               face = "bold",
                                               margin = margin(t = 5, r = 0, b = 0, l = 0)),
                       axis.title.y = element_text(color = "black", size = 8, 
                                               face = "bold"),
                       axis.text.x = element_text(color = "black", size = 8, 
                                              margin = margin(t = 2.5, r = 5, b = 0, l = 5)),
                       axis.text.y = element_text(color = "black", size = 8),
                       axis.ticks = element_line(linewidth = 0.25),
                       panel.grid.minor = element_blank(), #gets rid of grey and lines in the middle
                       panel.grid.major = element_blank())


## ------------------------------------------------------------------------------------------------------------
# read in raw data
allcells_prefilt_counts <- Read10X_h5("/home/labs/cbrooke_lab/LizT/sequencingData_forETpaper/scRNAseq_et02b04/processeddata/et02b04_filtered_feature_bc_matrix.h5")

# make the Seurat object
allcells_prefilt <- CreateSeuratObject(counts = allcells_prefilt_counts, min.cells = 1, 
                                       min.features = 400, project = "pIC_timecourse", names.field = 2,
                                       names.delim = "-")
# what does it look like?
allcells_prefilt@meta.data


## ------------------------------------------------------------------------------------------------------------
# organize sample order
sample_order <- read.csv("/home/labs/cbrooke_lab/LizT/sequencingData_forETpaper/scRNAseq_et02b04/supplemental/aggr_timepoints.csv")

sample_order$sample_id[1:10] <- c("sample1a", "sample1b", "sampleaa", "sample2b", "sample3a", "sample3b", "sample4a", "sample4b", "sample5a", "sample5b")

sample_order$number <- as.character(1:nrow(sample_order))
levels(allcells_prefilt$orig.ident) <- sample_order$sample_id[order(sample_order$number)]
sample_order$time <- factor(sample_order$time, levels = c("0hr", "4hr", "8hr", "12hr", "16hr"))
sample_order <- sample_order %>% arrange(time)

# what does it look like now?
allcells_prefilt@meta.data


## ------------------------------------------------------------------------------------------------------------
# Filtering and normalization steps

# Remove mitochondria genes to get rid of dead cells
allcells_prefilt[["percent.mt"]] <- PercentageFeatureSet(allcells_prefilt, pattern = "^MT-")
allcells_prefilt <- subset(allcells_prefilt, subset = percent.mt < 25)

# define cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# normalize data
allcells_postnorm <- NormalizeData(allcells_prefilt, normalization.method = "LogNormalize", scale.factor = 10000)

# find top most variable genes
allcells_postvar <- FindVariableFeatures(allcells_postnorm, selection.method = "vst", nfeatures = 2000)

# from this data, we can classify what part of the cell cycle each cell was in
allcells_postvar_cycle <- CellCycleScoring(allcells_postvar, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# **longest part** tranform the data
allcells_posttrans <- SCTransform(allcells_postvar_cycle,
                        method = "glmGamPoi",
                        vars.to.regress = c("S.Score", "G2M.Score"),
                        return.only.var.genes = TRUE)


## ------------------------------------------------------------------------------------------------------------
#Run pca

allcells_posttrans <- RunPCA(allcells_posttrans, verbose = FALSE)

print(allcells_posttrans[["pca"]], dims = 1:5, nfeatures = 5)

allcells_posttrans <- FindNeighbors(allcells_posttrans, dims = 1:15)
allcells_posttrans <- FindClusters(allcells_posttrans, resolution = 0.3)


#Do UMAP showing 0.3 res clusters

allcells_posttrans_UMAP <- RunUMAP(allcells_posttrans,
                    dims = 1:15,
                    verbose = FALSE)
                    
DimPlot(allcells_posttrans_UMAP, reduction = "umap")


## ------------------------------------------------------------------------------------------------------------
# bypass all the initial processing and load the final data
allcells_posttrans <- readRDS("data_seurat_posttrans_80525.rds")


## ------------------------------------------------------------------------------------------------------------
# bypass even more and only load the counts of each gene, all the barcodes, and the timepoint column
all4000genes <- read.csv("all4000genes_12.8.25.csv")
# just want a list of the genes, too
genes_vector1 <- colnames(all4000genes)
genes_vector <- genes_vector1[-c(1,2)]


## ------------------------------------------------------------------------------------------------------------
#read in all of the probs and barcodes from tarun
probsandbarcodes <- read.csv("barcodes_umap_probabilities_info.csv")
# rename cells to barcodes to match the 4000genes file
names(probsandbarcodes)[names(probsandbarcodes) == "cells"] <- "barcodes"
# makes the barcodes column the same format as in the counts data
probsandbarcodes$barcodes <- sub("\\.(\\d+)$", "-\\1", probsandbarcodes$barcodes)



## ------------------------------------------------------------------------------------------------------------
# take just the barcodes and terminal one columns
probsandbarcodesterm1 <- probsandbarcodes[,c(1,9)]
# join the counts with the probability data
probsandcounts_all <- full_join(all4000genes, probsandbarcodesterm1, by = "barcodes")
probsandcounts_all <- subset(probsandcounts_all, !is.na(probability_terminal1))
# get rid of the random X column
probsandcounts_all <- probsandcounts_all %>%
  select(-X)


## ------------------------------------------------------------------------------------------------------------
# what's a good cutoff for "high" positive IFNL1 expressors?
probsandcounts_all%>% filter(orig.ident == "16hr")%>% group_by(IFNL1) %>% summarise(n = n()) %>%
    ggplot(aes(x = IFNL1, y = n)) +
    geom_line() +
    geom_vline(xintercept = 5)+
    geom_point(size = 2) +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    boxplot_theme +
    ggtitle("Shit happens", "for reals")



## ------------------------------------------------------------------------------------------------------------
# using the "high" cutoff, plotting the fraction of cells that are "high"
fracthigh <- fractbytime %>%
    filter(orig.ident == "16hr") %>%
    ggplot(aes(x = factor(orig.ident, levels = c("0hr", "4hr", "8hr", "12hr", "16hr")),
               y = fraction,
               fill = factor(orig.ident, levels = c("0hr", "4hr", "8hr", "12hr", "16hr"))))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#909082", "#D4C156", "#48D89A", "#52B5CF", "#C366EB"))+
        scale_y_continuous(
            limits = c(0.0, 0.5),
            breaks = seq(0.00, 0.5, 0.1),
            minor_breaks = seq(0.00, 1.00, 0.125))+
        boxplot_theme+
        theme(legend.position = "none")

           
fracthigh
#ggsave(filename="fractover0.3_bytime_8.26.25.pdf", plot = fracthigh, device="pdf", units="in", height=2, width=1.75)


## ------------------------------------------------------------------------------------------------------------
# looking at the probability distribution for each timepoint
boxplot_probdist <- ggplot(probsandcounts_all, aes(y = X16hr_1, x = factor(orig.ident, levels = c("0hr", "4hr", "8hr", "12hr", "16hr")), 
                               fill = factor(orig.ident, levels = c("0hr", "4hr", "8hr", "12hr", "16hr"))))+
    geom_hline(yintercept = 0.306, linetype = "dashed", linewidth = 0.1)+
    geom_boxplot(outlier.size = 0.005, size = 0.25, outlier.alpha = 0.1)+
    scale_fill_manual(values = c("#909082", "#D4C156", "#48D89A", "#52B5CF", "#C366EB"))+
    #scale_y_log10()+
    labs(x = "Time (Hours)", 
        y = "Probability of Transitioning to Terminal State 1")+
    boxplot_theme+
    theme(legend.position = "none")

boxplot_probdist

#ggsave(filename="boxplot_probdist_bytime_8.26.25.pdf", plot = boxplot_probdist, device="pdf", units="in", height=2, width=1.75)


## ------------------------------------------------------------------------------------------------------------
# same thing but with a histogram
hist_probdist <- ggplot(probsandcounts_all, aes(x = X16hr_1,
                               fill = factor(orig.ident, levels = c("16hr","12hr","8hr","4hr","0hr"))))+
    geom_histogram(aes(y = after_stat(ifelse(count > 0, count, NA))), stat = "bin", position = "identity", binwidth = 0.004)+
    scale_fill_manual(values = c("#C366EB","#52B5CF","#48D89A","#D4C156","#909082"))+
    geom_vline(xintercept = 0.31, linewidth = 0.4)+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(trans = "log10", expand = c(0.005, 0.005))+
    coord_cartesian(ylim = c(NA, 1e4))+
    boxplot_theme+
    theme(legend.position = "none")

hist_probdist

#ggsave(filename="hist_probdist_bytime_10.28.25.pdf", plot = hist_probdist, device="pdf", units="in", height=2, width=2)


## ------------------------------------------------------------------------------------------------------------
# subset to have just 0hr
probsandcounts_0hr <- subset(probsandcounts_all, orig.ident == "0hr")

box0hr <- ggplot(probsandcounts_0hr, aes(y = X16hr_1, x = orig.ident, fill = orig.ident))+
    geom_boxplot(outlier.size = 0.005, size = 0.25)+
    scale_fill_manual(values = "#909082")+
    geom_hline(yintercept = 0.31, linewidth = 0.2)+
    boxplot_theme+
    theme(legend.position = "none")

box0hr

#ggsave(filename="probdist_0hr_box_10.21.25.pdf", plot = box0hr, device="pdf", units="in", height=2, width=1.5)


## ------------------------------------------------------------------------------------------------------------
# formatting and stuff to make the genes columns just one column
probsandcounts_all$probability_terminal1 <- as.numeric(probsandcounts_all$probability_terminal1) 

probsandcounts_all$high_prob <- ifelse(!is.na(probsandcounts_all$probability_terminal1) & probsandcounts_all$probability_terminal1 > 0.31, "High", "Low")

# loooong file
probs_long <- probsandcounts_all %>%
  pivot_longer(
    cols = where(is.numeric) & !all_of("probability_terminal1"),
    names_to = "gene",
    values_to = "counts"
  )

# subset for just 0hr
probslong_0hr <- subset(probs_long, orig.ident == "0hr")


## ------------------------------------------------------------------------------------------------------------
# fix anything in the dataframe
fix_atom <- function(x) {
  if (is.list(x)) purrr::map_chr(x, ~ as.character(.x)[1]) else as.character(x)
}

df_clean <- probslong_0hr %>%
  ungroup() %>%                         # drop any prior grouping
  mutate(
    gene       = fix_atom(gene),
    # include this only if you actually have orig.ident; remove otherwise
    orig.ident = if ("orig.ident" %in% names(.)) fix_atom(orig.ident) else NA_character_,
    prob       = as.numeric(probability_terminal1),
    counts     = as.numeric(counts)
  )

df_clean$high_prob <- ifelse(df_clean$high_prob == "High", 1, 0)



## ------------------------------------------------------------------------------------------------------------
library(data.table)

# assume df has columns: gene, prob, counts (plus others)
DTw <- as.data.table(df_clean)
cutoff <- 0.31

# 1) define groups once (0 = low_prob, 1 = high_prob)
DTw[, grp := as.integer(prob > cutoff)]

# 2) tie correction term per gene: sum(t^3 - t) over ties of counts
#    (needed for variance of U under ties)
tie_tab <- DTw[, .N, by = .(gene, orig.ident, counts)]
tie_corr <- tie_tab[, .(tie_sum = sum(N^3 - N)), by = .(gene, orig.ident)]

# 3) per-gene U statistic and effect size using fast ranks
stat <- DTw[, {
  N  <- .N
  n1 <- sum(grp == 1L)
  n0 <- N - n1
  if (n1 == 0L || n0 == 0L) {
    list(n1 = n1, n0 = n0, U = NA_real_, effect = NA_real_)
  } else {
    r  <- frank(counts, ties.method = "average")     # fast ranks within gene
    R1 <- sum(r[grp == 1L])
    U  <- R1 - n1 * (n1 + 1) / 2
    eff <- U / (n1 * n0)                             # Prob(high > low)
    glass_rb <- (2 * U) / (n1 * n0) - 1
    list(n1 = n1, n0 = n0, U = U, effect = eff, glass_rb = glass_rb)
  }
}, by = .(gene, orig.ident)]

# 4) merge tie correction, compute z and p-value with tie-adjusted sigma
stat[tie_corr, on = .(gene, orig.ident), tie_sum := i.tie_sum]
stat[is.na(tie_sum), tie_sum := 0]

stat[, N := n1 + n0]
stat[, sigma := {
  # Var(U) = n1*n0/12 * [(N+1) - (sum(t^3 - t))/(N*(N-1))]  (tie correction)
  tc <- ifelse(N > 1, tie_sum / (N * (N - 1)), 0)
  sqrt((n1 * n0) / 12 * ((N + 1) - tc))
}]

stat[, z := (U - (n1 * n0) / 2) / sigma]
stat[!is.finite(z) | sigma == 0, z := NA_real_]

stat[, p_value := ifelse(is.na(z), NA_real_, 2 * pnorm(abs(z), lower.tail = FALSE))]
stat[, p_adj := {
  pv <- p_value
  pv[!is.na(pv)] <- p.adjust(pv[!is.na(pv)], method = "fdr")
  pv
}, by = orig.ident]

res_mw_dt <- stat[, .(gene, orig.ident, n0, n1, U, effect, p_value, p_adj, glass_rb)]
res_mw_dt[]


## ------------------------------------------------------------------------------------------------------------
# filtering for just the top genes
res_mw_0hr_filt_glass <- res_mw_dt %>%
    subset(!is.na(effect))%>%
    filter(p_adj <= 0.05,
          glass_rb > 0)%>%
    arrange(-glass_rb)
# all of the genes
res_mw_0hr_nofilt <- res_mw_dt %>%
    subset(!is.na(glass_rb))%>%
    arrange(-glass_rb)
res_mw_0hr_filt_glass


## ------------------------------------------------------------------------------------------------------------
scatter_theme <- theme(legend.title = element_text(color = "black", size = 8, 
                                                  face = "bold"),
                       panel.background = element_rect(color ="black", linewidth = 0.5, fill = "white"),
                       legend.text = element_text(color = "black", size = 8),
                       legend.key = element_blank(),
                       axis.title.x = element_text(color = "black", size = 8, 
                                               face = "bold",
                                               margin = margin(t = 5, r = 0, b = 0, l = 0)),
                       axis.title.y = element_text(color = "black", size = 8, 
                                               face = "bold"),
                       axis.text.x = element_text(color = "black", size = 8, 
                                              margin = margin(t = 2.5, r = 5, b = 0, l = 5)),
                       axis.text.y = element_text(color = "black", size = 8),
                       axis.ticks = element_line(linewidth = 0.25, color = "black"),
                       panel.grid.minor = element_blank(), #gets rid of grey and lines in the middle
                       panel.grid.major = element_blank())


## ------------------------------------------------------------------------------------------------------------
# looking at the distribution of correlations
top_five <- c("IFIT1", "IFIT2", "ISG15", "MX1", "OASL")
mw_dist_plot <- res_mw_0hr_nofilt %>%
    ggplot(aes(x = effect, y = -log10(p_adj)))+
        geom_point(alpha = 0.2, size = 0.1)+
        geom_hline(yintercept = 1.3, linewidth = 0.1, linetype = "dashed")+
        geom_vline(xintercept = 0.5,  linewidth = 0.1, linetype = "dashed")+
        scale_y_continuous(
            breaks = seq(0.00, 300, 50))+
        scale_x_continuous(
            breaks = seq(0, 1, 0.1))+
       geom_point(data = res_mw_0hr_filt_glass,
                   color = "red", size = 0.1, alpha = 0.2)+
        geom_label_repel(data = subset(res_mw_0hr_nofilt, gene %in% top_five), aes(label = gene),
                size = 2, box.padding = 0.4) +
        scatter_theme
mw_dist_plot


## ------------------------------------------------------------------------------------------------------------
sessionInfo()
# R version 4.5.2 (2025-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Ventura 13.0
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#         [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Chicago
# tzcode source: internal
# 
# attached base packages:
#         [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#         [1] ggrepel_0.9.6               data.table_1.17.8           magrittr_2.0.4             
# [4] BiocSingular_1.26.0         scran_1.38.0                scater_1.38.0              
# [7] scuttle_1.20.0              simpleSingleCell_1.34.0     DropletUtils_1.30.0        
# [10] SingleCellExperiment_1.32.0 SummarizedExperiment_1.40.0 Biobase_2.70.0             
# [13] GenomicRanges_1.62.1        Seqinfo_1.0.0               IRanges_2.44.0             
# [16] S4Vectors_0.48.0            BiocGenerics_0.56.0         generics_0.1.4             
# [19] MatrixGenerics_1.22.0       matrixStats_1.5.0           glmGamPoi_1.22.0           
# [22] sctransform_0.4.2           lubridate_1.9.4             forcats_1.0.1              
# [25] stringr_1.6.0               dplyr_1.1.4                 purrr_1.2.0                
# [28] readr_2.1.6                 tidyr_1.3.1                 tibble_3.3.0               
# [31] tidyverse_2.0.0             dittoSeq_1.22.0             ggplot2_4.0.1              
# [34] Seurat_5.3.1                SeuratObject_5.2.0          sp_2.2-0                   
# [37] BiocManager_1.30.27        
# 
# loaded via a namespace (and not attached):
#         [1] RcppAnnoy_0.0.22          splines_4.5.2             later_1.4.4              
# [4] CodeDepends_0.6.6         R.oo_1.27.1               polyclip_1.10-7          
# [7] graph_1.88.1              XML_3.99-0.20             fastDummies_1.7.5        
# [10] lifecycle_1.0.4           edgeR_4.8.1               processx_3.8.6           
# [13] globals_0.18.0            lattice_0.22-7            MASS_7.3-65              
# [16] rmarkdown_2.30            limma_3.66.0              plotly_4.11.0            
# [19] yaml_2.3.12               metapod_1.18.0            httpuv_1.6.16            
# [22] otel_0.2.0                spam_2.11-1               spatstat.sparse_3.1-0    
# [25] reticulate_1.44.1         cowplot_1.2.0             pbapply_1.7-4            
# [28] RColorBrewer_1.1-3        abind_1.4-8               Rtsne_0.17               
# [31] R.utils_2.13.0            irlba_2.3.5.1             listenv_0.10.0           
# [34] spatstat.utils_3.2-0      pheatmap_1.0.13           goftest_1.2-3            
# [37] RSpectra_0.16-2           spatstat.random_3.4-3     dqrng_0.4.1              
# [40] fitdistrplus_1.2-4        parallelly_1.46.0         DelayedMatrixStats_1.32.0
# [43] codetools_0.2-20          DelayedArray_0.36.0       tidyselect_1.2.1         
# [46] farver_2.1.2              viridis_0.6.5             ScaledMatrix_1.18.0      
# [49] spatstat.explore_3.6-0    jsonlite_2.0.0            BiocNeighbors_2.4.0      
# [52] progressr_0.18.0          ggridges_0.5.7            survival_3.8-3           
# [55] tools_4.5.2               ica_1.0-3                 Rcpp_1.1.0               
# [58] glue_1.8.0                gridExtra_2.3             SparseArray_1.10.6       
# [61] xfun_0.54                 HDF5Array_1.38.0          withr_3.0.2              
# [64] fastmap_1.2.0             bluster_1.20.0            rhdf5filters_1.22.0      
# [67] rsvd_1.0.5                callr_3.7.6               digest_0.6.39            
# [70] timechange_0.3.0          R6_2.6.1                  mime_0.13                
# [73] colorspace_2.1-2          scattermore_1.2           tensor_1.5.1             
# [76] spatstat.data_3.1-9       R.methodsS3_1.8.2         h5mread_1.2.1            
# [79] httr_1.4.7                htmlwidgets_1.6.4         S4Arrays_1.10.1          
# [82] uwot_0.2.4                pkgconfig_2.0.3           gtable_0.3.6             
# [85] lmtest_0.9-40             S7_0.2.1                  XVector_0.50.0           
# [88] htmltools_0.5.9           dotCall64_1.2             scales_1.4.0             
# [91] png_0.1-8                 spatstat.univar_3.1-5     knitr_1.50               
# [94] rstudioapi_0.17.1         tzdb_0.5.0                reshape2_1.4.5           
# [97] nlme_3.1-168              zoo_1.8-14                rhdf5_2.54.1             
# [100] KernSmooth_2.23-26        vipor_0.4.7               parallel_4.5.2           
# [103] miniUI_0.1.2              pillar_1.11.1             grid_4.5.2               
# [106] vctrs_0.6.5               RANN_2.6.2                promises_1.5.0           
# [109] beachmat_2.26.0           xtable_1.8-4              cluster_2.1.8.1          
# [112] beeswarm_0.4.0            evaluate_1.0.5            locfit_1.5-9.12          
# [115] cli_3.6.5                 compiler_4.5.2            rlang_1.1.6              
# [118] future.apply_1.20.1       ps_1.9.1                  ggbeeswarm_0.7.3         
# [121] plyr_1.8.9                stringi_1.8.7             viridisLite_0.4.2        
# [124] deldir_2.0-4              BiocParallel_1.44.0       lazyeval_0.2.2           
# [127] spatstat.geom_3.6-1       Matrix_1.7-4              RcppHNSW_0.6.0           
# [130] hms_1.1.4                 patchwork_1.3.2           sparseMatrixStats_1.22.0 
# [133] future_1.68.0             Rhdf5lib_1.32.0           statmod_1.5.1            
# [136] shiny_1.12.1              ROCR_1.0-11               igraph_2.2.1             

