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

# bypass even more and only load the counts of each gene, all the barcodes, and the timepoint column
all4000genes <- read.csv("all4000genes_12.8.25.csv")
# just want a list of the genes, too
genes_vector1 <- colnames(all4000genes)
genes_vector <- genes_vector1[-c(1,2)]

#read in all of the probs and barcodes from tarun
probsandbarcodes <- read.csv("barcodes_umap_probabilities_info.csv")

names(probsandbarcodes)[names(probsandbarcodes) == "cells"] <- "barcodes"

# makes the barocdes column the same format as in the counts data
probsandbarcodes$barcodes <- sub("\\.(\\d+)$", "-\\1", probsandbarcodes$barcodes)


# take just the barcodes and terminal columns
probsandbarcodesallterm <- probsandbarcodes[,c(1,9,10,11,12,13)]
# join the counts with the probability data
probsandcounts_all <- full_join(all4000genes, probsandbarcodesallterm, by = "barcodes")
probsandcounts_all <- subset(probsandcounts_all, !is.na(probability_terminal1))
probsandcounts_all <- probsandcounts_all %>%
  select(-X)

#extracting the cells that have a probability of transitioning to terminal at 100%
terminal1 <- probsandcounts_all %>%
filter(orig.ident == "16hr") %>%
filter(probability_terminal1 == 1)

terminal1$terminal <- "one"

terminal2 <- probsandcounts_all %>%
filter(orig.ident == "16hr") %>%
filter(probability_terminal2 == 1)

terminal2$terminal <- "two"

terminal3 <- probsandcounts_all %>%
filter(orig.ident == "16hr") %>%
filter(probability_terminal3 == 1)

terminal3$terminal <- "three"

terminal4 <- probsandcounts_all %>%
filter(orig.ident == "16hr") %>%
filter(probability_terminal4 == 1)

terminal4$terminal <- "four"

terminal5 <- probsandcounts_all %>%
filter(orig.ident == "16hr") %>%
filter(probability_terminal5 == 1)

terminal5$terminal <- "five"

# Combine all of the data into one dataframe
probscountsterminals <- rbind(terminal1, terminal2, terminal3, terminal4, terminal5) 
probscountsterminals$terminal <- factor(probscountsterminals$terminal , levels = c("one","two","three","four","five"))

# Check which terminal has the most cells that are IFNL+
percentterminal <- probscountsterminals %>%
group_by(terminal)%>%
summarise(fraction = mean(IFNL1 >= 5)*100)%>%
ggplot(aes(x = terminal, y = fraction, fill = terminal))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c("#3D4171","#5B4F90","#8460B0","#EE82EF","#B571D0"))+
    boxplot_theme+
    theme(legend.position = "none")
percentterminal
#ggsave(filename="percentterminalIFNL1_9.11.25.pdf", plot = percentterminal, device="pdf", units="in", height=1.5, width=1.5)

probsandcounts_all %>%
group_by(orig.ident)%>%
filter(probability_terminal1 > 0.31) %>%
summarise(fraction = mean(IFNL1 >1)*100)#%>%
#ggplot(aes(x = terminal, y = fraction))+
    #geom_bar(stat = "identity")

# getting the high probabability cells from terminal 1
probscounts_long <- probsandcounts_all %>%
group_by(orig.ident)%>%
filter(probability_terminal1 > 0.31)%>%
pivot_longer(cols = where(is.numeric) & !all_of(c("probability_terminal1", "probability_terminal2",
                                                  "probability_terminal3", "probability_terminal4", 
                                                  "probability_terminal5")),
    names_to = "gene",
    values_to = "counts")


probscounts_long

# getting the low probability cells from terminal1
probscounts_long_nothigh <- probsandcounts_all %>%
group_by(orig.ident)%>%
filter(probability_terminal1 <= 0.31)%>%
pivot_longer(cols = where(is.numeric) & !all_of(c("probability_terminal1", "probability_terminal2",
                                                  "probability_terminal3", "probability_terminal4", 
                                                  "probability_terminal5")),
    names_to = "gene",
    values_to = "counts")

probscounts_long_nothigh

# top genes from the Mann-Whitney analysis
top <- c("IFIT1", "IFIT2", "ISG15", "MX1", "OASL", "IFIT3", "DDX58", "IFIH1", "IFI6", "OAS2",
          "IRF7", "HELZ2", "CMPK2", "APOL2", "PARP9", "TNFSF10", "USP18", "IFI27", "PMAIP1", "HERC5")
top5 <- c("IFIT1", "IFIT2", "ISG15", "MX1", "OASL")

# genes known for the RNA sensing to IFNL1 pathway
interest <- c("DDX58","MAVS","TBK1","IRF3","IFNL1")

probscounts_long %>% 
filter(gene %in% top) %>% 
filter(orig.ident =="0hr")

probscounts_long_nothigh %>% 
filter(gene %in% top) %>% 
filter(orig.ident =="0hr")

# just looking at the expression of the top genes at 0hr
probscounts_long %>% 
filter(gene %in% top) %>% 
filter(orig.ident =="0hr") %>%
ggplot(aes(x = factor(gene, levels = top), y = counts))+
    geom_boxplot(outlier.size = 0.1, size = 0.25) + 
    scale_y_log10()+
    annotation_logticks(sides = "l")+
    boxplot_theme+
    theme(legend.position = "none")

# seeing what fraction of top highly terminal1 probable cells are producing the top genes at 0hr
probscounts_long %>% 
filter(gene %in% top) %>% 
filter(orig.ident =="0hr")%>%
group_by(gene)%>%
summarise(fraction = mean(counts >= 1)*100)

# seeing what fraction of lowly terminal1 probable cells are producing the top genes at 0hr
probscounts_long_nothigh %>% 
filter(gene %in% top) %>% 
filter(orig.ident =="0hr")%>%
group_by(gene)%>%
summarise(fraction = mean(counts >= 1)*100)

topgenespercent <- probscounts_long %>% 
filter(gene %in% top) %>% 
filter(orig.ident =="0hr")%>%
group_by(gene)%>%
summarise(fraction = mean(counts >= 1)*100)%>%
ggplot(aes(x = factor(gene, levels = top), y = fraction))+
    geom_bar(stat = "identity") +
    scale_y_continuous(limits = c(0,100), expand = c(0.025,0.025))+
    boxplot_theme
topgenespercent

#ggsave(filename="topmwpercent0hr_9.11.25.pdf", plot = topgenespercent, device="pdf", units="in", height=2, width=7)
lowhighprobs_topgenes <- full_join(probscounts_long_nothigh, probscounts_long, by = "gene")
# i plotted these values in Prism

# looking at the counts of the top interesting genes at 0hrs
probscounts_long %>% filter(gene %in% topinterest) %>% filter(orig.ident =="0hr")%>%
ggplot(aes(x = gene, y = counts))+
    geom_boxplot(outlier.size = 0.1, size = 0.25) + 
    scale_fill_manual(values = c("#909082", "#D4C156", "#48D89A", "#52B5CF", "#C366EB"))+
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), expand = c(0,0))+
    annotation_logticks(sides = "l")+
    boxplot_theme+
    theme(legend.position = "none")

# checking the expression distribution of OASL at 16hrs
probsandcounts_all%>% filter(orig.ident == "16hr")%>% group_by(OASL) %>% summarise(n = n()) %>%
    ggplot(aes(x = OASL, y = n)) +
    geom_line() +
    geom_vline(xintercept = 5)+
    geom_point(size = 2) +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    boxplot_theme +
    ggtitle("Shit happens", "for reals")

sessionInfo()
