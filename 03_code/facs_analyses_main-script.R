Flow cytometry 

library(BiocManager)
BiocManager::install(c("flowCore", "ggcyto"))

setwd("/Users/sd21/Documents/WORK/02_PROJECTS/single_cell-nuclei_transcriptomics/snRNAseq_optimisation/Cytometry/")

#libraries
library(flowCore)
library(ggcyto)
library(tidyverse)
library(patchwork)
library(mclust)


setwd("/Users/sd21/Documents/WORK/02_PROJECTS/single_cell-nuclei_transcriptomics/snRNAseq_optimisation/Paper_Cytometry")
source("nuclei_sybr_peak_analysis.R")

setwd("/Users/sd21/Documents/WORK/02_PROJECTS/single_cell-nuclei_transcriptomics/snRNAseq_optimisation/Paper_Cytometry")
source("facs_plots_all.R")