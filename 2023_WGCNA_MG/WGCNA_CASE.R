# "Weighted Gene Correlation Network Analysis for Larval CASE RNAseq"

# LOAD PACKAGES
library(WGCNA) # installed this package with `BiocManager::install("WGCNA")`
library(dplyr)
library(zoo)
library(DESeq2)
# for heatmap
# library(devtools)
# install_github("jokergoo/ComplexHeatmap") first run these - commented out to avoid running a second time...
library(ComplexHeatmap)
library(circlize)
library(reshape)
library(ggplot2)
library(hrbrthemes)

# SET WORKING DIRECTORY 
setwd("/home/mguidry/repos/Larval-Oyster-CASE-RNA/2022_WGCNA_MG")
getwd()

# LOAD DATA - Count matrix file and treatment/metadata file
#this data is a count matrix with rownames=sample ID and colnames=gene ID
count_data <- read.csv(file ="", sep =',', header=TRUE)
#here's a snapshot of the matrix...
