#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
PATH <- "PATH/TO/DIRECTORY" # PATH <- "C:/Users/jarno/OneDrive/Documents/GitHub/HDAC6Network"
setwd(paste0(PATH, "/CaseStudies"))

# Load packages
library(tidyverse)
library(biomaRt)

# Get pathway genes
pathwayGenes <- read.delim(paste0(PATH, "/ComparisonWithDatabasePPI/Pathways_nodes_v5.txt"))


#******************************************************************************#
# AD data set
#******************************************************************************#

# Load probe annotations
load("AD/AD_data.RData")

# Get significant and non-significant genes
sigGenes_AD <- unique(AD_data$ensembl_gene_id[(AD_data$P.Value < 0.05)])
nonsigGenes_AD <- unique(setdiff(AD_data$ensembl_gene_id,sigGenes_AD))

# Get pathways and non-pathways genes
pathGenes <- unique(intersect(pathwayGenes$Identifier, AD_data$ensembl_gene_id))
nopathGenes <- unique(setdiff(AD_data$ensembl_gene_id, pathGenes))

# Prepare data for ORA
i11 <- length(intersect(pathGenes, sigGenes_AD))
i12 <- length(intersect(pathGenes,nonsigGenes_AD)) 
i21 <- length(intersect(nopathGenes,sigGenes_AD))
i22 <- length(intersect(nopathGenes,nonsigGenes_AD))


dat <- data.frame(
  "DE" = c(i11, i21),
  "nonDE" = c(i12, i22),
  row.names = c("Path", "nonPath"),
  stringsAsFactors = FALSE
)

# Perform Chi-square test: expected frequencies is higher than 5 so no Fisher's exact test
chisq.test(dat)$expected
chisq.test(dat)$observed
chisq.test(dat)

# Make table for visualization
AD_vis <- AD_data[AD_data$ensembl_gene_id %in% pathwayGenes$Identifier,]
sum(AD_vis$P.Value < 0.05)
write.table(AD_vis, file = "AD_vis.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#******************************************************************************#
# PD data set
#******************************************************************************#
# Load probe annotations
load("PD/PD_data.RData")

# Get significant and non-significant genes
sigGenes_PD <- unique(PD_data$ensembl_gene_id[(PD_data$P.Value < 0.05)])
nonsigGenes_PD <- unique(setdiff(PD_data$ensembl_gene_id,sigGenes_PD))

# Get pathways and non-pathways genes
pathGenes <- unique(intersect(pathwayGenes$Identifier, PD_data$ensembl_gene_id))
nopathGenes <- unique(setdiff(PD_data$ensembl_gene_id, pathGenes))

# Prepare data for ORA
i11 <- length(intersect(pathGenes, sigGenes_PD))
i12 <- length(intersect(pathGenes,nonsigGenes_PD)) 
i21 <- length(intersect(nopathGenes,sigGenes_PD))
i22 <- length(intersect(nopathGenes,nonsigGenes_PD))


dat <- data.frame(
  "DE" = c(i11, i21),
  "nonDE" = c(i12, i22),
  row.names = c("Path", "nonPath"),
  stringsAsFactors = FALSE
)

# Perform Chi-square test: expected frequencies is higher than 5 so no Fisher's exact test
chisq.test(dat)$expected
chisq.test(dat)$observed
chisq.test(dat)

# Make table for visualization
PD_vis <- PD_data[PD_data$ensembl_gene_id %in% pathwayGenes$Identifier,]
sum(PD_vis$P.Value < 0.05)
write.table(PD_vis, file = "PD_vis.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#******************************************************************************#
# ALS data set
#******************************************************************************#

# Load probe annotations
load("ALS/ALS_data.RData")

# Get significant and non-significant genes
sigGenes_ALS <- unique(ALS_data$ensembl_gene_id[(ALS_data$PValue < 0.05)])
nonsigGenes_ALS <- unique(setdiff(ALS_data$ensembl_gene_id,sigGenes_ALS))

# Get pathways and non-pathways genes
pathGenes <- unique(intersect(pathwayGenes$Identifier, ALS_data$ensembl_gene_id))
nopathGenes <- unique(setdiff(ALS_data$ensembl_gene_id, pathGenes))

# Prepare data for ORA
i11 <- length(intersect(pathGenes, sigGenes_ALS))
i12 <- length(intersect(pathGenes,nonsigGenes_ALS)) 
i21 <- length(intersect(nopathGenes,sigGenes_ALS))
i22 <- length(intersect(nopathGenes,nonsigGenes_ALS))


dat <- data.frame(
  "DE" = c(i11, i21),
  "nonDE" = c(i12, i22),
  row.names = c("Path", "nonPath"),
  stringsAsFactors = FALSE
)

# Perform Chi-square test: expected frequencies is higher than 5 so no Fisher's exact test
chisq.test(dat)$expected
chisq.test(dat)$observed
chisq.test(dat)

# Make table for visualization
ALS_vis <- ALS_data[ALS_data$ensembl_gene_id %in% pathwayGenes$Identifier,]
sum(ALS_vis$PValue < 0.05)
write.table(ALS_vis, file = "ALS_vis.txt", sep = "\t", quote = FALSE, row.names = FALSE)


