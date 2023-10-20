#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
PATH <- "PATH/TO/DIRECTORY"
setwd(paste0(PATH, "/CaseStudies/AD"))

# Load packages
library(tidyverse)
library(biomaRt)

#==============================================================================#
# Prepare data for analysis
#==============================================================================#

# load data
DE_data <- read.delim("table_AD_PD_vs_contr_frontal_rma.tab")

# Get pathway genes
pathwayGenes <- read.delim(paste0(PATH, "/ComparisonWithDatabasePPI/Pathways_nodes_v5.txt"))

# Annotate probes with HGNC symbol and ENSEMBL ID
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
annotations <- getBM(attributes=c("ensembl_gene_id",
                                  "hgnc_symbol",
                                  "affy_hugene_1_0_st_v1"), 
                     filters = "affy_hugene_1_0_st_v1",
                     values = DE_data$X,
                     mart = ensembl)

# Save annotations
save(annotations, file = "ProbeAnnotations.RData")


#******************************************************************************#
# AD data set
#******************************************************************************#

# Load probe annotations
load("ProbeAnnotations.RData")

# Get annotations for each probe
AD_data <- inner_join(DE_data[,1:5], annotations, by = c("X" = "affy_hugene_1_0_st_v1"))

# Remove probes without associated ensembl ID
AD_data <- AD_data[!is.na(AD_data$ensembl_gene_id),]

# Select the most significant probe for each gene
AD_data <- arrange(AD_data, by = P.Value_AD)
AD_data <- AD_data[!duplicated(AD_data$ensembl_gene_id),]

# Get significant and non-significant genes
sigGenes_AD <- unique(AD_data$ensembl_gene_id[(AD_data$P.Value_AD < 0.05)])
nonsigGenes_AD <- unique(setdiff(AD_data$ensembl_gene_id,sigGenes_AD))

# Get pathways and non-pathways genes
pathGenes <- pathwayGenes[pathwayGenes$Identifier %in% annotations$ensembl_gene_id,]
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
sum(AD_vis$P.Value_AD < 0.05)
write.table(AD_vis, file = "AD_vis.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#******************************************************************************#
# PD data set
#******************************************************************************#

# Get annotations for each probe
PD_data <- inner_join(DE_data[,c(1,6:9)], annotations, by = c("X" = "affy_hugene_1_0_st_v1"))

# Remove probes without associated ensembl ID
PD_data <- PD_data[!is.na(PD_data$ensembl_gene_id),]

# Select the most significant probe for each gene
PD_data <- arrange(PD_data, by = P.Value_PD)
PD_data <- PD_data[!duplicated(PD_data$ensembl_gene_id),]

# Get significant and non-significant genes
sigGenes_PD <- unique(PD_data$ensembl_gene_id[(PD_data$P.Value_PD < 0.05)])
nonsigGenes_PD <- unique(setdiff(PD_data$ensembl_gene_id,sigGenes_PD))

# Get pathways and non-pathways genes
pathGenes <- pathwayGenes[pathwayGenes$Identifier %in% annotations$ensembl_gene_id,]
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
sum(PD_vis$P.Value_PD < 0.05)
write.table(PD_vis, file = "PD_vis.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Combine AD and PD visualization
all_vis <- inner_join(AD_vis, PD_vis, by = c("ensembl_gene_id" = "ensembl_gene_id"))
write.table(all_vis, file = "all_vis.txt", sep = "\t", quote = FALSE, row.names = FALSE)

