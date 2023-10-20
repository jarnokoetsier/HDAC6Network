
#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
PATH <- "PATH/TO/DIRECTORY"
setwd(paste0(PATH, "/ComparisonWithDatabasePPI"))

# Load packages
library(tidyverse)
library(biomaRt)

#==============================================================================#
# Filter protein expression data
#==============================================================================#

# Human protein atlas data, version 23.0 (Ensembl version 109)
normal_tissue <- read.delim("normal_tissue.tsv")

# Selected tissues
selected <- c("caudate", "cerebellum", "cerebral cortex", "dorsal raphe", "hippocampus", "hypothalamus", "substantia nigra")

# Filter data for selected tissues
filtered_data <- normal_tissue[normal_tissue$Tissue %in% selected, ]
# Filter rows with not detected
filtered_data <- filtered_data[filtered_data$Level != "Not detected", ]
# Filter rows with uncertain reliability
filtered_data <- filtered_data[filtered_data$Reliability != "Uncertain", ]

# make protein levels sortable low < medium < high
filtered_data$Level <- factor(filtered_data$Level, ordered=TRUE, levels = c("Low", "Medium", "High"))

# sort data by gene id and decreasing protein levels (highest expression level on top)
filtered_data <- filtered_data[order(filtered_data$Gene, filtered_data$Level, decreasing = TRUE), ] #sort by id and reverse of abs(value)

# remove all duplicates (always keeps the first --> highest entry)
data <- filtered_data[ !duplicated(filtered_data$Gene), ]

#compute number of measurements (genes) per tissue
table(filtered_data$Tissue)

# remove tissue, cell type, reliability (not important anymore)
data <- data[,-c(3,4,6)]

# saves resulting file
write.table(data, file="brain-protein-expression.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

#clear the unzipped file
file.remove("normal_tissue.tsv")
