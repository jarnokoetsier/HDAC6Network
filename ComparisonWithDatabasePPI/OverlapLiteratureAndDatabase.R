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

# Read expression data from the ProteinAtlas
brainExpr <- read.delim("brain-protein-expression.txt")

#==============================================================================#
# Create database-driven network
#==============================================================================#

# Load networks (these were downloaded from the corresponding Cytoscapes apps)
biogridNetwork <- read.csv("BioGrid.csv")
intactNetwork <- read.csv("IntAct.csv")

# Get the union of genes (Entrez IDs)
unionGenes <- union(as.character(biogridNetwork$ncbi_gene_id),
                    as.character(intactNetwork$Entrez.Gene))

# Remove NA genes
unionGenes <- unionGenes[!is.na(unionGenes)]

# Use biomart to get HGNC symbol and ENSEMBL IDs for each gene
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
annotations <- getBM(attributes=c("entrezgene_id",
                                  "ensembl_gene_id", 
                                  "hgnc_symbol"), 
                     filters = 'entrezgene_id',
                     values = unionGenes,
                     mart = ensembl)



# Combine with brain protein expression data
databaseNetwork1 <- left_join(annotations, brainExpr, by = c("ensembl_gene_id" = "Gene"))

# Select Entrez IDs that were in the union of the two networks
databaseNetwork1 <- databaseNetwork1[databaseNetwork1$entrezgene_id %in% unionGenes,]

# Change expression level to an ordered factor
databaseNetwork1$Level <- factor(databaseNetwork1$Level, levels = c("High", "Medium", "Low"))

# Select highest expression
databaseNetwork1 <- arrange(databaseNetwork1, by = Level)
databaseNetwork1 <- databaseNetwork1[!duplicated(databaseNetwork1$entrezgene_id),]

# Save the network
save(databaseNetwork1, file = "databaseNetwork.RData")


#==============================================================================#
# Create literature-based network
#==============================================================================#

# Read data nodes from PathVisio
literatureNetwork <- read.delim("Pathways_nodes_v5.txt")

# Combine with brain protein expression data
literatureNetwork1 <- left_join(literatureNetwork, brainExpr, by = c("Identifier" = "Gene"))

#Change expression level to an ordered factor
literatureNetwork1$Level <- factor(literatureNetwork1$Level, levels = c("High", "Medium", "Low"))

# Select highest expression
literatureNetwork1 <- arrange(literatureNetwork1, by = Level)
literatureNetwork1 <- literatureNetwork1[!duplicated(literatureNetwork1$Identifier),]

# Save the network
save(literatureNetwork1, file = "literatureNetwork.RData")


#==============================================================================#
# Compare both networks
#==============================================================================#

# Intersect
i1 <- intersect(literatureNetwork1$Identifier,databaseNetwork1$ensembl_gene_id)
table(literatureNetwork1$Level[literatureNetwork1$Identifier %in% i1])
sum(is.na(literatureNetwork1$Level[literatureNetwork1$Identifier %in% i1]))

# Unique to literature-based network
i2 <- setdiff(literatureNetwork1$Identifier,databaseNetwork1$ensembl_gene_id)
table(literatureNetwork1$Level[literatureNetwork1$Identifier %in% i2])
sum(is.na(literatureNetwork1$Level[literatureNetwork1$Identifier %in% i2]))

# Unique to database-driven network
i3 <- setdiff(databaseNetwork1$ensembl_gene_id,literatureNetwork1$Identifier)
table(databaseNetwork1$Level[databaseNetwork1$ensembl_gene_id %in% i3])
sum(is.na(databaseNetwork1$Level[databaseNetwork1$ensembl_gene_id %in% i3]))



