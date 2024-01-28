#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
# e.g. PATH <- "C:/Users/jarno/OneDrive/Documents/GitHub/HDAC6Network"
PATH <- "PATH/TO/DIRECTORY"
setwd(paste0(PATH, "/CaseStudies"))

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


network_types <- c("Full network",
                   "AggresomeSubnetwork",
                   "ComplexSubnetwork",
                   "DeacetylationSubnetwork",
                   "PhosphorylationSubnetwork",
                   "RegulatorySubnetwork")

network_names <- c("Full network",
                   "Aggresome-autophagy ",
                   "Protein complexes",
                   "(De)acetylation",
                   "Phosphorylation",
                   "Regulatory")

output_AD <- rep(NA, length(network_types))
names(output_AD) <- network_names

for (i in 1:length(network_types)){
    
    # Read network
    network <- read.csv(paste0(network_types[i],".csv"))
    genes_network <- unique(network$Ensembl)
    genes_network <- genes_network[genes_network != ""]
    
    # Network genes
    genes_network <- intersect(AD_data$ensembl_gene_id, genes_network)
    genes_nonetwork <- setdiff(AD_data$ensembl_gene_id, genes_network)
    
    # Significant genes
    sigGenes <- AD_data$ensembl_gene_id[AD_data$P.Value_AD < 0.05]
    nosigGenes <- setdiff(AD_data$ensembl_gene_id, sigGenes)
    
    # Prepare data for ORA
    i11 <- length(intersect(genes_network, sigGenes))
    i12 <- length(intersect(genes_network,nosigGenes)) 
    i21 <- length(intersect(genes_nonetwork,sigGenes))
    i22 <- length(intersect(genes_nonetwork,nosigGenes))
    
    
    dat <- data.frame(
      "DE" = c(i11, i21),
      "nonDE" = c(i12, i22),
      row.names = c("Path", "nonPath"),
      stringsAsFactors = FALSE
    )
    
    output_AD[i] <- fisher.test(dat)$p.value
  
}



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


network_types <- c("Full network",
                   "AggresomeSubnetwork",
                   "ComplexSubnetwork",
                   "DeacetylationSubnetwork",
                   "PhosphorylationSubnetwork",
                   "RegulatorySubnetwork")

network_names <- c("Full network",
                   "Aggresome-autophagy ",
                   "Protein complexes",
                   "(De)acetylation",
                   "Phosphorylation",
                   "Regulatory")

output_PD <- rep(NA, length(network_types))
names(output_PD) <- network_names

for (i in 1:length(network_types)){
  
  # Read network
  network <- read.csv(paste0(network_types[i],".csv"))
  genes_network <- unique(network$Ensembl)
  genes_network <- genes_network[genes_network != ""]
  
  # Network genes
  genes_network <- intersect(PD_data$ensembl_gene_id, genes_network)
  genes_nonetwork <- setdiff(PD_data$ensembl_gene_id, genes_network)
  
  # Significant genes
  sigGenes <- PD_data$ensembl_gene_id[PD_data$P.Value_PD < 0.05]
  nosigGenes <- setdiff(PD_data$ensembl_gene_id, sigGenes)
  
  # Prepare data for ORA
  i11 <- length(intersect(genes_network, sigGenes))
  i12 <- length(intersect(genes_network,nosigGenes)) 
  i21 <- length(intersect(genes_nonetwork,sigGenes))
  i22 <- length(intersect(genes_nonetwork,nosigGenes))
  
  
  dat <- data.frame(
    "DE" = c(i11, i21),
    "nonDE" = c(i12, i22),
    row.names = c("Path", "nonPath"),
    stringsAsFactors = FALSE
  )
  
  output_PD[i] <- fisher.test(dat)$p.value
  
}


#******************************************************************************#
# Make plot
#******************************************************************************#

# Prepare data for plotting
output_AD_fil <- output_AD[names(output_AD) != "Full network"]
output_PD_fil <- output_PD[names(output_PD) != "Full network"]

plotDF <- data.frame(Pvalue = c(output_AD_fil, output_PD_fil),
                     Subnetwork = c(names(output_AD_fil), names(output_PD_fil)),
                     Group = c(rep("Alzheimer's disease", length(output_AD_fil)),
                               rep("Parkinson's disease", length(output_PD_fil))),
                     FullP = c(rep(output_AD["Full network"], length(output_AD_fil)),
                               rep(output_PD["Full network"], length(output_PD_fil))))

# Set colors
colors <- setNames(c("black","#3A61FF","#8EB342","#FF6767","#FF9966","#FF1773"),
                   c("Full network",
                     "Aggresome-autophagy ",
                     "Protein complexes",
                     "(De)acetylation",
                     "Phosphorylation",
                     "Regulatory"))

# Make plot
p <- ggplot(plotDF) +
  geom_hline(yintercept = -log10(0.05), 
             linewidth = 1, linetype = "dashed")+
  geom_segment(aes(x = Subnetwork, xend = Subnetwork, 
                   y = -log10(Pvalue), yend = -log10(FullP),
                   color = Subnetwork),
               linewidth = 1.5) +
  geom_point(aes(x = Subnetwork, 
                 y = -log10(Pvalue),
                 color = Subnetwork),
             size = 5) +
  geom_hline(data = unique(plotDF[,c("FullP", "Group")]),
             aes(yintercept = -log10(FullP)), 
             linewidth = 1) +
  facet_grid(rows = vars(Group)) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression(-log[10]~"P value")) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = colors)

# Save plot
ggsave(p, file = "SubnetworkAnalysis.png", width = 6, height = 4)