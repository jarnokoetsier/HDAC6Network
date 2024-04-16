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

# Get pathway genes
pathwayGenes <- read.delim(paste0(PATH, "/ComparisonWithDatabasePPI/Pathways_nodes_v5.txt"))


#******************************************************************************#
# PD data set
#******************************************************************************#

# Load data
load("PD/PD_data.RData")


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
n_subnetwork_PD <- rep(NA, length(network_types))
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
  sigGenes <- PD_data$ensembl_gene_id[PD_data$P.Value < 0.05]
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
  n_subnetwork_PD[i] <- i11
  
}

#******************************************************************************#
# AD data set
#******************************************************************************#

# Load data
load("AD/AD_data.RData")

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
n_subnetwork_AD <- rep(NA, length(network_types))

for (i in 1:length(network_types)){
  
  # Read network
  network <- read.csv(paste0(network_types[i],".csv"))
  genes_network <- unique(network$Ensembl)
  genes_network <- genes_network[genes_network != ""]
  
  # Network genes
  genes_network <- intersect(AD_data$ensembl_gene_id, genes_network)
  genes_nonetwork <- setdiff(AD_data$ensembl_gene_id, genes_network)
  
  # Significant genes
  sigGenes <- AD_data$ensembl_gene_id[AD_data$P.Value < 0.05]
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
  n_subnetwork_AD[i] <- i11
  
}

#******************************************************************************#
# ALS data set
#******************************************************************************#

# Load ALS data
load("ALS/ALS_data.RData")

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

output_ALS <- rep(NA, length(network_types))
n_subnetwork <- rep(NA, length(network_types))
n_subnetwork_ALS <- rep(NA, length(network_types))

names(output_ALS) <- network_names

for (i in 1:length(network_types)){
  
  # Read network
  network <- read.csv(paste0(network_types[i],".csv"))
  genes_network <- unique(network$Ensembl)
  genes_network <- genes_network[genes_network != ""]
  
  # Network genes
  genes_network <- intersect(ALS_data$ensembl_gene_id, genes_network)
  genes_nonetwork <- setdiff(ALS_data$ensembl_gene_id, genes_network)
  
  # Significant genes
  sigGenes <- ALS_data$ensembl_gene_id[ALS_data$PValue < 0.05]
  nosigGenes <- setdiff(ALS_data$ensembl_gene_id, sigGenes)
  
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
  
  output_ALS[i] <- fisher.test(dat)$p.value
  n_subnetwork[i] <- length(genes_network)
  n_subnetwork_ALS[i] <- i11
}



#******************************************************************************#
# Make plot
#******************************************************************************#

# Prepare data for plotting
output_AD_fil <- output_AD[names(output_AD) != "Full network"]
output_PD_fil <- output_PD[names(output_PD) != "Full network"]
output_ALS_fil <- output_ALS[names(output_ALS) != "Full network"]

plotDF <- data.frame(Pvalue = c(output_AD_fil, output_PD_fil, output_ALS_fil),
                     Subnetwork = c(names(output_AD_fil), names(output_PD_fil),
                                    names(output_ALS_fil)),
                     Group = c(rep("Alzheimer's disease", length(output_AD_fil)),
                               rep("Parkinson's disease", length(output_PD_fil)),
                               rep("Amyotrophic lateral sclerosis", length(output_ALS_fil))),
                     FullP = c(rep(output_AD["Full network"], length(output_AD_fil)),
                               rep(output_PD["Full network"], length(output_PD_fil)),
                               rep(output_ALS["Full network"], length(output_ALS_fil))),
                     Nsig = c(n_subnetwork_AD[-1], n_subnetwork_PD[-1],
                              n_subnetwork_ALS[-1]),
                     Ntotal = c(rep(n_subnetwork[-1],3)),
                     PercT = rep(c(n_subnetwork_AD[1], n_subnetwork_PD[1],
                                     n_subnetwork_ALS[1]), each = 5)/100
                     )

plotDF$Perc <- plotDF$Nsig/plotDF$Ntotal

# Set colors
colors <- setNames(c("black","#3A61FF","#8EB342","#FF6767","#FF9966","#FF1773"),
                   c("Full network",
                     "Aggresome-autophagy ",
                     "Protein complexes",
                     "(De)acetylation",
                     "Phosphorylation",
                     "Regulatory"))
plotDF$Sig <- ifelse(plotDF$Pvalue < 0.05, "Yes", "No")


p <- ggplot(plotDF) +
  geom_rect(data = unique(plotDF[,c("PercT", "Group")]),
            aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = PercT),
            fill  ="lightgrey", alpha = 0.5) +

  geom_hline(data = unique(plotDF[,c("PercT", "Group")]), aes(yintercept = PercT), 
             linewidth = 1, linetype = "dashed")+
  geom_segment(aes(x = Subnetwork, xend = Subnetwork, 
                   y = 0, yend = Perc,
                   color = -log10(Pvalue)), linewidth = 1) +
  geom_point(aes(x = Subnetwork, 
                 y = Perc,
                 color = -log10(Pvalue),
                 size = -log10(Pvalue))) +
  facet_grid(rows = vars(Group)) +
  guides(size = "none") +
  labs(color = expression(-log[10]~"P value ")) +
  coord_flip() +
  xlab(NULL) +
  ylab("Proportion of DEGs") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_viridis_c(option = "viridis")


# Save plot
ggsave(p, file = "SubnetworkAnalysis.png", width = 6, height = 7)
