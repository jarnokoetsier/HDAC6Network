#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
PATH <- "PATH/TO/DIRECTORY"
setwd(paste0(PATH, "/ComparisonWithDatabasePPI"))

# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer) 

#==============================================================================#
# Make plot
#==============================================================================#

# Load database-derived network
load("databaseNetwork.RData")

# Load literature-based network
load("literatureNetwork.RData")

# Get genes that are unique to the database-derived network
uniqueGenes <- setdiff(databaseNetwork1$ensembl_gene_id,literatureNetwork1$Identifier)
databaseNetwork_fil <- databaseNetwork1[databaseNetwork1$ensembl_gene_id %in% uniqueGenes,]

# Set levels for plotting
databaseNetwork_fil$Level1 <- as.character(databaseNetwork_fil$Level)
databaseNetwork_fil$Level1[is.na(databaseNetwork_fil$Level1)] <- "NA"
databaseNetwork_fil$LevelName <- "1"
databaseNetwork_fil$LevelName[databaseNetwork_fil$Level1 == "Low"] <- "2"
databaseNetwork_fil$LevelName[databaseNetwork_fil$Level1 == "Medium"] <- "3"
databaseNetwork_fil$LevelName[databaseNetwork_fil$Level1 == "High"] <- "4"

write.csv(databaseNetwork_fil, file = "uniqueDBGenes.csv")

# create a data frame for making the network

# 1. From HDAC6 to the different groups
d1 <-  data.frame(from = "HDAC6", to = c("HDAC6",
                                       "High",
                                       "group1",
                                       "Medium",
                                       "group2",
                                       "Low",
                                       "group3",
                                       "NA",
                                       "group4"))

# 2. From different groups to genes
d2 = data.frame(from = c(databaseNetwork_fil$Level1[1:(which(databaseNetwork_fil$Level1 == "Medium")[1]-1)],
                         "group1", 
                         databaseNetwork_fil$Level1[(which(databaseNetwork_fil$Level1 == "Medium")[1]):(which(databaseNetwork_fil$Level1 == "Low")[1]-1)],
                         "group2", 
                         databaseNetwork_fil$Level1[(which(databaseNetwork_fil$Level1 == "Low")[1]):(which(databaseNetwork_fil$Level1 == "NA")[1]-1)],
                         "group3", 
                         databaseNetwork_fil$Level1[(which(databaseNetwork_fil$Level1 == "NA")[1]):nrow(databaseNetwork_fil)],
                         "group4"), 
                to = c(databaseNetwork_fil$hgnc_symbol[1:(which(databaseNetwork_fil$Level1 == "Medium")[1]-1)],
                       "gene1", 
                       databaseNetwork_fil$hgnc_symbol[(which(databaseNetwork_fil$Level1 == "Medium")[1]):(which(databaseNetwork_fil$Level1 == "Low")[1]-1)],
                       "gene2", 
                       databaseNetwork_fil$hgnc_symbol[(which(databaseNetwork_fil$Level1 == "Low")[1]):(which(databaseNetwork_fil$Level1 == "NA")[1]-1)],
                       "gene3", 
                       databaseNetwork_fil$hgnc_symbol[(which(databaseNetwork_fil$Level1 == "NA")[1]):nrow(databaseNetwork_fil)],
                       "gene4")
                )

# Combine to create edges
edges <- rbind(d1, d2)

# Set as factor
edges$from <- factor(edges$from, 
                     levels = c("HDAC6",
                                "High",
                                "group1",
                                "Medium",
                                "group2",
                                "Low",
                                "group3",
                                "NA",
                                "group4"))

# create a vertices data.frame. One line per object of our hierarchy
vertices = data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to))) , 
  value = 0.2
) 
# Let's add a column with the group of each name. It will be useful later to color points
vertices$group = edges$from[ match( vertices$name, edges$to ) ]

#Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
#calculate the ANGLE of the labels
vertices$id=NA
myleaves=which(is.na( match(vertices$name, edges$from) ))
nleaves=length(myleaves)
vertices$id[ myleaves ] = seq(1:nleaves)
vertices$angle= 90 - 360 * vertices$id / nleaves

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

# Create a graph object
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# Set colors
colors_edges <- c("black","#99000D", "white","#EF3B2C","white","#FCBBA1", "white","grey","white")
colors_nodes <- c(rep("white",4),"#99000D","#FCBBA1", "#EF3B2C", "grey")

# Set alpha
alpha_text <- rep(1, length(unique(edges$to)))
alpha_text[str_detect(unique(edges$to), "gene")] <- 0
names(alpha_text) <- unique(edges$to)


# Make the plot:
# Note that the plot is not finished yet. Some final adjustments to it's layout
# are still needed
p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(aes(colour=as.character(from))) +
  scale_edge_colour_manual(values = colors_edges) +
  #scale_edge_alpha_manual(values = alpha_edges) +
  geom_node_text(aes(x = x*1.05, y=y*1.05, filter = leaf, 
                     label=name, angle = angle, hjust=hjust,
                     alpha = name), size=2.7) +
  geom_node_point(aes(filter = leaf, x = x, y=y, color=group), size=0.8, alpha=1) +
  scale_size_continuous( range = c(0.1,10) ) +
  scale_color_manual(values = colors_nodes) +
  scale_alpha_manual(values = alpha_text) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))

# Save plot
ggsave(p, file = "DatabaseNetwork2.png", width = 8, height = 8)