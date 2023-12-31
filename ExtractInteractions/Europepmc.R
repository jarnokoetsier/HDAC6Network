#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load required libraries
library(europepmc)

# set working directory to annotations directory
# Example: setwd("C:/Users/jarno/OneDrive/Documents/GitHub/HDAC6Network/ExtractInteractions/ANNOTATIONS")
setwd("PATH/TO/ANNOTATIONS")

#==============================================================================#
# Extract HDAC6 interactions
#==============================================================================#

# 1) load the corpus file names
corpus <- read.delim("corpus.txt", as.is=TRUE, header=FALSE)
corpus <- corpus[,1]

# 2) transform the names into unique Pubmed and PMC identifiers
corpusIDs <- unique(unlist(lapply(strsplit(corpus,".",fixed=TRUE),function(x) x[1])))

# 3) add databases Xrefs to the entries, for use in the epmc functions:

#initially MED: to all, for Pubmed
corpusIDs <- paste0("MED:",corpusIDs)

#the replace the MED: for the PMC identifiers to PMC:
corpusIDs <- sub("MED:PMC","PMC:PMC",corpusIDs)

# 4) run epmc on the corpus
epmcann <- epmc_annotations_by_id(corpusIDs)

#how many papers do have hits in other parts than abstract or title?
rowSums(table(epmcann$ext_id,epmcann$section)[,-c(2,14)])

# 5) save the results to a tab delimited text file
write.table(epmcann, file="epmc_annotations_HDAC6_corpus.txt", sep="\t", col.names=NA)

