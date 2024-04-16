#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
PATH <- "PATH/TO/DIRECTORY" # PATH <- "C:/Users/jarno/OneDrive/Documents/GitHub/HDAC6Network"
setwd(paste0(PATH, "/CaseStudies/AD"))

# Load needed libraries
library(oligo)
library(limma)
library(bioDist)
library(gplots)
library(gcrma)
library(biomaRt)
library(GEOquery)
library(R2HTML)
library(tidyverse)

# Include some functions from ArrayAnalysis.org
source(paste0(PATH,"/CaseStudies/functions_ArrayAnalysis_v2.R"))

#==============================================================================#
# read and set up data variables
#==============================================================================#

#read the dataset
dat <- read.celfiles(list.celfiles("GSE36980_RAW", full.names=TRUE))

#adjust sample names to remove 'CEL'
sampleNames(dat) <- sub(".CEL.gz","",sampleNames(dat),fixed=TRUE)

#create a backup for the QC phase
dat.b <- dat
#extract the probe expression data and set to 2log scale
range(exprs(dat))
# 20 54766
dat.b <- log(exprs(dat.b),2)

#create metadata information
#based on study description as given in the record
geoFile <- getGEO("GSE36980")
desc <- geoFile$GSE36980_series_matrix.txt.gz@phenoData@data
desc <- desc[,c("geo_accession","source_name_ch1", "Sex:ch1", "age:ch1", "tissue:ch1")]
colnames(desc) <- c("SampleID", "Disease", "Sex", "Age", "Tissue")
desc$Disease <- factor(ifelse(str_detect(desc$Disease, "non-AD"), "CTR", "AD"),
                       levels = c("CTR", "AD"))
desc$Sex <- factor(ifelse(desc$Sex == "male", "M", "F"),
                   levels = c("M", "F"))
desc$Age <- as.numeric(desc$Age)

# Only include hippocampus
desc <- desc[desc$Tissue == "Hippocampus",]

# Remove outlier
desc <- desc[desc$SampleID != "GSM4764672",]

#reorder rows of desc to match columns of dat and dat.b
dat.b <- dat.b[,desc$SampleID]
dat <- dat[,desc$SampleID]
all(rownames(desc) == colnames(dat.b))
all(rownames(desc) == colnames(dat))
# -> TRUE, samples are now indeed in same order in dat.b and desc

#==============================================================================#
# pre-processing and QC
#==============================================================================#

#create QC plots for raw data
factors <- c("Sex","Disease","Age")
if(!dir.exists("QC_raw")) dir.create("QC_raw")
setwd("QC_raw")
createQCPlots(dat.b, factors, Table=desc, normMeth="", postfix="")
setwd("..")

#normalisation#
sum(sampleNames(dat)==rownames(desc))==length(sampleNames(dat))
# -> TRUE, samples are in same order in dat and desc
norm.rma <- oligo::rma(dat)
norm.rma <- exprs(norm.rma)


#create QC plots for norm data
if(!dir.exists("QC_normRMA")) dir.create("QC_normRMA")
setwd("QC_normRMA")
createQCPlots(norm.rma, factors, Table=desc, normMeth="", postfix="RMA")
setwd("..")

#save the normalised data tables
write.table(norm.rma,file="AD_norm_RMA.txt",sep="\t",quote=FALSE,col.names=NA)


#==============================================================================#
# Prepare for statistics
#==============================================================================#

#use norm.rma, for convenience rename this to norm
norm <- norm.rma

#are there NAs in dat?
sum(is.na(norm)) #0 -> no

#what is the range of the data?
range(norm)
#  1.058543 14.697479 -> these indeed are 2log scaled data; there are no negative values

##note: don't do any filtering here, as we want the data set mainly for pathway mapping
##so also keep non/lowly expressed genes in, as being non/lowly expressed also conveys information


#==============================================================================#
# Statistical modelling
#==============================================================================#

#model design: take group_ as a fixed effect for an intercept model
design <- model.matrix(~Disease,data=desc)
colnames(design) <- gsub("[()]","",colnames(design))

#fit model
fit <- lmFit(norm,design)
fit <- eBayes(fit)

#build the (contrast) matrix to compute some group differences of interest based on the model parameters
cont.matrix <- makeContrasts(
  AD_vs_contr = DiseaseAD,
  levels = colnames(design)
)
#compute the contrast fits
contrast.fit <- contrasts.fit(fit, cont.matrix)
contrast.fit <- eBayes(contrast.fit)

#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
files <- saveStatOutput(cont.matrix,contrast.fit,postfix="rma",annotation=NULL)

#create summary table of the contrast results
createPvalTab(files,postfix="AD_rma",html=TRUE)


#==============================================================================#
# Add gene annotations
#==============================================================================#

# Read statistics table
top.table <- read.delim("table_AD_vs_contr_rma.tab")

# Add gene annotations using biomaRt
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
annotations <- getBM(attributes=c("ensembl_gene_id",
                                  "hgnc_symbol",
                                  "affy_hugene_1_0_st_v1"), 
                     filters = "affy_hugene_1_0_st_v1",
                     values = as.character(top.table$X),
                     mart = ensembl)

# Get annotations for each probe
AD_data <- inner_join(top.table, annotations, by = c("X" = "affy_hugene_1_0_st_v1"))

# Remove probes without associated ensembl ID
AD_data <- AD_data[!is.na(AD_data$ensembl_gene_id),]

# Select the most significant probe for each gene
AD_data <- arrange(AD_data, by = P.Value)
AD_data <- AD_data[!duplicated(AD_data$ensembl_gene_id),]

# Save data
save(AD_data, file = "AD_data.RData")
