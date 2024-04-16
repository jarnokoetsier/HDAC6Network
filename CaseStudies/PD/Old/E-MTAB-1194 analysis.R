#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 
 
# Set working directory
PATH <- "PATH/TO/DIRECTORY"
setwd(paste0(PATH, "/CaseStudies/PD"))

# Load needed libraries
library(oligo)
library(limma)
library(bioDist)
library(gplots)
library(gcrma)
library(biomaRt)
library(R2HTML)

# Include some functions from ArrayAnalysis.org
source(paste0(PATH,"/CaseStudies/functions_ArrayAnalysis_v2.R"))

#==============================================================================#
# read and set up data variables
#==============================================================================#

#read the dataset
dat <- read.celfiles(list.celfiles("E-MTAB-1194", full.names=TRUE))
#adjust sample names to remove 'CEL'
sampleNames(dat) <- sub(".CEL","",sampleNames(dat),fixed=TRUE)

#create a backup for the QC phase
dat.b <- dat
#extract the probe expression data and set to 2log scale
range(exprs(dat))
# 14.0 15741.2
dat.b <- log(exprs(dat.b),2)

#create metadata information
#based on study description as given in the record
desc <- read.delim("sample_table.txt",as.is=TRUE,row.names=6)
desc$sex <- factor(desc$sex,levels=c("m","f"))
desc$tissue <- factor(desc$tissue,levels=c("frontal"))
desc$disease <- factor(desc$disease,levels=c("control","PD"))

#reorder rows of desc to match columns of dat and dat.b
desc <- desc[match(colnames(dat.b),rownames(desc)),]
#check order of samples in data and metadata
sum(colnames(dat.b)==rownames(desc)) == dim(dat.b)[2]
# -> TRUE, samples are now indeed in same order in dat.b and desc

#==============================================================================#
# pre-processing and QC
#==============================================================================#

#create QC plots for raw data
factors <- c("sex","disease","age")
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


#based on the QC plots, sample 21_G05_P2 is an outlier
# -> remove it
norm.rma.S <- norm.rma[,-match("21_G05_P2",colnames(norm.rma))]
desc.S <- desc[-match("21_G05_P2",rownames(desc)),]
#check order of samples in data and metadata
sum(colnames(norm.rma.S)==rownames(desc.S)) == dim(norm.rma.S)[2]
# -> TRUE, samples are in same order in norm.rma.S and desc.S


#create QC plots for norm.S (selected) data
if(!dir.exists("QC_normRMA_sel")) dir.create("QC_normRMA_sel")
setwd("QC_normRMA_sel")
createQCPlots(norm.rma.S, factors, Table=desc.S, normMeth="", postfix="RMA_sel")
setwd("..")


#save the normalised data tables
write.table(norm.rma.S,file="E-MTAB-1194_norm_RMA.txt",sep="\t",quote=FALSE,col.names=NA)


#==============================================================================#
# Prepare for statistics
#==============================================================================#

#use norm.rma, for convenience rename this to norm
norm <- norm.rma.S

#are there NAs in dat?
sum(is.na(norm)) #0 -> no

#what is the range of the data?
range(norm)
# 1.252467 13.106817 -> these indeed are 2log scaled data; there are no negative values

##note: don't do any filtering here, as we want the data set mainly for pathway mapping
##so also keep non/lowly expressed genes in, as being non/lowly expressed also conveys information


#create annotation variables, having HGNC gene symbols, gene types, and gene descriptions
#downloaded from Ensembl BioMart GRCh38.p13 on July 08, 2023
#note this file has the HuGene 1.0 array annotations (BioMart doesn't offer HuGene 1.1 array)
ann <- read.delim("../E-GEOD-36980/Ensembl_Biomart_export_GRCh38_p13_08072023.txt",as.is=TRUE)
ann <- ann[,c(1,5,4,3,2)]
#order ann to have Ensembl IDs in increasing order
ann <- ann[order(ann$Gene.stable.ID),]
#by this ordering, for each probeset ID, the gene with the lowest ENSG identifier will be retrieved
#remove all rows that do not have a AFFY HuGene 1.0st_v1 probeset entry
ann <- ann[!is.na(ann$AFFY.HuGene.1.0.st.v1.probe),]
#rownames(ann) <- ann$ensembl_gene_id
#ann <- ann[,-1]


#==============================================================================#
# Statistical modelling
#==============================================================================#

#model design: take group_ as a fixed effect for an intercept model
design <- model.matrix(~disease,data=desc.S)
colnames(design) <- gsub("[()]","",colnames(design))

#fit model
fit <- lmFit(norm,design)
fit <- eBayes(fit)

#build the (contrast) matrix to compute some group differences of interest based on the model parameters
cont.matrix <- makeContrasts(
  PD_vs_contr_frontal = diseasePD,
  levels = colnames(design)
)
#compute the contrast fits
contrast.fit <- contrasts.fit(fit, cont.matrix)
contrast.fit <- eBayes(contrast.fit)

#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
files <- saveStatOutput(cont.matrix,contrast.fit,postfix="rma",annotation=NULL)

#create summary table of the contrast results
createPvalTab(files,postfix="frontal_rma",html=TRUE)

#add gene annotation to the statistics table
statsF <- read.delim(files,as.is=TRUE,row.names=1)
colnames(statsF) <- paste(colnames(statsF),"frontal",sep="_")
statsF <-  statsF[,c(1,3,5,6)]
statsF <- cbind(statsF,ann[match(rownames(statsF),ann$AFFY.HuGene.1.0.st.v1.probe),])
statsF <- statsF[,-match("AFFY.HuGene.1.0.st.v1.probe",colnames(statsF))]

#save the combined statistics table to file
write.table(statsF,file="E-MTAB-1194_stats.txt",sep="\t",quote=FALSE,col.names=NA)
