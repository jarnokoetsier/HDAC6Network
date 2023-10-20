#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 
 
# Set working directory
PATH <- "PATH/TO/DIRECTORY"
setwd(paste0(PATH, "/CaseStudies/AD"))

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
# Read and set up data variables
#==============================================================================#

#read the dataset
dat <- read.celfiles(list.celfiles("E-GEOD-36980", full.names=TRUE))

#adjust sample names to remove 'CEL'
sampleNames(dat) <- sub(".CEL","",sampleNames(dat),fixed=TRUE)

#create a backup for the QC phase
dat.b <- dat
#extract the probe expression data and set to 2log scale
range(exprs(dat))
# 20 54766
dat.b <- log(exprs(dat.b),2)

#create metadata information
#based on study description as given in the record
desc <- read.delim("sample_table.txt",as.is=TRUE,row.names=3)
desc$sex <- factor(desc$sex,levels=c("m","f"))
desc$tissue <- factor(desc$tissue,levels=c("hippocampus","frontal","temporal"))
desc$disease <- factor(desc$disease,levels=c("control","AD"))

#check order of samples in data and metadata
sum(colnames(dat.b)==desc$assay_name) == dim(dat.b)[2]
# -> TRUE, samples are in same order in dat.b and desc
#rename columns names to more informative names
colnames(dat.b) <- rownames(desc)

#==============================================================================#
# pre-processing and QC
#==============================================================================#

#create QC plots for raw data
factors <- c("sex","tissue","disease","age")
if(!dir.exists("QC_raw")) dir.create("QC_raw")
setwd("QC_raw")
createQCPlots(dat.b, factors, Table=desc, normMeth="", postfix="")
setwd("..")

#normalisation
sum(sampleNames(dat)==desc$assay_name)==length(sampleNames(dat))
# -> TRUE, samples are in same order in dat and desc
sampleNames(dat) <- rownames(desc)
norm.rma <- oligo::rma(dat)
norm.rma <- exprs(norm.rma)

#also load RMA normalised data by the authors
files <- dir("E-GEOD-36980_RMA")
norm.auth <- read.delim(paste("E-GEOD-36980_RMA",files[1],sep="/"),as.is=TRUE,row.names=1)
colnames(norm.auth) <- substr(files[1],1,9)
for (i in (2:length(files))) {
	norm.tmp <- read.delim(paste("E-GEOD-36980_RMA",files[i],sep="/"),as.is=TRUE,row.names=1)
	colnames(norm.tmp) <- substr(files[i],1,9)
	if(identical(rownames(norm.tmp),rownames(norm.auth))) {
		norm.auth <- cbind(norm.auth,norm.tmp)
	} else {
		stop("sample tables don't have equal gene names")
	}
}
sum(colnames(norm.auth)==desc$assay_name)==length(colnames(norm.auth))
# -> TRUE, samples are in same order in norm.auth and desc
colnames(norm.auth) <- rownames(desc)


#create QC plots for norm data
if(!dir.exists("QC_normRMA")) dir.create("QC_normRMA")
setwd("QC_normRMA")
createQCPlots(norm.rma, factors, Table=desc, normMeth="", postfix="RMA")
setwd("..")

if(!dir.exists("QC_normAuth")) dir.create("QC_normAuth")
setwd("QC_normAuth")
createQCPlots(norm.auth, factors, Table=desc, normMeth="", postfix="Auth")
setwd("..")

#save the normalised data tables
write.table(norm.rma,file="E-GEOD-36980_norm_RMA.txt",sep="\t",quote=FALSE,col.names=NA)
write.table(norm.auth,file="E-GEOD-36980_norm_authors.txt",sep="\t",quote=FALSE,col.names=NA)

#==============================================================================#
# prepare for statistics
#==============================================================================#

#use norm.auth, for convenience rename this to norm
norm <- norm.auth

#are there NAs in dat?
sum(is.na(norm)) #0 -> no

#what is the range of the data?
range(norm)
# 1.06455 14.74830 -> these indeed are 2log scaled data; there are no negative values

##note: don't do any filtering here, as we want the data set mainly for pathway mapping
##so also keep non/lowly expressed genes in, as being non/lowly expressed also conveys information


#create annotation variables, having HGNC gene symbols, gene types, and gene descriptions
#downloaded from Ensembl BioMart GRCh38.p13 on July 08, 2023
ann <- read.delim("Ensembl_Biomart_export_GRCh38_p13_08072023.txt",as.is=TRUE)
ann <- ann[,c(1,5,4,3,2)]
#order ann to have Ensembl IDs in increasing order
ann <- ann[order(ann$Gene.stable.ID),]
#by this ordering, for each probeset ID, the gene with the lowest ENSG identifier will be retrieved
#remove all rows that do not have a AFFY HuGene 1.0st_v1 probeset entry
ann <- ann[!is.na(ann$AFFY.HuGene.1.0.st.v1.probe),]
#rownames(ann) <- ann$ensembl_gene_id
#ann <- ann[,-1]

#==============================================================================#
# statistical modelling
#==============================================================================#

#model design: take group_ as a fixed effect for an intercept model
design <- model.matrix(~disease*tissue,data=desc)
colnames(design) <- gsub("[()]","",colnames(design))
colnames(design) <- gsub(":","_x_",colnames(design))

#we cannot take pairing into account as we don't know which samples belong to the same person
#also, between different tissues the extra correlation likely isn't that high

#fit model
fit <- lmFit(norm,design)
fit <- eBayes(fit)

#build the (contrast) matrix to compute some group differences of interest based on the model parameters
cont.matrix <- makeContrasts(
  AD_vs_contr_hippocampus = diseaseAD,
  AD_vs_contr_frontal = diseaseAD+diseaseAD_x_tissuefrontal,
  AD_vs_contr_temporal = diseaseAD+diseaseAD_x_tissuetemporal,
  levels = colnames(design)
)
#compute the contrast fits
contrast.fit <- contrasts.fit(fit, cont.matrix)
contrast.fit <- eBayes(contrast.fit)

#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
files <- saveStatOutput(cont.matrix,contrast.fit,postfix="rma",annotation=NULL)

#create summary table of the contrast results
createPvalTab(files,postfix="rma",html=TRUE)


#==============================================================================#
# per tissue modeling
#==============================================================================#

#==============================================================================#
# prepare for statistics and QC
#==============================================================================#

#split the dataset and the metadata descriptions per tissue
sum(colnames(norm.auth)==rownames(desc))==length(colnames(norm.auth))
# -> TRUE, samples are in same order in norm.auth and desc
norm.hippocampus <- norm[,desc$tissue=="hippocampus"]
norm.frontal <- norm[,desc$tissue=="frontal"]
norm.temporal <- norm[,desc$tissue=="temporal"] 
desc.hippocampus <- desc[desc$tissue=="hippocampus",]
desc.frontal <- desc[desc$tissue=="frontal",]
desc.temporal <- desc[desc$tissue=="temporal",]
factors2 <- factors[factors!="tissue"]

#create QC plots for norm data per tissue
if(!dir.exists("QC_normAuth_hippocampus")) dir.create("QC_normAuth_hippocampus")
setwd("QC_normAuth_hippocampus")
createQCPlots(norm.hippocampus, factors2, Table=desc.hippocampus, normMeth="", postfix="AuthHip")
setwd("..")

#create QC plots for norm data per tissue
if(!dir.exists("QC_normAuth_frontal")) dir.create("QC_normAuth_frontal")
setwd("QC_normAuth_frontal")
createQCPlots(norm.frontal, factors2, Table=desc.frontal, normMeth="", postfix="AuthFro")
setwd("..")

#create QC plots for norm data per tissue
if(!dir.exists("QC_normAuth_temporal")) dir.create("QC_normAuth_temporal")
setwd("QC_normAuth_temporal")
createQCPlots(norm.temporal, factors2, Table=desc.temporal, normMeth="", postfix="AuthTem")
setwd("..")

#==============================================================================#
# statistical modelling
#==============================================================================#

#model design: take group_ as a fixed effect for an intercept model
designH <- model.matrix(~disease,data=desc.hippocampus)
colnames(designH) <- gsub("[()]","",colnames(designH))
designF <- model.matrix(~disease,data=desc.frontal)
colnames(designF) <- gsub("[()]","",colnames(designF))
designT <- model.matrix(~disease,data=desc.temporal)
colnames(designT) <- gsub("[()]","",colnames(designT))

#fit model
fitH <- lmFit(norm.hippocampus,designH)
fitH <- eBayes(fitH)
fitF <- lmFit(norm.frontal,designF)
fitF <- eBayes(fitF)
fitT <- lmFit(norm.temporal,designT)
fitT <- eBayes(fitT)

#build the (contrast) matrix to compute some group differences of interest based on the model parameters
cont.matrixH <- makeContrasts(
  AD_vs_contr_hippocampus = diseaseAD,
  levels = colnames(designH)
)
cont.matrixF <- makeContrasts(
  AD_vs_contr_frontal = diseaseAD,
  levels = colnames(designF)
)
cont.matrixT <- makeContrasts(
  AD_vs_contr_temporal = diseaseAD,
  levels = colnames(designT)
)

#compute the contrast fits
contrast.fitH <- contrasts.fit(fitH, cont.matrixH)
contrast.fitH <- eBayes(contrast.fitH)
contrast.fitF <- contrasts.fit(fitF, cont.matrixF)
contrast.fitF <- eBayes(contrast.fitF)
contrast.fitT <- contrasts.fit(fitT, cont.matrixT)
contrast.fitT <- eBayes(contrast.fitT)

#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
filesH <- saveStatOutput(cont.matrixH,contrast.fitH,postfix="rma",annotation=NULL)
filesF <- saveStatOutput(cont.matrixF,contrast.fitF,postfix="rma",annotation=NULL)
filesT <- saveStatOutput(cont.matrixT,contrast.fitT,postfix="rma",annotation=NULL)

#create summary table of the contrast results
createPvalTab(filesH,postfix="hippocampus_rma",html=TRUE)
createPvalTab(filesF,postfix="frontal_rma",html=TRUE)
createPvalTab(filesT,postfix="temporal_rma",html=TRUE)

#combine the results into one file for further processing(mapping to the HDAC6 pathway)
statsH <- read.delim(filesH,as.is=TRUE,row.names=1)
statsH <- statsH[order(rownames(statsH)),]
colnames(statsH) <- paste(colnames(statsH),"hippocampus",sep="_")
statsF <- read.delim(filesF,as.is=TRUE,row.names=1)
statsF <- statsF[order(rownames(statsF)),]
colnames(statsF) <- paste(colnames(statsF),"frontal",sep="_")
statsT <- read.delim(filesT,as.is=TRUE,row.names=1)
statsT <- statsT[order(rownames(statsT)),]
colnames(statsT) <- paste(colnames(statsT),"temporal",sep="_")
sum(rownames(statsH)==rownames(statsF))==dim(statsH)[1]
sum(rownames(statsH)==rownames(statsT))==dim(statsH)[1]

#create the combined statistics table
statsHFT <-  cbind(statsH[,c(1,3,5,6)],statsF[,c(1,3,5,6)],statsT[,c(1,3,5,6)])
statsHFT <- cbind(statsHFT,ann[match(rownames(statsHFT),ann$AFFY.HuGene.1.0.st.v1.probe),])
statsHFT <- statsHFT[,-match("AFFY.HuGene.1.0.st.v1.probe",colnames(statsHFT))]

#save the combined statistics table to file
write.table(statsHFT,file="E-GEOD-36980_stats.txt",sep="\t",quote=FALSE,col.names=NA)
