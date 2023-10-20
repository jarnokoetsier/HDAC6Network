#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
PATH <- "PATH/TO/DIRECTORY"
setwd(paste0(PATH, "/CaseStudies/AD"))

#==============================================================================#
# Merge AD and PD statistics
#==============================================================================#

#read the AD and PD statistics tables
AD.dat <- read.delim("AD/table_AD_vs_contr_frontal_rma.tab",as.is=TRUE,row.names=1)
PD.dat <- read.delim("PD/table_PD_vs_contr_frontal_rma.tab",as.is=TRUE,row.names=1)

#change the column headers to indicate which data set they are from
colnames(AD.dat) <- paste(colnames(AD.dat),"AD",sep="_")
colnames(PD.dat) <- paste(colnames(PD.dat),"PD",sep="_")

#reorder the data sets to increasing Affymetrix probeID order
AD.dat <- AD.dat[order(rownames(AD.dat)),]
PD.dat <- PD.dat[order(rownames(PD.dat)),]

#check whether both datasets contain exactly the same probeIDs and in the same order
identical(rownames(AD.dat),rownames(PD.dat))
#TRUE -> identical

#combine the two data tables into one merged table
dat <- cbind(AD.dat,PD.dat)

#keep only the columns we want to map later (remove Fold.Change, which is derived from logFC, and t)
dat <- dat[,-c(2,4,8,10)]

#save the combined table to file
write.table(dat,file="table_AD_PD_vs_contr_frontal_rma.tab",sep="\t",col.names=NA,quote=FALSE)
