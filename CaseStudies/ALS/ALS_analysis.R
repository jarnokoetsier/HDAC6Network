# Clear workspace and console
rm(list = ls())
cat("\014") 

library(GEOquery)
library(edgeR)
library(biomaRt)
library(tidyverse)

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE76220", "file=GSE76220_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)


# Load meta data

test <- GEOquery::getGEO("GSE76220")

meta <- test$GSE76220_series_matrix.txt.gz@phenoData@data
meta <- meta[colnames(tbl),]

all(colnames(tbl) == rownames(meta))

# Make DGE list object
d <- DGEList(counts = tbl,
             group = make.names(meta$description))

# remove lowly expressed miRNAs
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE]

# Normalize the data
d <- calcNormFactors(d)

group <- factor(make.names(meta$description))
design <- model.matrix(~0+group)

# Estimate dispersion
d <- estimateDisp(d, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(d, design)

# Compare RTT vs WT
contr <- makeContrasts(groupsporadic.amyotrophic.lateral.sclerosis - groupControl,
                       levels = colnames(coef(fit)))

# Test for significance
test <- glmQLFTest(fit, contrast = contr[,1])
top.table <- topTags(test, n = Inf)$table

top.table$GeneID <- as.character(rownames(top.table))

ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
annotations <- getBM(attributes=c("ensembl_gene_id",
                                  "hgnc_symbol",
                                  "entrezgene_id"), 
                     filters = "entrezgene_id",
                     values = top.table$GeneID,
                     mart = ensembl)
annotations$entrezgene_id <- as.character(annotations$entrezgene_id)

# Get annotations for each probe
ALS_data <- inner_join(top.table, annotations, by = c("GeneID" = "entrezgene_id"))

# Remove probes without associated ensembl ID
ALS_data <- ALS_data[!is.na(ALS_data$ensembl_gene_id),]

# Select the most significant probe for each gene
ALS_data <- arrange(ALS_data, by = PValue)
ALS_data <- ALS_data[!duplicated(ALS_data$ensembl_gene_id),]

save(ALS_data, file = "ALS/ALS_data.RData")
