
# NOTE:
# This version adds the requirement that for "Gene Disease Relationship", 
# the selected annotation should be from uniprot. This allows for only picking 
# the protein ('gene') and not the disease in addition it avoids picking part 
# of the false positives because of an error in epmc for some papers 
# (being a few tokens off) to avoid also the 'gene' false positives, 
# one could take the protein referred to from uniprot instead of the name 
# as given by epmc (not done for now, as maybe epmc sometimes has a 
# non-primary id, not sure, but to avoid missing any; also this should be fixed 
# at the epmc tool, not our code)

#==============================================================================#
# Prepare working environment
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load required library
library(stringr)
library(tokenizers)
library(quanteda)
library(stringi)

# set working directory
# Example: setwd("C:/Users/jarno/OneDrive/Documents/GitHub/HDAC6Network/ExtractInteractions/OUTPUT")
setwd("PATH/TO/OUTPUT")

# Set path to corpus paper files
# Example: corpus_path <- "C:/Users/jarno/OneDrive/Documents/GitHub/HDAC6Network/ExtractInteractions/CORPUS"
corpus_path <- "PATH/TO/CORPUS"


# Set path to epmc annotations
# Example: annotation_path <- "C:/Users/jarno/OneDrive/Documents/GitHub/HDAC6Network/ExtractInteractions/ANNOTATIONS"
annotation_path <- "PATH/TO/ANNOTATIONS"

#==============================================================================#
# Make a vector of PDF file names
# --> ONLY RUN ONCE
# NOTE: "xpdf-tools-win-4.04.zip" needs to be unzipped first!
#==============================================================================#

# List PDF files
myfiles <- list.files(path = corpus_path, pattern = "pdf",  full.names = TRUE)

# Convert each PDF file that is named in the vector into a text file 
# Text file is created in the same directory as the PDFs
lapply(myfiles, function(i) system(paste('"xpdf-tools-win-4.04/bin64/pdftotext.exe"', 
                                        paste0('"', i, '"')), wait = FALSE) )

#==============================================================================#
# Perform cleaning
#==============================================================================#

# List txt files
mytxtfiles <- list.files(path = corpus_path, pattern = "txt",  full.names = TRUE)

# Read epmc annotations
epmc_results <- read.delim(paste0(annotation_path,"/epmc_annotations_HDAC6_corpus.txt"), row.names=1, quote="\"", as.is=TRUE)

# Filter annotations for those of type "Gene_Proteins" and "Gene Disease Relationship"
epmc_results_genes <- epmc_results[(epmc_results$type=="Gene_Proteins") | ((epmc_results$type=="Gene Disease Relationship") & grepl("uniprot",epmc_results$uri)),]

# Get all unique gene names from the annotation
gene_names <- unique(epmc_results_genes$name)
labelsSel <- data.frame(labels=gene_names)

# Create log file
logFile <- "example_output.txt"
cat("Output file", file=logFile, append=FALSE, sep = "\n")

print("Start")
for(i in 1:length(mytxtfiles)) {
  
  # Read txt files
  result <- readLines(mytxtfiles[i],encoding = "AINSI")
  print (paste("paper",i,":",mytxtfiles[i]))
  cat(mytxtfiles[i], file=logFile, append=TRUE, sep = "\n")
  
  for (lines in result){
    
    # Get sentences from the paper
    sentences <- tokenize_sentences(stri_enc_toutf8(lines, is_unknown_8bit = TRUE), lowercase = FALSE)
    
    for (sentence in sentences[[1]] ){
      
      # Get words from each sentence
      words <- tokens(sentence) 
      flag <- FALSE
      for (word in types(words)){
        
        # Check if word is HDAC6
        if ((is.element("HDAC6" , types(words)) | (is.element("HDAC-6" , types(words))) | (is.element("Hdac6" , types(words))) |(is.element("hdac6" , types(words)))) ){
          flag <- TRUE; 
        }
        
        # Check if word is any other gene/protein name
        if (    is.element(word, labelsSel$labels)  &&  (!is.element(word, "HDAC6" ))& (!is.element(word, "HDAC-6" )) & (!is.element(word, "Hdac6" )) & (!is.element(word, "hdac6" ))){
          
          
          if (flag) {
            
            # Get gene/protein names
            words_in_labels <-  str_trim(gsub("\t\t","\t",gsub("NA","",paste(labelsSel[unique(match(types(words),labelsSel$labels)),,],collapse="\t"),fixed=TRUE),fixed=TRUE))
            
            # Save gene/protein names and the corresponding sentence
            cat(sentence, words_in_labels, "\n", file=logFile, append=TRUE, sep = "\t")
            break()
          }
        }
      }
    }
  }
}
print("Finish")
