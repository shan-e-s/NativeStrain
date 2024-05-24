# Script to set-up server for running bioinformatics scripts
# Install required packages

if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

#BiocManager::install(version = '3.16')

BiocManager::install("dada2")

BiocManager::install("phyloseq")

#BiocManager::install("Biostrings") #this is included in dada i think

print("Finished installing")