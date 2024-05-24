# Script that build all phyloseq related data products
## source this script when any upstream changes are made to data building scripts
## what about R_bioinformatics/04_? rmd files so not easy to source

rm(list = ls())

source("R_Scripts/01_DataSetup_bySample.R")

source("R_Scripts/02_DataFilter_samples.R")

source("R_Scripts/03_DataFilter_decontam.R")

source("R_Scripts/04_DataRarefy.R")

source("R_Scripts/05_recalculate_AlphaDiversity.R")

rm(list = ls())

# The phyloseq objects i have :
# phylo-summer: has all samples and all the taxa, output by 01_DataSetup_bySample.R
# phylo-experimental: has only bird samples with low read samples removed (<7500 reads), all lab controls and duplicates removed
# phylo-experimental-all: has only bird samples includes low read samples, all lab controls and duplicates removed
# phylo-experimental-decontam: only bird samples with decontaminant taxa removed. Alpha diversity recalculated. 
# phylo-experimental-decontamPrevalence: *main analysis object* only bird samples with decontaminant taxa removed using stricter 'Prevalence' method. Alpha diversity recalculated. 
# phylo-rarefied: only bird samples with decontaminants removed AND samples have been rarefied to even depth. Alpha diversity recalculated.
