# Set up phyloseq object for analyses, Feb 2023
## Match experiment data i.e. birdData.csv with 'LabSampleData.csv'

rm(list=ls())

#install.packages("car")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("metagenomeSeq")

library(tidyverse)
library(phyloseq)
library(microbiome)
#library(metagenomeSeq)

# read in experiment metadata, sample metadata and sequencing data
birdData <- read_csv("Data/Dataframes_raw_csv/birdData.csv") # alternatively could read in individual df and filter out birds i dont have samples for | maybe those samples dropped automatically by phyloseq?
LabSampleData <- read_csv("Data/Dataframes_R_out/LabSampleData.csv")
load("Data/dadaR_out/dada_outputCOMPLETE.rda") 

# remove unnecessary objects 
rm(list = c("seqtab.all", "seqtab.nochim", "taxa"))

# verify phyloseq object read in ok
ps

# change ps sample names to match metadata
## drop .fastq.gz
sample_names(ps) <- sub("_R1_001.fastq.gz","", sample_names(ps))

####################### MERGE LAB DATA AND BIRD DATA ##############################################
# rename sample_ID to sample.ID to aid merge
LabSampleData <- LabSampleData %>% dplyr::rename(Sample.ID = Sample_ID)

# make sample IDs (and nest IDs) consistent by dropping zeroes between site abbreviation and nest number
## if first occurrence of number == 0 then remove this number
### need to remove "\" in sample names (diet samples) first as this is upsetting regex
### no actually regex being upset by no matches i.e. diet samples have no digits so are returning "logical(0)" instead of FALSE
### so adding logical statement that says if there is no digits, and hence character(0) is returned by regmatch then logical statement becomes false
### i will still remove the "\"'s because they may cause trouble later
for (i in 1:length(LabSampleData$Sample.ID)) {
  LabSampleData$Sample.ID[i] <- gsub("\\/",".",LabSampleData$Sample.ID[i])
  LabSampleData$Sample.ID[i] <- ifelse(regmatches(LabSampleData$Sample.ID[i], regexpr("\\d", LabSampleData$Sample.ID[i])) == "0" && 
                                         !identical(regmatches(LabSampleData$Sample.ID[i], regexpr("\\d", LabSampleData$Sample.ID[i])), character(0)),
                                       sub("0","",LabSampleData$Sample.ID[i]),LabSampleData$Sample.ID[i])
}

# drop zeroes in birdData as well
for (i in 1:length(birdData$Sample.ID)) {
  birdData$Sample.ID[i] <- ifelse(regmatches(birdData$Sample.ID[i], regexpr("\\d", birdData$Sample.ID[i])) == "0" && 
                                         !identical(regmatches(birdData$Sample.ID[i], regexpr("\\d", birdData$Sample.ID[i])), character(0)),
                                       sub("0","", birdData$Sample.ID[i]), birdData$Sample.ID[i])
}

# combine experiment and sample metadata
metadata <- merge(birdData, LabSampleData, by = "Sample.ID", all.y = TRUE) # all.y necessary to keep control samples

# Verify ring numbers match from different dataframes
# Rename "Bird_ID" column i.e. ring numbers
#check <- cbind(metadata$Ring.Mark, metadata$Bird_ID)
metadata <- metadata %>% dplyr::rename("Ring.Mark.lab" ="Bird_ID")

# Consolidate the 2 ring/ID columns
## if ring.mark is NA replace withring.mark.lab
## note: in some cases ring.mark has only nestling mark e.g. "BLUE-RIGHT" or "UNR" and i am missing this case for now
metadata$Ring.Mark <- ifelse(is.na(metadata$Ring.Mark), metadata$Ring.Mark.lab, metadata$Ring.Mark)

# set rownames to match sequencing names ie. Azenta sample names
rownames(metadata) <- metadata$Azenta_noPeriod
#############################################################################################


######################## FILTERING ##############################################
#Note: removing chloroplasts and mitochondria BEFORE dropping low read samples otherwise ends up w/ low read samples in phylo.object
#Note: subset_taxa also removes NA's ie. "The OTUs retained in the dataset is equivalent to x[subset & !is.na(subset)]" unless is.na() given.
# remove chloroplasts
ps <- subset_taxa(ps, Order!="Chloroplast" | is.na(Order)) 
# remove mitochondria
ps <- subset_taxa(ps,Family!="Mitochondria"| is.na(Family))
# remove Eukaryota
ps <- subset_taxa(ps,Kingdom!="Eukaryota"| is.na(Kingdom)) 
# remove Archaea
ps <- subset_taxa(ps,Kingdom!="Archaea"| is.na(Kingdom)) # no Archaea

# add sequence read numbers to samples
seq_reads <- sample_sums(ps)
metadata <- merge(metadata, seq_reads, by = "row.names") %>% dplyr::rename("NumberOfReads"="y")
rownames(metadata) <- metadata$Azenta_noPeriod
metadata <- metadata[, !colnames(metadata) %in% c("Row.names")] # drop "Row.names" column introduced by merge

# remove low read samples <40,000
## maybe try 5000 or 10000, as there is natural separation there and it isnt so strict
#metadata <- metadata[metadata$NumberOfReads>10000,] # drops 258-247=11 samples

############################################################################################

######################## CREATE NEW VARIABLES ################################################ 
#Add distance to edge variable for each nest. Calculating in separate script and sourcing.
source("R_scripts/01a_DistanceToEdge.R")

# need to merge on nest.ID
metadata <- merge.data.frame(metadata, nestbox.edge[,c("nest","DistanceToEdge")], by.x = "Nest", by.y = "nest", all.x=TRUE)

######################### brood size variable ######################## 
# add brood size when sampled
nest.data <- read.csv("Data/Dataframes_raw_csv/NestRecords2021.csv", header = TRUE)
nest.data <- nest.data[,c("site.boxNumber", "BroodSize.D8", "BroodSize.D15", "firstEggLayDate", "lastEggLayDate")]

# only take nests for which we have birds
nest.data <- nest.data[nest.data$site.boxNumber %in% metadata$Nest,] %>% dplyr::rename("Nest" = "site.boxNumber")

metadata <- merge(metadata, nest.data, by="Nest", all.x = T)
########################  Prepare metadata to join phyloseq object ################################################ 

# merge w/ phylo object
phylo.metadata <- sample_data(metadata) # labelling as metadata
rownames(phylo.metadata) <- phylo.metadata$Azenta_noPeriod
phylo.summer <- merge_phyloseq(phylo.metadata, ps) 


################################################################################################  
########################  CALCULATE ALPHA DIVERSITY ########################################################################  
ps.diversity <- estimate_richness(phylo.summer, measures=c("Observed", "Shannon", "Chao1"))

# add phylogenetic diversity
# need to run script 01b_faithsPD.R to get this csv, takes a long time to run so not sourcing
faithsPD <- read.csv(file = "Data/R_savedObjects/faithsPD.csv") %>% column_to_rownames(., var = "X")
ps.diversity <- merge(ps.diversity, faithsPD, all.x = T, by = "row.names") %>% 
  column_to_rownames(., var = "Row.names")

# merge diversity w/ phylo object
rownames(metadata) <- metadata$Azenta_noPeriod
metadata <- merge(metadata, ps.diversity, by = "row.names")
metadata <- metadata[, !colnames(metadata) %in% c("Row.names")] # drop "Row.names" column introduced by merge


########################  Prepare metadata to re-join phyloseq object ################################################ 

# set as numeric
numeric.columns <- c( "Day","Wing_mm", "Tarsus_mm", "Weight_g", "Plate_sample_no.", "Qubit_prePCR", 
                      "Qubit_prePool", "NumberOfReads", "DistanceToEdge", "Observed", 
                      "Chao1", "se.chao1", "Shannon", "BroodSize.D8", "BroodSize.D15", "firstEggLayDate", "lastEggLayDate",
                      "PD", "SR")

# set columns as numeric vectors
metadata <- metadata %>% mutate_at( numeric.columns, as.numeric)

# convert date variable to date type
metadata$Date_extracted <- as.Date(metadata$Date_extracted,"%d/%m/%Y")

phylo.metadata <- sample_data(metadata) # labelling as metadata
rownames(phylo.metadata) <- phylo.metadata$Azenta_noPeriod
phylo.summer <- merge_phyloseq(phylo.metadata, ps) 

########################  Redefine 'Probiotic' as L. kimchicus ########################################################################  
phylo.summer@sam_data$Treatment <- ifelse(phylo.summer@sam_data$Treatment=="Probiotic", "L. kimchicus", phylo.summer@sam_data$Treatment)

################################################################################################
######################## SAVE OUTPUT ################################################
# save data as R object
saveRDS(phylo.summer, file = "Data/Phylo_objects/phylo-summer.rds")

