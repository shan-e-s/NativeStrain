# Filter phylo oject
## filter out low read samples
# remove controls and duplicate samples for main analysis
# decide which repeats/common samples to keep

rm(list=ls())

library(phyloseq)
library(microbiome)

ps <- readRDS("Data/Phylo_objects/phylo-summer.rds")

################################
# which repeats, commons, duplicates to keep

repeat.swapIn <- c("KB10xCxD8xP2", "DW15xAxD15xP2", "CB2xBxD15xP2", "KB916xAxD8xP2", 
                   "KB10xFxP2", "KB910xCxD8xP2", "KB10xMxP2", "KB29xBxD8xP3",
                   "DW49xAxD15xP3", "KB04xCxD8xP3", "CB26xFexP3", "IN25xAxD15xP3",
                   "KB15xFxD8xP3", "CB2xFexP3", "KB29xMx2xP3", "CB20xFexP3", 
                   "CB62xBxD15xP3", "KB22xCxD8xP3")

main.swapOut <- c("KB10xCxD8xP1", "DW15xAxD15xP1","CB2xBxD15xP1","KB916xAxD8xP1",
                  "KB10xFxP1", "KB910xCxD8xP1", "KB10xMxP1","KB29xBxD8xP2", 
                  "DW49xAxD15xP2", "KB04xCxD8xP2", "CB26xFexP2", "IN25xAxD15xP2",
                  "KB15xFxD8xP2", "CB2xFexP2", "KB29xMx2xP2", "CB20xFexP2",
                  "CB62xBxD15xP2", "KB22xCxD8xP2") # keep main: "CB2xMxP2"

contam.discard <- c("KB29xMx1xP3", "DW27xMxP3", "KB22xCxD8xerrorxP3", "IN1xFexP2",
                    "CB2xCxD15xP1", "DW15xFex2xP2", "DW49xFexP1", "KB151xAxD8xP1",
                    "KB151xDxD8xaltxP1")


# edit the following samples to be sample.type = repeat.swapIn | main.swapOut

## in cases where the ID is included in the repeat.swapIn vector list, edit the sampleType to be repeat.swapIn
#ps@sam_data$SampleType <- ifelse(ps@sam_data$Azenta_noPeriod %in% repeat.swapIn, "repeat.swapIn", NULL)
ps@sam_data[ps@sam_data$Azenta_noPeriod %in% repeat.swapIn, "SampleType"] <- "repeat.swapIn"

## in cases where the ID is included in the main.swapOut vector list, edit the sampleType to be main.swapOut
ps@sam_data[ps@sam_data$Azenta_noPeriod %in% main.swapOut, "SampleType"] <- "main.swapOut"

## discard these
ps@sam_data[ps@sam_data$Azenta_noPeriod %in% contam.discard, "SampleType"] <- "contam.discard"

# remove controls, duplicates, mealworm, repeat samples
ps.experimental <- prune_samples(ps@sam_data$SampleType %in% c("main", "repeat.swapIn"), ps) # define those i want to keep

# save phylo object where low read samples are included
saveRDS(ps.experimental, "Data/Phylo_objects/phylo-experimental-all.rds")

################################################################

# According to 05_DataExploreMicrobiome script, i think sample completeness occurs at around 7500 ASVs

#ps1 <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE) # prevalence filter, drop taxa which occur in one sample
#ps <- prune_taxa(taxa_sums(ps)>1, ps) #abundance filter, drop taxa that only occur once

ps.experimental <- prune_samples(sample_sums(ps.experimental)>7500, ps.experimental)

# remove ASV's which now have zero abundance because of removed samples
ps.experimental <- prune_taxa(taxa_sums(ps.experimental) > 1, ps.experimental) 

################################################################
# remove unnecessary columns in metadata
keep.metadata <- c("Nest", "Sample.ID", "Bird.ID", "Date", "Day", "Ring.Mark", 
  "Colour.Ring", "pit.Tag", "Site", "Chick.LetterID", "Age.code", 
  "Age.category", "Sex", "feather.Size", "Wing_mm", "Tarsus_mm", 
  "Weight_g", "Faecal.Sample", "Breeding", "newRing", "Treatment", 
  "Notes", "Main.sample", "Plate", "Plate_sample_no.", 
  "Plate_well", "Azenta_ID", "Azenta_noPeriod", 
  "Qubit_prePCR", "Qubit_prePool", "Date_extracted", 
  "SampleType", "Extraction_notes", "LibPrep_notes", 
  "Ring.Mark.lab", "Post_lab_notes", "NumberOfReads", "DistanceToEdge", 
  "Observed", "Chao1", "se.chao1", "Shannon", "BroodSize.D8", "BroodSize.D15",
  "firstEggLayDate", "lastEggLayDate", "PD", "SR")

ps.experimental@sam_data <- ps.experimental@sam_data[,colnames(ps.experimental@sam_data) %in% keep.metadata]

# save phylo object intended for downstream analysis
saveRDS(ps.experimental, "Data/Phylo_objects/phylo-experimental.rds")