# Script for recalcualting alpha diversity scores after removing contaminant taxa
## uses phyloseq's estimate_richness function
## could use Tjazi get_asympotic_richness

rm(list=ls())

library(phyloseq)
library(microbiome)

ps <- readRDS("Data/Phylo_objects/phylo-experimental-decontam.rds")

ps.diversity <- estimate_richness(ps, measures=c("Observed", "Shannon", "Chao1"))

# merge diversity w/ phylo object
## putting into if statement in case phyloseq ever changes and upsets order of samples
if (all(rownames(ps@sam_data)==rownames(ps.diversity))) {
ps@sam_data$Observed <- ps.diversity$Observed
ps@sam_data$Chao1 <- ps.diversity$Chao1
ps@sam_data$Shannon <- ps.diversity$Shannon
} else {print("Trying to match diversity scores to incorrect samples.")}


# rownames(metadata) <- metadata$Azenta_noPeriod
# metadata <- merge(metadata, ps.diversity, by = "row.names")
# metadata <- metadata[, !colnames(metadata) %in% c("Row.names")] # drop "Row.names" column introduced by merge

saveRDS(ps, "Data/Phylo_objects/phylo-experimental-decontam.rds")

################################################################################
# recalculate for rarefied data
ps <- readRDS("Data/Phylo_objects/phylo-rarefied.rds")

ps.diversity <- estimate_richness(ps, measures=c("Observed", "Shannon", "Chao1"))

# merge diversity w/ phylo object
## putting into if statement in case phyloseq ever changes and upsets order of samples
if (all(rownames(ps@sam_data)==rownames(ps.diversity))) {
  ps@sam_data$Observed <- ps.diversity$Observed
  ps@sam_data$Chao1 <- ps.diversity$Chao1
  ps@sam_data$Shannon <- ps.diversity$Shannon
} else {print("Trying to match diversity scores to incorrect samples.")}

saveRDS(ps, "Data/Phylo_objects/phylo-rarefied.rds")

################################################################################
# recalculate for prevalence filtered data
ps <- readRDS("Data/Phylo_objects/phylo-experimental-decontamPrevalence.rds")

ps.diversity <- estimate_richness(ps, measures=c("Observed", "Shannon", "Chao1"))

# merge diversity w/ phylo object
## putting into if statement in case phyloseq ever changes and upsets order of samples
if (all(rownames(ps@sam_data)==rownames(ps.diversity))) {
  ps@sam_data$Observed <- ps.diversity$Observed
  ps@sam_data$Chao1 <- ps.diversity$Chao1
  ps@sam_data$Shannon <- ps.diversity$Shannon
} else {print("Trying to match diversity scores to incorrect samples.")}

saveRDS(ps, "Data/Phylo_objects/phylo-experimental-decontamPrevalence.rds")

