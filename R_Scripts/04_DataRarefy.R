# Script to create rarefied version of phylo-experimental object

rm(list=ls())

library(phyloseq)

ps <- readRDS("Data/Phylo_objects/phylo-experimental-decontam.rds")

ps.rarefied <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)),
                                rngseed = 1189, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

saveRDS(ps.rarefied, file = "Data/Phylo_objects/phylo-rarefied.rds")
