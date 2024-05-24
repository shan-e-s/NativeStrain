# Decontam
## script to run decontam process on phylo-experimental i.e. main phylo object for downstream analysis
## following tutorial: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

rm(list=ls())

#BiocManager::install("decontam")

library(phyloseq)
library(decontam)
library(ggplot2)

ps <- readRDS("Data/Phylo_objects/phylo-summer.rds")
ps.experimental <- readRDS("Data/Phylo_objects/phylo-experimental.rds")

# remove unnecessary samples
ps <- prune_samples(ps@sam_data$SampleType %in% c("main", "control"), ps)

# inspect library sizes
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=SampleType)) + geom_point()

##################################################################################
# ID contaminants: frequency
contamdf.freq <- isContaminant(ps, method="frequency", conc="Qubit_prePool")
head(contamdf.freq)

table(contamdf.freq$contaminant)

head(which(contamdf.freq$contaminant))

# inspect a couple more of the ASVs that were classified as contaminants to ensure they look like what we expect:
set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),10)], conc="Qubit_prePool") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

# remove contaminants from the main experimental phyloseq object
## keep taxa that are not in the contaminants.ASV vector
contaminanats.ASV <- row.names(contamdf.freq[contamdf.freq$contaminant,])

ps.noncontam <- prune_taxa(!taxa_names(ps.experimental) %in% contaminanats.ASV, ps.experimental)

#ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps.experimental)

ps.noncontam

saveRDS(ps.noncontam, file = "Data/Phylo_objects/phylo-experimental-decontam.rds")


##################################################################################
# ID contaminants: prevalence

# summarize control vs main as logical variable
sample_data(ps)$is.neg <- sample_data(ps)$SampleType == "control"

# the default threshold for a contaminant is that it reaches a probability of 0.1 
#threshold=0.5, that will identify as contaminants all sequences thare are more prevalent in negative controls than in positive samples.
# batch option allows sequencing runs (or extraction or PCR runs) to be analysed separately, as they might be expected to have different contaminant profiles
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5, batch = sample_data(ps)$Plate)

table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# Make phyloseq object of presence-absence in negative controls and true samples
# ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
# ps.pa.neg <- prune_samples(sample_data(ps.pa)$SampleType == "control", ps.pa) #ps.pa.neg <- prune_samples(sample_sums(ps.pa.neg)>0, ps.pa.neg)
# ps.pa.pos <- prune_samples(sample_data(ps.pa)$SampleType == "main", ps.pa)
# # Make data.frame of prevalence in positive and negative samples
# df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
#                     contaminant=contamdf.prev$contaminant)
# ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
#   xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants from the main experimental phyloseq object
## keep taxa that are not in the contaminants.ASV vector
contaminanats.ASV <- row.names(contamdf.prev[contamdf.prev$contaminant,])

ps.noncontam <- prune_taxa(!taxa_names(ps.experimental) %in% contaminanats.ASV, ps.experimental)

ps.noncontam

saveRDS(ps.noncontam, file = "Data/Phylo_objects/phylo-experimental-decontamPrevalence.rds")
