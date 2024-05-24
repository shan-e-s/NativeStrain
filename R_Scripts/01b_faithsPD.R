# Script to calculate faiths PD for diversity analyses following PhD viva corrections, October 2023

rm(list = ls())

# Susan Holmes recommendeded picante::comm.pd 

#install.packages("btools") # doesnt work for this version of r
#install.packages("picante")
library(picante)
library(phyloseq)

########################################################################################
# build tree

set.seed(1189)

#BiocManager::install("DECIPHER", "phangorn")
library(phangorn)
library(dada2)
library(DECIPHER)

ps <- readRDS("Data/Phylo_objects/phylo-experimental-decontamPrevalence.rds")

seqs <- ps@refseq
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose=T, processors = 10)

# The phangorn R package is then used to construct a phylogenetic tree. 
# Here we first construct a neighbor-joining tree, and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) 
# maximum likelihood tree using the neighbor-joining tree as a starting point.

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)

# originally ran with stochastic arrangement and had to kill process after 122 hours
# ran NNI subsequently which was much faster and compared Faiths PD results with each tree
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "NNI", control = pml.control(trace = 0))

detach("package:phangorn", unload=TRUE)

# add tree to phylo object
tree <- phy_tree(fitGTR$tree)
ps.tree <- merge_phyloseq(ps, tree)

# not workiung because taxa names dont match
# probably need to correct names(seq) command above
# in the meantime i will simply rename as asvs assuming seqs in same order.. ## VERIFIED THIS IS TRUE
tree.ASV.labels <- tree
tree.ASV.labels$tip.label <-  taxa_names(ps)

# add tree to phylo object
ps.tree <- merge_phyloseq(ps, tree.ASV.labels)

#######################################
## Picante uses the phylo format implemented in the ape package to represent phylogenetic relationships
##among taxa.

## Picante uses the same community data format as the vegan package - a matrix
## or data.frame with sites/samples in the rows and taxa in the columns. The elements of this data frame should be numeric values indicating the abundance or
## presence/absence (0/1) of taxa in different samples

# comm <- otu_table(ps)
# phy <- phy_tree(fitGTR$tree)

comm <- otu_table(ps.tree)
phy <- tree.ASV.labels

# getting error that i need rooted tree but include root command not fixing it
# had to use include.root=F
comm.pd <- pd(comm, phy, include.root = F) #, include.root = T
head(comm.pd)

# save faiths PD
write.csv(comm.pd, "Data/R_savedObjects/faithsPD_NNI.csv")

# NNI and (arrested) stochastic trees yield similar PD values
