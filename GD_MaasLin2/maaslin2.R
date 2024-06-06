#Gabrielle Davidson's MaAsLin2 analyses for Native Strain microbiome paper 2024

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")

library(Maaslin2)
#MAASLIN

meta_data <-read.csv("metadata.csv", header=TRUE)
otu_data <- read.table("otu_table.tsv",  header=TRUE)
rownames(otu_data) <- otu_data[,1]
otu_data[,1] <- NULL

rownames(meta_data) <- meta_data[,1]
meta_data[,1] <- NULL


fit_data = Maaslin2(
  input_data = otu_data, 
  input_metadata = meta_data, 
  output = "Maaslin2_output_AllFixedTerms", 
  fixed_effects = c("Treatment", "Age.category"),
  random_effects = c("Site", "Nest"),
  reference = c("Age.category,D8"),
  plot_scatter = TRUE)

KO <- read.table(file='KO_pred_metagenome_unstrat.tsv', header=TRUE)
#Assigning row names from 1st column 
rownames(KO) <- KO[,1] 
KO = subset(KO, select = -function. )
KO[] <- lapply(KO, as.integer)
PA <- read.table(file='path_abun_unstrat.tsv', header=TRUE)
#Assigning row names from 1st column 
rownames(PA) <- PA[,1] 
##remove the duplicate column
PA = subset(PA, select = -pathway )
##make values intergers
PA[] <- lapply(PA, as.integer)


fit_data2 = Maaslin2(
  input_data = KO, 
  input_metadata = meta_data, 
  output = "Maaslin2_output_KO_AllFixedTerms", 
  fixed_effects = c("Treatment", "Age.category"),
  random_effects = c("Site", "Nest"),
  reference = c("Age.category,D8"),
  plot_scatter = FALSE)

fit_data3 = Maaslin2(
  input_data = PA, 
  input_metadata = meta_data, 
  output = "Maaslin2_output_PA_AllFixedTerms", 
  fixed_effects = c("Treatment", "Age.category"),
  random_effects = c("Site", "Nest"),
  reference = c("Age.category,D8"),
  plot_scatter = FALSE)