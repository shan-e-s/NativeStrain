dadaShane <- function(trimF, trimR, truncF, truncR){
  # run from directory which contains script directory
  # sets specific output and input folders
  # ensure path to fastq is set correctly
  # expects taxonomy reference database eg silva in working directory
  
  library(dada2)
  library(phyloseq)
  
  start.time <- format(Sys.time(), "%a %b %d %X %Y")
  
  in.path <- "fastq/" # CHANGE FILEPATH TO FASTQ DIRECTORY
  out.path <- "Data/dadaR_out/"
  referenceDBs <- "ReferenceDBs/" # e.g. Silva reference DB
  multithread <- 10 # CHANGE TO SUIT MACHINE, alternative to number of cores is to set to TRUE/FALSE
  
  
  # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
  fnFs <- sort(list.files(in.path, pattern="_R1_001.fastq.gz$", full.names = TRUE)) # fine to use .gz files and not unzip?
  fnRs <- sort(list.files(in.path, pattern="_R2_001.fastq.gz$", full.names = TRUE))
  
  # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
  sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
  
  # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(in.path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(in.path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  print(paste0("Start filterAndTrim()"))
  # Well use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. 
  # The maxEE parameter sets the maximum number of "expected errors" allowed in a read,
  # which is a better filter than simply averaging quality scores.
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(trimF,trimR), truncLen=c(truncF,truncR),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=multithread, verbose = TRUE) # On Windows set multithread=FALSE
  print("Finished filterAndTrim")
  head(out)
  ## The standard filtering parameters are starting points, not set in stone. If you want to speed up downstream computation, 
  ## consider tightening maxEE. If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the 
  ## reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails. 
  ## Remember though, when choosing  truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later.
  
  # estimate errors (needed for denoising machine learning)
  set.seed(2023)
  errF <- learnErrors(filtFs, multithread=multithread)
  errR <- learnErrors(filtRs, multithread=multithread)
  
  # Sanity check
  
  #dev.new()
  # plotErrors(errF, nominalQ=TRUE)
  
  # apply core sample inference to derep.ed data
  dadaFs <- dada(filtFs, err=errF, multithread=multithread)
  dadaRs <- dada(filtRs, err=errR, multithread=multithread)
  
  # Merge paired reads, to obtain full denoised sequences. Merging aligns forward and reversed reads to construct "contig" sequences
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
  
  # Inspect the merger data.frame from the first sample
  head(mergers[[1]])
  ## if most reads do not merge  check upstream parameters
  
  # Construct Sequence table, [samples,ASV's]
  seqtab <- makeSequenceTable(mergers)
  dim(seqtab)
  
  # save seqtab
  saveRDS(seqtab, file = paste0(out.path, "seqtab.rds"))
  # save others
  save(sample.names, out, errF, errR, dadaFs, dadaRs, mergers, file = paste0(out.path, "dada_out_pt1.rda"))
  
  
  
  #############################################################
  #end dada + create seqtab
  #############################################################
  print("Finished pt.I starting pt.II")
  
  
  # init seqtab.list
  seqtab.list <- character(0)
  
  # seqtab.all <- mergeSequenceTables(tables = seqtab.list) 
  seqtab.all <- seqtab # as only used 1 flow cell/run for this sequencing data i dont need to make a list of seqtabs
  
  
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=multithread, verbose=TRUE)
  sum(seqtab.nochim)/sum(seqtab.all)
  ## if too many seqs lost at this step, possible primers or adapters havent been removed, try trimming first x bases
  print(paste("Chimera's removed. sum(seqtab.nochim)/sum(seqtab.all) = ", sum(seqtab.nochim)/sum(seqtab.all)))
  #saveRDS(seqtab.nochim, file = "seqtab-nochim.rds")
  #seqtab.nochim <- readRDS("seqtab-nochim.rds")
  
  #Assign taxonomy
  ## species level assignment based on exact matching
  taxa <- assignTaxonomy(seqtab.nochim, paste0(referenceDBs, "silva_nr99_v138.1_train_set.fa.gz"), multithread=multithread)
  #saveRDS(taxa, file = paste0("R_ouput/","taxa-genus.rds"))
  taxa <- addSpecies(taxa,paste0(referenceDBs, "silva_species_assignment_v138.1.fa.gz")) # can allowMultiple=T, which returns multiple possibilities for species if there are multiple matches
  # check taxa assignations
  taxa.print <- taxa # Removing sequence rownames for display only
  rownames(taxa.print) <- NULL
  head(taxa.print)
  
  print("Taxonomy assigned")
  print(head(taxa.print))
  
  # create phyloseq object
  ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                 tax_table(taxa))
  
  # change sequence names to ASV### format and store sequences in refseq slot
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  ps
  
  # Track reads through pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  head(track)
  write.csv2(track, file = paste0(out.path, "trackReads.csv"))
  
  # save necessary files
  ## adding seqtab in case i cant merge seqtab.nochim
  save(seqtab.all, seqtab.nochim, taxa, ps, file = paste0(out.path, "dada_outputCOMPLETE.rda"))
  print("DONE")
  
  end.time <- format(Sys.time(), "%a %b %d %X %Y")
  
  print(paste0("Start time: ", start.time, " | End time: ", end.time))
}

