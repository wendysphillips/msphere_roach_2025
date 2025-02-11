
library(dada2)
library(phyloseq)
library(Biostrings)
library(tidyverse)
library(janitor)
library(phangorn)
library(DECIPHER)

# Set path to fastq files
path <- ""
# Set path to working directory
wd <- ""

# Get lists of filenames for processing
fnFs <- sort(list.files(path, pattern="R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2.fastq", full.names = TRUE))

# Simplify sample names, removing prefix that is common to all
sample_names <- fnFs
sample_names <- gsub("_..........-trimmed-R1.fastq", "",sample_names)

# Create file path for filtered files in filtered_uncompressed/ subdirectory
filtFs <- file.path(path, "filtered_uncompressed", paste0(sample_names, "_F.fastq"))
filtRs <- file.path(path, "filtered_uncompressed", paste0(sample_names, "_R.fastq"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# Filter and trim fastq files
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,180),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=FALSE, multithread=FALSE)

# Learn errors in reads
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

#Perform dada analysis on forward and reverse reads
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

# Merge read pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Make initial ASV table
asv_table <- makeSequenceTable(mergers)

# Remove ASV chimeras
asv_table_nochim <- removeBimeraDenovo(asv_table, method="consensus", multithread=FALSE, verbose=TRUE)

# Because assignment of ASV's to taxa is a computationally intensive step 
#   and low count ASV's more likely to be sequencing errors
# Remove ASV's that have overall low abundance
asv_filt <- temp_asv[(rowSums(temp_asv) >= 10) ,]

# Save the filtered ASV set with name final_filtered_asv_table_t.
final_filtered_asv_table_t <- t(asv_filt)

## Assign taxonomy using dada2 and Silva
taxa_silva_full <- assignTaxonomy(final_filtered_asv_table_t, paste0(wd,"data/silva_nr99_v138.1_wSpecies_train_set.fa.gz"), multithread=FALSE)
taxa_silva_sp <- addSpecies(taxa_silva_full, paste0(wd, "data/silva_species_assignment_v138.1.fa.gz"))
taxa_map_silva_sp <- data.frame(tax_table(taxa_silva_sp))
asv_taxa_silva_sp <- merge(taxa_map_silva_sp, asv_filt, by = 0, all.x = TRUE)
asv_taxa_silva_sp <- asv_taxa_silva_sp  |> column_to_rownames("Row.names")

# Read in metadata
meta <- read.delim("../combined_metadata.tsv")
meta <- meta |> column_to_rownames("Sample")
meta <- meta[meta$Experiment == "Experiment 1",]

# Subset to just samples for this analysis
final_filtered_asv_table_t <- final_filtered_asv_table_t[rownames(final_filtered_asv_table_t) %in% rownames(meta),]

### Make phylogenetic tree
# Rename ASVs with shorter names (ASV + number)
ASVs_filt_renamed = DNAStringSet(colnames(final_filtered_asv_table_t))
names(ASVs_filt_renamed) = paste0("ASV", 1:ncol(final_filtered_asv_table_t))

# Make otutable with new numbered ASV names
asv_df_filtered_renamed = final_filtered_asv_table_t
colnames(asv_df_filtered_renamed) = names(ASVs_filt_renamed)

# Make taxa to ASV map with numbered ASVs
taxa_silva_renamed = taxa_silva_full
rownames(taxa_silva_renamed) = names(ASVs_filt_renamed)

# Using DECIPHER, "perform profile-to-profile alignment of multiple unaligned sequences following a guide tree".
alignment2 = AlignSeqs(ASVs_filt_renamed, anchor=NA)

# Transform into phyDat format
phang_align2 <- phyDat(as(alignment2, "matrix"), type="DNA")

# Compute pairwise distances
dm2 <- dist.ml(phang_align2)

# Perform neighbor-joining tree estimation
treeNJ2 <- NJ(dm2)

# Compute the likeliness of the tree
fit2 = pml(treeNJ2, data=phang_align2)
fitGTR2 <- update(fit2, k=4, inv=0.2)

# Build phyloseq objects ----
# Build phyloseq object using phylogenetic tree without optimization
ps_exp1 = phyloseq(
  otu_table(asv_df_filtered_renamed, taxa_are_rows=FALSE),
  sample_data(meta),
  tax_table(taxa_silva_renamed),
  refseq(ASVs_filt_renamed),
  phy_tree(fitGTR2$tree))

saveRDS(object = ps_exp1, file = "ps_exp1c.Rds")
