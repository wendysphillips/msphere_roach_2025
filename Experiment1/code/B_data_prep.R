# General preparation of data for analysis

# Load required libraries ----
library(tidyverse)
library(phyloseq)
library(conflicted)
library(ggpubr)
library(microViz)
library(RColorBrewer)
library(Maaslin2)
library(ggVennDiagram)
library(magrittr)
library(vegan)
library(cowplot)

# Utility functions ----
rm_zero_rows <- function(x) {
  x[rowSums(x) > 0, ]
}

# Set plot defaults ----
theme_set(theme_pubr(base_size = 16))

# Define color palettes ----
# Basic color schemes for different plot types
colors3 <- c("#43C699", "#277EB0", "#614D56")
colors5 <- c("#2D7696", "#59B3D9", "#800C0A", "#BA2B29", "#E06967")

# Extended color palette for taxa plots
mycolors <- c(
  "#F5CECE", "#F58989", "#F54949", "#E09CB4F0", "#F02B91E7", "#7D590A", "#A80C0CE4", "#9E6E85",
  "#6B304C", "#E0C1E6", "#DA96E6", "#D55BEB", "#AF0CCC", "#A27EA8", "#7F348A", "#AD7ECC",
  "#8824C7", "#856CEB", "#5B45B5", "#2A0CA3", "#AEB8D6", "#8198DE", "#4E6DD4", "#0B36D1",
  "#7BAAE8", "#7BAAE8", "#317AE0", "#0D51B8", "#B6E1E3", "#62CBD1", "#22B7BF", "#0DD6B8",
  "#08A18A", "#8AEDB1", "#51D686", "#36A864", "#C6E3BF", "#6F9E3C", "#F9CAB2", "#C4CC68",
  "#98A31A", "#6D7511", "#EBD8AE", "#D6B058", "#F2B118", "#9C7206", "#DE35B4", "#A35710",
  "#347199", "#E68C43"
)
# Randomize colors for consistent but random assignment
set.seed(199)
my_colors_random <- sample(mycolors, length(mycolors), replace = FALSE)

# Create yellow-green-blue color gradient
cols <- brewer.pal(9, "YlGnBu")
ygb_pal <- colorRampPalette(cols)

# Load phyloseq object ----
ps_exp1 <- readRDS("ps_exp1.Rds")

# Define sample groupings ----
dig_tissues <- c("Small", "Cecum", "Large")

# Initial phyloseq object filtering ----
# Remove non-bacterial and organelle sequences
ps_exp1_filt <- subset_taxa(ps_exp1, Kingdom %in% c("Bacteria"))
ps_exp1_filt <- subset_taxa(ps_exp1_filt, Order != "Chloroplast")
ps_exp1_filt <- subset_taxa(ps_exp1_filt, Family != "Mitochondria")

# Remove low-depth samples
ps_exp1_filt <- prune_samples(sample_sums(ps_exp1_filt) > 1000, ps_exp1_filt)

# Metadata processing ----
meta_ps_exp1 <- sample_data(ps_exp1_filt)

# Add group column
meta_ps_exp1$group <- if_else(meta_ps_exp1$Housing == "Clean", "Clean", "Dirty")
meta_ps_exp1$group[meta_ps_exp1$Housing == "Dirty"] <- if_else(meta_ps_exp1$Cage[meta_ps_exp1$Housing  == "Dirty"] =="D1_B", "Dirty A", "Dirty B")
meta_ps_exp1$group[meta_ps_exp1$Housing == "Pet store"] <- if_else(meta_ps_exp1$Cage[meta_ps_exp1$Housing  == "Pet store"] =="D1_B", "Dirty A", "Dirty B")

# Set factor levels for consistent ordering
meta_ps_exp1$Housing_Infection <- factor(meta_ps_exp1$Housing_Infection,
  levels = c("Clean_Uninfected", "Clean_Infected", "Dirty_Uninfected", "Dirty_Infected")
)
meta_ps_exp1$Infection <- factor(meta_ps_exp1$Infection, levels = c("Uninfected", "Infected"))
meta_ps_exp1$Tissue <- factor(meta_ps_exp1$Tissue, levels = dig_tissues)

# Update phyloseq object with processed metadata
sample_data(ps_exp1_filt) <- meta_ps_exp1

# Clean taxonomy labels ----
ps_exp1_filt_fix <- ps_exp1_filt |>
  tax_fix(unknowns = c(NA, "NA", "unclassified", "Unclassified", "Incertae Sedis")) |>
  phyloseq_validate(remove_undetected = TRUE)

# Read depth analysis ----
# Convert to dataframe for analysis
ps_exp1_filt_fix_asv <- data.frame(t(ps_exp1_filt_fix@otu_table))

# Calculate read depths
sample_read_totals <- data.frame(colSums(ps_exp1_filt_fix_asv))
colnames(sample_read_totals)[1] <- "Total_counts"
sample_read_totals <- merge(sample_read_totals, meta_ps_exp1, by = 0)

# Check minimum read counts by tissue
min(sample_read_totals$Total_counts[sample_read_totals$Tissue == "Small"])
min(sample_read_totals$Total_counts[sample_read_totals$Tissue == "Cecum"])
min(sample_read_totals$Total_counts[sample_read_totals$Tissue == "Large"])

# Sample subset creation ----
# Add sample IDs to metadata
meta_ps_exp1$Sample <- rownames(meta_ps_exp1)
meta_ps_exp1 <- data.frame(meta_ps_exp1, check.names = FALSE)
dig_samples <- meta_ps_exp1$Sample

# Create non-pet store subset
meta_ps_exp1_np <- meta_ps_exp1[meta_ps_exp1$Housing != "Pet store", ]
meta_ps_exp1_np_samples <- meta_ps_exp1_np$Sample

# Create uninfected housing subsets
dig_housing <- meta_ps_exp1[meta_ps_exp1$Infection == "Uninfected", ]
# Uninfected, no pet store
dig_housing_np <- dig_housing[dig_housing$Housing != "Pet store", ]

# Taxonomic aggregation and subsetting ----
# Aggregate to genus level
phylo_gen_dig <- ps_exp1_filt_fix |> tax_agg("Genus")

# Create various subsets 
phylo_gen_dig_np <- microViz::ps_filter(phylo_gen_dig, Housing != "Pet store")
phylo_gen_dig_housing <- microViz::ps_filter(phylo_gen_dig, Infection == "Uninfected")
phylo_gen_dig_np_housing <- microViz::ps_filter(phylo_gen_dig_np, Infection == "Uninfected")
phylo_gen_dig_small <- microViz::ps_filter(phylo_gen_dig, Tissue %in% c("Small"))
phylo_gen_dig_cl <- microViz::ps_filter(phylo_gen_dig, Tissue %in% c("Cecum", "Large"))

# Create count and relative abundance tables
counts_gen_dig <- data.frame(phylo_gen_dig@otu_table, check.names = FALSE)
gen_dig_table <- data.frame(t(phylo_gen_dig@otu_table), check.names = FALSE)
gen_dig_rel <- sweep(gen_dig_table, 2, colSums(gen_dig_table), "/")

# Rarefaction ----
# Rarefy combined Tissues to lowest depth
phylo_gen_dig_rar <- rarefy_even_depth(phylo_gen_dig,
  rngseed = 999,
  sample.size = min(sample_sums(phylo_gen_dig)),
  replace = FALSE
)

# Separate rarefaction for small intestine
phylo_gen_dig_small_rar <- rarefy_even_depth(phylo_gen_dig_small,
  rngseed = 999,
  sample.size = min(sample_sums(phylo_gen_dig_small)),
  replace = FALSE
)

phylo_gen_dig_cl_rar <- rarefy_even_depth(phylo_gen_dig_cl,
  rngseed = 999,
  sample.size = 0.9 * min(sample_sums(phylo_gen_dig_cl)),
  replace = FALSE
)

phylo_gen_dig_housing_rar <- rarefy_even_depth(phylo_gen_dig_housing,
  rngseed = 999,
  sample.size = min(sample_sums(phylo_gen_dig_housing)),
  replace = FALSE
)


# Process rarefied data ----
# Convert rarefied phyloseq objects to dataframes
gen_dig_small_rar_df <- data.frame(t(phylo_gen_dig_small_rar@otu_table), check.names = FALSE)
gen_dig_housing_rar_df <- data.frame(t(phylo_gen_dig_housing_rar@otu_table), check.names = FALSE)
gen_dig_cl_rar_df <- data.frame(t(phylo_gen_dig_cl_rar@otu_table), check.names = FALSE)

small_samples <- dig_housing_np$Sample[dig_housing_np$Tissue == "Small"]
gen_dig_housing_rar_df_small <- gen_dig_housing_rar_df[, small_samples]
gen_dig_housing_rar_df_small <- rm_zero_rows(gen_dig_housing_rar_df_small)

# Get samples for cecum and large intestine
cecum_large_samples <- dig_housing_np$Sample[dig_housing_np$Tissue %in% c("Cecum", "Large")]
gen_dig_housing_rar_df_cecum_large <- gen_dig_cl_rar_df[, cecum_large_samples]

gen_dig_housing_rar_df_cecum_large <- rm_zero_rows(gen_dig_housing_rar_df_cecum_large)

# Combine rarefied data ----
# Merge small intestine with cecum/large intestine data
gen_dig_housing_rar_df_combined <- merge(gen_dig_housing_rar_df_small, 
  gen_dig_housing_rar_df_cecum_large, 
  by = 0, 
  all = TRUE) |>
  column_to_rownames("Row.names")
gen_dig_housing_rar_df_combined[is.na(gen_dig_housing_rar_df_combined)] <- 0

# Combine rarefied small and separately rarefied cecum & large
gen_dig_rar_df_combined <- merge(gen_dig_small_rar_df, gen_dig_cl_rar_df, by = 0, all = TRUE) |>
  column_to_rownames("Row.names")

# Change NA values to 0
gen_dig_rar_df_combined[is.na(gen_dig_rar_df_combined)] <- 0

# Initialize tissue-specific data structures ----
tissue_tables <- list()
for (tis in dig_tissues) {
  temp <- meta_ps_exp1[meta_ps_exp1$Tissue == tis, ]
  temp <- meta_ps_exp1[meta_ps_exp1$Tissue == tis, ]
  tissue_tables[[tis]][["meta"]] <- temp
  
  # Clean up temporary variables
  rm(list = ls(pattern = "^temp"))
}

# Prepare taxonomy tables for Maaslin2 compatibility
phylo_gen_dig_tax <- data.frame(phylo_gen_dig@tax_table, check.names = FALSE)
phylo_gen_dig_tax_new_names <- phylo_gen_dig_tax
temp_names <- data.frame(
  matrix(
    nrow = 1, 
    ncol = length(rownames(phylo_gen_dig_tax)), 
    dimnames = list(NULL, rownames(phylo_gen_dig_tax))
  )
)
rownames(phylo_gen_dig_tax_new_names) <- colnames(temp_names)
phylo_gen_dig_tax_new_names$feature <- rownames(phylo_gen_dig_tax_new_names)
name_map <- phylo_gen_dig_tax_new_names |> dplyr::select(Genus, feature)

