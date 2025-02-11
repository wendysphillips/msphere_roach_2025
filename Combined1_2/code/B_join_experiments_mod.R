# Read in data ----
# Load Experiment data from RDS files
exp1_path <- "Experiment1/"
exp2_path <- "Experiment2/"
ps_exp1 <- readRDS(paste0(exp1_path, "ps_exp1.Rds"))
ps_exp2 <- readRDS(paste0(exp2_path, "ps_exp2.Rds"))

# Extract out data tables -----
# Convert phyloseq objects components into data frames
map_exp1 <- data.frame(tax_table(ps_exp1))
otu_exp1 <- data.frame(t(otu_table(ps_exp1)))

map_exp2 <- data.frame(tax_table(ps_exp2))
otu_exp2 <- data.frame(t(otu_table(ps_exp2)))

# Create taxatables -----
# Merge OTU and taxonomy tables for both experiments
tax_exp1 <- merge( map_exp1,otu_exp1, by = 0)
tax_exp1$Row.names <- NULL
tax_exp1b <- tax_exp1 |> tidyr::unite(Kingdom:Genus, col = "taxa", sep = ";")
tax_exp1b$Species <- NULL
tax_exp1b <- tax_exp1b |>
  group_by(taxa) |>
  summarize_all(list(sum))

tax_exp2 <- merge(map_exp2, otu_exp2, by = 0)
tax_exp2$Row.names <- NULL
tax_exp2 <- tax_exp2 |> tidyr::unite(Kingdom:Genus, col = "taxa", sep = ";")
tax_exp2$Species <- NULL
tax_exp2 <- tax_exp2 |>
  group_by(taxa) |>
  summarize_all(list(sum))

# Join taxatables ----
# Merge the taxatables from both experiments and replace NA with 0
tax_all <- merge(tax_exp2, tax_exp1, by = "taxa", all = TRUE)
tax_all[is.na(tax_all)] <- 0

# Make phyloseq type tables ----
# Separate merged taxa into proper tables and prepare for phyloseq object creation
tax_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
otu_all <- tax_all |> separate(taxa, into = tax_names, fill = "right", sep = ";")
rownames(otu_all) <- paste0("ASV", 1:nrow(otu_all))
tax_phy <- otu_all[, 1:7]
otu_all <- otu_all[, (colnames(otu_all) %in% tax_names) == FALSE]

# Import combined metadata ----
# Load and format the combined metadata file, then subset by sample order of otu_all
meta_all <- read.delim("combined_metadata.tsv")
meta_all$Housing_Infection <- gsub("_", " ", meta_all$Housing_Infection)
meta_all$Housing_Infection[meta_all$Housing == "Pet store"] <- "Pet store"
HI_levels <- c("Clean Uninfected", "Clean Infected", "Dirty Uninfected", "Dirty Infected", "Pet store")
meta_all$Housing_Infection <- factor(meta_all$Housing_Infection, 
  levels = HI_levels)
meta_all$Infection <- factor(meta_all$Infection, levels = c("Uninfected", "Infected"))
meta_all <- meta_all |> remove_rownames() |> column_to_rownames("Sample")
meta_all <- meta_all[colnames(otu_all), ]

# Create phyloseq object ------
# Assemble the phyloseq object from OTU, metadata, and taxonomy tables
psall <- phyloseq(
  otu_table(as.matrix(otu_all), taxa_are_rows = TRUE),
  sample_data(meta_all),
  tax_table(as.matrix(tax_phy))
)

# Subset to include only small intestine
ps_si <- microViz::ps_filter(psall, Tissue == "Small")

# Initial phyloseq object filtering ----
# Remove non-bacterial sequences and low-depth samples
ps_si <- subset_taxa(ps_si, Kingdom %in% c("Bacteria"))
ps_si <- subset_taxa(ps_si, Order != "Chloroplast")
ps_si <- subset_taxa(ps_si, Family != "Mitochondria")
ps_si <- prune_samples(sample_sums(ps_si) > 1000, ps_si)

# Extract metadata back out
meta_si <- data.frame(sample_data(ps_si))

# Metadata without pet store mice
meta_no_pet <- meta_si |> dplyr::filter(Genotype == "B6")

# Aggregate to genus -----
# Summarize data at the genus level after taxonomic fixes and validation
ps_si_gen <- ps_si %>%
  tax_fix(unknowns = c(NA, "NA", "unclassified", "Unclassified", "Incertae Sedis")) |>
  tax_agg("Genus") %>%
  phyloseq_validate(remove_undetected = T)

map_si <- data.frame(tax_table(ps_si_gen))
otu_si <- data.frame(otu_table(ps_si_gen))
taxa_si <- merge(otu_si, map_si, by = 0, all = FALSE, all.x = TRUE)
taxa_si$Row.names <- NULL
taxa_si <- taxa_si |> select(Kingdom, Phylum, Class, Order, Family, Genus, everything())

# Make table of just names
taxa_si_names <- taxa_si[, 1:6]

# Now subset out all taxa columns but genus and set rownames accordingly
taxa_si <- taxa_si[c("Genus", rownames(meta_si))]
length(unique(taxa_si$Genus)) == nrow(taxa_si)
taxa_si <- taxa_si |> column_to_rownames("Genus")

# Rarefy then aggregate ----
# Perform rarefaction and aggregate at the genus level on the rarefied data
ps_si_rar <- rarefy_even_depth(ps_si, replace = F, rngseed = 999)

ps_si_rar_gen <- ps_si_rar %>%
  tax_fix(unknowns = c(NA, "NA", "unclassified", "Unclassified", "Incertae Sedis")) |>
  tax_agg("Genus") %>%
  phyloseq_validate(remove_undetected = T)

# Create a mapping from Genus to feature names (for MaAsLin2 outputs)
temp_names <- data.frame(matrix(
  nrow = 1, ncol = length(taxa_si_names$Genus),
  dimnames = list(NULL, taxa_si_names$Genus)
))
taxa_si_names$feature <- colnames(temp_names)

name_map <- taxa_si_names[, c("Genus", "feature")]

# Relative abundance tables ---
# Calculate relative abundances for downstream analyses
temp_df <- taxa_si
temp_df <- temp_df[, rownames(meta_no_pet)]
temp_rel <- sweep(temp_df, 2, colSums(temp_df), "/")
small_rel <- temp_rel
small_rel_log10 <- small_rel
small_rel_log10[small_rel_log10 == 0] <- NA
small_rel_log10 <- log10(small_rel_log10)

# Clean up some things no longer needed
rm(map_exp1, map_exp2, meta_all, otu_all, otu_exp1, otu_exp2, ps_exp1, ps_exp2, psall, tax_all, tax_exp1, tax_exp2, temp_names)
