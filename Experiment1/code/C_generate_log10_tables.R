# Generates count log10 count tables
#  so that a count of 0 becomes 1e-5 before transformation

# Control flag for taxa name formatting
# When TRUE, names will match Maaslin2 output format
fix_names <- FALSE

# Initialize list to store log10-transformed abundance tables for each tissue
tables_for_plots_all_log10 <- list()

# Process each tissue type separately
for (tis in dig_tissues) {
  # Filter metadata to get tissue-specific samples, excluding pet store mice
  temp_meta <- data.frame(meta_ps_exp1_np) |>
    dplyr::filter(Tissue == tis)  # updated from tissue == tis

  # Extract counts for selected samples from rarefied data
  temp_df <- gen_dig_rar_df_combined[, rownames(temp_meta)]

  # Convert counts to relative abundances
  temp_rel <- sweep(temp_df, 2, colSums(temp_df), "/")

  # Replace zeros with small value before log transformation
  temp_rel[temp_rel == 0] <- 1e-5
  temp_rel2 <- log10(temp_rel)

  # Optionally format taxon names to match Maaslin2 output
  if (fix_names == TRUE) {
    temp_names <- data.frame(matrix(
      nrow = 1, 
      ncol = length(rownames(temp_rel2)), 
      dimnames = list(NULL, rownames(temp_rel2))
    ))
    rownames(temp_rel2) <- colnames(temp_names)
  }

  # Convert to long format and add metadata
  temp_rel2_long <- reshape2::melt(temp_rel2 |>
    rownames_to_column("taxa")) |>
    rename(Sample = variable) |>
    left_join(y = temp_meta, by = "Sample")

  # Store processed data in list
  tables_for_plots_all_log10[[tis]] <- temp_rel2_long

  # Clean up temporary variables
  rm(list = ls(pattern = "^temp"))
}
