# Run Maaslin2 for each tissue and group ----
m_infection_results <- list()
for (tis in dig_tissues) {
  print(tis)
  # Get tissue metadata subset, exclude pet store mice
  temp_meta <- data.frame(meta_ps_exp1_np) |>
    dplyr::filter(Tissue == tis)

  # Run analysis for each group within tissue
  for (gp in unique(temp_meta$group)) {
    print(gp)
    # Filter metadata for current group
    temp_meta_group <- temp_meta[temp_meta$group == gp, ]
    # Get abundance data for current samples
    temp_df <- gen_dig_rel[, rownames(temp_meta_group)]
    # Filter low abundance features
    temp_df <- temp_df[, apply(temp_df, 2, (function(x) sum(x > 0.001) >= 2))]
    # Run Maaslin2 analysis
    temp_fit_data <- Maaslin2(temp_df, temp_meta_group, paste0("maaslin2_", tis, "_", gp),
      fixed_effect = c("Infection"),
      reference = c("Uninfected"),
      standardize = FALSE
    )
    m_infection_results[[tis]][[gp]] <- temp_fit_data
  }
  # Clean up temporary variables
  rm(list = ls(pattern = "^temp"))
}

group_v <- unique(meta_ps_exp1_np$group)

# Process Maaslin2 results ----
# Add original taxonomy name to results
m_list2 <- list()
for (tis in dig_tissues) {
  for (gp in group_v) {
    temp <- m_infection_results[[tis]][[gp]]$results
    temp <- left_join(temp, name_map)
    m_list2[[tis]][[gp]] <- temp
    rm(temp)
  }
}

# Combine results by tissue
m_small_infection <- data.table::rbindlist(m_list2$Small, use.names = TRUE, fill = TRUE, idcol = "Group")
m_cecum_infection <- data.table::rbindlist(m_list2$Cecum, use.names = TRUE, fill = TRUE, idcol = "Group")
m_large_infection <- data.table::rbindlist(m_list2$Large, use.names = TRUE, fill = TRUE, idcol = "Group")

# Create list of results by tissue
m_infection_results_list <- list(
  Small = m_small_infection,
  Cecum = m_cecum_infection,
  Large = m_large_infection
)

# Combine all results into single dataframe
m_infection_results_df <- data.table::rbindlist(m_infection_results_list,
  use.names = TRUE,
  fill = TRUE,
  idcol = "Tissue"
)

# Add full taxonomy information
m_infection_results_df <- left_join(m_infection_results_df,
  phylo_gen_dig_tax_new_names |>
    select(-Genus),
  by = "feature"
)

# Reorder columns for clarity
m_infection_results_df <- m_infection_results_df |>
  dplyr::select(Tissue, feature, Kingdom, Phylum, Class, Order, Family, Genus, everything())

# Extract significant results ----
# Filter for significant features
m_infection_results_sig <- lapply(m_infection_results_list, function(x) {
  x |> dplyr::filter(qval <= 0.05)
})

# Extract genus names of significant features
m_infection_results_sig_taxa <- lapply(m_infection_results_sig, function(x) {
  x |> dplyr::select(Genus)
})

# Create comprehensive results for significant taxa
m_sig_taxa_results_all <- list()
for (tis in dig_tissues) {
  # Get all results for significant genera
  temp <- m_infection_results_list[[tis]] |>
    dplyr::filter(Genus %in% m_infection_results_sig_taxa[[tis]]$Genus)
  m_sig_taxa_results_all[[tis]] <- temp
}
