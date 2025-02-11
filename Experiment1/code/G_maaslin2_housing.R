# Produces Figures 1G,H,I,J; Supplementary figures 1A,B,C
# Run Maaslin2 analysis by tissue ----
m_results_housing <- list()
for (tis in dig_tissues) {
  print(tis)
  # Filter metadata for uninfected samples only
  temp_meta <- data.frame(meta_ps_exp1_np) |>
    dplyr::filter(Tissue == tis) |>
    dplyr::filter(Infection == "Uninfected")

  # Prepare abundance data
  temp_df <- t(gen_dig_rel[, rownames(temp_meta)])

  # Filter low abundance features
  temp_df <- temp_df[, apply(temp_df, 2, (function(x) sum(x > 0.001) >= 2))]

  # Run Maaslin2
  temp_fit_data <- Maaslin2(temp_df, temp_meta, paste0("maaslin2_housing_", tis),
    fixed_effect = c("Housing"),
    reference = c("Clean"),
    standardize = FALSE
  )
  m_results_housing[[tis]] <- temp_fit_data

  # Clean up temporary variables
  rm(list = ls(pattern = "^temp"))
}

# Process Maaslin2 results ----
# Set significance threshold
p_lim <- 0.05

# Add taxonomy information to results
m_housing_list2 <- list()
for (tis in dig_tissues) {
  temp <- m_results_housing[[tis]]$results
  temp <- left_join(temp, name_map)
  m_housing_list2[[tis]] <- temp

  # Clean up temporary variables
  rm(list = ls(pattern = "^temp"))
}

# Combine results across tissues
m_housing_results_df <- data.table::rbindlist(m_housing_list2,
  use.names = TRUE,
  fill = TRUE,
  idcol = "Tissue"
)

# Add full taxonomy information
m_housing_with_tax_df <- left_join(m_housing_results_df,
  phylo_gen_dig_tax_new_names,
  by = "feature"
)
m_housing_with_tax_df <- m_housing_with_tax_df |>
  dplyr::select(Tissue, feature, Kingdom, Phylum, Class, Order, Family, everything())

# Extract significant results ----
# Filter for significant features
m_housing_results_sig <- lapply(
  m_housing_list2,
  function(x) x |> dplyr::filter(qval <= p_lim)
)

m_housing_results_sig_taxa <- lapply(
  m_housing_results_sig,
  function(x) x |> dplyr::select(Genus)
)

# Process significant taxa results
m_housing_sig_taxa_results <- list()
for (tis in dig_tissues) {
  temp <- m_housing_list2[[tis]] |>
    dplyr::filter(Genus %in% m_housing_results_sig_taxa[[tis]]$Genus)
  m_housing_sig_taxa_results[[tis]] <- temp
  # Clean up temporary variables
  rm(temp)
}

# Extract genus names for Venn diagram
for (tis in dig_tissues) {
  temp <- m_housing_results_sig_taxa[[tis]]$Genus
  m_housing_results_sig_taxa[[tis]] <- temp
  # Clean up temporary variables
  rm(temp)
}

# Generate Venn diagram of significant taxa ----
# Produces Figure 1J
ggVennDiagram(m_housing_results_sig_taxa, label = "count", show_intersect = F) +
  scale_fill_distiller(palette = "Spectral") +
  theme(
    legend.position = "none",
    text = element_text(family = "ArialMT")
  )

# Visualization preparation ----
# Define color schemes for annotations
annoCol <- list(Housing = c(Clean = colors5[2], Dirty = colors5[5], Pet_store = colors5[3]))
annoCol_np <- list(Housing = c(Clean = colors5[2], Dirty = colors5[5]))

# Generate visualizations for each tissue ----
for (tis in dig_tissues) {
  # Prepare metadata for plotting
  temp_meta <- data.frame(meta_ps_exp1) |>
    dplyr::filter(Tissue == tis) |>
    dplyr::filter(Infection == "Uninfected") |>
    dplyr::arrange(group)

  temp_meta$Housing <- gsub(" ", "_", temp_meta$Housing)

  # Prepare annotation data frame
  temp_ann_df <- data.frame(temp_meta[, "Housing"])
  colnames(temp_ann_df)[1] <- "Housing"
  rownames(temp_ann_df) <- rownames(temp_meta)

  # Order taxa by effect size
  temp_order <- m_housing_results_sig[[tis]] |>
    dplyr::arrange(coef)

  # Prepare abundance data for heatmap
  temp_df <- gen_dig_rel[temp_order$Genus, temp_meta$Sample]
  temp_df[temp_df == 0] <- NA
  temp_df_log10 <- log10(temp_df)

  # Generate heatmap
  # Produces Supplementary figures 1A,B,C
  temp_plot <- pheatmap::pheatmap(temp_df_log10,
    scale = "none",
    annotation_col = temp_ann_df,
    na_col = "grey88",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_colors = annoCol,
    cellwidth = 20,
    cellheight = 20,
    fontfamily = "ArialMT",
    color = ygb_pal(20),
    show_colnames = FALSE
  )

  # Generate tile plot of mean abundances
  temp <- tables_for_plots_all_log10[[tis]] |>
    dplyr::filter(taxa %in% temp_order$Genus)

  temp_sub <- temp |>
    dplyr::select(taxa, Sample, group, value, Housing)

  temp_sub_sum <- temp_sub |>
    group_by(taxa, group, Housing) |>
    summarize(mean_val = mean(value)) |>
    dplyr::arrange(Housing, mean_val)

  # Create and save tile plot
  # Produces Figures 1G,H,I
  temp_sub_sum$taxa <- factor(temp_sub_sum$taxa, levels = unique(temp_sub_sum$taxa))
  ggplot(temp_sub_sum, aes(y = taxa, x = group, fill = mean_val)) +
    geom_tile() +
    scale_fill_gradientn(colors = ygb_pal(20)) +
    coord_fixed() +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    labs(y = "", x = "", title = tis, fill = "Log10(RA)")

  # Clean up temporary variables
  rm(list = ls(pattern = "^temp"))
}
