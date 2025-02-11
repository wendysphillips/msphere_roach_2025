# For Supplementary figures 1D & 1E

# Define color scheme for housing and infection status visualization
annoCol <- list(
  Housing = c(
    "Clean" = colors5[2],
    "Dirty" = colors5[5]
  ),
  Infection = c(Infected = "grey33", Uninfected = "grey66")
)

# Loop through each tissue type
for (tis in dig_tissues) {
  # Filter significant taxa based on infection results
  temp <- tables_for_plots_all_log10[[tis]] |>
    dplyr::filter(taxa %in% m_infection_results_sig[[tis]]$Genus)

  # Select relevant columns for analysis
  temp_sub <- temp |>
    dplyr::select(taxa, Sample, group, value, Housing, Infection)

  # Calculate mean values for each taxa by group and infection status
  temp_sub_sum <- temp_sub |>
    group_by(taxa, group, Infection) |>
    summarize(mean_val = mean(value))

  # Sort data by group, infection status, and mean value
  temp_sub_sum <- temp_sub_sum |>
    dplyr::arrange(group, Infection, mean_val)

  # Join with infection results and process coefficients
  temp_sub_sum2 <- left_join(temp_sub_sum,
    m_infection_results_list[[tis]] |>
      rename(group = Group, taxa = feature),
    by = c("group", "taxa")
  )

  # Set coefficient to NA for uninfected samples
  temp_sub_sum2$coef[temp_sub_sum2$Infection == "Uninfected"] <- NA

  # Sort and format data for visualization
  temp_sub_sum2 <- temp_sub_sum2 |>
    dplyr::arrange(group, Infection, mean_val)
  temp_sub_sum2$taxa <- factor(temp_sub_sum2$taxa, levels = unique(temp_sub_sum2$taxa))

  # Generate detailed heatmap with pheatmap -----
  # Prepare data in wide format
  temp_wide <- temp |>
    pivot_wider(names_from = taxa, id_cols = Sample, values_from = value) |>
    column_to_rownames("Sample")

  # Filter and arrange metadata
  temp_meta <- meta_ps_exp1 |>
    dplyr::filter(Tissue == tis) |>
    dplyr::filter(Housing != "Pet store") |>
    dplyr::arrange(Housing, Infection)

  # Prepare data matrix and handle missing values
  temp_wide_plot <- temp_wide[temp_meta$Sample, rev(levels(temp_sub_sum2$taxa))]
  temp_wide_plot[temp_wide_plot == -5] <- NA
  
  # Create annotation data frame
  temp_ann_df <- data.frame(temp_meta[, "Infection"])
  colnames(temp_ann_df)[1] <- "Infection"
  temp_ann_df$Housing <- temp_meta$Housing
  rownames(temp_ann_df) <- rownames(temp_meta)

  # Generate pheatmap visualization
  # This produces Supplementary figures 1D & 1E
  p <- pheatmap::pheatmap(t(temp_wide_plot),
    scale = "none",
    na_col = "grey88",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    cellwidth = 20,
    cellheight = 20,
    fontfamily = "ArialMT",
    color = ygb_pal(20),
    show_colnames = FALSE,
    annotation_col = temp_ann_df,
    annotation_colors = annoCol
  )

  # Clean up temporary variables
  rm(list = ls(pattern = "^temp"))
}
print(p)
