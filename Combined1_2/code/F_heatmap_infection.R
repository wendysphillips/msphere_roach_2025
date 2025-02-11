# Produces Supplemental figure 1G

# Plot all significant taxa----
feat_to_plot <- sig_infection_taxa_vector

# Arrange samples and remove pet store
meta_si_sort <- meta_si |>
  dplyr::arrange(Housing, Infection, Experiment) |>
  filter(Genotype == "B6")

# Create annotation information df
temp_ann_df <- data.frame(meta_si_sort[, "Infection"])
colnames(temp_ann_df)[1] <- "Infection"
rownames(temp_ann_df) <- rownames(meta_si_sort)
temp_ann_df$Housing <- meta_si_sort$Housing
temp_ann_df$Experiment <- meta_si_sort$Experiment
temp_ann_df$Infection <- factor(temp_ann_df$Infection, levels = c("Uninfected", "Infected"))
# Set annotation colors
annoCol <- list(
  Infection = c("Infected" = "grey33", "Uninfected" = "grey66"),
  Housing = c("Clean" = colors5[2], "Dirty" = colors5[4]),
  Experiment = c("Experiment 1" = exp_cols[1], "Experiment 2" = exp_cols[2])
)

# Log10 transform counts -----
# Relative abundance no pet
temp_df <- taxa_si
temp_df <- temp_df[, rownames(meta_no_pet)]
small_rel <- sweep(temp_df, 2, colSums(temp_df), "/")
small_rel_log10 <- small_rel
small_rel_log10[small_rel_log10 == 0] <- NA
small_rel_log10 <- log10(small_rel_log10)

# Subset data to plot
temp_to_plot <- small_rel_log10[feat_to_plot, rownames(meta_si_sort)]

# Produces Supplemental figure 1G ----
temp_plot <- pheatmap::pheatmap(temp_to_plot,
  scale = "none",
  na_col = "grey92",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  cellwidth = 20,
  cellheight = 20,
  fontfamily = "ArialMT",
  color = ygb_pal(20),
  annotation_col = temp_ann_df,
  annotation_colors = annoCol,
  show_colnames = FALSE
)