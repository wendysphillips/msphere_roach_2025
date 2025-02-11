# Produces Figures 1C,D

# Create count data frame from phyloseq object
gen_dig_rar_table <- data.frame(t(phylo_gen_dig_rar@otu_table))

# Calculate alpha diversity metrics
alpha <- data.frame(Sample = colnames(gen_dig_rar_table))
alpha$Shannon <- vegan::diversity(gen_dig_rar_table, "shannon", MARGIN = 2)
alpha$Simpson <- vegan::diversity(gen_dig_rar_table, "simpson", MARGIN = 2)

# Join with metadata
alpha_dig <- left_join(alpha, meta_ps_exp1, by = "Sample")

# Filter datasets for different analyses
alpha_dig_np <- alpha_dig |> dplyr::filter(Housing != "Pet store")
alpha_housing <- alpha_dig |> dplyr::filter(Infection == "Uninfected")
alpha_housing_np <- alpha_housing |> dplyr::filter(Housing != "Pet store")

# Housing Effect Analysis ----

# Filter data for housing analysis
alpha_plot_data <- alpha_housing_np |>
  dplyr::filter(Infection == "Uninfected") |>
  dplyr::filter(Housing != "Pet store")

# Analyze Simpson diversity by housing
stat_test <- alpha_plot_data |>
  group_by(Tissue) %>%  
  rstatix::wilcox_test(Simpson ~ Housing)

# Adjust and round p-values
stat_test$padj <- p.adjust(stat_test$p, method = "fdr")
stat_test$padj <- round(stat_test$padj, 3)

# Generate plots for housing effect on Simpson diversity
# Produces Figure 1D
ggplot(alpha_plot_data, aes(x = Housing, y = Simpson, color = Tissue)) +
  geom_boxplot(outlier.colour = NA, fatten = 1, linewidth = 0.25) +
  geom_jitter(
    height = 0, width = 0.2, alpha = 0.8, stroke = 2, size = 2,
    aes(color = Tissue, shape = Housing)
  ) +
  facet_grid(~Tissue) +  # updated from ~tissue
  theme_pubr(base_size = 14) +
  stat_pvalue_manual(stat_test, y.position = rep(0.93, 3), label = "padj") +
  scale_color_manual(values = colors3) +
  scale_shape_manual(values = c(1, 0)) +
  labs(shape = "") +
  theme(legend.position = "none")

# Analyze Shannon diversity by housing
stat_test <- alpha_plot_data |>
  group_by(Tissue) %>%  
  rstatix::wilcox_test(Shannon ~ Housing)

# Adjust and round p-values
stat_test$padj <- p.adjust(stat_test$p, method = "fdr")
stat_test$padj <- round(stat_test$padj, 3)

# Generate plots for housing effect on Shannon diversity
# Produces Figure 1C
ggplot(alpha_plot_data, aes(x = Housing, y = Shannon, color = Tissue)) +
  geom_boxplot(outlier.colour = NA, fatten = 1, linewidth = 0.25) +
  geom_jitter(
    height = 0, width = 0.2, alpha = 0.8, stroke = 2, size = 2,
    aes(color = Tissue, shape = Housing)
  ) +
  facet_grid(~Tissue) +  # updated from ~tissue
  theme_pubr(base_size = 14) +
  stat_pvalue_manual(stat_test, y.position = rep(2.95, 3), label = "padj") +
  scale_color_manual(values = colors3) +
  scale_shape_manual(values = c(1, 0)) +
  labs(shape = "") +
  scale_y_continuous(limits = c(0.9, 3.1)) +
  theme(legend.position = "none")

rm(stat_test)
