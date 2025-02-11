# Produces Figures 2E,F

# Calculate alpha diversity metrics using Shannon and Simpson indices
alpha_si <- estimate_richness(ps_si_rar_gen, measures = c("Shannon", "Simpson"))

# Merge metadata with alpha diversity metrics using sample identifiers
info_si <- left_join(
  meta_si |> rownames_to_column("sample"),
  alpha_si |> rownames_to_column("sample"),
  by = "sample"
)

# Filter the data for Genotype 'B6'
info_np <- info_si |> dplyr::filter(Genotype == "B6")

# Shannon housing_infection analysis:
# Run Dunn's test to compare Shannon diversity across Housing_Infection groups
stat_test <- info_np |>
  rstatix::dunn_test(Shannon ~ Housing_Infection)

# Filter the test comparisons to retain only those of interest
stat_test_shannon <- stat_test |> dplyr::filter(
  ((stat_test$group1 == "Clean Uninfected") & (stat_test$group2 == "Clean Infected")) |
    ((stat_test$group1 == "Clean Uninfected") & (stat_test$group2 == "Dirty Uninfected")) |
    ((stat_test$group1 == "Dirty Uninfected") & (stat_test$group2 == "Dirty Infected")) |
    ((stat_test$group1 == "Clean Infected") & (stat_test$group2 == "Dirty Infected"))
)

# Adjust the p-values using the Benjamini-Hochberg method and round q-values
stat_test_shannon$q <- p.adjust(stat_test_shannon$p, method = "BH")
stat_test_shannon$q <- round(stat_test_shannon$q, 3)

# Figure 2E:
# Plot Shannon diversity comparing Housing_Infection groups with jitter for individual points
ggplot(info_np, aes(y = Shannon, x = Housing_Infection)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(
    aes(color = Experiment, shape = Housing_Infection),
    width = 0.2,
    height = 0,
    alpha = 0.8,
    size = 2.5,
    stroke = 2
  ) +
  scale_color_manual(values = exp_cols) +
  scale_shape_manual(values = c(1, 16, 0, 15, 2)) +
  labs(x = "Infection status", color = "") +
  theme_pubr(base_size = 14) +
  scale_x_discrete(labels = c("CU", "CI", "DU", "DI")) +
  guides(shape = "none") +
  stat_pvalue_manual(
    stat_test_shannon,
    y.position = seq(2, 3.2, 0.4),
    label = "q"
  )

# Simpson housing_infection analysis:
# Run Dunn's test for Simpson diversity across Housing_Infection groups
stat_test <- info_np |>
  rstatix::dunn_test(Simpson ~ Housing_Infection)

# Filter the test results to include only selected pairwise comparisons
stat_test_simpson <- stat_test |> dplyr::filter(
  ((stat_test$group1 == "Clean Uninfected") & (stat_test$group2 == "Clean Infected")) |
    ((stat_test$group1 == "Clean Uninfected") & (stat_test$group2 == "Dirty Uninfected")) |
    ((stat_test$group1 == "Dirty Uninfected") & (stat_test$group2 == "Dirty Infected")) |
    ((stat_test$group1 == "Clean Infected") & (stat_test$group2 == "Dirty Infected"))
)

# Adjust p-values and round q-values for Simpson tests
stat_test_simpson$q <- p.adjust(stat_test_simpson$p, method = "BH")
stat_test_simpson$q <- round(stat_test_simpson$q, 3)

# Figure 2F:
# Plot Simpson diversity with customized aesthetics and group labels
ggplot(info_np, aes(y = Simpson, x = Housing_Infection)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(
    aes(color = Experiment, shape = Housing_Infection),
    width = 0.2,
    height = 0,
    alpha = 0.8,
    size = 2.5,
    stroke = 2
  ) +
  scale_color_manual(values = exp_cols) +
  scale_shape_manual(values = c(1, 16, 0, 15, 2)) +
  labs(x = "Infection status", color = "") +
  theme_pubr(base_size = 14) +
  scale_x_discrete(labels = c("CU", "CI", "DU", "DI")) +
  guides(shape = "none") +
  stat_pvalue_manual(
    stat_test_simpson,
    y.position = c(0.75, 0.85, 0.95, 1.05),
    label = "q"
  )
