# Produces Figures 2H,I

# Data set up ------
# Initialize and clean taxa and metadata
taxa_with_pet <- taxa_si
taxa_no_pet <- taxa_with_pet[, rownames(meta_no_pet)]

meta_no_pet$Housing_Infection <- as.character(meta_no_pet$Housing_Infection)
meta_no_pet$Housing_Infection <- gsub(" ", "_", meta_no_pet$Housing_Infection)

# Filter out low abundance, low frequency taxa ------
# Compute relative abundances; retain taxa with sufficient prevalence
taxa_no_pet_rel <- sweep(taxa_no_pet, 2, colSums(taxa_no_pet), "/")
taxa_no_pet_rel2 <- taxa_no_pet_rel[apply(taxa_no_pet_rel, 1, (function(x) sum(x > 0.001) >= 4)), ]

# Remove empty taxa ------
# Keep only rows with nonzero total abundance
taxa_no_pet_rel_filt <- taxa_no_pet_rel2[rowSums(taxa_no_pet_rel2) > 0, ]

# Maaslin infection analysis -------
# Loop over housing status (Clean and Dirty) and run Maaslin2 for infection-differential analysis
expt2_results_all_list <- list()
expt2_results_all_full <- list()

for (st in c("Clean", "Dirty")) {
  # Subset to samples for given housing status
  temp_df <- taxa_no_pet
  temp_meta <- meta_no_pet[meta_no_pet$Housing == st, ]
  temp_df <- temp_df[, rownames(temp_meta)]

  # Filter out low abundance, low frequency taxa for current group
  temp_rel <- sweep(temp_df, 2, colSums(temp_df), "/")
  temp_rel2 <- temp_rel[apply(temp_rel, 1, (function(x) sum(x > 0.001) >= 4)), ]

  # Set reference values depending on housing status
  if (st == "Clean") {
    rv <- "Cage,C1_U"
  } else {
    rv <- "Cage,D1_B"
  }

  # Run Maaslin2 analysis using infection and cage as fixed effects
  temp_fit_data <- Maaslin2(
    t(temp_rel2),
    temp_meta,
    paste0("maaslin2_expt2_", st),
    fixed_effect = c("Infection", "Cage"),
    reference = c("Infection,Uninfected", rv),
    standardize = F,
    min_prevalence = 0,
    normalization = "NONE"
  )

  temp_results <- temp_fit_data[["results"]]
  temp_results <- left_join(temp_results, name_map)
  expt2_results_all_list[[st]] <- temp_results
  expt2_results_all_full[[st]] <- temp_fit_data
}

# Compile all results -------
# Combine and reorder results from both housing conditions
expt2_results_all <- data.table::rbindlist(
  expt2_results_all_list,
  use.names = TRUE,
  fill = TRUE,
  idcol = "Housing"
)
expt2_results_all_tax <- left_join(
  expt2_results_all,
  taxa_si_names |> select(-Genus),
  by = "feature"
)
expt2_results_all_tax <- expt2_results_all_tax |>
  dplyr::select(
    Housing, metadata, feature, Kingdom, Phylum,
    Class, Order, Family, Genus, everything()
  ) |>
  dplyr::arrange(
    Housing, metadata, Kingdom, Phylum,
    Class, Order, Family, Genus
  )

# Analyze infection parameter ------
# Subset results to only those associated with Infection
expt2_results_infection <- expt2_results_all |>
  dplyr::filter(metadata == "Infection")

# Set significance thresholds ------
q_lim <- 0.05
abs_lim <- 1

# Identify significant features across all parameters ------
expt2_results_all_sig <- expt2_results_all |>
  dplyr::filter(qval <= q_lim) |>
  dplyr::filter(abs(coef) > abs_lim)

# Summary counts by parameter ------
expt2_results_all_sig |>
  group_by(metadata) |>
  count()

# Focus on significant infection features ------
expt2_results_infection_sig <- expt2_results_all_sig |>
  dplyr::filter(metadata == "Infection")

# Count taxa with lower/higher coefficients and total ------
expt2_results_infection_sig |>
  group_by(Housing) |>
  filter(coef < 0) |>
  count()
expt2_results_infection_sig |>
  group_by(Housing) |>
  filter(coef > 0) |>
  count()
expt2_results_infection_sig |>
  group_by(Housing) |>
  count()

# Identify taxa significant in both housing conditions ------
in_both <- expt2_results_infection_sig |>
  group_by(Genus) |>
  filter(n() > 1) |>
  ungroup() |>
  arrange(Housing, coef)

# Identify taxa significant in only one housing condition ------
in_one <- expt2_results_infection_sig |>
  group_by(Genus) |>
  filter(n() == 1) |>
  ungroup() |>
  arrange(Housing, coef)

# Create vector for plotting order ------
sig_infection_taxa_vector <- c(unique(in_both$Genus), in_one$Genus)

# Assemble all test results for infection taxa ------
expt2_results_infection_all_sig_taxa <- expt2_results_infection |>
  dplyr::filter(Genus %in% sig_infection_taxa_vector) |>
  dplyr::filter(metadata == "Infection")

# Fill missing taxa with NA and set non-significant q-values ------
expt2_results_infection_all_sig_taxa <- expt2_results_infection_all_sig_taxa |> complete(Housing, Genus)
expt2_results_infection_all_sig_taxa$qval[is.na(expt2_results_infection_all_sig_taxa$qval)] <- 2
expt2_results_infection_all_sig_taxa$sig <- ifelse(
  expt2_results_infection_all_sig_taxa$qval < 0.05,
  "yes",
  "no"
)

# Plotting preparations ------
# Reorder factors for plotting
expt2_results_infection_all_sig_taxa$Genus <- factor(
  expt2_results_infection_all_sig_taxa$Genus,
  levels = rev(sig_infection_taxa_vector)
)
expt2_results_infection_all_sig_taxa$Housing <- factor(
  expt2_results_infection_all_sig_taxa$Housing,
  levels = c("Clean", "Dirty")
)

# Produces Figure 2H ----
# Generate bar plot for infection coefficients by taxa and housing
ggplot(
  expt2_results_infection_all_sig_taxa,
  aes(x = Genus, y = coef, group = rev(Housing), fill = Housing, color = sig)
) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.75) +
  scale_fill_manual(values = colors5[c(2, 4)], name = "") +
  scale_color_manual(values = c("white", "black")) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "", y = "Maaslin2 coefficient") +
  guides(color = "none") +
  scale_x_discrete(limits = rev)

# Interaction plot preparation -----
# Subset taxa for interaction plot based on significant taxa in one condition
temp_int_df <- taxa_no_pet_rel_filt[rownames(taxa_no_pet_rel_filt) %in% in_one$Genus, ]

# Replace zero values with a low number for log transformation
temp_int_df[temp_int_df == 0] <- 0.00001

# Log10 transform the abundance data for plotting
temp_sig_df <- log10(temp_int_df)

# Convert data from wide to long format for ggplot
temp_sig_df <- temp_sig_df |> rownames_to_column("taxa")
temp_sig_df_long <- reshape2::melt(temp_sig_df)
temp_sig_df_long <- temp_sig_df_long |> dplyr::rename(sample = variable)

# Merge with metadata and create grouping variable ------
temp_sig_df_long <- left_join(
  temp_sig_df_long,
  meta_no_pet |> rownames_to_column("sample"),
  by = "sample"
)
temp_sig_df_long$grouping <- paste0(
  temp_sig_df_long$Experiment, "_",
  temp_sig_df_long$Cage, "_",
  temp_sig_df_long$Housing
)
temp_sig_df_long$Housing_Infection <- factor(
  temp_sig_df_long$Housing_Infection,
  levels = c("Clean_Uninfected", "Clean_Infected", "Dirty_Uninfected", "Dirty_Infected")
)

# Add stats -----
# Conduct Wilcoxon tests for each taxa and record results
stat_results <- data.frame(
  matrix(
    nrow = 0,
    ncol = 10,
    dimnames = list(
      NULL,
      c(
        ".y.", "group1", "group2", "n1", "n2",
        "statistic", "p", "p.adj", "taxa", "qval"
      )
    )
  )
)

for (st in c("Clean", "Dirty")) {
  print(st)
  temp <- expt2_results_infection[expt2_results_infection$Housing == st, ] |>
    filter(metadata == "Infection")
  temp <- temp[temp$Genus %in% in_one$Genus[in_one$Housing == st], ]
  temp_test <- temp_sig_df_long |>
    dplyr::filter(Housing == st) |>
    dplyr::filter(taxa %in% temp$Genus)
  for (tx in unique(temp$Genus)) {
    print(tx)
    temp_test2 <- temp_test |>
      dplyr::filter(taxa == tx)
    try(stat_test <- temp_test2 %>%
      rstatix::wilcox_test(value ~ Infection), silent = T)
    stat_test$`p.adj` <- NA
    stat_test$group1 <- paste0(st, "_", stat_test$group1)
    stat_test$group2 <- paste0(st, "_", stat_test$group2)
    stat_test$taxa <- tx

    stat_test$qval <- expt2_results_infection$qval[
      (expt2_results_infection$Genus == tx) & (expt2_results_infection$Housing == st)
    ]
    stat_results <- rbind(stat_results, stat_test)
  }
}

# Assign significance symbols based on q-value thresholds
stat_results$qsymbol <- NA
stat_results$qsymbol[stat_results$qval > 0.05] <- "ns"
stat_results$qsymbol[stat_results$qval < 0.05] <- "*"
stat_results$qsymbol[stat_results$qval < 0.01] <- "**"
stat_results$qsymbol[stat_results$qval < 0.001] <- "***"

# Shorten long taxa names for better display in plots
stat_results$taxa <- gsub(
  "Burkholderia-Caballeronia-Paraburkholderia",
  "Burk-Cab-Paraburk",
  stat_results$taxa
)
temp_sig_df_long$taxa <- gsub(
  "Burkholderia-Caballeronia-Paraburkholderia",
  "Burk-Cab-Paraburk",
  temp_sig_df_long$taxa
)
sig_infection_taxa_vector_sub <- gsub(
  "Burkholderia-Caballeronia-Paraburkholderia",
  "Burk-Cab-Paraburk",
  sig_infection_taxa_vector
)
temp_sig_df_long$taxa <- factor(
  temp_sig_df_long$taxa,
  levels = sig_infection_taxa_vector_sub
)

# Produces Figure 2I ----
# Generate boxplot with jitter and annotated significance symbols for the interaction plot
ggplot(temp_sig_df_long, aes(x = Housing_Infection, y = value)) +
  geom_boxplot(outlier.colour = NA, color = "grey", fatten = 2) +
  geom_jitter(aes(color = Experiment, shape = Housing_Infection),
    width = 0.2, height = 0,
    alpha = 0.8, size = 2.5, stroke = 2
  ) +
  scale_color_manual(values = exp_cols) +
  scale_shape_manual(values = c(1, 16, 0, 15, 2)) +
  guides(shape = "none") +
  labs(y = "Log10(Relative read abundance)", x = "", color = "") +
  scale_y_continuous(limits = c(-5.5, 0.5), breaks = c(-5, -4, -3, -2, -1)) +
  scale_x_discrete(labels = c("CU", "CI", "DU", "DI")) +
  stat_pvalue_manual(stat_results, label = "qsymbol", y.position = 0.15) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text(family = "ArialMT", color = "black"),
    strip.text = element_text(family = "ArialMT", color = "black", size = 9),
    text = element_text(family = "ArialMT"),
    legend.position = "top"
  ) +
  facet_wrap(~taxa, nrow = 3)
