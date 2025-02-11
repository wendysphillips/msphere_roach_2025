# Produces Figure 1E

# Process phyloseq object for combined tissue analysis ----
# Clean taxonomy and validate
pseq_all <- phylo_gen_dig_housing_rar %>%
  tax_fix(unknowns = "NA") %>%
  phyloseq_validate(remove_undetected = TRUE)

# Generate combined tissues PCoA plot ----
# Produces Figure 1E
p <- pseq_all |>
  tax_agg("Genus") |>
  dist_calc("bray") |>
  ord_calc("PCoA") |>
  ord_plot(
    color = "Tissue", 
    shape = "Housing",
    size = 3,
    alpha = 0.8,
    stroke = 2,
    axes = c(1, 3)
  ) +
  theme_pubr(base_size = 14) +
  theme(
    legend.text = element_text(size = 14),
    legend.position = "right"
  ) +
  scale_color_manual(values = c(colors3), name = "") +
  scale_shape_manual(values = c(1, 0, 2), name = "")
print(p)

# Initialize lists for statistical results
permanova_housing <- list()

# Process each tissue separately ----
for (tis in dig_tissues) {
  # Create tissue-specific subsets
  # Small intestine requires different filtering than cecum/large intestine
  if (tis == "Small") {
    # Filter housing samples for small intestine
    # Additional filtering for non-pet store mice
    pseq_np <- phylo_gen_dig_housing_rar |>
      ps_filter(Infection == "Uninfected") |>
      ps_filter(Housing != "Pet store") |>
      ps_filter(Tissue == tis)  
  } else {
    # Filter housing samples for cecum/large intestine
    # Additional filtering for non-pet store mice
    pseq_np <- phylo_gen_dig_cl_rar |>
      ps_filter(Infection == "Uninfected") |>
      ps_filter(Housing != "Pet store") |>
      ps_filter(Tissue == tis) 
  }

  # Calculate Bray-Curtis distances
  temp_d <- pseq_np |>
    tax_agg("Genus") |>
    dist_calc("bray")

  # Perform PERMANOVA test
  temp_p <- temp_d |>
    dist_permanova(
      seed = 999,
      variables = c("Housing"),
      n_processes = 1,
      n_perms = 9999
    )

  # Store PERMANOVA results
  permanova_housing[[tis]] <- temp_p@permanova

  # Clean up temporary variables
  rm(list = ls(pattern = "^temp"))
}

# Save statistical results ----
# Combine PERMANOVA results
housing_permanova_results <- data.table::rbindlist(permanova_housing,
  use.names = TRUE,
  fill = TRUE,
  idcol = "Tissue"
)