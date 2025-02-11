# Produces Figures 2A,B,C

# Set taxonomic level for analysis
tax <- "Genus"

# Process beta diversity by tissue ----
# Filter out pet store samples
pseq <- phylo_gen_dig |>
  microViz::ps_filter(Housing != "Pet store")

# Generate plots for each tissue type
for (i in 1:3) {
  tis <- dig_tissues[i]

  # Calculate Bray-Curtis distances for each tissue
  temp_d <- pseq |>
    ps_filter(Tissue == tis) |>
    rarefy_even_depth(rngseed = 999, replace = F) |>
    tax_agg(tax) |>
    dist_calc("bray")

  # Produces Figures 2C
  # Generate basic PCoA plot for all tissues
  temp_d |>
    ord_calc("PCoA") |>
    ord_plot(color = colors3[i], shape = "Housing_Infection", size = 3, stroke = 2, alpha = 0.8) +
    scale_shape_manual(values = c(1, 16, 0, 15, 5), name = "Status") +
    theme_pubr(base_size = 14) +
    theme(legend.position = "right") +
    ggtitle(tis) +
    guides(shape = guide_legend(override.aes = list(linetype = 0), title = ""))

  # For only small and large, plot with ellipses
  if (i != 1) {
    # Produces Figures 2A,B
    # Generate PCoA plot with ellipses and legend
    temp_d |>
      ord_calc("PCoA") |>
      ord_plot(color = colors3[i], shape = "Housing_Infection", size = 3, stroke = 2, alpha = 0.8) +
      scale_shape_manual(values = c(1, 16, 0, 15, 5), name = "Status") +
      stat_ellipse(aes(linetype = Housing_Infection),
        color = colors3[i],
        show.legend = FALSE,
        level = 0.95
      ) +
      theme_pubr(base_size = 14) +
      theme(legend.position = "right") +
      ggtitle(tis)
  }
  rm(temp_d)
}
