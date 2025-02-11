# Produces Figure 2D and Supplemental figure 1F

# Prepare the phyloseq object:
# - Fix taxonomic names
# - Validate and remove undetected taxa
pseq <- ps_si %>% 
  tax_fix(unknowns = c(NA, "NA", "unclassified", "Unclassified", "Incertae Sedis")) %>% 
  phyloseq_validate(remove_undetected = TRUE)

# PCoA calculation:
# - Rarefy the dataset for even depth
# - Validate after rarefaction
# - Aggregate taxa at the Genus level
# - Compute Bray-Curtis distance and perform PCoA
p_pcoa <- pseq |>
  rarefy_even_depth(rngseed = 99, replace = FALSE) |>
  phyloseq_validate(remove_undetected = TRUE) |>
  tax_agg("Genus") |>
  dist_calc("bray") |>
  ord_calc("PCoA")

# Figure 2D:
# Create a PCoA plot using axes 1 and 2. Customize point aesthetics and legend.
p_plot_1_2 <- p_pcoa |>
  ord_plot(
    axes = c(1, 2),
    shape = "Housing_Infection",
    color = "Experiment",
    size = 3,
    stroke = 2,
    alpha = 0.8
  ) +
  theme_pubr(base_size = 16) +
  theme(
    legend.position = "right", 
    legend.box = "vertical", 
    legend.margin = margin(), 
    legend.text = element_text(size = 10)
  ) +
  guides(
    shape = guide_legend(ncol = 1), 
    color = guide_legend(override.aes = list(shape = 16))
  ) +
  scale_shape_manual(
    values = c(1, 16, 0, 15, 2),
    name = ""
  ) +
  scale_colour_manual(
    values = exp_cols,
    name = ""
  )

# Print the Figure 2D plot
print(p_plot_1_2)

# Supplemental Figure 1F:
# Create a PCoA plot using axes 1 and 3 with similar styling as Figure 2D.
p_plot <- p_pcoa |>
  ord_plot(
    axes = c(1, 3),
    shape = "Housing_Infection",
    color = "Experiment",
    size = 3,
    stroke = 2,
    alpha = 0.8
  ) +
  theme_pubr(base_size = 16) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(),
    legend.text = element_text(size = 10)
  ) +
  guides(
    shape = guide_legend(ncol = 1),
    color = guide_legend(override.aes = list(shape = 16))
  ) +
  scale_shape_manual(values = c(1, 16, 0, 15, 2), name = "") +
  scale_colour_manual(values = exp_cols, name = "")

# Print the Supplemental Figure 1F plot
print(p_plot)

# PERMANOVA Analysis:
# Filter the phyloseq object to include only Genotype 'B6', then process it similarly
temp_p <- pseq |>
  microViz::ps_filter(Genotype == "B6") |>
  rarefy_even_depth(rngseed = 99, replace = FALSE) |>
  phyloseq_validate(remove_undetected = TRUE) |>
  tax_agg("Genus")

# Compute Bray-Curtis distance from the aggregated OTU table
temp_d <- temp_p |>
  dist_calc("bray")

# Extract OTU table and sample metadata for PERMANOVA analysis
si_otu <- as.matrix(data.frame(t(temp_d@otu_table)))
si_meta <- data.frame(temp_d@sam_data)
si_BC <- vegdist(si_otu, method = "bray")

# Perform overall PERMANOVA for the Housing_Infection factor
adonis_housing_infection <- vegan::adonis2(si_BC ~ Housing_Infection, data = si_meta)

# Run pairwise PERMANOVA for Housing_Infection groups
adonis_housing_infection_pairwise <- pairwiseAdonis::pairwise.adonis2(si_BC ~ Housing_Infection, data = si_meta)

# Reformat the pairwise PERMANOVA results:
# Add a 'parameter' column to each result and combine them into one data frame.
adonis_housing_infection_pairwise <- lapply(
  adonis_housing_infection_pairwise[2:7],
  function(x) x |> rownames_to_column("parameter")
)
adonis_housing_infection_pairwise_df <- data.table::rbindlist(
  adonis_housing_infection_pairwise,
  use.names = TRUE,
  fill = TRUE,
  idcol = "comparison"
)
