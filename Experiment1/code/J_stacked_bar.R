# Produces Figure 1F

# Function to prepare data for stacked bar plot
# Filters taxa based on abundance threshold and restructures data
# Orders taxa alphabetically
# Also handles the "Other" category differently
prep_for_stacked_alph <- function(table = gen_dig_rel, p = 0.1) {
  # Filter rows where maximum abundance exceeds threshold p
  rmx <- matrixStats::rowMaxs(base::as.matrix(table, rownames.force = TRUE)) > p
  # Subset table to keep only taxa passing threshold
  x <- table[rmx, ]
  # Sort taxa by their total abundance across samples (descending)
  x <- x |> dplyr::arrange(desc(rowSums(x)))
  # Store taxa names
  z <- rownames(x)
  # Transpose matrix and convert to data frame (samples as rows, taxa as columns)
  df <- data.frame(base::t(x))
  # Set column names to taxa names
  colnames(df) <- z
  # Calculate row sums (total abundance per sample)
  RS <- rowSums(df)
  # Initialize "Other" category for remaining abundance
  df$Other <- 1 - RS
  # If total abundance is less than 100%
  if (sum(RS) < (length(RS) * 1)) {
    # Calculate remaining abundance for "Other" category
    df$Other <- 1 - RS
    # Convert row names to a column named "Sample"
    df <- df |> tibble::rownames_to_column("Sample")
    # Reshape data from wide to long format
    df_melt <- reshape2::melt(df, by = c(Sample))
    # Rename columns for clarity
    colnames(df_melt) <- c("Sample", "Taxa", "Abundance")
    # Sort taxa alphabetically
    zs <- sort(z)
    # Set factor levels with "Other" at the end
    df_melt$Taxa <- factor(df_melt$Taxa, levels = c(zs, "Other"))
  } else {
    # If total abundance is 100%, process without "Other" category
    df <- df |> tibble::rownames_to_column("Sample")
    df_melt <- reshape2::melt(df, by = c(Sample))
    colnames(df_melt) <- c("Sample", "Taxa", "Abundance")
    zs <- sort(z)
    # Set factor levels without "Other"
    df_melt$Taxa <- factor(df_melt$Taxa, levels = zs)
  }
  # Return processed data frame
  df_melt
}

# Create color mapping for taxa
# Generate temporary dataset to get qualifying taxa
temp <- prep_for_stacked_alph(table = gen_dig_rel, p = 0.05)

# Assign colors to taxa names
# Grey color is assigned to "Other" category
colors_sb <- my_colors_random
names(colors_sb) <- sort(unique(temp$Taxa))
colors_sb[["Other"]] <- "grey75"

rm(temp)

# Define plot appearance parameters
font_size <- 10
legend_position <- "right"
label <- "Taxa"

# Generate stacked bar plot
temp1 <- data.frame(meta_ps_exp1) |>
  dplyr::filter(Infection == "Uninfected")

# Prepare data for plotting
temp2 <- prep_for_stacked_alph(gen_dig_rel[, temp1$Sample], p = 0.05)
temp2 <- left_join(temp2, temp1, by = "Sample")
temp2$Housing <- factor(temp2$Housing, levels = c("Clean", "Pet store", "Dirty"))
temp2 <- temp2 |> dplyr:::arrange(Housing, group, )

# Create stacked bar plot
g <- ggplot2::ggplot(temp2, aes(x = Sample, y = Abundance, fill = forcats::fct_rev(Taxa))) +
  theme_pubr() +
  geom_bar(stat = "identity", position = "stack", color = "grey33", linewidth = 0.25) +
  # Set theme elements for consistent formatting
  theme(
    legend.key.size = unit(0.1, "cm"),
    legend.key.width = unit(.5, "cm"),
    legend.text = element_text(size = font_size),
    legend.direction = "vertical",
    legend.position = legend_position,
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16)
  ) +
  # Configure plot elements
  guides(fill = guide_legend(ncol = 1, title = label)) +
  scale_fill_manual(values = colors_sb) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
  labs(y = "Relative abundance", x = "") +
  facet_grid(~ Tissue + group, scales = "free", space = "free")
print(g)
