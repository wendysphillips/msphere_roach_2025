# Produces Figure 2G

# Function to prepare data for stacked bar plot
prep_for_stacked_alph <- function(table = small_rel, p = 0.05) {
  # Find rows where maximum value exceeds threshold p
  rmx <- matrixStats::rowMaxs(base::as.matrix(table, rownames.force = TRUE)) > p
  # Subset table to rows exceeding threshold
  x <- table[rmx, ]
  # Arrange rows by descending row sums
  x <- x %>% dplyr::arrange(desc(rowSums(x)))
  # Store row names (taxa names)
  z <- rownames(x)
  # Transpose matrix and convert to data frame
  df <- data.frame(base::t(x))
  colnames(df) <- z
  # Calculate row sums for abundance
  RS <- rowSums(df)
  # Add "Other" category for remaining abundance
  df$Other <- 1 - RS
  # If total abundance is less than 100%
  if (sum(RS) < (length(RS) * 1)) {
    # Add "Other" category
    df$Other <- 1 - RS
    # Convert row names to Sample column
    df <- df %>% tibble::rownames_to_column("Sample")
    # Reshape data from wide to long format
    df_melt <- reshape2::melt(df, by = c(Sample))
    # Rename columns
    colnames(df_melt) <- c("Sample", "Taxa", "Abundance")
    # Sort taxa names
    zs <- sort(z)
    # Set factor levels including "Other"
    df_melt$Taxa <- factor(df_melt$Taxa, levels = c(zs, "Other"))
  } else {
    # If abundance equals 100%, follow similar steps without "Other"
    df <- df %>% tibble::rownames_to_column("Sample")
    df_melt <- reshape2::melt(df, by = c(Sample))
    colnames(df_melt) <- c("Sample", "Taxa", "Abundance")
    zs <- sort(z)
    df_melt$Taxa <- factor(df_melt$Taxa, levels = zs)
  }
  # Return the melted data frame
  df_melt
}

# Subset to taxa > 3% in a sample
temp <- prep_for_stacked_alph(table = small_rel, p = 0.0305)

# Modify Housing_Infection name
meta_no_pet$Housing_Infection <- as.character(meta_no_pet$Housing_Infection)
meta_no_pet$Housing_Infection <- gsub(" ", "_", meta_no_pet$Housing_Infection)

# Join abundances with metadata
temp <- left_join(temp, meta_no_pet |> rownames_to_column("Sample"))
temp$Housing_Infection <- factor(temp$Housing_Infection,
  levels = c("Clean_Uninfected", "Dirty_Uninfected","Clean_Infected", "Dirty_Infected")
)

# Assign colors to taxa in the data subset
stack_colors <- my_colors_random
names(stack_colors) <- sort(unique(temp$Taxa))
stack_colors[["Other"]] <- "grey55"

# Produces Figure 2G
g <- ggplot2::ggplot(
  temp,
  aes(
    x = Sample,
    y = Abundance,
    fill = forcats::fct_rev(Taxa)
  )
) +
  # Base plot elements
  geom_bar(
    stat = "identity",
    position = "stack",
    color = "grey33",
    linewidth = 0.25
  ) +
  # Scales
  scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
  scale_fill_manual(values = stack_colors) +
  # Faceting
  facet_grid(~ Housing_Infection + Experiment,
    scales = "free_x",
    space = "free_x"
  ) +
  # Labels
  labs(
    y = "Relative abundance",
    x = ""
  ) +
  # Theme customization
  theme_pubr() +
  theme(
    # Legend formatting
    legend.direction = "vertical",
    legend.position = "right",
    legend.key.size = unit(0.1, "cm"),
    legend.key.width = unit(.5, "cm"),
    legend.text = element_text(size = 10),
    # Axis formatting
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16)
  ) +
  # Legend guides
  guides(fill = guide_legend(ncol = 1, title = "Taxa"))
print(g)
