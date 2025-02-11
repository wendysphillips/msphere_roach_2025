library(tidyverse)
library(RColorBrewer)
library(openxlsx2)
library(phyloseq)
library(ggpubr)
library(microViz)
library(Maaslin2)
library(vegan)

colors5 <- c("#2D7696", "#59B3D9", "#800C0A", "#BA2B29" ,"#E06967" )
colorblind_safe <- c("#F99634", "#BD1300", "#720494", "#15768A", "#4AC368")
exp_cols <- c("#3455a3","#FFB90F")

my_colors_random <- c("#317AE0", "#E09CB4F0", "#A80C0CE4", "#98A31A", "#A27EA8",  "#0D51B8","#E0C1E6", "#8198DE", "#51D686", "#4E6DD4", "#7F348A", "#F02B91E7", "#A35710", "#5B45B5", "#856CEB", "#F2B118", "#C6E3BF", "#D55BEB", "#6B304C", "#22B7BF", "#AEB8D6", "#9C7206", "#F9FFB2",  "#7BAAE8", "#347199", "#0B36D1", "#C4CC68", "#AD7ECC", "#DA96E6", "#F58989", "#0DD6B8", "#6D7511", "#9E6E85", "#D6B058", "#E68C43", "#08A18A", "#B6E1E3", "#8AEDB1", "#6F9E3C", "#2A0CA3",  "#EBD8AE", "#7D590A", "#36A864", "#8824C7", "#62CBD1", "#AF0CCC", "#F54949")
cage_cols <- c("#CC6677", "#332288", "#117733", "#88CCEE", "#882255", "#44AA99", "#999933", "#AA4499")
cage_cols2 <- c("#F99634", "#BD1300", "#720494", "#15768A", "#4AC368",  "#A27EA8", "#8198DE", "#A35710")
cols <- brewer.pal(9, "YlGnBu")
ygb_pal <- colorRampPalette(cols)

theme_set(theme_bw() + theme(
  legend.position = "right",
  axis.text = element_text(family = "ArialMT", color = "black"),
  legend.title = element_text(family = "ArialMT", color = "black"),
  legend.text = element_text(family = "ArialMT", color = "black"),
  text = element_text(family = "ArialMT")
))