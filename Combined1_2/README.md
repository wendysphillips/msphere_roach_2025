# README

This code joins phyloseq objects from Experiment 1 and Experiment 2. The merged phyloseq object is then subset to be only small intestine samples. Microbiome analytics are preformed on the combined small intestinal data.

These two scripts must be run first in this order.  

1. A_Start_here_analysis.R
2. B_join_experiments.R: requires as input `Experiment1/ps_exp1.Rds` and `Experiment2/ps_exp2.Rds`

After running the above, the following scripts can be run. They should run in any order independent of each other unless noted.   

- C_alpha.R
- D_beta_plots_bray.R
- E_differential_features.R
- F_heatmap_infection.R (run after E_differential_features.R)
- G_stacked_bar.R

## Figures Produced

- Figure 2D - [code/D_beta_plots_bray.R](./code/D_beta_plots_bray.R)
- Figure 2E - [code/C_alpha.R](./code/C_alpha.R)
- Figure 2F - [code/C_alpha.R](./code/C_alpha.R)
- Figure 2G - [code/G_stacked_bar.R](./code/G_stacked_bar.R)
- Figure 2H - [code/E_differential_features.R](./code/E_differential_features.R)
- Figure 2I - [code/E_differential_features.R](./code/E_differential_features.R)
- Supplemental Figure 1F - [code/D_beta_plots_bray.R](./code/D_beta_plots_bray.R)
- Supplemental Figure 1G - [code/F_heatmap_infection.R](./code/F_heatmap_infection.R)