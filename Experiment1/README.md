# README
## Analysis for Experiment 1 data

This code processes 16S amplicon sequencing reads, creates taxa tables from QC'd reads, and analyzes the taxa tables to characterize bacterial communities in mouse experimental samples.

Code used for the analyses and figures in the MSphere manuscript are in the directory `code/`. 

These four scripts must be run first in this order:

1. cutadapt.txt 
    - Requires 
        - the location of the directory containing the raw fastq files to be processed
2. A_dada2_to_phyloseq.R : produces **ps_exp1.Rds** file used for all downstream analysis
    - Requires 
        - the location of the directory containing the processed fastq files
        - `data/experiment1_meta.tsv` file containing metadata for the experiment
        - Silva 138.1 databases downloaded from https://zenodo.org/records/4587955      
3. B_data_prep.R
4. C_generate_log10_tables.R

The following scripts should run subsequent to the above steps without being interdependent unless noted.   

- D_alpha.R : produces Figures 1C,D
- E_beta_microviz_housing.R : produces Figure 1E
- F_beta_microviz_infection.R : produces Figures 2A,B,C
- G_maaslin2_housing.R : produces Figures 1G-J
- H_maaslin2_infection.R
- I_maaslin2_infection_plots.R (requires G_maaslin2_infection.R to be run first) : produces Supplementary figures 1D,E
- J_stacked_bar.R : : produces Figure 1F