# README

This code processes 16S amplicon sequencing reads, creates taxa tables from QC'd reads, and analyzes the taxa tables to characterize bacterial communities in mouse experimental samples.

The first three files to run are coded in R Markdown.   

1. Reads_to_ASVs.Rmd
2. ASVs_to_taxa.Rmd
3. Phylogenies.Rmd


**Reads_to_ASVs.Rmd**   

*Inputs* 

- Working directory should be set to the directory containing the raw fastq files to be processed

*Outputs* 

- Reads with sequencing adapter trimmed off in fastq format, located in a directory named `trimmed/`
- Reads "denoised" and filtered by dada2 in the directory `trimmed/uncompressed_fastq/`
- An RData file with objects created during processing: `project16_cutadapt_dada2_out.RData`

**ASVs_to_taxa.Rmd**   

*Inputs* 

- `project16_cutadapt_dada2_out.RData` from above
- Silva 138.1 databases downloaded from https://zenodo.org/records/4587955

*Outputs* 

- `outputs/project16_taxa_dada_silva.csv`: a table with ASVs and their Silva taxonomic assignments, with read counts for each sample
- `outputs/qc_tracking_project16.csv`: a table with the number of reads surviving each QC step
- `RData_files/project16_cutadapt_dada2_taxa.RData`: contains all objects in the R environment

**Phylogenies.Rmd**  

*Inputs*

- `RData_files/project16_cutadapt_dada2_taxa.RData` from above
- `project16_meta.csv` found in `data/` directory

*Outputs*   

- `RData_files/project16_phyloseq_objects.RData`: an RData file with objects created during processing, most importantly the phyloseq objects, modified metadata, and taxatables used in downstream analysis


**`A_Process_ASVs.R`**

*Inputs*

- `RData_files/project16_phyloseq_objects.RData`: from above

*Outputs*   

- `project16_phyloseq.rds`: contains the phyloseq object created by the code