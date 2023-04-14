# rRNA_amplicons
These are scripts used to analyze metabarcoding data (16S and 18S rRNA gene amplicon sequencing) using cutadapt and dada2. There is also R codes to process the metabarcoding data and relate ASV and taxonomic tables with metadata including greenhouse gas concentrations.

## individual scripts
### bacteria_cutadapt_Svalbard_Finse.slurm 
Slurm script to run cutadapt on the raw demultiplexed fastq files from bacterial 16S rRNA gene amplicons.

### bacteria_dada2_Svalbard.R 
This is an R script that runs dada2 on bacterial 16S rRNA gene amplicons. Outputs are ASV tables and taxonomy files.

### bacteria_run_dada2_Svalbard.slurm
Slurm script to run dada2 on the 16S rRNA gene amplicons.

### eukarya_cutadapt_Svalbard_Finse.slurm 
Slurm script to run cutadapt on the raw demultiplexed fastq files from eukaryotic 18S rRNA gene amplicons.

### eukarya_dada2_Svalbard.R 
This is an R script that runs dada2 on bacterial 18S rRNA gene amplicons. Outputs are ASV tables and taxonomy files.

### stats_final.R
This contains the statistical analysis and data visualizatio on the outputs from dada2 (ASV tables and taxonomy filess) both from 16S and 18S gene amplicons. Here we analyzed alpha and beta diversity across the chronosequences as well as individual ASVs and taxonomic groups.


