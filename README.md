# chronosequences

Using genomic and chemical data from freshwater chronosequences in Arctic Svalbard and Alpine Norway to study the functional and taxonomic succession of the microbiome upon glacial retreat. This repository includes slurm scripts to run Squeezemeta https://github.com/jtamames/SqueezeMeta on our system (https://documentation.sigma2.no/hpc_machines/saga.html) for metagenomic analyses (raw data processing, assmebly and annotation, and binning) as well as amplicon analysis using dada2 (https://benjjneb.github.io/dada2). Additional scripts are to perform statistics and visualization of the data in R.  

## pull repository

```
cd path to repositories
git clone https://gitlab.com/eiler_lab/chronosequences.git
```

## Authors and acknowledgment
Scripts were written by Alexander Eiler, Laurent Fontaine and Jing Wei.

## License
This Code is subject to the terms of the MIT License. 

## Project status
Results from this project have been submitted to a peer-reviewed scientifc journal.

## Folders and code
The analyses code is divided into two folders "metags" and "rRNA_amplicons" representing the code to analyze whole shotgun metagenomic and SSU rRNA gene amplicon data, respectively.

### metadata.R
This represents the R code to perform the statistical analysis and data visualization on the metadata such as nutrient and gas concentrations along the chronosequences.

## released version 1.0.0
https://zenodo.org/badge/latestdoi/588562152
