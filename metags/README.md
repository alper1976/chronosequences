# metags
These are scripts used to run the co-assembly and sequentially merged assembly modes of squeezemeta (slurm scripts) and R codes to process the outputs of squeezemeta for both the gene (stats_metags.R) and genome (stats_mags.R) centric approaches 

## individual scripts
### SqueezeMeta_conf.pl 
This contains the configurations of the pipeline.

### extract_metags.R 
This is an R script that merges outputs from the individual sample assemblies and annotation

### parameters.pl 
This containes the parameter used for running SqueezeMeta

### run_squeezemeta_coassembly.slurm
This is a slurm script specifying the metagenomic analysis to run SqueezeMeta in co-assembly mode.

### run_squeezemeta_seqmerge.slurm
This is a slurm script specifying the metagenomic analysis to run SqueezeMeta in indiviudal assembly mode and later merging the assemblies.

### stats_mags.R 
This contains the code for statistical analysis ofthe metagenomes assembled genomes (MAGs) from the co-assembly as well as data visualization. We refer to it as the genome centric analysis as annotations and coverage of high quality MAGs were analyzed.

### stats_metags.R
This contains the code for statistical analysis of the individual assembled metagenomes and data visualization. Here we analyzed annotations from the contigs. We refer to it as the gene centric analysis as gene annotations from MAGs were used in the analysis.


