# Statistical analysis of alpha and beta diversity in relation to meta data



################## Work in progress ################################

## saga
#module load R/4.0.0-foss-2020a
#module load GDAL/3.0.4-foss-2020a-Python-3.8.2
#module load MariaDB-connector-c/3.1.7-GCCcore-9.3.0

# load(file=file.path(figsPath, "stats.RData"))


.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_0_0", .libPaths()))
rootPath  = file.path("/cluster/projects/nn9745k/jing/02_results/svalbard") #  =--------- Change to your path
metadataPath  = file.path("/cluster/projects/nn9745k/jing/scripts/metadata")

bacterialPath  = file.path(rootPath, "bacteria")
eukaryaPath  = file.path(rootPath, "eukarya")

figsPath  = file.path(rootPath, "figs")

##load packages

if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2)
   }
if (!require("reshape2")) {
   install.packages("reshape2", dependencies = TRUE)
   library(reshape2)
   }
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
if (!require("phyloseq")) {
   BiocManager::install("phyloseq")
   library(phyloseq)
   }
if (!require("devtools")) {
   install.packages("devtools", dependencies = TRUE)
   library(devtools)
   }
if (!require("mice")){
  install.packages("mice")
  library(mice)
  }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }
if (!require("mgcv")) {
  install.packages("mgcv", dependencies = TRUE)
  library(mgcv)
}
if (!require("nlme")) {
  install.packages("nlme", dependencies = TRUE)
  library(nlme)
}
if (!require("sp")) {
  install.packages("sp", dependencies = TRUE)
  library(sp)
}
if (!require("gratia")) {
  install.packages("gratia", dependencies = TRUE)
  library(gratia)
}

metadata

summary(metadata)
summary(metadata[which(metadata$gl_dist <= 350),])
summary(metadata[which(metadata$gl_dist >= 3350),])
apply(metadata[which(metadata$gl_dist <= 350),], 2, sd, na.rm = TRUE)
apply(metadata[which(metadata$gl_dist >= 3350),], 2, sd, na.rm = TRUE)

