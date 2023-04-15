### This is an R script to analyse to infer greenhouse gas emissions and metabarcoding
### data from Svalbard


## saga
#module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
#module load GDAL/3.0.4-foss-2020a-Python-3.8.2
#module load MariaDB-connector-c/3.1.7-GCCcore-9.3.0
path_scripts = "/cluster/projects/nn9745k/jing/scripts"
path_out = "/cluster/projects/nn9745k/jing/02_results/svalbard"
path_euk = file.path(path_out, "eukarya/dada2")
path_bac = file.path(path_out, "bacteria/dada2optimized")

figs_dir <- file.path(path_out, "figs")


# load(file=file.path(figs_dir, "stats.RData"))
# save.image(file=file.path(figs_dir, "stats.RData"))

#In R
.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_0_0", .libPaths()))

## LOAD PACKAGES

if (!require("phyloseq", quietly = TRUE)){
	BiocManager::install("phyloseq")
	library(phyloseq)
	}
if (!require("stringr", quietly = TRUE)){
	install.packages("stringr")
	library(stringr)
	}
if (!require("Hmisc")) {
   install.packages("Hmisc", dependencies = TRUE)
   library(Hmisc)
   }
if (!require("corrplot")) {
  install.packages("corrplot", dependencies = TRUE)
  library(corrplot)
  }
if (!require("R.utils")) {
   install.packages("R.utils", dependencies = TRUE)
   library(R.utils)
   }
if (!require("plsdepot")) {
   install.packages("plsdepot", dependencies = TRUE)
   library(plsdepot)
   }
if (!require("vegan")) {
   install.packages("vegan", dependencies = TRUE)
   library(vegan)
   }
if (!require("mgcv")) {
  install.packages("mgcv", dependencies = TRUE)
  library(mgcv)
}
# if (!require("mgcViz")) {
#   install.packages("mgcViz", dependencies = TRUE)
#   library(mgcViz)
# }
if (!require("nlme")) {
  install.packages("nlme", dependencies = TRUE)
  library(nlme)
}
if (!require("sp")) {
  install.packages("sp", dependencies = TRUE)
  library(sp)
}
if (!require("visreg")) {
  install.packages("visreg", dependencies = TRUE)
  library(visreg)
}
if (!require("mvabund")) {
  install.packages("mvabund", dependencies = TRUE)
  library(mvabund)
}
if (!require("gllvm")) {
  install.packages("gllvm", dependencies = TRUE)
  library(gllvm)
}
if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}
if (!require("pheatmap")) {
  install.packages("pheatmap", dependencies = TRUE)
  library(pheatmap)
}
if (!require("minerva")) {
   install.packages("minerva", dependencies = TRUE)
   library(minerva)
   }

library(biomformat)

## LOAD local scripts
sourceDirectory('/cluster/projects/nn9745k/scripts/amplicon_analysis_package')

## specific functions

extract.xyz <- function(obj) {
    xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
    xyz <- cbind(xy, c(obj$grid$z))
    names(xyz) <- c("x", "y", "z")
    return(xyz)
}

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_palette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

color_palette_euk = c("#D1BBD7", "#AE76A3", "#882E72", "#1965B0",
                        "#5289C7", "#7BAFDE", "#4EB265", "#90C987",
                        "#CAE0AB", "#F7F056", "#F6C141", "#F1932D",
                        "#E8601C", "#DC050C", "#72190E")
color_palette_bac = c("#114477", "#4477AA", "#77AADD", "#117755",
                          "#44AA88", "#99CCBB", "#777711", "#AAAA44",
                          "#DDDD77", "#771111", "#AA4444", "#DD7777",
                          "#771144", "#AA4477", "#DD77AA")




## LOAD VARIABLES
# path_arch
ASV_euk = as.matrix(read.table(file.path(path_euk,"ASV_table.tsv"), header=T, check.names=FALSE))
tax_euk = as.matrix(read.table(file.path(path_euk,"Taxonomy_table_pr2.tsv"), header=T, check.names=FALSE))

ASV_bac = as.matrix(read.table(file.path(path_bac,"ASV_table.tsv"), header=T, check.names=FALSE))
tax_bac = as.matrix(read.table(file.path(path_bac,"Taxonomy_table_silva.tsv"), header=T, check.names=FALSE))
ASV_no_taxo = readDNAStringSet(file.path(path_bac,"dada2/ASV_no_taxo.fasta"))


# load metadata ## NEEDS FIXING
cls <- c("sample_id" = "factor",
		 "sample_name" = "factor",
		 "sample_date" = "factor",
		 "locality" = "factor",
		 "altitude" = "numeric",
		 "north" = "numeric",
		 "east" = "numeric",
     "gl_dist" = "numeric",
     "Ccounts" = "numeric",
		 "Temp" = "numeric",
		 "Cond" = "numeric",
		 "pH" = "numeric",
		 "TOC" = "numeric",
		 "DOC" = "numeric",
		 "TN" = "numeric",
		 "DN" = "numeric",
		 "TP" = "numeric",
		 "DP" = "numeric",
     "O2" = "numeric",
     "N2" = "numeric",
     "CO2" = "numeric",
     "CH4" = "numeric",
     "N2O" = "numeric",
     "bird" = "numeric",
     "N2_sat" = "numeric",
     "O2_sat" = "numeric",
     "CO2_sat" = "numeric",
     "CH4_sat" = "numeric",
     "N2O_sat" = "numeric")

metadata = read.csv(file.path(path_scripts, "metadata", "svfi19.tsv"), header=T, check.names=FALSE, sep = "\t", stringsAsFactors=FALSE,colClasses=cls)

## Parse data
rownames(metadata) = metadata$sample_id
rownames(ASV_euk) = str_split(rownames(ASV_euk), pattern = "_", simplify = TRUE)[, 1]
rownames(ASV_bac) = str_split(rownames(ASV_bac), pattern = "_", simplify = TRUE)[, 1]


## MERGE INTO PHYLOSEQ OBJECT

# Eukaryota
dim(ASV_euk)
euk_otu <- otu_table(t(ASV_euk), taxa_are_rows = TRUE)
dim(euk_otu)
euk_tax <- tax_table(tax_euk)
euk_physeq <- phyloseq(euk_otu, euk_tax, sample_data(data.frame(metadata)))
euk_physeq

euk_physeq_clean <- subset_taxa(euk_physeq, Kingdom == "Eukaryota")
## remove metazoa??
metazoa_physeq_clean <- subset_taxa(euk_physeq, Division == "Metazoa")
eukn_physeq_clean <- subset_taxa(euk_physeq_clean, Division != "Metazoa")
arthropoda_physeq_clean <- subset_taxa(euk_physeq, Class == "Arthropoda")

# Bacteria
dim(ASV_bac)
bac_otu <- otu_table(t(ASV_bac), taxa_are_rows = TRUE)
dim(bac_otu)
bac_tax <- tax_table(tax_bac)
bac_physeq <- phyloseq(bac_otu, bac_tax, sample_data(data.frame(metadata)))
bac_physeq

bac_physeq_clean  = subset_taxa(bac_physeq, Kingdom == "Bacteria")
bac_physeq_clean  = subset_taxa(bac_physeq_clean, Order != "Chloroplast")
bac_physeq_clean  = subset_taxa(bac_physeq_clean, Family != "Mitochondria")


## RUN statistics and visualization
# summary statistics of metadata

sample_metadata = sample_data(euk_physeq_clean)[,c(7, 10, 12:15, 17:33)]

scaled_metadata = scale(sample_metadata )

corr_spearman = rcorr(as.matrix(scaled_metadata), type = "spearman")

# get Holm's adjusted p-value
corr_spearman_adjusted = print_rcorr_adjust(rcorr.adjust(scaled_metadata, type = "spearman",
                                              use = "pairwise.complete.obs"))

df_1 <- dplyr::mutate_all(as.data.frame(corr_spearman_adjusted, stringsAsFactors = FALSE), function(x) as.numeric(as.character(x)))
df_1[is.na(df_1)] <- 0.0001
df_2 <- as.matrix(df_1)

postscript(file.path(figs_dir, "correlation_matrix_metadata_adjusted.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE)
  corrplot::corrplot(corr_spearman$r, type="upper", p.mat = df_2, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()

postscript(file.path(figs_dir, "correlation_matrix_metadata.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE)
  corrplot::corrplot(corr_spearman$r, type="upper", p.mat = corr_spearman$P, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()

# rarefaction curves to check where to set cutoff for sequence reads
cairo_ps(file.path(figs_dir, "rarefaction_curves_euk.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  rarecurve(t(phyloseq::otu_table(euk_physeq_clean)),
                      step = 100, xlab = "Sample Size",
                      ylab = "ASVs",
                      label = FALSE,
                      cex.lab = 2,
                      cex.axis = 1.5)
dev.off()

cairo_ps(file.path(figs_dir, "rarefaction_curves_metaz.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  rarecurve(t(phyloseq::otu_table(metazoa_physeq_clean)),
                      step = 100, xlab = "Sample Size",
                      ylab = "ASVs",
                      label = FALSE,
                      cex.lab = 2,
                      cex.axis = 1.5)
dev.off()

cairo_ps(file.path(figs_dir, "rarefaction_curves_bac.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  rarecurve(t(phyloseq::otu_table(bac_physeq_clean)),
                      step = 100, xlab = "Sample Size",
                      ylab = "ASVs",
                      label = FALSE,
                      cex.lab = 2,
                      cex.axis = 1.5)
dev.off()

colSums(phyloseq::otu_table(euk_physeq_clean))
colSums(phyloseq::otu_table(metazoa_physeq_clean))
colSums(phyloseq::otu_table(bac_physeq_clean))


## based on rarefaction curves cutoffs for number of sequences were set.

# alpha diversity - richness and evenness
# Eukaryota
euk_physeq_20000 <- prune_samples(sample_sums(euk_physeq_clean)>=20000, euk_physeq_clean) # removes two samples SV011, SV021
rarefied_euk <- rarefy_even_depth(euk_physeq_20000, sample.size = min(sample_sums(euk_physeq_20000)),
          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) #Normalizing species data
richness_euk <- estimate_richness(rarefied_euk)
min(sample_sums(euk_physeq_20000))
summary(richness_euk)

# Bacteria
bac_physeq_9000 <- prune_samples(sample_sums(bac_physeq_clean)>=9000, bac_physeq_clean) # removes two samples SV011, SV021
rarefied_bac <- rarefy_even_depth(bac_physeq_9000, sample.size = min(sample_sums(bac_physeq_9000)),
          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) #Normalizing species data
richness_bac <- estimate_richness(rarefied_bac)
min(sample_sums(bac_physeq_9000))
summary(richness_bac)

colnames(richness_euk) = c("observed_euk",
                           "Chao1_euk",
                           "se.chao1_euk",
                           "ACE_euk",
                           "se.ACE_euk",
                           "Shannon_euk",
                           "Simpson_euk",
                           "InvSimpson_euk",
                           "Fisher_euk")
colnames(richness_bac) = c("observed_bac",
                           "Chao1_bac",
                           "se.chao1_bac",
                           "ACE_bac",
                           "se.ACE_bac",
                           "Shannon_bac",
                           "Simpson_bac",
                           "InvSimpson_bac",
                           "Fisher_bac")
# relations to metadata
# correlation matrix
sample_metadata = sample_data(rarefied_euk)[,c(7, 10, 12:15, 17:33)]
scaled_metadata = scale(sample_metadata)

diversity_metadata = cbind(richness_bac, richness_euk, scaled_metadata)

corr_spearman <- rcorr(as.matrix(diversity_metadata), type = "spearman")

# get Holm's adjusted p-value
corr_spearman_adjusted = print_rcorr_adjust(rcorr.adjust(diversity_metadata, type = "spearman",
                                              use = "pairwise.complete.obs"))

df_1 = dplyr::mutate_all(as.data.frame(corr_spearman_adjusted, stringsAsFactors = FALSE),
        function(x) as.numeric(as.character(x)))
df_1[is.na(df_1)] = 0.0001
df_2 = as.matrix(df_1)

postscript(file.path(figs_dir, "correlation_matrix_diversity_adjusted.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE)
  corrplot::corrplot(corr_spearman$r, type="upper", p.mat = df_2, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()

postscript(file.path(figs_dir, "correlation_matrix_diversity.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE)
  corrplot::corrplot(corr_spearman$r, type="upper", p.mat = corr_spearman$P, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()


# PLS/PCA
# Eukaryota
independent_variables = scale(sample_metadata[c(1:5,10:29),c(1:3, 5, 7:9, 11,13, 15:18,21:23)])
dependent_variables = richness_euk[c(1:5,10:29),c("ACE_euk", "Simpson_euk")]

pls2_reg_euk = plsreg2(independent_variables[5:nrow(independent_variables),], dependent_variables[5:nrow(dependent_variables),], comps=5, crosval=TRUE)
summary(pls2_reg_euk)

pls2_reg_euk$expvar
pls2_reg_euk$Q2cum
pls2_reg_euk$VIP
cairo_ps(file.path(figs_dir, "PLS_diversity_euk.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)

  plot(pls2_reg_euk, xlab = paste0("PC1, R2X = ",
                                   round(pls2_reg_euk$expvar[1], 2),
                                   ", R2Y = ",
                                   round(pls2_reg_euk$expvar[11],2)),
                     ylab = paste0("PC2, R2X = ",
                                   round(pls2_reg_euk$expvar[2], 2),
                                   ", R2Y = ",
                                   round(pls2_reg_euk$expvar[12],2)),
                 cex = 0.7,
                 cex.axis = 2,
                 cex.lab = 2,
                 main = NA)
dev.off()

# Bacteria
dependent_variables = richness_bac[c(1:5,10:29),c("ACE_bac", "Simpson_bac")]

pls2_reg_bac = plsreg2(independent_variables[5:nrow(independent_variables),], dependent_variables[5:nrow(dependent_variables),], comps=5, crosval=TRUE)
summary(pls2_reg_bac)

pls2_reg_bac$expvar
pls2_reg_bac$Q2cum
pls2_reg_bac$VIP

cairo_ps(file.path(figs_dir, "PLS_diversity_bac.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)

  plot(pls2_reg_bac, xlab = paste0("PC1, R2X = ",
                                   round(pls2_reg_bac$expvar[1], 2),
                                   ", R2Y = ",
                                   round(pls2_reg_bac$expvar[11],2)),
                     ylab = paste0("PC2, R2X = ",
                                   round(pls2_reg_bac$expvar[2], 2),
                                   ", R2Y = ",
                                   round(pls2_reg_bac$expvar[12],2)),
                 cex = 0.7,
                 cex.axis = 2,
                 cex.lab = 2,
                 main = NA)
dev.off()

## MAKE figure 2

library(pls)
library(tidyverse)

y_bac <- as.matrix(diversity_metadata[c(1:5,10:29),c(4)])
x_bac <- as.matrix(sample_metadata[c(1:5,10:29),c(1:3, 5, 7:9, 11,13, 15:18,21:23)])
df_plsr <- mvr(y_bac  ~ x_bac , ncomp = 2, method = "oscorespls" , scale = T)

df_coef <- as.data.frame(coef(df_plsr, ncomp =1:2))
df_coef <- df_coef %>%
  dplyr::mutate(variables = rownames(.)) %>%
  dplyr::mutate(variables = factor(variables,
                                   levels = c(
                                      "altitude", "gl_dist",  "Ccounts",  "Cond",
                                      "TOC", "DOC", "TN", "TP", "O2", "CO2", "CH4",
                                      "N2O", "bird", "CO2_sat", "CH4_sat", "N2O_sat")
                                      ))

colnames(df_coef)[1] <- "regression_coefficients"
df_coef  <- df_coef[, -c(2)]

# VIP

VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")

  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


vip_bac <- as.data.frame(VIP(df_plsr))
vip_bac_comp1 <- vip_bac[1, ]
row.names(vip_bac_comp1)<-NULL

combine_vip_coef <- vip_bac_comp1 %>%
  tidyr::gather("variables", "VIP") %>%
  dplyr::full_join(., df_coef, by = c("variables")) %>%
  dplyr::mutate(variables = factor(variables,
                                   levels = c(
                                      "altitude", "gl_dist",  "Ccounts",  "Cond",
                                      "TOC", "DOC", "TN", "TP", "O2", "CO2", "CH4",
                                      "N2O", "bird", "CO2_sat", "CH4_sat", "N2O_sat")
                                      ))

p1_vip <-
  ggplot(combine_vip_coef, aes(x = variables,  y = VIP, group =1))+
  geom_bar(stat="identity",  fill = "black") +
  geom_hline(yintercept = 1, size = 0.55, linetype = 3) +
  theme_bw()+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 6),
        axis.title.y = element_text(size = 10))+
  labs(x= "")

p2_coef <-
  ggplot(df_coef, aes(x = variables, y = regression_coefficients, group = 1))+  geom_bar(stat = "identity",  fill = "black")+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 8),
        axis.title.y = element_text(size = 2))+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  labs(x="")

cairo_ps(file.path(figs_dir, "PLS_diversity_bac_coeff.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)
p2_coef
dev.off()

cairo_ps(file.path(figs_dir, "PLS_diversity_bac_vip.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)
p1_vip
dev.off()

y_euk <- as.matrix(diversity_metadata[c(1:5,10:29),c(13)])
x_euk <- as.matrix(sample_metadata[c(1:5,10:29),c(1:3, 5, 7:9, 11,13, 15:18,21:23)])
df_plsr <- mvr(y_euk  ~ x_euk , ncomp = 2, method = "oscorespls" , scale = T)

df_coef <- as.data.frame(coef(df_plsr, ncomp =1:2))
df_coef <- df_coef %>%
  dplyr::mutate(variables = rownames(.)) %>%
  dplyr::mutate(variables = factor(variables,
                                   levels = c(
                                      "altitude", "gl_dist",  "Ccounts",  "Cond",
                                      "TOC", "DOC", "TN", "TP", "O2", "CO2", "CH4",
                                      "N2O", "bird", "CO2_sat", "CH4_sat", "N2O_sat")
                                      ))

colnames(df_coef)[1] <- "regression_coefficients"
df_coef  <- df_coef[, -c(2)]

# VIP
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")

  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

vip_euk <- as.data.frame(VIP(df_plsr))
vip_euk_comp1 <- vip_euk[1, ]
row.names(vip_euk_comp1)<-NULL

combine_vip_coef <- vip_euk_comp1 %>%
  tidyr::gather("variables", "VIP") %>%
  dplyr::full_join(., df_coef, by = c("variables")) %>%
  dplyr::mutate(variables = factor(variables,
                                   levels = c(
                                      "altitude", "gl_dist",  "Ccounts",  "Cond",
                                      "TOC", "DOC", "TN", "TP", "O2", "CO2", "CH4",
                                      "N2O", "bird", "CO2_sat", "CH4_sat", "N2O_sat")
                                      ))


p1_vip <-
  ggplot(combine_vip_coef, aes(x = variables,  y = VIP, group =1))+
  geom_bar(stat="identity",  fill = "black") +
  geom_hline(yintercept = 1, size = 0.55, linetype = 3) +
  theme_bw()+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 6),
        axis.title.y = element_text(size = 10))+
  labs(x= "")

p2_coef <-
  ggplot(df_coef, aes(x = variables, y = regression_coefficients, group = 1))+  geom_bar(stat = "identity",  fill = "black")+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 8),
        axis.title.y = element_text(size = 2))+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  labs(x="")

cairo_ps(file.path(figs_dir, "PLS_diversity_euk_coeff.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)
p2_coef
dev.off()

cairo_ps(file.path(figs_dir, "PLS_diversity_euk_vip.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)
p1_vip
dev.off()

## Make Figure 1
## Boxplots gas saturation
diversity_metadata = cbind(richness_bac, richness_euk, sample_metadata)

diversity_metadata$dist_level[diversity_metadata$gl_dist <= 550] = "closer than 550 m"
diversity_metadata$dist_level[diversity_metadata$gl_dist > 550] = "between 550 and 3300 m"
diversity_metadata$dist_level[diversity_metadata$gl_dist > 3300] = "over 3300 m"

gas_conc = c(diversity_metadata$CO2_sat * 100,
               diversity_metadata$CH4_sat * 100,
               diversity_metadata$N2O_sat * 100)
gas = rep(c("CO2", "CH4", "N2O"),
              each = length(gas_conc)/3)
dist_level = rep(diversity_metadata$dist_level, 3)

gasses <- data.frame(log(gas_conc), gas, dist_level)

cairo_ps(file.path(figs_dir, "boxplot_gasses.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
ggplot(gasses, aes(x = gas, y = gas_conc,
               colour = gas,
               shape = dist_level)) +
  geom_boxplot(aes(group = gas), outlier.shape = NA) +
  geom_jitter() +
  scale_y_continuous(trans='log10') +
  theme_bw()
dev.off()

cairo_ps(file.path(figs_dir, "boxplot_gasses_dist.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
ggplot(gasses, aes(x = gas, y = gas_conc,
               colour = gas,
               shape = dist_level)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_y_continuous(trans='log10') +
  theme_bw() + 
  theme(legend.position="none")
dev.off()

## GAMs
# CO2
summary(co2m1 <- gam(CO2 ~ s(TN) + s(TOC),
                    select=TRUE, data=data.frame(diversity_metadata[c(1:5, 10:29),])))
summary(gam(CO2 ~ TN + TOC, select=TRUE, data=data.frame(diversity_metadata[c(1:5, 10:29),])))

summary(co2sm1 <- gam(CO2_sat ~ s(TN) + s(DOC),
                    select=TRUE, data=data.frame(diversity_metadata[c(1:5, 10:29),])))
summary(gam(CO2_sat ~ TN + DOC, select=TRUE, data=data.frame(diversity_metadata[c(1:5, 10:29),])))

summary(co2sm2 <- gam(CO2_sat ~ s(gl_dist) ,
                    select=TRUE, data=data.frame(diversity_metadata[c(1:5, 10:29),])))


summary(co2m2 <- gam(CO2 ~ s(TN) + s(ACE_bac),
                    select=TRUE, data=data.frame(diversity_metadata[c(1:5, 10:29),])))

summary(co2sm1 <- gam(CO2_sat ~ s(Cond) + s(Simpson_bac),
                    select=TRUE, data=data.frame(diversity_metadata[c(1:5, 10:29),])))

# evaluation plots
cairo_ps(file.path(figs_dir, "co2m1_appraise.eps"), width=240/25.4, height=240/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  appraise(co2m1)
dev.off()
cairo_ps(file.path(figs_dir, "co2sm1_appraise.eps"), width=240/25.4, height=240/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  appraise(co2m2)
dev.off()
cairo_ps(file.path(figs_dir, "co2m2_appraise.eps"), width=240/25.4, height=240/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  appraise(co2m2)
dev.off()

# model plots
cairo_ps(file.path(figs_dir, "co2m1_visreg.eps"), width=240/25.4, height=240/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  visreg(co2m1, "TN", gg=TRUE) + theme_bw()
dev.off()

co2sm1 <- getViz(co2sm1)

o1 <- plot( sm(co2sm1, 1) ) +
    l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw()
o2 <- plot( sm(co2sm1, 2) ) +
    l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw()

cairo_ps(file.path(figs_dir, "co2sm1_plotsTN.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
o1
dev.off()
cairo_ps(file.path(figs_dir, "co2sm1_plotsDOC.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
o2
dev.off()
cairo_ps(file.path(figs_dir, "co2sm1_appraise.eps"), width=240/25.4, height=240/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  appraise(co2sm1)
dev.off()

# CH4
summary(ch4m1 <- gam(CH4 ~ s(TN) + s(DOC),
                    select=TRUE, data=diversity_metadata))
summary(ch4sm1 <- gam(CH4_sat ~ s(TN) + s(DOC),
                    select=TRUE, data=diversity_metadata))
summary(ch4sm3 <- gam(CH4_sat ~ TN + DOC,
                    select=TRUE, data=diversity_metadata))
summary(ch4m2 <- gam(CH4 ~ s(gl_dist) + s(DOC),
                    select=TRUE, data=diversity_metadata))
summary(ch4sm2 <- gam(CH4_sat ~ s(gl_dist) + s(TOC),
                    select=TRUE, data=diversity_metadata))

ch4sm1 <- getViz(ch4sm1)

o1 <- plot( sm(ch4sm1, 1) ) +
    l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw()
o2 <- plot( sm(ch4sm1, 2) ) +
    l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw()

cairo_ps(file.path(figs_dir, "ch4sm1_plotsTN.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
o1
dev.off()
cairo_ps(file.path(figs_dir, "ch4sm1_plotsDOC.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
o2
dev.off()
cairo_ps(file.path(figs_dir, "ch4sm1_appraise.eps"), width=240/25.4, height=240/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  appraise(ch4sm1)
dev.off()

summary(ch4sm1 <- gam(CH4_sat ~ s(TN) + s(DOC),
                    select=TRUE, data=diversity_metadata))
summary(ch4m2 <- gam(CH4 ~ s(gl_dist),
                    select=TRUE, data=diversity_metadata))

# plot
ch4m2 <- getViz(ch4m2)

o1 <- plot( sm(ch4m2, 1) ) +
    l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw()

cairo_ps(file.path(figs_dir, "ch4m2_plotsgl.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
o1
dev.off()
cairo_ps(file.path(figs_dir, "ch4m2_appraise.eps"), width=240/25.4, height=240/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  appraise(ch4m1)
dev.off()

summary(ch4sm3 <- gam(CH4 ~ s(gl_dist) + s(DOC),
                    select=TRUE, data=diversity_metadata))

summary(ch4sm3 <- gam(CH4 ~ s(ACE_bac) + s(TN),
                    select=TRUE, data=diversity_metadata))
summary(ch4m3 <- gam(CH4 ~ s(TN),
                    select=TRUE, data=diversity_metadata))
summary(ch4m4 <- gam(CH4 ~ s(ACE_bac),
                    select=TRUE, data=diversity_metadata))

summary(ch4m5 <- lm(CH4 ~ scale(gl_dist),
                    select=TRUE, data=diversity_metadata))
ch4m4 <- getViz(ch4m4)
o1 <- plot( sm(ch4m4, 1) ) +
    l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw()

# Plot
cairo_ps(file.path(figs_dir, "ch4m4_plotACE.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
o1
dev.off()

# N2O
summary(n2o1 <- gam(N2O_sat ~ s(ACE_bac),
                    select=TRUE, data=diversity_metadata))
summary(n2o2 <- gam(N2O_sat ~ s(TN),
                    select=TRUE, data=diversity_metadata))
summary(n2o3 <- gam(N2O_sat ~ s(TN) + s(DOC),
                    select=TRUE, data=diversity_metadata))
summary(n2o4 <- gam(N2O_sat ~ s(ACE_bac) + s(TN),
                    select=TRUE, data=diversity_metadata))

cairo_ps(file.path(figs_dir, "n2o3_appraise.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  appraise(n2o3)
dev.off()

o1 <- plot( sm(n2o3, 1) ) +
    l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw()

cairo_ps(file.path(figs_dir, "n2o3_draw.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  o1
dev.off()

## Bacterial cell counts and GHG saturation stats
summary(metadata)
summary(metadata[which(metadata$gl_dist <= 550),])
summary(metadata[which(metadata$gl_dist >= 3350),])
apply(metadata[which(metadata$gl_dist <= 550),], 2, sd, na.rm = TRUE)
apply(metadata[which(metadata$gl_dist >= 3350),], 2, sd, na.rm = TRUE)

## Alpha diversity
# ACE
summary(m2b <- gam(ACE_bac ~ s(TN) + s(gl_dist),
                    select=TRUE, data=diversity_metadata))

# Simpson
summary(m2b <- gam(Simpson_euk ~ s(TN) + as.factor(bird),
                    select=TRUE, data=diversity_metadata))
summary(m2b <- gam(Simpson_euk ~ s(TN),
                    select=TRUE, data=diversity_metadata))
summary(m2b <- gam(Simpson_euk ~ s(DOC),
                    select=TRUE, data=diversity_metadata))

## beta diversity
# NMDS - grouping of samples
# Eukaryota
otu_euk <- otu_table(rarefied_euk)
otu_euk <- otu_euk[rowSums(otu_euk) != 0,]
otu_euk_stand = decostand(otu_euk, MARGIN = 1, method="hellinger")
saveRDS(otu_euk_stand, file = file.path(figs_dir, "otu_euk_stand.rds"))

otu_euk_dist = vegdist(t(otu_euk_stand), method = "bray")
# saveRDS(otu_euk_dist, file = file.path(figs_dir, "otu_euk_dist.rds"))
euk_meta <- sample_data(rarefied_euk)

mds_euk = vegan::metaMDS(t(otu_euk_stand), distance = "bray", autotransform = FALSE, try = 1000)
# saveRDS(mds_euk, file = file.path(figs_dir, "mds_euk.rds"))
mds_euk_scores = data.frame(vegan::scores(mds_euk, "sites"))
mds_euk_scores$sites <- rownames(mds_euk_scores)
names(mds_euk_scores)[c(1, 2)] <- c("x", "y")
mds_euk_scores$z <- NA

euk_ordisurf <- ordisurf(mds_euk ~ euk_meta$"gl_dist", plot = FALSE, scaling = 3)
head(euk_ordisurf)

contour_vals <- extract.xyz(obj = euk_ordisurf)
head(contour_vals)

p <- ggplot(data = contour_vals, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + theme_bw()
nmds_euk <- p + geom_text(data = mds_euk_scores, aes(x = x, y = y, label = rownames(mds_euk_scores)),
    colour = "red") + coord_equal() + theme_bw() + labs(x = "NMDS1", y = "NMDS2") +
    theme(panel.border = element_rect(fill = NA), axis.text.x = element_blank(),
        axis.text.y = element_blank(), legend.position = "none")

# Bacteria
otu_bac <- otu_table(rarefied_bac)
otu_bac <- otu_bac[rowSums(otu_bac) != 0,]
otu_bac_stand = decostand(otu_bac, MARGIN = 1, method="hellinger")
saveRDS(otu_bac_stand, file = file.path(figs_dir, "otu_bac_stand.rds"))

otu_bac_dist = vegdist(t(otu_bac_stand), method = "bray")
# saveRDS(otu_bac_dist, file = file.path(figs_dir, "otu_bac_dist.rds"))
bac_meta <- sample_data(rarefied_bac)

mds_bac = vegan::metaMDS(t(otu_bac_stand), distance = "bray", autotransform = FALSE, try = 1000)
# saveRDS(mds_bac, file = file.path(figs_dir, "mds_bac.rds"))
mds_bac_scores = data.frame(vegan::scores(mds_bac))
mds_bac_scores$sites <- rownames(mds_bac_scores)
names(mds_bac_scores)[c(1, 2)] <- c("x", "y")
mds_bac_scores$z <- NA

bac_ordisurf <- ordisurf(mds_bac ~ bac_meta$"gl_dist", plot = FALSE, scaling = 3)
head(bac_ordisurf)

contour_vals <- extract.xyz(obj = euk_ordisurf)
head(contour_vals)

p <- ggplot(data = contour_vals, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(-1, 1.5)) + theme_bw()


nmds_bac <- p + geom_text(data = mds_bac_scores, aes(x = x, y = y, label = rownames(mds_bac_scores)),
    colour = "red") + coord_equal() + theme_bw() + labs(x = "NMDS1", y = "NMDS2") +
    theme(panel.border = element_rect(fill = NA), axis.text.x = element_blank(),
        axis.text.y = element_blank(), legend.position = "none")

# NMDS plot - plot distance from glacier in color and the 4 localities with symbols
cairo_ps(file.path(figs_dir, "nmds_euk.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  nmds_euk
dev.off()

cairo_ps(file.path(figs_dir, "nmds_bac.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  nmds_bac
dev.off()

## RDA and envfit
# Eukaryota
otu_euk_stand_subset = otu_euk_stand[,5:ncol(otu_euk_stand)]
mds_euk_subset = vegan::metaMDS(t(otu_euk_stand_subset), distance = "bray", autotransform = FALSE, try = 1000)

envfit_results <- vegan::envfit(mds_euk_subset, independent_variables, na.rm = TRUE, permu = 999)
envfit_euk_bray <- data.frame(round((envfit_results$vectors)$arrows, 3), round((envfit_results$vectors)$r, 3), round((envfit_results$vectors)$pvals, 3))
colnames(envfit_euk_bray) <- c("DIM1", "DIM2", "R", "p")
write.csv(envfit_euk_bray, file.path(figs_dir, "Table_envfit_euk.csv"))

# Bacteria
otu_bac_stand_subset = otu_bac_stand[,5:ncol(otu_bac_stand)]
mds_bac_subset = vegan::metaMDS(t(otu_bac_stand_subset), distance = "bray", autotransform = FALSE, try = 1000)

envfit_results <- vegan::envfit(mds_bac_subset, independent_variables, na.rm = TRUE, permu = 999)
envfit_bac_bray <- data.frame(round((envfit_results$vectors)$arrows, 3), round((envfit_results$vectors)$r, 3), round((envfit_results$vectors)$pvals, 3))
colnames(envfit_bac_bray) <- c("DIM1", "DIM2", "R", "p")
write.csv(envfit_bac_bray, file.path(figs_dir, "Table_envfit_bac.csv"))

# Make Figure envfit variables
tmp <- data.frame(R_value = c(envfit_euk_bray$R, envfit_bac_bray$R),
                  p_value = c(envfit_euk_bray$p, envfit_bac_bray$p),
                  Domain = rep(c("Eukarya", "Bacteria"),
                               each = length(envfit_euk_bray$p)),
                  env_var = rownames(envfit_euk_bray))
tmp = tmp[c(1:8,13:24,29:32),]
cairo_ps(file.path(figs_dir, "envfit_plot.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
ggplot(tmp, aes(x = Domain, y = env_var)) +
  geom_point(aes(size = 1/p_value, color = R_value), alpha = 0.7) +
  scale_color_gradient(low = "#E69F00", high = "#56B4E9") +
  scale_size(breaks = c(20,100,1000)) +
  theme_bw()
dev.off()

################## mvabund + gllvm #######################
# analyis of single ASVs and their distribution along the environmental gradients
			 
otu_euk_mva <- mvabund(otu_euk)
rownames(otu_bac) <-
otu_bac_mva <- mvabund(otu_bac)

# check mean-variance relationship
cairo_ps(file.path(figs_dir, "meanvar_plot_euk.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  meanvar.plot(otu_euk_mva)
dev.off()
cairo_ps(file.path(figs_dir, "meanvar_plot_bac.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  meanvar.plot(otu_bac_mva)
dev.off()

# prep data for gllvm
#euk
colnames(sample_data(rarefied_euk))
metadata_matrix_euk = as.matrix(sample_data(rarefied_euk)[,c(10,18,20,21,28)])

independent_variables_euk = scale(metadata_matrix_euk)
most_abundant_euk = otu_euk[rowSums(otu_euk)>100,]
most_abundant_euk_biom = make_biom(most_abundant_euk, sample_metadata = NULL,
                            observation_metadata = NULL,
                            id = NULL, matrix_element_type = "int")
write_biom(most_abundant_euk_biom, biom_file=file.path(figs_dir, "most_abundant_euk.biom"))

header <- taxa_names(rarefied_euk)
fasta <- dplyr::data_frame(name = header,
                                seq = header)
write_fasta(fasta, file.path(figs_dir, "euk-seqs.fna"))

#bac
metadata_matrix_bac = as.matrix(sample_data(rarefied_bac)[,c(10,18,20,21,28)])

independent_variables_bac = scale(metadata_matrix_bac)
most_abundant_bac = otu_bac[rowSums(otu_bac)>100,]
most_abundant_bac_biom = make_biom(most_abundant_bac, sample_metadata = NULL,
                            observation_metadata = NULL,
                            id = NULL, matrix_element_type = "int")
write_biom(most_abundant_bac_biom, biom_file=file.path(figs_dir, "most_abundant_bac.biom"))
ASV_no_taxo = readDNAStringSet(file.path(path_bac,"dada2/ASV_no_taxo.fasta"))

# fix names
names(ASV_no_taxo) = seq(ASV_no_taxo)

# subset most abundant sequences
ASV_no_taxo_most_abund = ASV_no_taxo[names(ASV_no_taxo) %in% rownames(most_abundant_bac)]

writeXStringSet(ASV_no_taxo_most_abund, file.path(figs_dir,"seq_bac_most_abund.fasta"), format = "fasta")

# use subset of ASVs - most abundant

# use genera - most abundant
seed_num = 1234
gllvm_euk1 = gllvm(t(most_abundant_euk), family = "negative.binomial")

criterias <- NULL
for(i in 1:5){
  fiti_bac <- gllvm(t(most_abundant_euk), independent_variables_euk,
                family = "negative.binomial", num.lv = i, sd.errors = FALSE,
                formula = ~ gl_dist + DN + DOC, seed = 1234)
  criterias[i + 1] <- summary(fiti_bac)$AICc
  names(criterias)[i + 1] = i
}

gllvm_euk2 = gllvm(t(most_abundant_euk), independent_variables_euk, family = "negative.binomial", num.lv = 4)

gllvm_bac2 = gllvm(t(most_abundant_bac), independent_variables_bac, family = "negative.binomial", num.lv = 4)

cairo_ps(file.path(figs_dir, "gllvm_euk2_eval.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  plot(gllvm_euk2, cex.ylab = 0.7)
dev.off()

cairo_ps(file.path(figs_dir, "gllvm_bac2_eval.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  plot(gllvm_bac2, cex.ylab = 0.7)
dev.off()

cairo_ps(file.path(figs_dir, "gllvm_euk2.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  coefplot(gllvm_euk2, cex.ylab = 0.3)
dev.off()

cairo_ps(file.path(figs_dir, "gllvm_bac2.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  coefplot(gllvm_bac2, cex.ylab = 0.3)
dev.off()

criterias <- NULL
for(i in 1:5){
  fiti_bac <- gllvm(t(most_abundant_bac), independent_variables_bac,
                family = "negative.binomial", num.lv = i, sd.errors = FALSE,
                formula = ~ gl_dist + DN + DOC, seed = 1234)
  criterias[i + 1] <- summary(fiti_bac)$AICc
  names(criterias)[i + 1] = i
}

######## procrustes test ##########
# relationship between bacterial and eukaryotic beta diversity			
pro <- procrustes(X = mds_euk_subset, Y = mds_bac_subset, symmetric = FALSE)
pro

pro <- procrustes(X = mds_euk_subset, Y = mds_bac_subset, symmetric = TRUE)
pro

protest(X = mds_euk_subset, Y = mds_bac_subset, scores = "sites", permutations = 999)

dist_bac_subset = vegan::vegdist(t(otu_bac), method = "bray")
dist_euk_subset = vegan::vegdist(t(otu_euk), method = "bray")
mantel(dist_euk_subset, dist_bac_subset, method = "spearman", permutations = 9999, na.rm = TRUE)


############## taxonomic composition ##################
# Eukaryotes
# BARPLOT - taxonomic composition

# RDA with taxonomic groups
taxonomy = "Class"

classGlom = tax_glom(rarefied_euk, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix')
rownames(taxon_table) = tax_matrix[,taxonomy]
tax_table = prop.table(taxon_table, margin = 2)*100
tax_table <- tax_table[order(rowSums(-taxon_table)),]
rowSums(taxon_table)
dim(taxon_table)
tax_table <- tax_table[1:15,]
rownames(tax_table)

# calculate proportions based on phylum data
melted_tax_table = melt(tax_table)
melted_tax_table <- dplyr::arrange(melted_tax_table, Var1, desc(value))
melted_tax_table$Var2 <- factor(melted_tax_table$Var2 , levels = unique(melted_tax_table$Var2))

class_colors <- setNames(color_palette_euk, levels(melted_tax_table$Var1))

svg(file.path(figs_dir, "boxplot_class_euk.svg"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE)
  ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
     geom_bar(stat = "identity", position = "stack", alpha = .5) +
     guides(fill = guide_legend(title = taxonomy)) +
     coord_flip() +
     theme(axis.text = element_text(size=6),
       axis.title = element_text(size=10, face="bold"),
       legend.text = element_text(size=8),
       plot.background = element_rect(fill = "white"),
       panel.background = element_rect(fill = "white"),
       axis.line.x = element_line(color = "grey")) +
     xlab("lake systems") +
     ylab("proportion of reads [%]") +
     scale_fill_manual(values = class_colors) +
     scale_x_discrete(limits = rev(levels(melted_tax_table$Var2))
     )
dev.off()

taxonomy = "Family"

classGlom = tax_glom(rarefied_euk, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix')
rownames(taxon_table) = tax_matrix[,taxonomy]
tax_table = prop.table(taxon_table, margin = 2)*100
tax_table <- tax_table[order(rowSums(-taxon_table)),]
rowSums(taxon_table)
dim(taxon_table)
tax_table <- tax_table[1:15,]
rownames(tax_table)

# calculate proportions based on phylum data
melted_tax_table = melt(tax_table)
melted_tax_table <- dplyr::arrange(melted_tax_table, Var1, desc(value))
melted_tax_table$Var2 <- factor(melted_tax_table$Var2 , levels = unique(melted_tax_table$Var2))

class_colors <- setNames(color_palette_euk, levels(melted_tax_table$Var1))

svg(file.path(figs_dir, "boxplot_family_euk.svg"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE)
  ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
     geom_bar(stat = "identity", position = "stack", alpha = .5) +
     guides(fill = guide_legend(title = taxonomy)) +
     coord_flip() +
     theme(axis.text = element_text(size=6),
       axis.title = element_text(size=10, face="bold"),
       legend.text = element_text(size=8),
       plot.background = element_rect(fill = "white"),
       panel.background = element_rect(fill = "white"),
       axis.line.x = element_line(color = "grey")) +
     xlab("lake systems") +
     ylab("proportion of reads [%]") +
     scale_fill_manual(values = class_colors) +
     scale_x_discrete(limits = rev(levels(melted_tax_table$Var2))
     )
dev.off()

taxonomy = "Genus"

classGlom = tax_glom(rarefied_euk, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix')
rownames(taxon_table) = tax_matrix[,taxonomy]
tax_table = prop.table(taxon_table, margin = 2)*100



# calculate proportions based on phylum data
melted_tax_table = melt(tax_table)
melted_tax_table <- dplyr::arrange(melted_tax_table, Var1, desc(value))
melted_tax_table$Var2 <- factor(melted_tax_table$Var2 , levels = unique(melted_tax_table$Var2))

class_colors <- setNames(color_palette_euk, levels(melted_tax_table$Var1))

svg(file.path(figs_dir, "boxplot_genera_euk.svg"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE)
  ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
     geom_bar(stat = "identity", position = "stack", alpha = .5) +
     guides(fill = guide_legend(title = taxonomy)) +
     coord_flip() +
     theme(axis.text = element_text(size=6),
       axis.title = element_text(size=10, face="bold"),
       legend.text = element_text(size=8),
       plot.background = element_rect(fill = "white"),
       panel.background = element_rect(fill = "white"),
       axis.line.x = element_line(color = "grey")) +
     xlab("lake systems") +
     ylab("proportion of reads [%]") +
     scale_fill_manual(values = class_colors) +
     scale_x_discrete(limits = rev(levels(melted_tax_table$Var2))
     )
dev.off()

## RDA on taxonomic diversity

metadata4 = data.frame(sample_data(classGlom))
metadata4 = metadata4[,c("gl_dist", "Temp", "TOC", "DN", "bird")]

rda_div_meta <- rda(decostand(t(tax_table), method="hellinger"), as.matrix(metadata4))

rda_div_meta
summary_rda_div_meta <- summary(rda_div_meta)
RsquareAdj(rda_div_meta)

metadata_class = cbind(metadata4, t(taxon_table))
metadata_class_spearman <- rcorr(as.matrix(metadata_class), type = "spearman")

cairo_ps(file.path(figs_dir, "metadata_eukgenera_spearman.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  corrplot::corrplot(metadata_class_spearman$r, type="upper", p.mat = metadata_class_spearman$P, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()

p_rda <- autoplot(rda_div_meta, layers = c("species", "biplot"), legend.position = "none")

rda_scores_env = vegan::scores(rda_div_meta, display = "bp")
rda_scores_species = vegan::scores(rda_div_meta, display = "sp")

names = rownames(rda_scores_species)
names
melted_tax_table$Var1


rda_plot_species <- ggplot(data.frame(rda_scores_species), aes(x = RDA1, y = RDA2, color = names)) +
                geom_point(size = 5, alpha = .5) +
                scale_color_manual(values = class_colors)

mult = 0.2

rda_biplot_class <- rda_plot_species +
  geom_segment(data = data.frame(rda_scores_env), aes(x = 0, xend = mult * RDA1,
                    y = 0, yend = mult * RDA2),
                arrow = arrow(length = unit(0.25, "cm")), colour = "red", alpha = .4) +
  geom_text(data = data.frame(rda_scores_env),
            aes(x = mult * RDA1, y = mult * RDA2, label = rownames(rda_scores_env),
                hjust = 0.5 * (1-sign(RDA1)), vjust = 0.5 * (1-sign(RDA2))),
                color = "red", size = 3, alpha = .4) +
  coord_cartesian(xlim = c(-0.50, 0.30), ylim = c(-0.3, 0.45)) +
  theme(plot.background = element_rect(fill = "white"),
                panel.background = element_rect(fill = "white"),
                panel.grid = element_line(color = "grey"),
                legend.position = "None")

svg(file.path(figs_dir, "rda_genera_euk.svg"), width=80/25.4, height=80/25.4, pointsize = 6, bg = FALSE)
  rda_biplot_class
dev.off()

# bac
taxonomy = "Genus"
classGlom = tax_glom(rarefied_bac, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix')
rownames(taxon_table) = tax_matrix[,taxonomy]
tax_table = prop.table(taxon_table, margin = 2)*100

# methanotrophs
methanotrophs = c("Methylibium",
                  "Methylacidimicrobium",
                  "Methyloacida",
                  "Methylobacter",
                  "Methylocaldum",
                  "Methylococcus",
                  "Methylocystis",
                  "Methylocella",
                  "Methylocapsa",
                  "Methyloferula",
                  "Methyloglobulus",
                  "Methylohalobius",
                  "Methylomagnum",
                  "Methylomarinum",
                  "Methylomarinovum",
                  "Methyloparacoccus",
                  "Methylosphaera",
                  "Methyloprofundus",
                  "Methylophilus",
                  "Methylosarcina",
                  "Methylosoma",
                  "Methylothermus",
                  "Methylogaea",
                  "Methylovulum",
                  "Crenothrix",
                  "Candidatus_Methylopumilus",
                  "Candidatus_Methylacidiphilum",
                  "Candidatus_Anammoximicrobium")

genus_methanotrophs = taxon_table[rownames(taxon_table) %in% methanotrophs,]

genus_methanotrophs_prop = genus_methanotrophs/colSums(taxon_table)*100

melted_tax_table = melt(genus_methanotrophs_prop)
melted_tax_table <- dplyr::arrange(melted_tax_table, Var1, desc(value))
melted_tax_table$Var2 <- factor(melted_tax_table$Var2 , levels = unique(melted_tax_table$Var2))

class_colors <- setNames(color_palette_bac, levels(melted_tax_table$Var1))

svg(file.path(figs_dir, "boxplot_methanotrophs.svg"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE)
  ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
     geom_bar(stat = "identity", position = "stack", alpha = .5) +
     guides(fill = guide_legend(title = taxonomy)) +
     coord_flip() +
     theme(axis.text = element_text(size=6),
       axis.title = element_text(size=10, face="bold"),
       legend.text = element_text(size=8),
       plot.background = element_rect(fill = "white"),
       panel.background = element_rect(fill = "white"),
       axis.line.x = element_line(color = "grey")) +
     xlab("lake systems") +
     ylab("proportion of reads [%]") +
     scale_fill_manual(values = class_colors) +
     scale_x_discrete(limits = rev(levels(melted_tax_table$Var2))
     )
dev.off()

#genera
tax_table <- tax_table[order(rowSums(-taxon_table)),]
rowSums(taxon_table)
dim(taxon_table)
tax_table <- tax_table[1:15,]
rownames(tax_table)

# calculate proportions based on phylum data
melted_tax_table = melt(tax_table)
melted_tax_table <- dplyr::arrange(melted_tax_table, Var1, desc(value))
melted_tax_table$Var2 <- factor(melted_tax_table$Var2 , levels = unique(melted_tax_table$Var2))

class_colors <- setNames(color_palette_bac, levels(melted_tax_table$Var1))

svg(file.path(figs_dir, "boxplot_genera_bac.svg"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE)
  ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
     geom_bar(stat = "identity", position = "stack", alpha = .5) +
     guides(fill = guide_legend(title = taxonomy)) +
     coord_flip() +
     theme(axis.text = element_text(size=6),
       axis.title = element_text(size=10, face="bold"),
       legend.text = element_text(size=8),
       plot.background = element_rect(fill = "white"),
       panel.background = element_rect(fill = "white"),
       axis.line.x = element_line(color = "grey")) +
     xlab("lake systems") +
     ylab("proportion of reads [%]") +
     scale_fill_manual(values = class_colors) +
     scale_x_discrete(limits = rev(levels(melted_tax_table$Var2))
     )
dev.off()


## fecal indicators
taxonomy = "Family"
classGlom = tax_glom(rarefied_bac, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix')
rownames(taxon_table) = tax_matrix[,taxonomy]
tax_table = prop.table(taxon_table, margin = 2)*100
dim(taxon_table)

fecal_family = c("Fusobacteriaceae", "Enterococcaceae", "Enterobacteriaceae",
                 "Clostridiaceae", "Bacteriodaceae")

family_fecal = tax_table[rownames(tax_table) %in% fecal_family,]

sample_data(classGlom)

fecal_meta = cbind(sample_data(classGlom), t(family_fecal))
fecal_meta = fecal_meta[,c(28, 34, 35)]

fecal_spearman = rcorr.adjust(as.matrix(fecal_meta), type = "spearman")


# identify ASVs that show trends with gl_dist - distance to glacier
# extract ASVs with confidence interval > 2 or < -2

fit_sign_split_gl_dist1 = fit_sign_split[fit_sign_split[,1]>=2,]
fit_sign_split_gl_dist2 = fit_sign_split[fit_sign_split[,2]<=-2,]
fit_sign_split_gl_dist = rbind(fit_sign_split_gl_dist1,fit_sign_split_gl_dist2)

# check taxonomic assignment
tax_bac_gl_dist1 = tax_bac[rownames(tax_bac) %in% rownames(fit_sign_split_gl_dist1),]
tax_bac_gl_dist2 = tax_bac[rownames(tax_bac) %in% rownames(fit_sign_split_gl_dist2),]

# identify ASVs that show trends with DN - dissolved nitrogen or more generl with nutrient status
fit_sign_bac_DN1 = fit_sign_split[fit_sign_split[,3]>=2,]
fit_sign_bac_DN2 = fit_sign_split[fit_sign_split[,4]<=-2,]

tax_bac_DN1 = tax_bac[rownames(tax_bac) %in% rownames(fit_sign_bac_DN1),]
tax_bac_DN2 = tax_bac[rownames(tax_bac) %in% rownames(fit_sign_bac_DN2),]

fit_sign_bac_DOC1 = fit_sign_split[fit_sign_split[,5]>=2,]
fit_sign_bac_DOC2 = fit_sign_split[fit_sign_split[,6]<=-2,]

tax_bac_DOC1 = tax_bac[rownames(tax_bac) %in% rownames(fit_sign_bac_DOC1),]
tax_bac_DOC2 = tax_bac[rownames(tax_bac) %in% rownames(fit_sign_bac_DOC2),]


## Eukaryotes
fit_euk_confin = confint(fit_env_euk, level = 0.95, parm = "Xcoef")


colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
col.names = c("gl_dist_low", "gl_dist_up", "DN_low", "DN_up", "DOC_low", "DOC_up")
fit_sign_split_euk <- read.table(text = "",
                 colClasses = colClasses,
                 col.names = col.names)

fit_sign_split_euk = cbind(fit_euk_confin[1:437,],fit_euk_confin[438:874,],fit_euk_confin[875:1311,])
colnames(fit_sign_split_euk) = col.names
rownames(fit_sign_split_euk) = rownames(fit_env_euk$param$Xcoef)

fit_sign_euk_gl_dist1 = fit_sign_split_euk[fit_sign_split_euk[,1]>=2,]
fit_sign_euk_gl_dist2 = fit_sign_split_euk[fit_sign_split_euk[,2]<=-2,]
fit_sign_euk_gl_dist = rbind(fit_sign_euk_gl_dist1,fit_sign_euk_gl_dist2)

tax_euk_gl_dist1 = tax_euk[rownames(tax_euk) %in% rownames(fit_sign_euk_gl_dist1),]
tax_euk_gl_dist2 = tax_euk[rownames(tax_euk) %in% rownames(fit_sign_euk_gl_dist2),]


fit_sign_euk_DN1 = fit_sign_split_euk[fit_sign_split_euk[,3]>=2,]
fit_sign_euk_DN2 = fit_sign_split_euk[fit_sign_split_euk[,4]<=-2,]

tax_euk_DN1 = tax_euk[rownames(tax_euk) %in% rownames(fit_sign_euk_DN1),]
tax_euk_DN2 = tax_euk[rownames(tax_euk) %in% rownames(fit_sign_euk_DN2),]

fit_sign_euk_DOC1 = fit_sign_split_euk[fit_sign_split_euk[,5]>=2,]
fit_sign_euk_DOC2 = fit_sign_split_euk[fit_sign_split_euk[,6]<=-2,]

tax_euk_DOC1 = tax_euk[rownames(tax_euk) %in% rownames(fit_sign_euk_DOC1),]
tax_euk_DOC2 = tax_euk[rownames(tax_euk) %in% rownames(fit_sign_euk_DOC2),]

########## MINE analysis #############
otu_bac <- otu_table(rarefied_bac)
otu_bac <- otu_bac[rowSums(otu_bac) != 0,]

otu_euk <- otu_table(rarefied_euk)
otu_euk <- otu_euk[rowSums(otu_euk) != 0,]

dim(otu_euk)
dim(otu_bac)

rownames(otu_bac) = paste(rownames(otu_bac), "B", sep = "")
rownames(otu_euk) = paste(rownames(otu_euk), "E", sep = "")

otu_bac_tax = tax_table(rarefied_bac)
otu_euk_tax = tax_table(rarefied_euk)
rownames(otu_bac_tax) = paste(rownames(otu_bac_tax), "B", sep = "")
rownames(otu_euk_tax) = paste(rownames(otu_euk_tax), "E", sep = "")

otu_tot_tax = rbind(otu_bac_tax[,c(1,6)], otu_euk_tax[,c(1,7)])

otu_tot = t(rbind(otu_bac, otu_euk))
otu_tot = otu_tot[, apply(otu_tot == 0, 2, sum) <= 25 ]

ticenull = minerva::mictools(otu_tot, nperm=1000, p.adjust.method = "BH")
ticenull_pvalues = as.vector(ticenull$pval$pval)

corr_res = Hmisc::rcorr(as.matrix(otu_tot), type = "spearman")
corr_r = reshape2::melt(corr_res$r)
corr_res_adj = rcorr_adjust(as.matrix(otu_tot), type = "spearman")

n_value = corr_res_adj$R$n
p_value_adj = corr_res_adj$P
r_value = round(corr_res_adj$R$r, 4)

mine_res = minerva::mine(otu_tot, alpha = 1)
mine_res$MAS[upper.tri(mine_res$MAS)] <- NA
mine_mas = reshape2::melt(mine_res$MAS, na.rm = TRUE)
mine_res$MIC[upper.tri(mine_res$MIC)] <- NA
mine_mic = reshape2::melt(mine_res$MIC, na.rm = TRUE)
mine_res$MEV[upper.tri(mine_res$MEV)] <- NA
mine_mev = reshape2::melt(mine_res$MEV, na.rm = TRUE)
mine_res$MCN[upper.tri(mine_res$MCN)] <- NA
mine_mcn = reshape2::melt(mine_res$MCN, na.rm = TRUE)

colnames(ticenull$pval) = c("pval", "I1", "I2", "Var2", "Var1", "adj.P.Val")

tmp = merge(ticenull$pval, mine_mas, by = c("Var1", "Var2"))
dim(tmp)
colnames(tmp)[ncol(tmp)] <- "MAS"
tmp = merge(tmp, mine_mic, by = c("Var1", "Var2"))
colnames(tmp)[ncol(tmp)] <- "MIC"
tmp = merge(tmp, mine_mev, by = c("Var1", "Var2"))
colnames(tmp)[ncol(tmp)] <- "MEV"
tmp = merge(tmp, mine_mcn, by = c("Var1", "Var2"))
colnames(tmp)[ncol(tmp)] <- "MCN"
tmp = merge(tmp, corr_r, by = c("Var1", "Var2"))
colnames(tmp)[ncol(tmp)] <- "r_corr"

comp_bac = tmp[tmp$Var1 %like% "B",]
bac_comp = nrow(comp_bac)
tot_comp = nrow(tmp)
comp_euk = tmp[tmp$Var1 %like% "E",]
comp_euk = comp_euk[comp_euk$Var2 %like% "E",]
euk_comp = nrow(comp_euk)

comp_sig_euk = comp_euk[comp_euk$adj.P.Val < 0.05, ]
comp_sig_bac = comp_bac[comp_bac$adj.P.Val < 0.05, ]

bac_sig_comp = nrow(comp_sig_bac)
euk_sig_comp = nrow(comp_sig_euk)

cont_bac_sig_comp = bac_sig_comp/bac_comp # prevalence of significant co-occurrences within the bacterial kingdom 2.3%
cont_euk_sig_comp = euk_sig_comp/euk_comp # prevalence of significant co-occurrences within the eukaryotic kingdom 3.2%

comp_comp1 = comp_bac[comp_bac$Var2 %like% "E",]
comp_sig_comp1 = comp_comp1[comp_comp1$adj.P.Val < 0.05, ]
comp_comp2 = tmp[tmp$Var1 %like% "E",]
comp_comp2 = comp_comp2[comp_comp2$Var2 %like% "B",]
comp_sig_comp2 = comp_comp2[comp_comp2$adj.P.Val < 0.05, ]

cont_comp_sig_comp = (nrow(comp_sig_comp1)+nrow(comp_sig_comp2))/(nrow(comp_comp1)+nrow(comp_comp2))


