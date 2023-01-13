### This is an R script to analyse to infer greenhouse gas emissions and metagenomic
### data from Svalbard


## saga
# module load R/4.0.0-foss-2020a
# path to metags

path_out = "/cluster/projects/nn9745k/jing/02_results/svalbard"
path_metags = file.path("/cluster/work/users/alexaei/02_results/13_svalbard_metaGs")

figs_dir <- file.path(path_out, "figs")

# load(file.path(figs_dir, "metags_stats.Rdata"))


#In R
## LOAD PACKAGES
.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_0_0", .libPaths()))

if (!require("vegan", quietly = TRUE)){
  install.packages("vegan")
  library(vegan)
}
if (!require("reshape2", quietly = TRUE)){
  install.packages("reshape2")
  library(reshape2)
}

if (!require("caret", quietly = TRUE)){
  install.packages("caret")
  library(caret)
}

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

#if (!require("Deseq2")) {
#  BiocManager::install("DESeq2", lib="/cluster/projects/nn9745k/Rpackages_4_0_0")
   library(DESeq2)
#}
#if (!require("Deseq2")) {
#  BiocManager::install("phyloseq", lib="/cluster/projects/nn9745k/Rpackages_4_0_0")
   library(phyloseq)
#}

# if (!require("gratia")) {
#   install.packages("gratia", dependencies = TRUE)
#   library(gratia)
# }
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
if (!require("mgcv")) {
  install.packages("mgcv", dependencies = TRUE)
  library(mgcv)
}
if (!require("mgcViz")) {
  install.packages("mgcViz", dependencies = TRUE)
  library(mgcViz)
}
if (!require("nlme")) {
  install.packages("nlme", dependencies = TRUE)
  library(nlme)
}
if (!require("gllvm")) {
  install.packages("gllvm", dependencies = TRUE)
  library(gllvm)
}
if (!require("pheatmap")) {
  install.packages("pheatmap", dependencies = TRUE)
  library(pheatmap)
}



## LOAD local scripts
extract.xyz <- function(obj) {
    xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
    xyz <- cbind(xy, c(obj$grid$z))
    names(xyz) <- c("x", "y", "z")
    return(xyz)
}

color_palette_euk = c("#D1BBD7", "#AE76A3", "#882E72", "#1965B0",
                        "#5289C7", "#7BAFDE", "#4EB265", "#90C987",
                        "#CAE0AB", "#F7F056", "#F6C141", "#F1932D",
                        "#E8601C", "#DC050C", "#72190E")

## LOAD VARIABLES
#######################################################
## Load metagenome read stats and summarize




#######################################################
# load metagenome tables (KEGG, pfams, COGs)


tpm_kegg_table = read.table(file = file.path(figs_dir, "tpm_kegg_table.csv"),
  sep = "\t")
tpm_kegg_table = tpm_kegg_table[,2:30]
colnames(tpm_kegg_table)[[2]] = "FI003"
tpm_pfam_table = read.table(file = file.path(figs_dir, "tpm_pfam_table.csv"),
  sep = "\t")
tpm_pfam_table = tpm_pfam_table[,2:30]
colnames(tpm_pfam_table)[[2]] = "FI003"
cno_kegg_table = read.table(file = file.path(figs_dir, "cno_kegg_table.csv"),
  sep = "\t")
cno_kegg_table = cno_kegg_table[,2:30]
colnames(cno_kegg_table)[[2]] = "FI003"
cno_pfam_table = read.table(file = file.path(figs_dir, "cno_pfam_table.csv"),
  sep = "\t")
cno_pfam_table = cno_pfam_table[,2:30]
colnames(cno_pfam_table)[[2]] = "FI003"

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

metadata = read.csv(file.path("/cluster/projects/nn9745k/jing/scripts/metadata/svfi19.tsv"), header=T, check.names=FALSE, sep = "\t", stringsAsFactors=FALSE,colClasses=cls)

metadata_1 = metadata[,c(7, 10, 12:15, 17:33)]

scaled_metadata = scale(metadata_1)


## Parse data
rownames(scaled_metadata) = metadata$sample_id


## merge metagenome tables

## get metadata for metagenomes




# rarefaction curves to check where to set cutoff for sequence reads
tpm_kegg_table_int = round(tpm_kegg_table)
tpm_kegg_table_int[is.na(tpm_kegg_table_int)] <- 0
cairo_ps(file.path(figs_dir, "rarefaction_curves_kegg_tpm.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  rarecurve(t(tpm_kegg_table_int),
                      step = 100, xlab = "Sample Size",
                      ylab = "ASVs",
                      label = FALSE,
                      cex.lab = 2,
                      cex.axis = 1.5)
dev.off()

tpm_pfam_table_int = round(tpm_pfam_table)
tpm_pfam_table_int[is.na(tpm_pfam_table_int)] <- 0
cairo_ps(file.path(figs_dir, "rarefaction_curves_pfam_tpm.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  rarecurve(t(tpm_pfam_table_int),
                      step = 100, xlab = "Sample Size",
                      ylab = "ASVs",
                      label = FALSE,
                      cex.lab = 2,
                      cex.axis = 1.5)
dev.off()

colSums(tpm_kegg_table_int)
colSums(tpm_pfam_table_int)


## based on rarefaction curves cutoffs for tpms were set.

## Make phyloseq objects

# PFAM
row.names(metadata) = metadata$sample_id
dim(tpm_pfam_table_int)
pfam_otu <- otu_table(tpm_pfam_table_int, taxa_are_rows = TRUE)
dim(pfam_otu)
pfam_physeq <- phyloseq(pfam_otu, sample_data(data.frame(metadata)))
pfam_physeq

# KEGG
dim(tpm_kegg_table_int)
kegg_otu <- otu_table(tpm_kegg_table_int, taxa_are_rows = TRUE)
dim(kegg_otu)
kegg_physeq <- phyloseq(kegg_otu, sample_data(data.frame(metadata)))
kegg_physeq

# alpha diversity
min(colSums(tpm_pfam_table_int))
pfam_physeq_20000 <- prune_samples(sample_sums(pfam_physeq)>=min(colSums(tpm_pfam_table_int)),
                                     pfam_physeq)

rarefied_pfam <- rarefy_even_depth(pfam_physeq_20000, sample.size = min(sample_sums(pfam_physeq_20000)),
          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) #Normalizing species data
richness_pfam <- estimate_richness(rarefied_pfam)
min(sample_sums(pfam_physeq_20000))
summary(richness_pfam)

# kegg
min(colSums(tpm_kegg_table_int))
kegg_physeq_1M <- prune_samples(sample_sums(kegg_physeq)>=min(colSums(tpm_kegg_table_int)),
                                    kegg_physeq)

rarefied_kegg <- rarefy_even_depth(kegg_physeq_1M, sample.size = min(sample_sums(kegg_physeq_1M)),
          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) #Normalizing species data
richness_kegg <- estimate_richness(rarefied_kegg)
min(sample_sums(kegg_physeq_1M))
summary(richness_kegg)

# beta diversity
# NMDS - grouping of samples
# KEGGs
otu_kegg <- otu_table(rarefied_kegg)
otu_kegg <- otu_kegg[rowSums(otu_kegg) != 0,]
otu_kegg_stand = decostand(otu_kegg, MARGIN = 1, method="hellinger")

mds_kegg = vegan::metaMDS(t(otu_kegg_stand), distance = "bray", autotransform = FALSE, try = 1000)
mds_kegg_scores = data.frame(vegan::scores(mds_kegg))

# Pfams
otu_pfam <- otu_table(rarefied_pfam)
otu_pfam <- otu_pfam[rowSums(otu_pfam) != 0,]
otu_pfam_stand = decostand(otu_pfam, MARGIN = 1, method="hellinger")

mds_pfam = vegan::metaMDS(t(otu_pfam_stand), distance = "bray", autotransform = FALSE, try = 1000)
mds_pfam_scores = data.frame(vegan::scores(mds_pfam))




# envfit
# kegg
scaled_metadata_kegg = scaled_metadata[rownames(scaled_metadata) %in% colnames(otu_kegg_stand), ]

envfit_results <- vegan::envfit(mds_kegg, data.frame(scaled_metadata_kegg), na.rm = TRUE, permu = 999)
envfit_table_bray <- data.frame(round((envfit_results$vectors)$arrows, 3), round((envfit_results$vectors)$r, 3), round((envfit_results$vectors)$pvals, 3))
colnames(envfit_table_bray) <- c("DIM1", "DIM2", "R", "p")
write.csv(envfit_table_bray, file.path(figs_dir, "Table_envfit_kegg.csv"))

# pfam
scaled_metadata_pfam = scaled_metadata[rownames(scaled_metadata) %in% colnames(otu_pfam_stand), ]

envfit_results <- vegan::envfit(mds_pfam, data.frame(scaled_metadata_pfam), na.rm = TRUE, permu = 999)
envfit_table_bray <- data.frame(round((envfit_results$vectors)$arrows, 3), round((envfit_results$vectors)$r, 3), round((envfit_results$vectors)$pvals, 3))
colnames(envfit_table_bray) <- c("DIM1", "DIM2", "R", "p")
write.csv(envfit_table_bray, file.path(figs_dir, "Table_envfit_pfam.csv"))

# NMDS plot - plot distance from glacier in color and the 4 localities with symbols

# kegg
names(mds_kegg_scores)[c(1, 2)] <- c("x", "y")
mds_kegg_scores$z <- NA


kegg_ordisurf <- ordisurf(mds_kegg ~ scaled_metadata_kegg[,"gl_dist"], plot = FALSE, scaling = 3)
head(kegg_ordisurf)

contour_vals <- extract.xyz(obj = kegg_ordisurf)
head(contour_vals)

p <- ggplot(data = contour_vals, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + theme_bw()


nmds_kegg <- p + geom_text(data = mds_kegg_scores, aes(x = x, y = y, label = rownames(mds_kegg_scores)),
    colour = "red") + coord_equal() + theme_bw() + labs(x = "NMDS1", y = "NMDS2") +
    theme(panel.border = element_rect(fill = NA), axis.text.x = element_blank(),
        axis.text.y = element_blank(), legend.position = "none")


# kegg_heatmap_row <- pheatmap(data.matrix(otu_kegg),
#                             #dendrogram = "row",
#                             xlab = "", ylab = "",
#                             clustering_distance_col = "correlation",
#                             clustering_distance_row = "canberra",
#                             main = "",
#                             scale = "none",
#                             cutree_rows = 7
#                             )
#
# cairo_ps(file.path(figs_dir, "heatmap_kegg.eps"), width=80/25.4, height=80/25.4, # pointsize = 4, bg = FALSE, fallback_resolution = 300)
#   kegg_heatmap_row
# dev.off()

# pfam
names(mds_pfam_scores)[c(1, 2)] <- c("x", "y")
mds_pfam_scores$z <- NA

pfam_ordisurf <- ordisurf(mds_pfam ~ scaled_metadata_pfam[,"gl_dist"], plot = FALSE, scaling = 3)
head(pfam_ordisurf)

contour_vals <- extract.xyz(obj = pfam_ordisurf)
head(contour_vals)

p <- ggplot(data = contour_vals, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + theme_bw()


nmds_pfam <- p + geom_text(data = mds_pfam_scores, aes(x = x, y = y, label = rownames(mds_pfam_scores)),
    colour = "red") + coord_equal() + theme_bw() + labs(x = "NMDS1", y = "NMDS2") +
    theme(panel.border = element_rect(fill = NA), axis.text.x = element_blank(),
        axis.text.y = element_blank(), legend.position = "none")

# NMDS plot - plot distance from glacier in color and the 4 localities with symbols

cairo_ps(file.path(figs_dir, "nmds_kegg.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  nmds_kegg
dev.off()

cairo_ps(file.path(figs_dir, "nmds_pfam.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  nmds_pfam
dev.off()


# RDA on pathways



# procrustes test between pfam and kegg tables

pro <- procrustes(X = mds_kegg, Y = mds_pfam, symmetric = FALSE)
pro

protest(X = mds_kegg, Y = mds_pfam, scores = "sites", permutations = 999)

################## mvabund + gllvm #######################

otu_kegg_mva <- mvabund(otu_kegg)
# otu_pfam_mva <- mvabund(otu_pfam)

# check mean-variance relationship
cairo_ps(file.path(figs_dir, "meanvar_plot_kegg.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  meanvar.plot(otu_kegg_mva)
dev.off()
# cairo_ps(file.path(figs_dir, "meanvar_plot_pfam.eps"), width=80/25.4, height=80/25.4,
#   pointsize = 4, bg = FALSE, fallback_resolution = 300)
#   meanvar.plot(otu_pfam_mva)
# dev.off()

# prep data for gllvm
#kegg
colnames(sample_data(rarefied_kegg))
metadata_matrix_kegg = as.matrix(sample_data(rarefied_kegg)[,c(10,18,20,21,28)])

independent_variables_kegg = scale(metadata_matrix_kegg)
most_abundant_kegg = otu_kegg[rowSums(otu_kegg)>100,]

#pfam
# metadata_matrix_pfam = as.matrix(sample_data(rarefied_pfam)[,c(10,12,13,18, 20,21)])
#
# independent_variables_pfam = scale(metadata_matrix_pfam)
# most_abundant_pfam = otu_pfam[rowSums(otu_pfam)>100,]
# rownames(most_abundant_pfam) = substr(rownames(most_abundant_pfam), 0,7)

# use subset of ASVs - most abundant

# use genera - most abundant
seed_num = 1234
# gllvm_kegg1 = gllvm(t(most_abundant_kegg), family = "negative.binomial")
# gllvm_kegg2 = gllvm(t(most_abundant_kegg), independent_variables_kegg, family = "negative.binomial")
#
# gllvm_pfam1 = gllvm(t(most_abundant_pfam), family = "negative.binomial")
# gllvm_pfam2 = gllvm(t(most_abundant_pfam), independent_variables_pfam, family = "negative.binomial")
#
# cairo_ps(file.path(figs_dir, "gllvm_kegg2_eval.eps"), width=80/25.4, height=80/25.4,
#   pointsize = 4, bg = FALSE, fallback_resolution = 300)
#   par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
#   plot(gllvm_kegg2, var.colors = 1)
# dev.off()
#
# cairo_ps(file.path(figs_dir, "gllvm_pfam2_eval.eps"), width=80/25.4, height=80/25.4,
#   pointsize = 4, bg = FALSE, fallback_resolution = 300)
#   par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
#   plot(gllvm_pfam2, var.colors = 1)
# dev.off()

# criterias <- NULL
# for(i in 1:5){
#   fiti_kegg <- gllvm(t(most_abundant_kegg), independent_variables_kegg,
#                 family = "negative.binomial", num.lv = i, sd.errors = FALSE,
#                 formula = ~ gl_dist + DN + DOC, seed = 1234)
#   criterias[i + 1] <- summary(fiti_kegg)$AICc
#   names(criterias)[i + 1] = i
# }

# fit_env_kegg <- gllvm(t(most_abundant_kegg), independent_variables_kegg,
#                       family = "negative.binomial", num.lv = 4,
#                       formula = ~ gl_dist + DN + DOC,
#                       seed = 1234)
#
# fit_env_pfam <- gllvm(t(most_abundant_pfam), independent_variables_pfam,
#                       family = "negative.binomial", num.lv = 4,
#                       formula = ~ gl_dist + DN + DOC,
#                       seed = 1234)

# cairo_ps(file.path(figs_dir, "gllvm_kegg_coeff.eps"), width=80/25.4, height=80/25.4,
#   pointsize = 4, bg = FALSE, fallback_resolution = 300)
#   coefplot(fit_env_kegg, cex.ylab = 0.7, mar = c(4, 9, 2, 1),
#            xlim.list = list(NULL, NULL, c(-8, 8)), mfrow=c(3,1))
# dev.off()

# cairo_ps(file.path(figs_dir, "gllvm_pfam_coeff.eps"), width=80/25.4, height=80/25.4,
#   pointsize = 4, bg = FALSE, fallback_resolution = 300)
#   coefplot(fit_env_pfam, cex.ylab = 0.7, mar = c(4, 9, 2, 1),
#            xlim.list = list(NULL, NULL, c(-4, 4)), mfrow=c(1,1))
# dev.off()

# KEGG

most_interesting_kegg = c("K00855", # phosphoribulokinase                           Carbon fixation
                          "K01602", # RuBisCO small chain                           Carbon fixation
                          "K08684", # methane monooxygenase                         Aerobic methane oxidation
                          "K02256", # cytochrome c oxidase subunit I (coxI)         Aerobic respiration
                          "K02262", # cytochrome c oxidase subunit III (coxIII)     Aerobic respiration
                          "K02274", # cytochrome c oxidase subunit I (coxA)         Aerobic respiration
                          "K02276", # cytochrome c oxidase subunit III (coxC)"      Aerobic respiration
                          "K01187", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K01199", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K01210", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K01188", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K05349", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K01179", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K19357", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K01225", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K19668", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K01196", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K01176", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K01182", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K01194", # glycoside hydrolase family 31 protein         Glycoside hydrolases
                          "K00174", # 2-oxoglutarate:ferredoxin oxidoreductase subunit alpha  Anaerobic C-fixation
                          "K001752", # 2-oxoglutarate:ferredoxin oxidoreductase subunit beta  Anaerobic C-fixation
                          "K002442", # frdA; fumarate reductase flavoprotein subunit          Anaerobic C-fixation
                          "K016482", # adenosinetriphosphate (ATP) citrate lyase              Anaerobic C-fixation
                          "K001942", # CO dehydrogenase subunit delta                         Anaerobic C-fixation
                          "K001972", # CO dehydrogenase subunit gamma                         Anaerobic C-fixation
                          "K03518",  # CO dehydrogenase small subunit (coxS)                    CO oxidation
                          "K03519",  # cutM, coxM; carbon-monoxide dehydrogenase medium subunit CO oxidation
                          "K03520",  # cutL, coxL; carbon-monoxide dehydrogenase large subunit  CO oxidation
                          "K00016",  # L-lactate dehydrogenase                      Fermentation
                          "K00399",   # mcrA; Methanogenesis
                          "K00400",  # coenzyme M methyl reductase beta subunit (mcrB)  Methanogenesis
                          "K00401",  # methyl coenzyme M reductase system, component A2 Methanogenesis

                          "K03385",  # formate-dependent nitrite reductase periplasmic cytochrome c552 (nrfA) Nitrite Ammonification
                          "K05904",  # cytochrome c nitrite reductase (nrfA)                                  Nitrite Ammonification
                          "K15876",  # nrfH; cytochrome c nitrite reductase small subunit                     Nitrite Ammonification

                          "K10535",  # hydroxylamine oxidoreductase/hydrazine oxidoreducatse (hao/hzo) ANNAMOX

                          "K17877",  # NIT-6; nitrite reductase (NAD(P)H)                   Nitrite reduction
                          "K00362",  # nirB; nitrite reductase (NADH) large subunit         Nitrite reduction
                          "K00363",  # nirD; nitrite reductase (NADH) small subunit         Nitrite reduction
                          "K00366",  # nirA; ferredoxin-nitrite reductase                   Nitrite reduction
                          "K00368",  # nirK; nitrite reductase (NO-forming)                 Nitrite reduction
                          "K15864",  # nirS; nitrite reductase (NO-forming)/hydroxylamine reductase Nitrite reduction

                          "K04561",  # norB; nitric oxide reductase subunit B     Nitric oxide reduction
                          "K02305",  # norC; nitric oxide reductase subunit C     Nitric oxide reduction
                          "K15877",  # CYP55; fungal nitric oxide reductase       Nitric oxide reduction

                          "K00376",  # nosZ; nitrous-oxide reductase              Nitrous oxide reduction

                          "K02567",  # periplasmic nitrate reductase (napA)   Nitrate reduction
                          "K02568",  # cytochrome c-type protein (napB)       Nitrate reduction

                          "K10944",  # ammonia monooxygenase subunit A (amoA) Ammonia oxidation
                          "K10945",  # ammonia monooxygenase subunit B (amoB) Ammonia oxidation
                          "K10946",  # ammonia monooxygenase subunit C (amoC) Ammonia oxidation

                          "K00370",  # narG, narZ, nxrA; nitrate reductase / nitrite oxidoreductase, alpha subunit  nitrate reductase/nitrite oxidation
                          "K00371",  # narH, narY, nxrB; nitrate reductase / nitrite oxidoreductase, beta subunit nitrate reductase/nitrite oxidation
                          "K00367",  # narB; ferredoxin-nitrate reductase (REMOVE FROM narG and narH)

                          "K00531",  # nitrogenase                                            Nitrogen fixation
                          "K02586",  # nitrogenase molybdenum-iron protein alpha chain (nifD) Nitrogen fixation
                          "K02588",  # nitrogenase iron protein (nifH)                        Nitrogen fixation
                          "K02591",  # nitrogenase molybdenum-iron protein beta chain (nifK)  Nitrogen fixation

                          "K00265",  # glutamate synthase (NADPH/NADH) large chain (gltB) Nitrogen assimilation
                          "K00284",  # glutamate synthase (ferredoxin-dependent) (gltS)   Nitrogen assimilation
                          "K00360",  # assimilatory nitrate reductase                     Nitrogen assimilation
                          "K00367",  # assimilatory nitrate reductase                     Nitrogen assimilation
                          "K01915",  # glutamine synthetase (glnA)                        Nitrogen assimilation

                          "K00260",  # glutamate dehydrogenase                Nitrogen Mineralization
                          "K00261",  # glutamate dehydrogenase                Nitrogen Mineralization
                          "K00262",  # glutamate dehydrogenase                Nitrogen Mineralization
# Sulfite oxidation
                          "K00387",  # SUOX; sulfite oxidase
                          "K05301",  # sorA; sulfite dehydrogenase (cytochrome) subunit A
                          "K00386",  # sorB; sulfite dehydrogenase (cytochrome) subunit B
                          "K21307",  # soeA; sulfite dehydrogenase (quinone) subunit SoeA
                          "K21308",  # soeB; sulfite dehydrogenase (quinone) subunit SoeB
                          "K21309",  # soeC; sulfite dehydrogenase (quinone) subunit SoeC

                          "K00860",  # adenylylsulfate kinase (cysC)                  Assimilatory sulfate reduction
                          "K00956",  # sulfate adenylyltransferase subunit 1 (cysN)   Assimilatory sulfate reduction
                          "K00957",  # sulfate adenylyltransferase subunit 2 (cysD)   Assimilatory sulfate reduction
                          "K00394",  # adenylylsulfate reductase subunit A (aprA)     Dissimilatory sulfate reduction/sulfite oxidation
                          "K00395",  # adenylylsulfate reductase subunit B (aprB)     Dissimilatory sulfate reduction/sulfite oxidation
                          "K11180",  # sulfite reductase (dsrA)                       Dissimilatory sulfate reduction/sulfite oxidation
                          "K11181",  # sulfite reductase (dsrB)                       Dissimilatory sulfate reduction/sulfite oxidation
# Sulfite reduction
                          "K16950",  # asrA; anaerobic sulfite reductase subunit A
                          "K16951",  # asrB; anaerobic sulfite reductase subunit B
                          "K00385",  # asrC; anaerobic sulfite reductase subunit C
# Sulfide reduction
                          "K17229",  # fccB; sulfide dehydrogenase [flavocytochrome c] flavoprotein chain
                          "K17230",  # fccA; cytochrome subunit of sulfide dehydrogenase
# Sulfur reduction
                          "K17219",  # sreA; sulfur reductase molybdopterin subunit
                          "K17220",  # sreB; sulfur reductase FeS subunit
                          "K17221",  # sreC; sulfur reductase membrane anchor

                          "K08352",  # phsA, psrA; thiosulfate reductase / polysulfide reductase chain A
                          "K08353",  # phsB; thiosulfate reductase electron transport protein
                          "K08354",  # phsC; thiosulfate reductase cytochrome b subunit
# Sulfur oxidation
                          "K16952",  # sor; sulfur oxygenase/reductase
# Thiosulfate oxidation
                          "K17222",  # soxA; L-cysteine S-thiosulfotransferase
                          "K17223",  # soxX; L-cysteine S-thiosulfotransferase
                          "K17226",  # soxY; sulfur-oxidizing protein SoxY
                          "K17227",  # soxZ; sulfur-oxidizing protein SoxZ
                          "K17224",  # soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase
                          "K17225",  # soxC; sulfane dehydrogenase subunit SoxC
                          "K22622",  # soxD; S-disulfanyl-L-cysteine oxidoreductase SoxD

# Synthesis of oxygenic-photosynthetic reaction center
                          "K02703",  # psbA; photosystem II P680 reaction center D1 protein [EC:1.10.3.9]
                          "K02704",  # psbB; photosystem II CP47 chlorophyll apoprotein
                          "K02705",  # psbC; photosystem II CP43 chlorophyll apoprotein
                          "K02706",  # psbD; photosystem II P680 reaction center D2 protein [EC:1.10.3.9]
                          "K02707",  # psbE; photosystem II cytochrome b559 subunit alpha
                          "K02708",  # psbF; photosystem II cytochrome b559 subunit beta
                          "K02689",  # psaA; photosystem I P700 chlorophyll a apoprotein A1
                          "K02690",  # psaB; photosystem I P700 chlorophyll a apoprotein A2
                          "K02691",  # psaC; photosystem I subunit VII
                          "K02692",  # psaD; photosystem I subunit II
                          "K02693",  # psaE; photosystem I subunit IV
                          "K02694",  # psaF; photosystem I subunit III
# Synthesis of anoxygenic-photosynthetic reaction center
                          # "K08940",  # pscA; photosystem P840 reaction center large subunit
                          # "K08941",  # pscB; photosystem P840 reaction center iron-sulfur protein
                          # "K08942",  # pscC; photosystem P840 reaction center cytochrome c551
                          # "K08943",  # pscD; photosystem P840 reaction center protein PscD
                          "K08928",  # pufL; photosynthetic reaction center L subunit
                          "K08929",  # pufM; photosynthetic reaction center M subunit
# Bacteriorhodopsin
                          "K04641",   # bop; bacteriorhodopsin
# Methanotrophy
                          "K10944",   # pmoA; ammonia/methane monooxygenase component A
                          "K10945",   # pmoB; ammonia/methane monooxygenase component B
                          "K10946",   # pmoC; ammonia/methane monooxygenase component C
                          "K16157",   # mmoX; methane monooxygenase component A alpha chain
                          "K16256",   # mxaA; mxaA protein
                          "K23995",   # xoxF; lanthanide-dependent methanol dehydrogenase
                          "K22516",   # fdhA; formate dehydrogenase (coenzyme F420) alpha subunit
                          "K05299",   # fdhA; formate dehydrogenase (NADP+) alpha subunit
                          "K00148",   # fdhA; glutathione-independent formaldehyde dehydrogenase
                          "K10713",   # fae; 5,6,7,8-tetrahydromethanopterin hydrolyase
                          "K13812",   # fae; methanotrophy
                          "K00300",   # mtdA; methylenetetrahydrofolate/methylenetetrahydromethanopterin dehydrogenase (NADP+)
                          "K10714",   # mtdB; methylene-tetrahydromethanopterin dehydrogenase
                          "K01499",   # mch; methenyltetrahydromethanopterin cyclohydrolase
                          "K00672",   # ftr; formylmethanofuran--tetrahydromethanopterin N-formyltransferase
                          "K00200",   # fwdA, fmdA; formylmethanofuran dehydrogenase subunit A
                          "K00201",   # fwdB, fmdB; formylmethanofuran dehydrogenase subunit B
                          "K00202",   # fwdC, fmdC; formylmethanofuran dehydrogenase subunit C
                          "K00203",   # fwdD, fmdD; formylmethanofuran dehydrogenase subunit D
                          "K00204",   # fwdH; 4Fe-4S ferredoxin
                          "K00205",   # fwdF, fmdF; 4Fe-4S ferredoxin
                          "K11260",   # fwdG; 4Fe-4S ferredoxin
                          "K11261",   # fwdE, fmdE; formylmethanofuran dehydrogenase subunit E
                          "K01500",   # fchA; methenyltetrahydrofolate cyclohydrolase
                          "K01938",   # fhs; formate--tetrahydrofolate ligase
                          "K03396",   # gfa; S-(hydroxymethyl)glutathione synthase
                          "K00121",   # frmA, ADH5, adhC; S-(hydroxymethyl)glutathione dehydrogenase/alcohol dehydrogenase
                          "K01070",   # frmB, ESD, fghA; S-formylglutathione hydrolase
                          "K00122",   # FDH; formate dehydrogenase
                          "K00125",   # fdhB; formate dehydrogenase (coenzyme F420) beta subunit
                          "K15022",   # fdhB; formate dehydrogenase (NADP+) beta subunit
                          "K22015"    # fdhF; formate dehydrogenase (hydrogenase)

                          )
# most important KEGGs analysis
kegg_data = data.frame(otu_kegg)
kegg_selection = kegg_data[rownames(kegg_data) %in% most_interesting_kegg,]

fit_env_kegg_sel <- gllvm(t(kegg_selection), independent_variables_kegg,
                      family = "negative.binomial", num.lv = 4,
                      formula = ~ gl_dist + DN + DOC,
                      seed = 1234)

cairo_ps(file.path(figs_dir, "gllvm_kegg_sel_eval.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  plot(fit_env_kegg_sel, cex.ylab = 0.7)
dev.off()

cairo_ps(file.path(figs_dir, "gllvm_kegg_sel.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  coefplot(fit_env_kegg_sel, cex.ylab = 0.3)
dev.off()

fit_env_kegg_sel2 <- gllvm(t(kegg_selection), independent_variables_kegg,
                      family = "negative.binomial", num.lv = 1,
                      formula = ~ gl_dist + DN + DOC + TP + bird,
                      seed = 1234)

cairo_ps(file.path(figs_dir, "gllvm_kegg_sel_eval2.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  plot(fit_env_kegg_sel2, cex.ylab = 0.7)
dev.off()

cairo_ps(file.path(figs_dir, "gllvm_kegg_sel2.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  coefplot(fit_env_kegg_sel2, cex.ylab = 0.3)
dev.off()


### Comparison between ASV tables and functional tables
otu_euk_stand = readRDS(file = file.path(figs_dir, "otu_euk_stand.rds"))
otu_bac_stand = readRDS(file = file.path(figs_dir, "otu_bac_stand.rds"))

# mds_bac = readRDS(file = file.path(figs_dir, "mds_bac.rds"))
# mds_euk = readRDS(file = file.path(figs_dir, "mds_euk.rds"))
# KEGG
otu_kegg_euk = otu_kegg_stand[,colnames(otu_kegg_stand) %in% colnames(otu_euk_stand)]
otu_euk_kegg = otu_euk_stand[,colnames(otu_euk_stand) %in% colnames(otu_kegg_stand)]

otu_kegg_bac = otu_kegg_stand[,colnames(otu_kegg_stand) %in% colnames(otu_bac_stand)]
otu_bac_kegg = otu_bac_stand[,colnames(otu_bac_stand) %in% colnames(otu_kegg_stand)]

# EUK

mds_kegg_euk = vegan::metaMDS(t(otu_kegg_euk), distance = "bray", autotransform = FALSE, try = 1000)
mds_euk_kegg = vegan::metaMDS(t(otu_euk_kegg), distance = "bray", autotransform = FALSE, try = 1000)
dist_kegg_euk = vegan::vegdist(t(otu_kegg_euk), method = "bray")
dist_euk_kegg = vegan::vegdist(t(otu_euk_kegg), method = "bray")

pro <- procrustes(X = mds_kegg_euk, Y = mds_euk_kegg, symmetric = FALSE)
pro

protest(X = mds_kegg_euk, Y = mds_euk_kegg, scores = "sites", permutations = 999)
mantel(dist_kegg_euk, dist_euk_kegg, method = "spearman", permutations = 9999, na.rm = TRUE)

# BAC
mds_kegg_bac = vegan::metaMDS(t(otu_kegg_bac), distance = "bray", autotransform = FALSE, try = 1000)
mds_bac_kegg = vegan::metaMDS(t(otu_bac_kegg), distance = "bray", autotransform = FALSE, try = 1000)
dist_kegg_bac = vegan::vegdist(t(otu_kegg_bac), method = "bray")
dist_bac_kegg = vegan::vegdist(t(otu_bac_kegg), method = "bray")


pro <- procrustes(X = mds_kegg_bac, Y = mds_bac_kegg, symmetric = FALSE)
pro

protest(X = mds_kegg_bac, Y = mds_bac_kegg, scores = "sites", permutations = 999)
mantel(dist_kegg_bac, dist_bac_kegg, method = "spearman", permutations = 9999, na.rm = TRUE)


# PFAM
otu_pfam_euk = otu_pfam_stand[,colnames(otu_pfam_stand) %in% colnames(otu_euk_stand)]
otu_euk_pfam = otu_euk_stand[,colnames(otu_euk_stand) %in% colnames(otu_pfam_stand)]

otu_pfam_bac = otu_pfam_stand[,colnames(otu_pfam_stand) %in% colnames(otu_bac_stand)]
otu_bac_pfam = otu_bac_stand[,colnames(otu_bac_stand) %in% colnames(otu_pfam_stand)]

# EUK
mds_pfam_euk = vegan::metaMDS(t(otu_pfam_euk), distance = "bray", autotransform = FALSE, try = 1000)
mds_euk_pfam = vegan::metaMDS(t(otu_euk_pfam), distance = "bray", autotransform = FALSE, try = 1000)
dist_pfam_euk = vegan::vegdist(t(otu_pfam_euk), method = "bray")
dist_euk_pfam = vegan::vegdist(t(otu_euk_pfam), method = "bray")


pro <- procrustes(X = mds_pfam_euk, Y = mds_euk_pfam, symmetric = FALSE)
pro

protest(X = mds_pfam_euk, Y = mds_euk_pfam, scores = "sites", permutations = 999)
mantel(dist_pfam_euk, dist_euk_pfam, method = "spearman", permutations = 9999, na.rm = TRUE)

#BAC
mds_pfam_bac = vegan::metaMDS(t(otu_pfam_bac), distance = "bray", autotransform = FALSE, try = 1000)
mds_bac_pfam = vegan::metaMDS(t(otu_bac_pfam), distance = "bray", autotransform = FALSE, try = 1000)
dist_pfam_bac = vegan::vegdist(t(otu_pfam_bac), method = "bray")
dist_bac_pfam = vegan::vegdist(t(otu_bac_pfam), method = "bray")


pro <- procrustes(X = mds_pfam_bac, Y = mds_bac_pfam, symmetric = FALSE)
pro

protest(X = mds_pfam_bac, Y = mds_bac_pfam, scores = "sites", permutations = 999)
mantel(dist_pfam_bac, dist_bac_pfam, method = "spearman", permutations = 9999, na.rm = TRUE)

## Look at variability in each chronosequence
# Plot PFAMs of importance in xy plot
# PFAM

# meta_pfam = data.frame(sample_data(rarefied_pfam))
# cq = unique(data.frame(sample_data(rarefied_pfam))$locality)[1]
# data_pfam = data.frame(otu_pfam)
# sub_metadata = data.frame(sample_data(rarefied_pfam))
#
# for (cq in unique(data.frame(sample_data(rarefied_pfam))$locality)) {
#
#   pfam_sel = data_pfam[,colnames(data_pfam) %in% rownames(sub_metadata[sub_metadata$locality==cq,]) ]
#   rownames(pfam_sel) = substr(rownames(pfam_sel), 0,7)
#   pfam_selection = pfam_sel[rownames(pfam_sel) %in% most_interesting_pfam, ]
#
# }

# reads per 21480 sample pfam

# KEGG

meta_kegg = data.frame(sample_data(rarefied_kegg))
cq = unique(data.frame(sample_data(rarefied_kegg))$locality)[1]
data_kegg = data.frame(otu_kegg)
sub_metadata = data.frame(sample_data(rarefied_kegg))

for (cq in unique(data.frame(sample_data(rarefied_kegg))$locality)) {

  kegg_sel = data_kegg[,colnames(data_kegg) %in% rownames(sub_metadata[sub_metadata$locality==cq,]) ]
  rownames(kegg_sel)[1]
  kegg_selection = kegg_sel[rownames(kegg_sel) %in% most_interesting_kegg, ]

}

# Write function that calculates pathway contribution from KEGG reads coverage
traits_name = c("Carbon_fixation",
                "Aerobic_methane_oxidation",
                "Aerobic_respiration",
                "Anaerobic_C_fixation",
                "Glycoside hydrolases",
                "CO_oxidation",
                "Fermentation",
                "Methanogenesis",
                "Methanotrophy",
                # Nitrogen cycle
                "Ammonification",
                "ANAMMOX",
                "Nitrite_red",
                "Nitrous_oxide_red",
                "Nitric_oxide_red",
                "Nitrate_reduction",
                "Ammonia_oxid",
                "N_fixation",
                "Nitrite_oxid",
                # Sulfur cycle
                "Sulfite_oxid",
                "Sulfite_red",
                "Sulfur_red",
                "Sulfide_oxid",
                "Sulfur_oxid",
                "Sulfate_red",
                "Thiosulfate_oxid",
                "Thiosulfate_disp1",
                "Thiosulfate_disp2",
                "Polysulfide_reduction",
                # Photosynthesis
                "Oxy_photosynthesis",
                "Anoxy_photosynthesis",
                "Bacteriorhodopsin"


                )

# Function to sum up traits from KEGGs

summarize_kegg_traits = function(kegg_selection){
  trait_table = data.frame(matrix(ncol=ncol(kegg_selection),
                          nrow=length(traits_name),
                          dimnames=list(traits_name, colnames(kegg_selection))))

  # Calculate
  Carbon_fixation = colSums(kegg_selection["K00855",] + kegg_selection["K01602",], na.rm = TRUE)/2
  trait_table["Carbon_fixation",] = Carbon_fixation

  Aerobic_methane_oxidation = kegg_selection["K08684",]
  trait_table["Aerobic_methane_oxidation",] = Aerobic_methane_oxidation

  Aerobic_respiration = colSums(kegg_selection["K02256",] + kegg_selection["K02262",], na.rm = TRUE)/2 +
                        colSums(kegg_selection["K02274",] + kegg_selection["K02276",], na.rm = TRUE)/2
  trait_table["Aerobic_respiration",] = Aerobic_respiration

  Anaerobic_C_fixation = colSums(kegg_selection["K00174",] + kegg_selection["K00175",] +
                                 kegg_selection["K00244",] + kegg_selection["K01648",], na.rm = TRUE)/4 +
                        colSums(kegg_selection["K00194",] + kegg_selection["K00197",], na.rm = TRUE)/2
  trait_table["Anaerobic_C_fixation",] = Anaerobic_C_fixation

  Glycoside_hydrolases = colSums(kegg_selection["K01187",] + kegg_selection["K01199",] +
                                 kegg_selection["K01210",] + kegg_selection["K01188",] +
                                 kegg_selection["K05349",] + kegg_selection["K01179",] +
                                 kegg_selection["K19357",] + kegg_selection["K01225",] +
                                 kegg_selection["K19668",] + kegg_selection["K01196",] +
                                 kegg_selection["K01176",] + kegg_selection["K01182",] +
                                 kegg_selection["K01194",], na.rm = TRUE)
  trait_table["Glycoside_hydrolases",] = Glycoside_hydrolases

  CO_oxidation = colSums(kegg_selection["K03518",] + kegg_selection["K03519",] +
                                 kegg_selection["K03520",] , na.rm = TRUE)/3
  trait_table["CO_oxidation",] = CO_oxidation

  Fermentation = kegg_selection["K00016",]
  trait_table["Fermentation",] = Fermentation

  Methanogenesis = colSums(kegg_selection["K00399",] + kegg_selection["K00400",] +
                           kegg_selection["K00401",], na.rm = TRUE)/3
  trait_table["Methanogenesis",] = Methanogenesis

  Methanotrophy = colSums(kegg_selection["K10944",] + kegg_selection["K10945",] +
                          kegg_selection["K10946",] + kegg_selection["K16157",] +
                          kegg_selection["K16256",] + kegg_selection["K23995",] +
                          kegg_selection["K22516",] + kegg_selection["K05299",] +
                          kegg_selection["K00148",] + kegg_selection["K10713",] +
                          kegg_selection["K13812",] + kegg_selection["K00300",] +
                          kegg_selection["K10714",] + kegg_selection["K01499",] +
                          kegg_selection["K00672",] + kegg_selection["K00200",] +
                          kegg_selection["K00201",] + kegg_selection["K00202",] +
                          kegg_selection["K00203",] + kegg_selection["K00204",] +
                          kegg_selection["K00205",] + kegg_selection["K11260",] +
                          kegg_selection["K11261",] + kegg_selection["K01500",] +
                          kegg_selection["K01938",] + kegg_selection["K03396",] +
                          kegg_selection["K00121",] + kegg_selection["K01070",] +
                          kegg_selection["K00122",] + kegg_selection["K00125",] +
                          kegg_selection["K15022",] + kegg_selection["K22015",], na.rm = TRUE)/32
  trait_table["Methanotrophy",] = Methanotrophy

  Ammonification = colSums(kegg_selection["K03385",] + kegg_selection["K05904",] +
                           kegg_selection["K04561",], na.rm = TRUE)/3
  trait_table["Ammonification",] = Ammonification

  ANAMMOX = kegg_selection["K10535",]
  trait_table["ANAMMOX",] = ANAMMOX

  Nitrite_red = colSums(kegg_selection["K17877",], na.rm = TRUE) +
                colSums(kegg_selection["K00362",] + kegg_selection["K00363",] +
                        kegg_selection["K00366",] + kegg_selection["K00368",] +
                        kegg_selection["K15864",], na.rm = TRUE)/5
  trait_table["Nitrite_red",] = Nitrite_red

  Nitrous_oxide_red = colSums(kegg_selection["K00376",], na.rm = TRUE)
  trait_table["Nitrous_oxide_red",] = Nitrous_oxide_red

  Nitric_oxide_red = colSums(kegg_selection["K04561",] + kegg_selection["K02305",], na.rm = TRUE)/2 +
                     colSums(kegg_selection["K15877",] , na.rm = TRUE)
  trait_table["Nitric_oxide_red",] = Nitric_oxide_red

  Nitrate_red = colSums(kegg_selection["K02567",] + kegg_selection["K02568",], na.rm = TRUE)/2 +
                colSums(kegg_selection["K00367",])
  trait_table["Nitrate_red",] = Nitrate_red

  Ammonia_oxid = colSums(kegg_selection["K10944",] + kegg_selection["K10945",] +
                         kegg_selection["K10946",], na.rm = TRUE)/3
  trait_table["Ammonia_oxid",] = Ammonia_oxid

  N_fixation = colSums(kegg_selection["K00531",] + kegg_selection["K02586",] +
                       kegg_selection["K02588",] + kegg_selection["K02591",], na.rm = TRUE)/4
  trait_table["N_fixation",] = N_fixation

  Nitrite_oxid = colSums(kegg_selection["K00370",] + kegg_selection["K00371",], na.rm = TRUE)/2
  trait_table["Nitrite_oxid",] = Nitrite_oxid

  Sulfite_oxid = colSums(kegg_selection["K00387",], na.rm = TRUE) +
                 colSums(kegg_selection["K05301",] + kegg_selection["K00386",], na.rm = TRUE)/2 +
                 colSums(kegg_selection["K21307",] + kegg_selection["K21308",] +
                         kegg_selection["K21309",], na.rm = TRUE)/3
  trait_table["Nitrite_oxid",] = Nitrite_oxid

  Sulfite_red = colSums(kegg_selection["K16950",] + kegg_selection["K16951",] +
                        kegg_selection["K00385",], na.rm = TRUE)/3
  trait_table["Sulfite_red",] = Sulfite_red

  Sulfur_red = colSums(kegg_selection["K16952",] + kegg_selection["K17219",] +
                       kegg_selection["K17220",] + kegg_selection["K17221",], na.rm = TRUE)/4
  trait_table["Sulfur_red",] = Sulfur_red

  Sulfide_oxid = colSums(kegg_selection["K17229",] + kegg_selection["K17230",], na.rm = TRUE)/2
  trait_table["Sulfide_pxid",] = Sulfide_oxid

  Sulfur_oxid = colSums(kegg_selection["K16952",], na.rm = TRUE)
  trait_table["Sulfur_oxid",] = Sulfur_oxid

  Sulfate_red = colSums(kegg_selection["K00860",] + kegg_selection["K00956",] +
                        kegg_selection["K00957",], na.rm = TRUE)/3 +
                colSums(kegg_selection["K00394",] + kegg_selection["K00395",] +
                        kegg_selection["K11180",], na.rm = TRUE)/3 -
                Sulfite_oxid
  trait_table["Sulfate_red",] = Sulfate_red

  Thiosulfate_oxid = colSums(kegg_selection["K17222",] + kegg_selection["K17223",] +
                             kegg_selection["K17226",] + kegg_selection["K17227",] +
                             kegg_selection["K17224",] + kegg_selection["K17225",] +
                             kegg_selection["K22622",], na.rm = TRUE)/7
  trait_table["Thiosulfate_oxid",] = Thiosulfate_oxid

#  Thiosulfate_disp1 =
#  trait_table["Thiosulfate_disp1",] = Thiosulfate_disp1
#
#  Thiosulfate_disp2 =
#  trait_table["Thiosulfate_disp2",] = Thiosulfate_disp2

  Polysulfide_reduction = colSums(kegg_selection["K08352",], na.rm = TRUE)
  trait_table["Polysulfide_reduction",] = Polysulfide_reduction

  Oxy_photosynthesis = colSums(kegg_selection["K02703",] + kegg_selection["K02704",] +
                               kegg_selection["K02705",] + kegg_selection["K02706",] +
                               kegg_selection["K02707",] + kegg_selection["K02708",], na.rm = TRUE)/12 +
                       colSums(kegg_selection["K02689",] + kegg_selection["K02690",] +
                               kegg_selection["K02691",] + kegg_selection["K02692",] +
                               kegg_selection["K02693",] + kegg_selection["K02694",], na.rm = TRUE)/12
  trait_table["Oxy_photosynthesis",] = Oxy_photosynthesis

  Anoxy_photosynthesis = colSums(kegg_selection["K08940",] +
                                 kegg_selection["K08941",] +
                                 kegg_selection["K08942",] +
                                 kegg_selection["K08943",], na.rm = TRUE)/6 +
                         colSums(kegg_selection["K08928",] +
                                 kegg_selection["K08929",], na.rm = TRUE)/2
  trait_table["Anoxy_photosynthesis",] = Anoxy_photosynthesis

  Bacteriorhodopsin = kegg_selection["K04641",]
  trait_table["Bacteriorhodopsin",] = Bacteriorhodopsin
  return(trait_table)

}

traits_table = summarize_kegg_traits(tpm_kegg_table)
traits_table = traits_table[rowSums(is.na(traits_table)) != ncol(traits_table), ]
traits_table[is.na(traits_table)] <- 0


fit_env_traits_full <- gllvm(t(traits_table), independent_variables_kegg,
                            family = "negative.binomial", num.lv = 1,
                            formula = ~ gl_dist + DN + DOC,
                            seed = 1234)

fit_traits_full_confin = confint(fit_env_traits_full, level = 0.95, parm = "Xcoef")
fit_split_traits_confin = cbind(fit_traits_full_confin[1:27,],
                                fit_traits_full_confin[28:54,],
                                fit_traits_full_confin[55:81,])
rownames(fit_split_traits_confin) = rownames(traits_table)
col.names = c("gl_dist_cilow", "gl_dist_ciup","DN_cilow", "DN_ciup", "TOC_cilow", "TOC_ciup")
colnames(fit_split_traits_confin) = col.names

cairo_ps(file.path(figs_dir, "gllvm_traits_full_eval.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  plot(fit_env_traits_full, cex.ylab = 0.7)
dev.off()

cairo_ps(file.path(figs_dir, "gllvm_traits_full_coeff.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  coefplot(fit_env_traits_full, cex.ylab = 0.8)
dev.off()

fit_env_traits_full2 <- gllvm(t(traits_table), independent_variables_kegg,
                            family = "negative.binomial", num.lv = 1,
                            formula = ~ gl_dist + DN + DOC + TP + bird,
                            seed = 1234)

cairo_ps(file.path(figs_dir, "gllvm_traits_full_eval2.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  plot(fit_env_traits_full2, cex.ylab = 0.7)
dev.off()

cairo_ps(file.path(figs_dir, "gllvm_traits_full_coeff2.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  coefplot(fit_env_traits_full2, cex.ylab = 0.8)
dev.off()




## heatmaps

traits_table1 = summarize_kegg_traits(tpm_kegg_table)
traits_table1[is.na(traits_table1)] <- 0


traits_table1 = traits_table1[which(rowMeans(!traits_table1 == 0) > 1/29),]
traits_table1 = scale(t(traits_table1))

traits_heatmap_row <- pheatmap(data.matrix(traits_table1),
                            #dendrogram = "row",
                            xlab = "", ylab = "",
                            clustering_distance_col = "correlation",
                            clustering_distance_row = "canberra",
                            main = "",
                            scale = "none",
                            cutree_rows = 7
                            )

cairo_ps(file.path(figs_dir, "traits_heatmap_row.eps"))
  traits_heatmap_row

dev.off()

kegg_table1 = tpm_kegg_table
kegg_table1[is.na(kegg_table1)] <- 0


kegg_table1 = kegg_table1[which(rowMeans(!kegg_table1 == 0) > 1/29),]

kegg_table1 = scale(t(kegg_table1))

kegg_heatmap_row <- pheatmap(data.matrix(kegg_table1),
                            treeheight_col = 0,
                            xlab = "", ylab = "",
                            clustering_distance_col = "correlation",
                            clustering_distance_row = "canberra",
                            main = "",
                            scale = "none",
                            cutree_rows = 7
                            )

cairo_ps(file.path(figs_dir, "kegg_heatmap_row.eps"))
 kegg_heatmap_row

dev.off()


######### Plot taxonomy based on metagenomic reads ###############

genus_abund_table = read.table(file = file.path(figs_dir, "genus_abund_table.tsv"), sep = "\t")

# Plot 3 domains

arch_abund_table  = genus_abund_table[grep("k_Archaea", rownames(genus_abund_table)),]
bac_abund_table  = genus_abund_table[grep("k_Bacteria", rownames(genus_abund_table)),]
euk_abund_table  = genus_abund_table[grep("k_Eukaryota", rownames(genus_abund_table)),]
vir_abund_table  = genus_abund_table[grep("k_Viruses", rownames(genus_abund_table)),]
colSums(vir_abund_table, na.rm = TRUE)

domain_table = rbind(colSums(arch_abund_table, na.rm = TRUE),
                     colSums(bac_abund_table, na.rm = TRUE),
                     colSums(euk_abund_table, na.rm = TRUE))
rownames(domain_table) = c("Archaea", "Bacteria", "Eukaryota")
domain_table_prop = prop.table(domain_table, margin = 2)*100

# calculate proportions based on phylum data
melted_domain_table = melt(domain_table_prop)
melted_domain_table <- dplyr::arrange(melted_domain_table, Var2, desc(value))
melted_domain_table$Var2 <- factor(melted_domain_table$Var2 , levels = unique(melted_domain_table$Var2))

class_colors <- setNames(color_palette_euk, levels(melted_domain_table$Var1))
taxonomy = rownames(domain_table)

svg(file.path(figs_dir, "boxplot_meta_domain.svg"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE)
  ggplot(melted_domain_table, aes(x = Var2, y = value, fill = Var1))+
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
     scale_x_discrete(limits = rev(levels(melted_domain_table$Var2))
     )
dev.off()

summary(t(domain_table_prop))



save.image(file.path(figs_dir, "metags_stats.Rdata"))



