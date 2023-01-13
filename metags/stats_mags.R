### This is an R script to analyse to infer greenhouse gas emissions and metagenomic
### data from Svalbard


## saga
# module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
# path to metags

path_in = file.path("/cluster/projects/nn9744k/02_results/13_svalbard_metaGs")

figs_dir <- file.path(path_in, "figs")

# load(file.path(figs_dir, "metags_stats.Rdata"))


#In R
## LOAD PACKAGES
.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_0_0", .libPaths()))

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


library(plyr)
library(ComplexHeatmap)
library(phyloseq)


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
# load data


checkm = read.table(file = file.path(path_in, "17.coassembled.checkM"),
  sep = "\t", header = FALSE, fill = TRUE)
toDelete = seq(1, nrow(checkm), 2)
checkm = checkm[-toDelete, ]
checkm_colnames = c("Bin_Id", "Marker_lineage", "no_genomes",
                    "markers", "marker_sets", "0", "1", "2", "3",
                    "4", "5+", "Completeness", "Contamination",
                    "Strain_heterogeneity")
colnames(checkm) = checkm_colnames

bincov = read.table(file = file.path(path_in, "18.coassembled.bincov"),
  sep = "\t")
bintable = read.table(file = file.path(path_in, "18.coassembled.bintable"),
  sep = "\t", header = TRUE)
contigsinbins = read.table(file = file.path(path_in, "18.coassembled.contigsinbins"),
  sep = "\t")
contigtable = read.table(file = file.path(path_in, "19.coassembled.contigtable"),
  sep = "\t", header = TRUE)
kegg_pathways = read.table(file = file.path(path_in, "20.coassembled.kegg.pathways"), check.names = FALSE,
  sep = "\t")
kegg_tax = read.table(file = file.path(path_in, "kegg_tax.txt"), check.names = FALSE, header = FALSE,
  sep = "\n")

metacyc_pathways = read.table(file = file.path(path_in, "20.coassembled.metacyc.pathways"), check.names = FALSE, header = FALSE,
  sep = "\t")

metacyc_tax = read.table(file = file.path(path_in, "metacyc_tax.txt"), check.names = FALSE, header = FALSE,
  sep = "\n")
intres_metacyc_tax = read.table(file = file.path(path_in, "interes_metacyc.txt"), check.names = FALSE, header = FALSE,
  sep = "\n")


# metacyc_pathways = read.table(file = file.path(path_in, "20.coassembled.metacyc.pathways"), check.names = FALSE,
#  sep = "\t")
# coassembled_stats = read.table(file = file.path(path_in, "21.coassembled.stats"),
#  sep = "\t")

kegg_contigs = read.table(file = file.path(path_in, "results/07.coassembled.fun3.kegg"), sep = "\t")
colnames(kegg_contigs) = c("Contig.ID", "besthit", "bestaver")
kegg_contigs$Contig.ID = sub("^([^_]*_[^_]*).*", "\\1", kegg_contigs$Contig.ID)


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

metadata = read.csv(file.path("/cluster/home/alexaei/figs", "svfi19.tsv"),
  header=T, check.names=FALSE, sep = "\t", stringsAsFactors=FALSE,colClasses=cls)
rownames(metadata) = metadata$sample_id

metadata_1 = metadata[,c(7, 10, 12:15, 17:33)]

scaled_metadata = scale(metadata_1)


## Parse data
rownames(scaled_metadata) = metadata$sample_id

## get some MAG stats

nrow(checkm[checkm$Completeness > 70,])
nrow(checkm[checkm$Contamination < 5,])
nrow(checkm[checkm$Strain_heterogeneity < 1,])

# good_bins = checkm[checkm$Completeness > 70 & checkm$Contamination < 5,]

compl_bintable = bintable[complete.cases(bintable),]


good_bintable = compl_bintable[compl_bintable$Completeness > 70 & compl_bintable$Contamination < 5,]

good_contigtable = contigtable[contigtable$"Bin.ID" %in% good_bintable$Bin.ID,]

good_kegg_contigs = kegg_contigs[kegg_contigs$Contig.ID %in% good_contigtable$Contig.ID,]

# intres_bins = checkm[checkm$Completeness > 30 & checkm$Contamination < 5,]

intres_bintable = compl_bintable[compl_bintable$Completeness > 50 & compl_bintable$Contamination < 5,]

intres_contigtable = contigtable[contigtable$"Bin.ID" %in% intres_bintable$Bin.ID,]

intres_kegg_contigs = kegg_contigs[kegg_contigs$Contig.ID %in% intres_contigtable$Contig.ID,]

# bin stats
selection = seq(12, 68, by = 2)

bins_coverage_per_sample = colSums(intres_bintable[,1+selection])/10000
bin_coverage = rowSums(intres_bintable[,1+selection])/10000/length(selection)
names(bin_coverage) = intres_bintable$Bin.ID

intres_bintable_stats = cbind(intres_bintable[,c(1,3:11)],bin_coverage)
write.table(intres_bintable_stats, file.path("/cluster/home/alexaei/figs", "intres_bintable_stats.csv"), sep = "\t", row.names = FALSE)

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
kegg_selection = good_kegg_contigs[good_kegg_contigs$besthit %in% most_interesting_kegg,]

contigtable_kegg_selection = good_contigtable[good_contigtable$Contig.ID %in% kegg_selection$Contig.ID, ]

bin_kegg_selection = join(contigtable_kegg_selection, kegg_selection, type="full", by = "Contig.ID")

bin_kegg_selection$avTPM = rowMeans(bin_kegg_selection[,seq(6, 93, by = 3)])

write.table(bin_kegg_selection, file.path("/cluster/projects/nn9745k/jing/02_results/svalbard", "bin_kegg_selection.csv"), sep = "\t", row.names = FALSE)

# interesting part

intres_kegg_selection = intres_kegg_contigs[intres_kegg_contigs$besthit %in% most_interesting_kegg,]

intres_contigtable_kegg_selection = intres_contigtable[intres_contigtable$Contig.ID %in% intres_kegg_selection$Contig.ID, ]

intres_bin_kegg_selection = join(intres_contigtable_kegg_selection, intres_kegg_selection, type="full", by = "Contig.ID")

intres_bin_kegg_selection$avTPM = rowMeans(intres_bin_kegg_selection[,seq(6, 93, by = 3)])
colSums(intres_bin_kegg_selection[,seq(6, 93, by = 3)])



hydrolases =  c("K01187", # glycoside hydrolase family 31 protein         Glycoside hydrolases
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
                "K01194"  # glycoside hydrolase family 31 protein         Glycoside hydrolases
                )
intres_kegg_hydrolases = intres_kegg_contigs[intres_kegg_contigs$besthit %in% hydrolases,]

intres_contigtable_hydrolases = intres_contigtable[intres_contigtable$Contig.ID %in% intres_kegg_hydrolases$Contig.ID, ]

intres_bin_hydrolases = join(intres_contigtable_hydrolases, intres_kegg_hydrolases, type="full", by = "Contig.ID")

intres_bin_hydrolases$avTPM = rowMeans(intres_bin_hydrolases[,seq(6, 93, by = 3)])


methanotroph =  c("K10944",   # pmoA; ammonia/methane monooxygenase component A
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

intres_kegg_methanotroph = intres_kegg_contigs[intres_kegg_contigs$besthit %in% methanotroph,]

intres_contigtable_methanotroph = intres_contigtable[intres_contigtable$Contig.ID %in% intres_kegg_methanotroph$Contig.ID, ]

intres_bin_methanotroph = join(intres_contigtable_methanotroph, intres_kegg_methanotroph, type="full", by = "Contig.ID")

intres_bin_methanotroph$avTPM = rowMeans(intres_bin_methanotroph[,seq(6, 93, by = 3)])



# PLOT KEGG SELECTION PER BIN

# cairo_ps(file.path("/cluster/projects/nn9745k/jing/02_results/svalbard", "bin_kegg.eps"), width=80/25.4, height=80/25.4,
cairo_ps(file.path("/cluster/home/alexaei/figs", "bin_kegg.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  ggplot(bin_kegg_selection, aes(Bin.ID, besthit, fill= avTPM)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text = element_text(size=6),
    axis.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=8),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.line.x = element_line(color = "grey")) +
  scale_fill_distiller(palette = "RdPu") +
  geom_tile()


dev.off()

# KEGG pathways

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

  Methanogenesis = colSums(kegg_selection["K00400",] + kegg_selection["K00401",], na.rm = TRUE)/2
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

  # Thiosulfate_disp1 =
  # trait_table["Thiosulfate_disp1",] = Thiosulfate_disp1
#
  # Thiosulfate_disp2 =
  # trait_table["Thiosulfate_disp2",] = Thiosulfate_disp2

  Polysulfide_reduction = colSums(kegg_selection["K08352",], na.rm = TRUE)
  trait_table["Polysulfide_reduction",] = Polysulfide_reduction

  Oxy_photosynthesis = colSums(kegg_selection["K02703",] + kegg_selection["K02704",] +
                               kegg_selection["K02705",] + kegg_selection["K02706",] +
                               kegg_selection["K02707",] + kegg_selection["K02708",], na.rm = TRUE)/12 +
                       colSums(kegg_selection["K02689",] + kegg_selection["K02690",] +
                               kegg_selection["K02691",] + kegg_selection["K02692",] +
                               kegg_selection["K02693",] + kegg_selection["K02694",], na.rm = TRUE)/12
  trait_table["Oxy_photosynthesis",] = Oxy_photosynthesis

  Anoxy_photosynthesis = # colSums(kegg_selection["K08940",] + kegg_selection["K08941",] +
                         #         kegg_selection["K08942",] + kegg_selection["K08943",], na.rm = TRUE)/6 +
                         colSums(kegg_selection["K08928",] + kegg_selection["K08929",], na.rm = TRUE)/2
  trait_table["Anoxy_photosynthesis",] = Anoxy_photosynthesis

  Bacteriorhodopsin = kegg_selection["K04641",]
  trait_table["Bacteriorhodopsin",] = Bacteriorhodopsin
  return(trait_table)

}

# good bins
colnames(bin_kegg_selection)
bin_kegg_avTPM = bin_kegg_selection[,c("avTPM", "besthit", "Bin.ID")]
test = reshape(bin_kegg_avTPM, idvar = "besthit", timevar = "Bin.ID", direction = "wide")
rownames(test) = test$besthit
colnames(test) = sub("avTPM.", "", colnames(test))
bin_kegg_avTPM = test[,2:ncol(test)]

traits_table = summarize_kegg_traits(bin_kegg_avTPM)
traits_table[is.na(traits_table)] <- 0

traits_table = traits_table[rowSums(traits_table) > 0, ]
traits_table$besthit = rownames(traits_table)

traits_table = melt(traits_table, id.vars = c("besthit"), variable.name = "Bin.ID")
# PLOT KEGG Pathways PER BIN

# cairo_ps(file.path("/cluster/projects/nn9745k/jing/02_results/svalbard", "bin_kegg_pathways.eps"), width=80/25.4, height=80/25.4,
cairo_ps(file.path("/cluster/home/alexaei/figs", "bin_kegg_pathways.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  ggplot(traits_table, aes(besthit, Bin.ID, fill= -value)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
    axis.text = element_text(size=6),
    axis.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=8),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.line.x = element_line(color = "grey")) +
  scale_fill_distiller(palette = "RdPu") +
  geom_tile()


dev.off()

# interesting bins
colnames(intres_bin_kegg_selection)
intres_bin_kegg_avTPM = intres_bin_kegg_selection[,c("avTPM", "besthit", "Bin.ID")]
test = reshape(intres_bin_kegg_avTPM, idvar = "besthit", timevar = "Bin.ID", direction = "wide")
rownames(test) = test$besthit
colnames(test) = sub("avTPM.", "", colnames(test))
intres_bin_kegg_avTPM = test[,2:ncol(test)]

intres_traits_table = summarize_kegg_traits(intres_bin_kegg_avTPM)
intres_traits_table[is.na(intres_traits_table)] <- 0

intres_traits_table = intres_traits_table[rowSums(intres_traits_table) > 0, ]
intres_traits_table = intres_traits_table[,colSums(intres_traits_table) > 0]
intres_traits_table_scaled = scale(intres_traits_table)
intres_traits_table_binary = intres_traits_table
intres_traits_table_binary[intres_traits_table_binary>0] <-1
intres_traits_table_binary[intres_traits_table_binary<0] <-0


intres_traits_table$besthit = rownames(intres_traits_table)

intres_traits_table = melt(intres_traits_table, id.vars = c("besthit"), variable.name = "Bin.ID")
# PLOT KEGG Pathways PER BIN

cairo_ps(file.path("/cluster/home/alexaei/figs", "intres_bin_kegg_pathways.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  ggplot(intres_traits_table, aes(besthit, Bin.ID, fill= -value)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
    axis.text = element_text(size=6),
    axis.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=8),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.line.x = element_line(color = "grey")) +
  scale_fill_distiller(palette = "RdPu") +
  geom_tile()


dev.off()

inters_traits_heatmap_row <- pheatmap(data.matrix(intres_traits_table_binary),
                            #dendrogram = "row",
                            xlab = "", ylab = "",
                            clustering_distance_col = "correlation",
                            clustering_distance_row = "canberra",
                            main = "",
                            scale = "none",
                            cutree_rows = 5
                            )

cairo_ps(file.path("/cluster/home/alexaei/figs", "inters_traits_heatmap_row.eps"))
  inters_traits_heatmap_row

dev.off()

# plot bin dynamics
bintable_tpm = bintable[,seq(13, 69, by = 2)]
rownames(bintable_tpm) = bintable$Bin.ID
interes_bintable_tpm = bintable_tpm[rownames(bintable_tpm) %in% colnames(intres_traits_table_binary),]
interes_bintable_tpm = scale(interes_bintable_tpm)

RdPu.pal <- RColorBrewer::brewer.pal(7,"RdPu")
BuPu.pal <- RColorBrewer::brewer.pal(7,"BuPu")

inters_bins_heatmap_row <- pheatmap(data.matrix(interes_bintable_tpm),
                            #dendrogram = "row",
                            xlab = "", ylab = "",
                            clustering_distance_col = "correlation",
                            clustering_distance_row = "canberra",
                            color = RdPu.pal,
                            main = "",
                            scale = "none",
                            cutree_rows = 5
                            )

cairo_ps(file.path("/cluster/home/alexaei/figs", "inters_bins_heatmap_row.eps"))
  inters_bins_heatmap_row

dev.off()

#taxonomy
bintable[bintable$Bin.ID %in% rownames(data_ht1 )]

 # add categories of lakes (< 350 m, > 350 and < 3350, > 3350 )
 # same order as inters_traits_heatmap_row.eps

 # Heatmap 1

data_ht1 = t(intres_traits_table_binary)
data_ht1 <- data_ht1[ order(row.names(data_ht1)), ]

test = bintable[bintable$Bin.ID %in% rownames(data_ht1), c("Bin.ID", "Tax", "Disparity")]

test$phylum = str_extract(test$Tax, "p_.*")
test$phylum = sapply(strsplit(test$phylum, "[/;]"), "[", 1)

test$genus = str_extract(test$Tax, "g_.*")
test$genus = sapply(strsplit(test$genus, "[/;]"), "[", 1)
rownames(test) = test$Bin.ID
test <- test[ order(row.names(test)), ]

# metadata
colnames(interes_bintable_tpm) = sub("TPM.", "", colnames(interes_bintable_tpm))

metadata4 = metadata[rownames(metadata) %in% colnames(interes_bintable_tpm),]
metadata4 = metadata4[order(rownames(metadata4)),]
interes_bintable_tpm <- interes_bintable_tpm[, order(colnames(interes_bintable_tpm))]
interes_bintable_tpm <- interes_bintable_tpm[, order(metadata4$gl_dist)]
rownames(data_ht1) = sub("metabat2.", "", rownames(data_ht1))

ha = HeatmapAnnotation(
  df = data.frame(gl_dist = metadata4$gl_dist),
   annotation_height = unit(4, "mm")
  )

ht1 = Heatmap(data_ht1, name = "traits", km = 3,
              col = RdPu.pal,
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              column_names_gp = gpar(fontsize = 9),
              row_names_gp = gpar(fontsize = 9)) +
      Heatmap(test$phylum, name = "phylum", width = unit(3, "mm")) +
      Heatmap(test$genus, name = "genus", width = unit(3, "mm"))

cairo_ps(file.path("/cluster/home/alexaei/figs", "inters_traits_bin_heatmap.eps"))
  ht1

dev.off()

# Heatmap 2
data_ht2 <- interes_bintable_tpm[ order(row.names(interes_bintable_tpm)), ]
rownames(data_ht2) = sub("metabat2.", "", rownames(data_ht2))

ht2 = Heatmap(data_ht2, name = "scaled RPM",
        top_annotation = ha,
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        col = BuPu.pal,
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 9))
# Combine the two heatmaps
cairo_ps(file.path("/cluster/home/alexaei/figs", "inters_bins_heatmap_full.eps"))
ht1 + ht2
dev.off()

# boxplots
data_ht2

lakes_550 = data_ht2[,metadata4$gl_dist<551]
lakes_3300 = data_ht2[,metadata4$gl_dist>551 & metadata4$gl_dist<3300]
lakes_10000 = data_ht2[,metadata4$gl_dist>3300]

.boxplot1 = anno_boxplot(lakes_550, which = "row")
.boxplot2 = anno_boxplot(lakes_3300, which = "row")
.boxplot3 = anno_boxplot(lakes_10000, which = "row")

ha_mix_right1 = HeatmapAnnotation(bxplt = .boxplot1,
                              which = "row", width = unit(2, "cm"))
ha_mix_right2 = HeatmapAnnotation(bxplt = .boxplot2,
                              which = "row", width = unit(2, "cm"))
ha_mix_right3 = HeatmapAnnotation(bxplt = .boxplot3,
                              which = "row", width = unit(2, "cm"))

ht_test = Heatmap(data_ht1, name = "traits", km = 3,
              col = RdPu.pal,
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              column_names_gp = gpar(fontsize = 9),
              row_names_gp = gpar(fontsize = 9)) +
      Heatmap(test$phylum, name = "phylum", width = unit(3, "mm"),
        column_names_gp = gpar(fontsize = 9), col = color_palette_euk[1:9]) +
      Heatmap(test$genus, name = "genus", width = unit(3, "mm"),
        column_names_gp = gpar(fontsize = 9)) +
      ha_mix_right1 +
      ha_mix_right2 +
      ha_mix_right3

cairo_ps(file.path("/cluster/home/alexaei/figs", "test_figure4.eps"))
ht_test
dev.off()

# Metacyc

intres_metacyc_bin = metacyc_pathways[metacyc_pathways$"V1" %in% intres_bintable$Bin.ID,]

rownames(intres_metacyc_bin) = intres_metacyc_bin$V1
intres_metacyc_bin = intres_metacyc_bin[,3:ncol(intres_metacyc_bin)]

intres_metacyc_bin[intres_metacyc_bin == "NF"] <- 0

# select interesting pathways

test = intres_metacyc_bin[,match(intres_metacyc_tax$V1, metacyc_tax$V1)]

names_intres_metacyc = c("Calvin-Benson-Bassham cycle",
                         "Carbon fixation in photosynthetic organisms",
                         "Methane metabolism",
                         "Nitrogen metabolism",
                         "Photosynthesis",
                         "Reductive carboxylate cycle",
                         "Sulfur metabolism",
                         "aerobic respiration II",
                         "aerobic respiration III",
                         "ammonia oxidation II (anaerobic)",
                         "catechol degradation I",
                         "catechol degradation II",
                         "catechol degradation to beta-ketoadipate",
                         "catechol degradation to 2-oxopent-4-enoate I",
                         "catechol degradation to 2-oxopent-4-enoate II",
                         "methanogenesis from CO2",
                         "methanogenesis from acetate",
                         "methanol oxidation to formaldehyde I",
                         "nitrate reduction I",
                         "nitrate reduction II",
                         "nitrate reduction IV",
                         "nitrate reduction V",
                         "nitrate reduction VI",
                         "nitrogen fixation",
                         "oxygenic photosynthesis",
                         "phosphate acquisition",
                         "photorespiration",
                         "photosynthesis light reactions",
                         "reductive TCA cycle I",
                         "reductive TCA cycle II",
                         "sulfate reduction I",
                         "sulfate reduction II",
                         "sulfate reduction IV",
                         "sulfate reduction V",
                         "sulfide oxidation III",
                         "sulfite oxidation I",
                         "sulfite oxidation II",
                         "sulfite oxidation III",
                         "sulfite oxidation IV",
                         "sulfur oxidation I",
                         "sulfur oxidation II",
                         "thiosulfate disproportionation III",
                         "thiosulfate oxidation II",
                         "bacteriochlorophyll biosynthesis",
                         "carotenoid biosynthesis",
                         "sulfide oxidation A ferrooxidans",
                         "sulfide oxidation S novella",
                         "sulfur amino acid biosynthesis S cerevisiae",
                         "sulfur metabolism D sulfoexigens",
                         "sulfur oxidation A ambivalens")
colnames(test) = names_intres_metacyc

data_ht3 = apply(test, 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
rownames(data_ht3) = rownames(test)

data_ht3 = data_ht3[,c("Calvin-Benson-Bassham cycle",
                       "Methane metabolism",
                       "Reductive carboxylate cycle",
                        "aerobic respiration II",
                        "aerobic respiration III",
                        "ammonia oxidation II (anaerobic)",
                        "methanogenesis from CO2",
                        "methanogenesis from acetate",
                        "methanol oxidation to formaldehyde I",
                         "nitrate reduction I",
                         "nitrate reduction II",
                         "nitrate reduction IV",
                         "nitrate reduction V",
                         "nitrate reduction VI",
                         "nitrogen fixation",
                         "oxygenic photosynthesis",
                        "reductive TCA cycle I",
                        "reductive TCA cycle II",
                         "sulfate reduction I",
                         "sulfate reduction II",
                         "sulfate reduction IV",
                         "sulfate reduction V",
                         "sulfide oxidation III",
                         "sulfite oxidation I",
                         "sulfite oxidation II",
                         "sulfite oxidation III",
                         "sulfite oxidation IV",
                         "sulfur oxidation I",
                         "sulfur oxidation II",
                         "thiosulfate disproportionation III",
                         "thiosulfate oxidation II",
                         "bacteriochlorophyll biosynthesis",
                         "carotenoid biosynthesis",
                         "sulfide oxidation A ferrooxidans",
                         "sulfide oxidation S novella",
                         "sulfur oxidation A ambivalens"



  )]

data_ht3 = data_ht3[rowSums(data_ht3) > 0, colSums(data_ht3) > 0]

ht3 = Heatmap(data_ht3, name = "ht1", km = 3,
              col = RdPu.pal,
              column_names_gp = gpar(fontsize = 9),
              row_names_gp = gpar(fontsize = 9))

interes_metacycbin_tpm = bintable_tpm[rownames(bintable_tpm) %in% rownames(data_ht3),]
data_ht4 = scale(interes_metacycbin_tpm)

data_ht4 <- data_ht4[ order(row.names(data_ht4)), ]
ht4 = Heatmap(data_ht4, name = "ht2",
        col = BuPu.pal,
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 9))

cairo_ps(file.path("/cluster/projects/nn9745k/jing/02_results/svalbard", "inters_bins_heatmap_metacyc.eps"))
ht3 + ht4
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

# Plot Bacterial phyla

# Plot Archaeal phyla

# Plot Eukaryotic phyla

# Plot based on 16S rRNA genes

# Comparison between Metabarcoding - Metagenomic

### PLS models comparison and VIPs

set.seed(1234)
test.id <- sample(1:nrow(X), size = .2*nrow(X)) ## Randomly choose 20% of rows for test set

test.data <- data[test.id,]  ## Subset to include rows designated to test set
train.data <- data[-test.id,]  ## Exclude rows designated to test set

# ASV bacteria

path_out = "/cluster/projects/nn9745k/jing/02_results/svalbard"
path_bac = file.path(path_out, "bacteria/dada2optimized")
path_scripts = "/cluster/projects/nn9745k/jing/scripts"

ASV_bac = as.matrix(read.table(file.path(path_bac,"ASV_table.tsv"), header=T, check.names=FALSE))
tax_bac = as.matrix(read.table(file.path(path_bac,"Taxonomy_table_silva.tsv"), header=T, check.names=FALSE))

metadata = read.csv(file.path(path_scripts, "metadata", "svfi19.tsv"), header=T, check.names=FALSE, sep = "\t", stringsAsFactors=FALSE,colClasses=cls)

## Parse data
rownames(metadata) = metadata$sample_id
rownames(ASV_bac) = str_split(rownames(ASV_bac), pattern = "_", simplify = TRUE)[, 1]



dim(ASV_bac)
bac_otu <- otu_table(t(ASV_bac), taxa_are_rows = TRUE)
dim(bac_otu)
bac_tax <- tax_table(tax_bac)
bac_physeq <- phyloseq(bac_otu, bac_tax, sample_data(data.frame(metadata)))
bac_physeq

bac_physeq_clean  = subset_taxa(bac_physeq, Kingdom == "Bacteria")
bac_physeq_clean  = subset_taxa(bac_physeq_clean, Order != "Chloroplast")
bac_physeq_clean  = subset_taxa(bac_physeq_clean, Family != "Mitochondria")

bac_physeq_9000 <- prune_samples(sample_sums(bac_physeq_clean)>=9000, bac_physeq_clean) # removes two samples SV011, SV021
rarefied_bac <- rarefy_even_depth(bac_physeq_9000, sample.size = min(sample_sums(bac_physeq_9000)),
          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) #Normalizing species data

metadata = sample_data(rarefied_bac)

asv = t(otu_table(rarefied_bac))

CO2 = as.matrix(metadata$CO2_sat)
rownames(CO2) = rownames(metadata)

asv_co2 = merge(asv, CO2, by = 0)

asv_co2 = na.omit(asv_co2)

asv_co2 = asv_co2[,2:ncol(asv_co2)]

pa_co2 = asv_co2
pa_co2[pa_co2>0] = 1

test = asv_co2[,colSums(pa_co2) > 3]




control <- trainControl("repeatedcv", number = 10, selectionFunction = "oneSE")

prePoc = preProcess(
  test,
  method = "scale",
  thresh = 0.95,
  pcaComp = NULL,
  na.remove = TRUE,
  k = 5,
  knnSummary = mean,
  outcome = NULL,
  fudge = 0.2,
  numUnique = 3,
  verbose = FALSE,
  freqCut = 25/5,
  uniqueCut = 3,
  cutoff = 0.9,
  rangeBounds = c(0, 1)
)

co2_sat_pls_rcv <- train(V1 ~ ., data = test,
 method = "pls",
 tuneLength = 20,
 trControl = control,
 preProc = c("center","scale"),
 ncomp = 10
 )

co2_sat_pls_rcv$bestTune

summary(co2_sat_pls_rcv$finalModel)

# genera



# KEGG table

tpm_kegg_table


# traits table

traits_table

# MAGs


save.image(file.path(figs_dir, "metags_stats.Rdata"))



