### This is an R script to analyse to infer greenhouse gas emissions and metagenomic
### data from Svalbard


## saga
# module purge   # clear any inherited modules
# source /cluster/software/Anaconda3/2019.03/etc/profile.d/conda.sh
# conda activate /cluster/projects/nn9745k/scripts/conda_envs/squeezemeta

# load(file.path(figs_dir, "metags_stats.Rdata"))


#In R
## LOAD PACKAGES
library(SQMtools)


## LOAD VARIABLES
# path_arch
path_out = "/cluster/projects/nn9745k/jing/02_results/svalbard"
path_metags = file.path("/cluster/work/users/alexaei/02_results/13_svalbard_metaGs")

figs_dir <- file.path(path_out, "figs")

# load metagenome tables (KEGG, pfams, COGs)

dir_list = list.dirs(path = path_metags, full.names = FALSE, recursive = FALSE)[-1]
# dir_list
dir = "FI001"
data_sqm = loadSQM(file.path(path_metags, dir))

tpm_kegg_table  = data.frame(data_sqm$functions$KEGG$tpm)
tpm_pfam_table  = data.frame(data_sqm$functions$PFAM$tpm)
cno_kegg_table  = data.frame(data_sqm$functions$KEGG$copy_number)
cno_pfam_table  = data.frame(data_sqm$functions$PFAM$copy_number)


for (dir in dir_list) {
  data_sqm = loadSQM(file.path(path_metags, dir))
  tpm_kegg_data = data.frame(data_sqm$functions$KEGG$tpm)
  tpm_pfam_data = data.frame(data_sqm$functions$PFAM$tpm)
  cno_kegg_data = data.frame(data_sqm$functions$KEGG$copy_number)
  cno_pfam_data = data.frame(data_sqm$functions$PFAM$copy_number)

  tpm_kegg_table = merge(tpm_kegg_table, tpm_kegg_data, by="row.names", all = TRUE)
  rownames(tpm_kegg_table) = tpm_kegg_table$Row.names
  tpm_kegg_table = subset(tpm_kegg_table, select=-c(Row.names))

  tpm_pfam_table = merge(tpm_pfam_table, tpm_pfam_data, by="row.names", all = TRUE)
  rownames(tpm_pfam_table) = tpm_pfam_table$Row.names
  tpm_pfam_table = subset(tpm_pfam_table, select=-c(Row.names))

  cno_kegg_table = merge(cno_kegg_table, cno_kegg_data, by="row.names", all = TRUE)
  rownames(cno_kegg_table) = cno_kegg_table$Row.names
  cno_kegg_table = subset(cno_kegg_table, select=-c(Row.names))

  cno_pfam_table = merge(cno_pfam_table, cno_pfam_data, by="row.names", all = TRUE)
  rownames(cno_pfam_table) = cno_pfam_table$Row.names
  cno_pfam_table = subset(cno_pfam_table, select=-c(Row.names))

  # extract sequencing statistics

}

write.table(tpm_kegg_table, file = file.path(figs_dir, "tpm_kegg_table.csv"),
  append = FALSE, quote = FALSE, sep = "\t")
write.table(tpm_pfam_table, file = file.path(figs_dir, "tpm_pfam_table.csv"),
 append = FALSE, quote = FALSE, sep = "\t")
write.table(cno_kegg_table, file = file.path(figs_dir,"cno_kegg_table.csv"),
 append = FALSE, quote = FALSE, sep = "\t")
write.table(cno_pfam_table, file = file.path(figs_dir, "cno_pfam_table.csv"),
 append = FALSE, quote = FALSE, sep = "\t")


dir = "FI001"

genus_abund_table  = read.table(file.path(path_metags, dir, "results/tables/FI001.genus.nofilter.abund.tsv"), sep = "\t", header = TRUE)
rownames(genus_abund_table) = genus_abund_table$X
genus_abund_table = subset(genus_abund_table, select=-c(X))

files_list = list.files(path = path_metags, pattern='*\\.genus.nofilter.abund.tsv', recursive=TRUE)[-1]

file = files_list[5]

for (file in files_list) {
  genus_abund_data = read.table(file.path(path_metags, file), sep = "\t", header = TRUE)
  rownames(genus_abund_data) = genus_abund_data$X
  genus_abund_data = subset(genus_abund_data, select=-c(X))

  genus_abund_table = merge(genus_abund_table, genus_abund_data, by="row.names", all = TRUE)
  rownames(genus_abund_table) = genus_abund_table$Row.names
  genus_abund_table = subset(genus_abund_table, select=-c(Row.names))


    # extract sequencing statistics

}

write.table(genus_abund_table, file = file.path(figs_dir, "genus_abund_table.tsv"),
  append = FALSE, quote = FALSE, sep = "\t")

genus_abund_table  = read.table(file.path(figs_dir, "genus_abund_table.csv"), sep = "\t", header = TRUE)



write.table(pathway_abund_table, file = file.path(figs_dir, "pathway_abund_table.tsv"),
  append = FALSE, quote = FALSE, sep = "\t")


save.image(file.path(figs_dir, "metags_stats.Rdata"))

