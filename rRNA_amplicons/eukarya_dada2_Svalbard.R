#eukarya

#In R
.libPaths(c("/cluster/projects/nn9745k/Rpackages", .libPaths()))
#install.packages("/cluster/projects/nn9745k/Rpackages/BiocManager_1.30.10.tar.gz", repos = NULL, type="source")
#BiocManager::install("dada2", version = "3.10")

# load libraries
if (!require("dada2")) {
   install.packages("dada2", dependencies = TRUE)
   library(dada2)
   }
if (!require("kableExtra")) {
   install.packages("kableExtra", dependencies = TRUE)
   library(kableExtra)
   }
if (!require("stringr")) {
   install.packages("stringr", dependencies = TRUE)
   library(stringr)
   }
if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2)
   }
if (!require("dplyr")) {
   install.packages("dplyr", dependencies = TRUE)
   library(dplyr)
   }
if (!require("tibble")) {
   install.packages("tibble", dependencies = TRUE)
   library(tibble)
   }
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

if (!require("dada2", quietly = TRUE)){
	BiocManager::install("dada2")
	library(dada2)
	}
if (!require("phyloseq", quietly = TRUE)){
	BiocManager::install("phyloseq")
	library(phyloseq)
	}
if (!require("Biostrings", quietly = TRUE)){
	BiocManager::install("Biostrings")
	library(Biostrings)
	}


#### Define the name of directories to use
fastq_dir <- "/cluster/projects/nn9745k/jing/02_results/svalbard/eukarya/AdaptersRemoved"  # fastq directory with the samples we are using
path_out = "/cluster/projects/nn9745k/jing/02_results/svalbard/eukarya/dada2"
database_dir <- "/cluster/projects/nn9745k/03_databases/silva138"  # folder with databases

filtered_dir <- file.path(path_out, "fastq_filtered")  # fastq filtered
qual_dir <- file.path(path_out, "qual_pdf")  # qual pdf
dada2_dir <- file.path(path_out, "dada2")  # dada2 results
blast_dir <- file.path(path_out, "blast")  # blast2 results
figs_dir <- file.path(path_out, "figs")

dir.create(filtered_dir)
dir.create(qual_dir)
dir.create(dada2_dir)
dir.create(blast_dir)
dir.create(figs_dir)


#### Examine fastq files
# get a list of all fastq files in the fastq" directory and separate R1 and R2
fns <- sort(list.files(fastq_dir, full.names = TRUE))
fns <- fns[str_detect(basename(fns), ".fastq")]
fns_R1 <- fns[str_detect(basename(fns), "R1")]
fns_R2 <- fns[str_detect(basename(fns), "R2")]

df <- data.frame()

for (i in 1:length(fns_R1)) {

  # use the dada2 function fastq.geometry
  geom <- fastq.geometry(fns_R1[i])

  # extract the information on number of sequences and file name
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(fns_R1[i]))

  # add one line to data frame
  df <- bind_rows(df, df_one_row)
}

# Make some output tables and figures
write.table(df, file = file.path(figs_dir, "n_seq.txt"), sep='\t', row.names = FALSE, na='',
            quote=FALSE)

pdf(file.path(figs_dir, "quantity_profiles.pdf"))
	ggplot(df, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity",
                                            binwidth = 1000)
	hist(df$n_seq)
dev.off()

pdf(file.path(figs_dir, "quality_profiles.pdf"))
plotQualityProfile(fns_R1[5:14])
plotQualityProfile(fns_R2[5:14])
dev.off()

filt_R1 <- file.path(filtered_dir, basename(fns_R1))
filt_R2 <- file.path(filtered_dir, basename(fns_R2))

#### Filter and Trim
out = filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2,
	  truncLen=c(270,190), maxN=0, maxEE=c(2,2),
	  truncQ=20, rm.phix=TRUE, compress=TRUE,
	  multithread=FALSE, matchIDs = TRUE)	#Use multithread=TRUE when running on Saga

filt <- sort(list.files(filtered_dir, full.names = TRUE))
filt <- filt[str_detect(basename(filt), ".fastq")]
filt_R1 <- filt[str_detect(basename(filt), "R1")]
filt_R2 <- filt[str_detect(basename(filt), "R2")]

#### DADA2 error models
errF = learnErrors(filt_R1, multithread=TRUE, nbases = 1e7)
errR = learnErrors(filt_R2, multithread=TRUE, nbases = 1e7)

pdf(file.path(figs_dir, "error_plots.pdf"))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

derepFs = derepFastq(filt_R1, verbose=TRUE)
derepRs = derepFastq(filt_R2, verbose=TRUE)
names(derepFs) = basename(filt_R1)
names(derepRs) = basename(filt_R2)

dadaFs = dada(derepFs, err=errF, multithread=TRUE)
dadaRs = dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab = makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## get length information on sequences
length_table = table(nchar(getSequences(seqtab)))
write.table(length_table, file.path(figs_dir, "length_table.tsv"))

pdf(file.path(figs_dir, "length_profile.pdf"))
	plot(table(nchar(getSequences(seqtab)))) #simple plot of length distribution
dev.off()

## summary information on sequences
getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- basename(fns_R1)
track

write.table(track, file.path(figs_dir, "track_table.tsv"))

#### Transforming and saving the ASVs sequences
seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% rownames_to_column(var = "sequence") %>% rowid_to_column(var = "OTUNumber") %>% mutate(OTUNumber = sprintf("otu%04d",OTUNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))

seq_out <- Biostrings::DNAStringSet(seqtab.nochim_trans$sequence)
names(seq_out) <- seqtab.nochim_trans$OTUNumber
seq_out

Biostrings::writeXStringSet(seq_out, file.path(dada2_dir, "ASV_no_taxo.fasta"),
                            compress = TRUE, width = 20000)

## make training data from mothur output
# dada2:::makeTaxonomyFasta_Silva(file.path("/cluster/projects/nn9745k/03_databases/silva138/silva.nr_v138.align"), file.path("/cluster/projects/nn9745k/03_databases/silva138/silva.nr_v138.tax"), "/cluster/projects/nn9745k/03_databases/silva138/silva_nr_v138_train_set.fa.gz")

## assign taxonomy
#taxa_silva = assignTaxonomy(seqtab.nochim, file.path(database_dir, "silva_nr_v138_train_set.fa.gz"), multithread=TRUE, tryRC = TRUE)
#taxa_silva.print = taxa_silva
#rownames(taxa_silva.print) = NULL
#write.table(taxa_silva.print, file.path(path_out,"Taxonomy_table_silva.tsv"))
#head(taxa_silva.print)
#taxa_silva[which.max(colSums(seqtab.nochim)),]

#pr2_file <- "/cluster/projects/nn9745k/03_databases/pr2_version_4.12.0_18S_dada2.fasta.gz"
#PR2_tax_levels <- c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family",
#                    "Genus", "Species")

#taxa_pr2 = assignTaxonomy(seqtab.nochim, refFasta = pr2_file, taxLevels = PR2_tax_levels, multithread=TRUE, tryRC = TRUE)
#taxa_pr2.print = taxa_pr2
#rownames(taxa_pr2.print) = NULL
#write.table(taxa_pr2.print, file.path(path_out,"Taxonomy_table_pr2.tsv"))
#head(taxa_pr2.print)
#taxa_pr2[which.max(colSums(seqtab.nochim)),]


ASV.print = seqtab.nochim
colnames(ASV.print) = seq(1, ncol(ASV.print), 1)
write.table(ASV.print,  file.path(path_out,"ASV_table.tsv"))
ASV.table = as.matrix(read.table(file.path(path_out,"ASV_table.tsv"), header=T, check.names=FALSE))

save.image(file.path(dada2_dir, "dada2.Rdata"))


