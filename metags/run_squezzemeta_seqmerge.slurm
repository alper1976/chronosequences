#!/bin/bash

# Job name: squezzemeta_svalbard
#SBATCH --job-name=sqmeta
#SBATCH --account=nn9745k
#SBATCH --time=7-00:00:00
#SBATCH --mem=47616MB
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1

##Setup job environment
module purge   # clear any inherited modules
set -o errexit # exit on errors

source /cluster/software/Anaconda3/2019.03/etc/profile.d/conda.sh

conda activate /cluster/projects/nn9745k/scripts/conda_envs/squeezemeta

SAMPLE_LIST="/cluster/work/users/alexaei/13_svalbard_metaGs/list_files.txt"
SAMPLE_FOLDER="/cluster/work/users/alexaei/13_svalbard_metaGs"
OUTPUT_FOLDER="/cluster/work/users/alexaei/02_results/13_svalbard_metaGs/seqmerge"


SqueezeMeta.pl -m seqmerge -s $SAMPLE_LIST -f $SAMPLE_FOLDER -p $OUTPUT_FOLDER --cleaning -cleaning_options "ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" -miniden 60 -a megahit
# cd $OUTPUT_FOLDER
# cd ..
# restart.pl seqmerge