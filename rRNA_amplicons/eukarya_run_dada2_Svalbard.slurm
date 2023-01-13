#!/bin/bash
# Job name: dada2 Eukarya svalbard
#SBATCH --job-name=dada2_e
#
# Project (change to your :
#SBATCH --account=nn9745k
#
# Wall clock limit:
#SBATCH --time=10:0:0
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=30720MB

#
## Number of tasks (nodes):
##SBATCH --ntasks-per-node=12
##

##Setup job environment
#source /cluster/bin/jobsetup
# module purge   # clear any inherited modules
set -o errexit # exit on errors

#Load R module

module load R/3.6.0-foss-2019a

#Set variable with number of processors and load mpi

mpirun -np 1 R CMD BATCH --no-save --no-restore /cluster/projects/nn9745k/jing/scripts/eukarya/eukarya_run_dada2_Svalbard.R

#End of script
