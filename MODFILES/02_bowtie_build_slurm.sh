#!/bin/bash
##PURPOSE: Use bowtie2 to generate a bowtie build based on a modtools pseudogenome fasta
#
# Job name:
#SBATCH --job-name=bowtie2-build_CZIIdownsample
#SBATCH --output=bowtie2-build_output.CZII.downsample.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=4000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
bowtie2-build ref/CC_CC72_2F_RG_realigned.downsample_.907.fa ref/CC_CC72_2F_RG_realigned.downsample
