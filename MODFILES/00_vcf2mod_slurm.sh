#!/bin/bash
##PURPOSE: Use modtools to generate a pseudogenome using a vcf and reference fasta
#
# Job name:
#SBATCH --job-name=vcf2mod_CZIIdownample
#SBATCH --output=vcf2mod_output.CZII.downsampled.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
vcf2mod -o CC_CC72_2F_RG_realigned.downsample_.907.mod Mus_musculus.GRCm38.dna.primary_assembly Mus_musculus.ref.meta CC_CC72_2F NEW_CZII_VCF/CC_CC72_2F_RG_realigned.downsample_.907.genotyped.noHets.vcf.gz
