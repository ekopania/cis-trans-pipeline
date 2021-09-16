#!/bin/bash
##PURPOSE: Use slurm to run suspenders for a single file
#
# Job name:
#SBATCH --job-name=suspenders
#SBATCH --output=suspenders_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
#SBATCH --partition=good_lab_cpu
#
## Command(s) to run:
#find . -name "LLLL125*_toWSB.lapels_output.bam" | while read file; do
#find . -name "CCPP232MLZ_toPWK.hapCallerRef.lapels_output.bam" | while read file; do
find . -name "CCPP273MRS_toPWK.PWKhapCallerRef.lapels_output.bam" | while read file; do
#	sample=$(echo "${fe}" | cut -d "_" -f 1-2)
        sample=$(echo "${file}" | cut -d "_" -f 1)
        echo "${sample}"
        pysuspenders "${sample}.PwkCzii.june2020.suspenders_output.bam" "${file}" "${sample}_toCZII.downsampledCZII.lapels_output.bam" &> "${sample}.PwkCzii.june2020.suspenders.log"
done

echo "Done!"
