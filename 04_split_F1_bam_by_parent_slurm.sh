#!/bin/bash
##PURPOSE: Use slurm to split suspenders output bam into two bams: reads that confidently map best to parent 1 and reads that confidently map best to parent 2
#
# Job name:
#SBATCH --job-name=split_by_parent
#SBATCH --output=split_by_parent_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
ls CCPP211M??.PwkCzii.june2020.suspenders_output.bam | while read file; do
#        name=$(echo "${file}" | cut -d "." -f 1-2)
	name=$(echo "${file}" | cut -d "." -f 1-3)
        echo "${name}"
        samtools view -h "${file}" | grep -v "po:i:1\|po:i:3" | samtools view -bS - > "${name}.CZIIallelesONLY.bam"
        samtools view -h "${file}" | grep -v "po:i:2\|po:i:3" | samtools view -bS - > "${name}.PWKallelesONLY.bam"
done
