#!/bin/bash
##PURPOSE: Wrapper to run 06_PEtoSEbam.R1only.py, converts all reads in bam to SE format, getting rid of redundant mate pairs
#
# Job name:
#SBATCH --job-name=PEtoSE
#SBATCH --output=PEtoSE_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=192000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
ls PPLL_SRR2060955.trimmed.PwkLewes.*allelesONLY.bam | while read file; do
	echo "${file}"
	python 06_PEtoSEbam.R1only.py "${file}"
done

echo "Done!"
