#!/bin/bash
##PURPOSE: Use slurm to run lapels
#
# Job name:
#SBATCH --job-name=lapels
#SBATCH --output=lapels_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=16 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=40G #Not sure if I should mess with these...
##SBATCH --mem=192000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:

#FOR SINGLE FILES (change sample name and to<strain>):
pylapels -p 16 -n -o CCPP273MRS_toPWK.PWKhapCallerRef.lapels_output.bam MODFILES/PWK_PhJ.hapCallerRef.mod CCPP273MRS_toPWK_tophat_out_trimmomatic_PWKhapCallerRef/accepted_hits.bam &> CCPP273MRS_toPWK.PWKhapCallerRef.lapels.log

#FOR MULTIPLE SAMPLES (loop through same genotype and mapped to same strain):
#find . -wholename "./CCPP211*_toCZII_tophat_out_trimmomatic_CZIIdownsample/accepted_hits.bam" | while read file; do
#        name=$(echo "${file}" | cut -d "_" -f 1-2)
#        echo "${name}"
#	if [ -f "${file}.bai" ]; then
#		echo "bam index exists; make sure it is current"
#	fi
#        pylapels -p 16 -n -o "${name}.downsampledCZII.lapels_output.bam" MODFILES/CC_CC72_2F_RG_realigned.downsample_.907.mod "${file}" &> "${name}.downsampledCZII.lapels.log"
#done
echo "Done!"
