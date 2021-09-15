#!/bin/bash
##PURPOSE: Use slurm to run lapels
#
# Job name:
#SBATCH --job-name=lapels
#SBATCH --output=lapels_output-%j.txt
##SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user= #Where to send mail
#SBATCH --cpus-per-task=16 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=40G #Not sure if I should mess with these...
##SBATCH --mem=192000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=<partition_name>
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:

#FOR SINGLE FILES (change input filenames):
outbam=""
modfile=""
inbam=""
logfile=""
pylapels -p 16 -n -o ${outbam} ${modfile} ${inbam} &> ${logfile}

#FOR MULTIPLE SAMPLES (loop through same genotype and mapped to same strain):
file_list=$(find . -wholename "/accepted_hits.bam")
modfile=""
for file in "${file_list[@]}"; do
	name=$(echo "${file}" | cut -d "_" -f 1-2)
	echo "${name}"
	if [ -f "${file}.bai" ]; then
		echo "bam index exists for input file; make sure it is current"
	fi
	pylapels -p 16 -n -o "${name}.lapels_output.bam" ${modfile} "${file}" &> "${name}.lapels.log"
done

echo "Done!"
