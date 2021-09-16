#!/bin/bash
#PURPOSE: Downsample CZII bam to similar number of mapped reads as PWK bam to generate similar quality vcf
#
# Job name:
#SBATCH --job-name=CZII_downsample
#SBATCH --output=CZII_downsample.log
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
pwk_reads=$(samtools view -c -F 260 PWKPhJ.bam)
czii_reads=$(samtools view -c -F 260 CC_CC72_2F_RG_realigned.bam)
prop=$(bc <<< "scale=3; $pwk_reads / $czii_reads")
echo "${pwk_reads} ${czii_reads} ${prop}"
samtools view -b -h -s $prop CC_CC72_2F_RG_realigned.bam > CC_CC72_2F_RG_realigned.downsample_${prop}.bam
