#!/bin/bash
##PURPOSE: run HTSeqCount to count reads mapping to features (genes) on final bams (lapels suspenders pipeline, split by parent, converted to SE, and parents downsampled to match allele-specifc F1 readcounts)
#
# Job name:
#SBATCH --job-name=HTSeqCount
#SBATCH --output=HTSeqCount_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=192000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:

#PARENTS (downsampled)
#ls *katyaData.PPLL.mapGT10.bam | while read file; do
#ls *SRR*ONLY*downsample*.PPLL.NEWdownsample.bam | while read file; do
#ls *SRR*SE_R1orMateUnmappedOnly*downsample*.PPLL.NEWdownsample2.bam | while read file; do
#ls *_??.PwkLewes*LLPP.NEWdownsample.bam | while read file; do
#ls *M_RS.PwkLewes*SE_R1only.LLPP.NEWdownsample.bam | while read file; do
#ls *M_RS.*.downsample_.*SE_R1only.PPLL*bam | while read file; do
#ls CCCC*RS*allelesONLY.downsample*SE_R1only.CCPP.NEWdownsample.bam | while read file; do
#ls CCCC*LZ.PwkCzii*downsample*bam | while read file; do
#ls WWWW*LZ.WsbLewes*WWLL.NEWdownsample.bam | while read file; do
#ls LLLL*.WsbLewes*WWLL.NEWdownsample.bam | while read file; do
#        name=$(echo "${file}" | cut -d "." -f 1-8)
#        echo "${name}"
#        htseq-count --format=bam --order=name --stranded=no --type=exon --mode=union "${file}" /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.78.gtf &> "${name}.htseq_count.out"
#done

#F1s or parents w/ low read counts (NOT downsampled)
#ls PPLL_SRR*ONLY.bam | while read file; do
#ls LLPP_SRR*trimmomatic*ONLY.SE.bam | while read file; do
ls PPLL_SRR2060955*allelesONLY.SE_R1orMateUnmappedOnly.bam | while read file; do
#ls PPLL???M??.PwkLewes.*allelesONLY.SE_R1orMateUnmappedOnly.bam | while read file; do
#ls CCPP21*RS.PwkCzii.june2020.*allelesONLY.SE_R1orMateUnmappedOnly.noNeg.bam | while read file; do
#ls PPPP*LZ.PwkCzii.*allelesONLY.bam | while read file; do
#ls WWLL*.WsbLewes.*allelesONLY.SE_R1orMateUnmappedOnly.bam | while read file; do
#ls WWWW*.WsbLewes.*allelesONLY.bam | while read file; do
#        name=$(echo "${file}" | cut -d "." -f 1-3)
        name=$(echo "${file}" | cut -d "." -f 1-4)
        echo "${name}"
#       htseq-count --format=bam --order=pos --stranded=no --type=exon --mode=union "${file}" /data/DBs/genomes/mouse_reference/Ensembl/GRCm38_release78/Mus_musculus.GRCm38.78.gtf &> "${name}.katyaData.LLPP.NEWdownsample.htseq_count.out"
#       htseq-count --format=bam --order=name --stranded=no --type=exon --mode=union "${file}" /data/DBs/genomes/mouse_reference/Ensembl/GRCm38_release78/Mus_musculus.GRCm38.78.gtf &> "${name}.SE_R1only.LLPP.NEWdownsample2.htseq_count.out"
#         htseq-count --format=bam --order=name --stranded=no --type=exon --mode=union "${file}" /data/DBs/genomes/mouse_reference/Ensembl/GRCm38_release78/Mus_musculus.GRCm38.78.gtf &> "${name}.SE_R1only.PPLL.NEWdownsample.htseq_count.out" #musculus version
	htseq-count --format=bam --order=name --stranded=no --type=exon --mode=union "${file}" /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.78.gtf &> "${name}.SE_R1only.PPLL.NEWdownsample.htseq_count.out" #Griz version
#       htseq-count --format=bam --order=name --stranded=no --type=exon --mode=union "${file}" /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.78.gtf &> "${name}.SE_R1only.CCPP.NEWdownsample.htseq_count.out"
#	htseq-count --format=bam --order=name --stranded=no --type=exon --mode=union "${file}" /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.78.gtf &> "${name}.SE_R1only.WWLL.NEWdownsample.htseq_count.out"
done

####FULL COUNT####
#NOT split or downsampled
#ls CCPP273MRS.PwkCzii.june2020.suspenders_output.bam | while read file; do
#ls WW*RS*.suspenders_output.bam | while read file; do
#	name=$(echo "${file}" | cut -d "." -f 1-3)
#	echo "${name}"
#	htseq-count --format=bam --order=name --stranded=no --type=exon --mode=union "${file}" /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.78.gtf &> "${name}.FULLsample.htseq_count.out"
#done

#Insert comment character (#) into htseq output files so that R knows to only read in the gene counts in downstream scripts
#ls WWWW*.WWLL.NEWdownsample.htseq_count.out | while read file; do
#ls *katyaData.PPLL.NEWdownsample.htseq_count.out | while read file; do
ls PPLL_SRR2060955.trimmed.PwkLewes.*PPLL.NEWdownsample.htseq_count.out | while read file; do
#ls PPPP*LZ*SE_R1only.CCPP.NEWdownsample.htseq_count.out | while read file; do
#ls *RS*.FULLsample.htseq_count.out | while read file; do 
       sed -i '/^E/! s/^/#/' "${file}"
done
