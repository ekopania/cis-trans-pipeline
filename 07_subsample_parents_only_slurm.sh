#!/bin/bash
##PURPOSE: Subsample parent bams to get a similar number of reads as map to a single parent in F1s. This is to avoid the power issue of having ~2X the number of reads mapping to a single parent than to a parent allele if you have similar sequencing coverage between parents and F1s.
#
# Job name:
#SBATCH --job-name=subsample_parents
#SBATCH --output=subsample_parents_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=96000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
#Overall method:
#Determine mean number reads from experimental F1s that map preferentially to one parent (either parent, so for suspenders this means either po:i:1 or po:i:2)
#Determine proportion of reads mapping to correct parent from parental samples that will approximately equal the mean # F1 single-parent reads (above)
#Subsample this proportion of reads from parental samples (ONLY from reads mapping to correct parent!)

#Calculating mean number of single-parent F1 reads
total_reads=0
n=0
while read file; do
        file_reads=$( samtools view -c -F 260 "${file}" )
        total_reads=$(( $total_reads + $file_reads ))
        n=$(( $n + 1 ))
        echo "${total_reads}"
        echo "${n}"
done < <(ls PPPP98-*LZ.PwkCzii.PWKallelesONLY.bam)
#Piping creates a "subshell" and gets rid of all the variable w/in it when it ends; do this instead
#Other options, depending on genotypes being analyzed:
#(ls WWWW*RS.WsbLewes*allelesONLY.bam) #(ls WWWW*LZ.WsbLewes*allelesONLY.bam) #(ls WWLL*LZ.WsbLewes*allelesONLY.SE_R1orMateUnmappedOnly.bam) #(ls WWLL*RS.WsbLewes*allelesONLY.SE_R1orMateUnmappedOnly.bam) #(ls CCPP*RS.PwkCzii.*allelesONLY.SE_R1orMateUnmappedOnly.bam) #(ls PPLL*LZ.PwkLewes.*allelesONLY.SE_R1orMateUnmappedOnly.bam) #(ls LLPP*RS.PwkLewes.*allelesONLY.SE_R1orMateUnmappedOnly.bam) #(ls PPLL*RS.PwkLewes.*allelesONLY.SE_R1orMateUnmappedOnly.bam) #(ls LLPP*LZ.PwkLewes.*allelesONLY.SE_R1orMateUnmappedOnly.bam) #(ls PPLL_SRR*ONLY.bam)  #(ls LLPP_SRR*ONLY.bam) #(ls PPLL_SRR*PwkLewes*allelesONLY.SE_R1orMateUnmappedOnly.bam)
echo "${total_reads}"
echo "${n}"
m=$(( $total_reads / $n))
echo "${m}"

#Calculate proportion and subsample from parents (reads mapping to correct parent ONLY)
#ls LLWW_SRR*DOMallelesONLY.bam.rmdup.bam | while read file; do #Katya's data - DOM parents
#ls LLWW_SRR*trimmomatic*DOMallelesONLY.SE.bam | while read file; do #Katya's data - DOM parents
#ls LLWW_SRR*trimmomatic*LEWESallelesONLY.bam | while read file; do #Katya's data - DOM parents
#ls LLWW_SRR*trimmomatic*SE_R1orMateUnmappedOnly.bam | while read file; do #Katya's data - DOM parents
#ls LLWW_SRR*trimmomatic*DOMallelesONLY.rmIdenticalReads.bam | while read file; do #Katya's data - DOM parents
#ls LLLL125-*LZ.PwkLewes.LEWESallelesONLY.bam | while read file; do #Our data - PURE LEWES parent, LZ
#ls LLLL125-*RS.PwkLewes.LEWESallelesONLY.bam | while read file; do #Our data - PURE LEWES parent, RS
#ls LLLL125-*LZ.WsbLewes.LEWESallelesONLY.bam | while read file; do #Our data - intraF1, LZ
#ls LLLL125-*RS.WsbLewes.LEWESallelesONLY.bam | while read file; do #Our data - intraF1, RS
ls CCCC153-*LZ.PwkCzii.CZIIallelesONLY.bam | while read file; do #Our data - PURE CZII parent, LZ
#ls CCCC153-*RS.PwkCzii.CZIIallelesONLY.bam | while read file; do #Our data - PURE CZII parent, RS
#ls WWLL*MLZ.PwkLewes.LEWESallelesONLY.SE_R1orMateUnmappedOnly.bam | while read file; do #Our data - intrasubspec dom parent, LZ
#ls WWLL*MRS.PwkLewes.LEWESallelesONLY.SE_R1orMateUnmappedOnly.bam | while read file; do #Our data - intrasubspec dom parent, RS
#       name=$(echo "${file}" | cut -d "." -f 1)
        name=$(echo "${file}" | cut -d "." -f 1-3)
        echo "${name}"
        parent_readcount=$(samtools view -c -F 260 "${file}")
        if [ ${parent_readcount} -gt ${m} ]; then
                prop=$(bc <<< "scale=3; $m / $parent_readcount")
                echo "${parent_readcount}"
                echo "${prop}"
#		samtools view -b -h -s $prop $file > "${name}.downsample_${prop}.SE_R1only.WWLL.NEWdownsample.bam" #WWLL
               samtools view -b -h -s $prop $file > "${name}.downsample_${prop}.SE_R1only.CCPP.NEWdownsample.bam" #CCPP
#               samtools view -b -h -s $prop $file > "${name}.downsample_${prop}.SE_R1only.LLPP.NEWdownsample2.bam" #LLPP
#                samtools view -b -h -s $prop $file > "${name}.downsample_${prop}.SE_R1only.PPLL.NEWdownsample.bam" #PPLL
        else
                echo "${name} has lower readcount than F1 allele-specific mean"
        fi
done

#ls PPCC_SRR*MUSallelesONLY.bam.rmdup.bam | while read file; do #Katya's data - MUS parents
#ls PPCC_SRR*trimmomatic*MUSallelesONLY.SE.bam | while read file; do #Katya's data - MUS parents
#ls PPCC_SRR*trimmomatic*PWKallelesONLY.bam | while read file; do #Katya's data - MUS parents
#ls PPCC_SRR*trimmomatic*SE_R1orMateUnmappedOnly.bam | while read file; do #Katya's data - MUS parents
#ls PPCC_SRR*trimmomatic*MUSallelesONLY.rmIdenticalReads.bam | while read file; do #Katya's data - MUS parents
#ls WWWW8*LZ.WsbLewes.WSBallelesONLY.bam | while read file; do #Our data - intraF1, LZ
#ls WWWW8*RS.WsbLewes.WSBallelesONLY.bam | while read file; do #Our data - intraF1, RS
#ls PPPP98-*LZ.PwkLewes.PWKallelesONLY.bam | while read file; do #Our data - PURE PWK parent, LZ
#ls PPPP98-*RS.PwkLewes.PWKallelesONLY.bam | while read file; do #Our data - PURE PWK parent, RS
#ls PPPP98-*LZ.PwkCzii.PWKallelesONLY.bam | while read file; do #Our data - PURE PWK parent, LZ
#ls PPPP98-*RS.PwkCzii.PWKallelesONLY.bam | while read file; do #Our data - PURE PWK parent, RS
#ls CCPP*MLZ.PWKallelesONLY.bam.SE_R1orMateUnmappedOnly.bam | while read file; do #Our data - intrasubspec mus parent, LZ
#ls CCPP*MRS.PwkLewes.PWKallelesONLY.SE_R1orMateUnmappedOnly.bam | while read file; do #Our data - intrasubspec mus parent, RS
#        name=$(echo "${file}" | cut -d "." -f 1)
#        name=$(echo "${file}" | cut -d "." -f 1-3)
#        echo "${name}"
#        parent_readcount=$(samtools view -c -F 260 "${file}")
#        if [ ${parent_readcount} -gt ${m} ]; then
#                prop=$(bc <<< "scale=3; $m / $parent_readcount")
#                echo "${parent_readcount}"
#                echo "${prop}"
##		samtools view -b -h -s $prop $file > "${name}.downsample_${prop}.SE_R1only.WWLL.NEWdownsample.bam" #WWLL
#                samtools view -b -h -s $prop $file > "${name}.downsample_${prop}.SE_R1only.CCPP.NEWdownsample.bam" #CCPP
###               samtools view -b -h -s $prop $file > "${name}.downsample_${prop}.SE_R1only.LLPP.NEWdownsample2.bam" #LLPP
##                samtools view -b -h -s $prop $file > "${name}.downsample_${prop}.SE_R1only.PPLL.NEWdownsample.bam" #PPLL
#        else
#                echo "${name} has lower readcount than F1 allele-specific mean"
#        fi
#done

echo "Done!"
