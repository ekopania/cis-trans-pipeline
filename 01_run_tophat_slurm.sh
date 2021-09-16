#!/bin/bash
##PURPOSE: Use slurm to run tophat
#
# Job name:
#SBATCH --job-name=tophat
#SBATCH --output=tophat_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=32 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=4000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:

#############F1s###############
#For a single file:
##ls DATA/F1s/121228_QB3_HS4B/WW_LL6_1M_RS_001.trimmed_1P.fq.gz | while read file; do #ROUND SPERMATIDS
##Loop through all samples of a genotype
#ls DATA/F1s/121228_QB3_HS4B/CC_PP21_2M_LZ_001.trimmed_1P.fq.gz | while read file; do
#        name1=$(echo "${file}" | cut -d "/" -f 4 | cut -d "_" -f 1)
#        name2=$(echo "${file}" | cut -d "/" -f 4 | cut -d "_" -f 2)
#        name3=$(echo "${file}" | cut -d "/" -f 4 | cut -d "_" -f 3)
#        name4=$(echo "${file}" | cut -d "/" -f 4 | cut -d "_" -f 4)
#        if [ "${name4}" == "SP" ]; then #To avoid wasting time on cell types we're not looking at
#		continue
#	elif [ "${name4}" == "DIP" ]; then
#		continue
#	fi
#	name=$(echo "${name1}${name2}${name3}${name4}")
#        echo "${name}"
#        #All R1 and SE
#        stuff1a=$(ls -m DATA/F1s/*/${name1}*${name2}*${name3}*${name4}*trimmed_1P.fq.gz | xargs echo | sed 's/[[:space:]]*//g')
#        stuff1b=$(ls -m DATA/F1s/*/${name1}*${name2}*${name3}*${name4}*trimmed.fq.gz | xargs echo | sed 's/[[:space:]]*//g')
##        stuff1=$(echo ${stuff1a},${stuff1b})
#	stuff1=$stuff1a
#        echo "${stuff1}"
#        #All R2
#        stuff2=$(ls -m DATA/F1s/*/${name1}*${name2}*${name3}*${name4}*trimmed_2P.fq.gz| xargs echo | sed 's/[[:space:]]*//g')
#        echo "${stuff2}"
#        tophat2 -p 16 -g 150 -o "${name}_toLEWES_tophat_out_trimmomatic" MODFILES/ref/LEWES_EiJ.pseudogenome "${stuff1}" "${stuff2}"
#        tophat2 -p 16 -g 150 -o "${name}_toPWK_tophat_out_trimmomatic" MODFILES/ref/PWKPhJ_b38_f "${stuff1}" "${stuff2}"
#        tophat2 -p 16 -g 150 -o "${name}_toWSB_tophat_out_trimmomatic" MODFILES/ref/WSBEiJ_b38_f "${stuff1}" "${stuff2}"
#       tophat2 -p 16 -g 150 -o "${name}_toWSB_tophat_out_trimmomatic" MODFILES/ref/WSB.newVariants "${stuff1}" "${stuff2}"
#	tophat2 -p 32 -g 150 -o "${name}_toPWK_tophat_out_trimmomatic_minQ50minDP10" MODFILES/ref/PWK_PhJ.minQ50minDP10 "${stuff1}" "${stuff2}"
#	tophat2 -p 32 -g 150 -o "${name}_toCZII_tophat_out_trimmomatic_CZIIdownsample" MODFILES/ref/CC_CC72_2F_RG_realigned.downsample  "${stuff1}" "${stuff2}"
#	tophat2 -p 32 -g 150 -o "${name}_toPWK_tophat_out_trimmomatic_PWKhapCallerRef" MODFILES/ref/PWK_PhJ.hapCallerRef "${stuff1}" "${stuff2}"
#done

#############F1s - other format###############
#For a single file:
ls DATA/F1s/uoregon_2046_JG6/CCPP273MRS.trimmed.fq.gz | while read file; do
	name=$(echo "${file}" | cut -d "/" -f 4 | cut -d "." -f 1)
        echo "${name}"
        stuff=$(ls -m DATA/F1s/*/${name}*.trimmed.fq.gz | xargs echo | sed 's/[[:space:]]*//g')
        echo "${stuff}"
	tophat2 -p 32 -g 150 -o "${name}_toCZII_tophat_out_trimmomatic_CZIIdownsample" MODFILES/ref/CC_CC72_2F_RG_realigned.downsample  "${stuff}"
	tophat2 -p 32 -g 150 -o "${name}_toPWK_tophat_out_trimmomatic_PWKhapCallerRef" MODFILES/ref/PWK_PhJ.hapCallerRef "${stuff}"
done


#############PARENTS###############
#ls DATA/PARENTS/PPPP98-?M-LZ_3642-JG-000?_L001.trimmed.fq.gz | while read file; do
#ls DATA/PARENTS/LLLL125-?M-??_3642-JG-000?_L001.trimmed.fq.gz | while read file; do
#ls DATA/PARENTS/CCCC153-?M-??_3642-JG-000?_L001.trimmed.fq.gz | while read file; do
#ls DATA/PARENTS/WWWW8?-?M-??_3642-JG-000?_L001.trimmed.fq.gz | while read file; do
#        begin=$(echo "${file}" | cut -d "_" -f 1)
#        echo "${begin}"
#        name=$(echo "${begin}" | cut -d "/" -f 3)
#        echo "${name}"
#        stuff=$(ls -m ${begin}*trimmed.fq.gz | xargs echo | sed 's/[[:space:]]*//g')
#        echo "${stuff}"
#        tophat2 -p 16 -g 150 -o "${name}_toLEWES_tophat_out_trimmomatic" MODFILES/ref/LEWES_EiJ.pseudogenome "${stuff}"
#       tophat2 -p 16 -g 150 -o "${name}_toPWK_tophat_out_trimmomatic" MODFILES/ref/PWKPhJ_b38_f "${stuff}"
#       tophat2 -p 16 -g 160 -o "${name}_toCZII_tophat_out_trimmomatic" MODFILES/ref/CZII "${stuff}"
#        tophat2 -p 16 -g 150 -o "${name}_toWSB_tophat_out_trimmomatic" MODFILES/ref/WSBEiJ_b38_f "${stuff}"
#done

#############MACKETALDATA###############
#ls DATA/Macketal2016/PPLL_SRR2060955.trimmed_1P.fq.gz | while read file; do
#        name=$(echo "${file}" | cut -d "_" -f 1-2)
#        echo ${name}
#        tophat2 -p 32 -g 150 -o "${name}_toLEWES_tophat_out_trimmomatic_2" MODFILES/ref/LEWES_EiJ.pseudogenome "${file}" "${name}_2P.fq.gz"
#        tophat2 -p 32 -g 150 -o "${name}_toPWK_tophat_out_trimmomatic_2" MODFILES/ref/PWKPhJ_b38_f "${file}" "${name}_2P.fq.gz"
#done
