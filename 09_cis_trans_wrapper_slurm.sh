#!/bin/bash
##PURPOSE: Wrapper to run several scripts for cis vs trans (and other reg category) analyses. These include:
	##10_generate_pca.r - generate PCA of downsampled parents and F1 samples split by allele
	##11_ASE_DE_edgeR.r - use EdgeR to determine genes with significant DE (between parents) and ASE (between alleles within F1)
	##12_fisher_and_logFCplot.r - use EdgeR results and Fisher's exact test to assign genes to a regulatory category; plot parent and F1 logFC colored by reg category
	##13_divergence_by_reg_category.r - plot expression divergence as a function of regulatory category; compare median expression divergence across regulatory categories (COMING SOON: also plot logFC, dN/dS by reg category; do each of these for induced genes)
#
# Job name:
#SBATCH --job-name=cis_trans_wrapper
#SBATCH --output=cis_trans_wrapper_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=96000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
#CHANGE THESE PARAMETERS depending on which samples and cell types you want to run this for
mom_geno="LLLL"
dad_geno="PPPP"
F1_geno="LLPP"
comparison="PwkLewes"
allele1="LEWES"
allele2="PWK"
min_rpkm=1
induced="induced"
echo "Here are the parameters:"
echo "Mom: ${mom_geno}, Dad: ${dad_geno}, F1: ${F1_geno}, Comparison: ${comparison}, Allele 1: ${allele1}, Allele 2: ${allele2}, Min RPKM: ${min_rpkm}, induced: ${induced}"

#PCA
echo "Generating PCAs..."
Rscript 10_generate_pca.r ${mom_geno} ${dad_geno} ${F1_geno} ${comparison} ${allele1} ${allele2} LZ
Rscript 10_generate_pca.r ${mom_geno} ${dad_geno} ${F1_geno} ${comparison} ${allele1} ${allele2} RS
#Rscript 10_generate_pca.r ${mom_geno} ${dad_geno} ${F1_geno} ${comparison} ${allele1} ${allele2} WT

#EdgeR for DE and ASE
echo "Running EdgeR..."
Rscript 11_ASE_DE_edgeR.r ${min_rpkm} ${mom_geno} ${dad_geno} ${F1_geno} ${allele1} ${allele2} LZ
Rscript 11_ASE_DE_edgeR.r ${min_rpkm} ${mom_geno} ${dad_geno} ${F1_geno} ${allele1} ${allele2} RS
#Rscript 11_ASE_DE_edgeR.r ${min_rpkm} ${mom_geno} ${dad_geno} ${F1_geno} ${allele1} ${allele2} WT

#Fisher and logFC
echo "Determining regulatory categories..."
Rscript 12_fisher_and_logFCplot.r ${mom_geno} ${dad_geno} ${F1_geno} ${comparison} ${allele1} ${allele2} LZ LZind #${induced}
mv reg_category_output.txt reg_category_output.${mom_geno}mom.${F1_geno}.LZ.txt
Rscript 12_fisher_and_logFCplot.r ${mom_geno} ${dad_geno} ${F1_geno} ${comparison} ${allele1} ${allele2} RS RSind #${induced}
mv reg_category_output.txt reg_category_output.${mom_geno}mom.${F1_geno}.RS.txt
#Rscript 12_fisher_and_logFCplot.r ${mom_geno} ${dad_geno} ${F1_geno} ${comparison} ${allele1} ${allele2} WT ${induced}
#mv reg_category_output.txt reg_category_output.${mom_geno}mom.${F1_geno}.WT.txt

#Divergence by reg cat
echo "Plotting divergence by regulatory category..."
Rscript 13_divergence_by_reg_category.r reg_category_output.${mom_geno}mom.${F1_geno}.LZ.txt reg_category_output.${mom_geno}mom.${F1_geno}.RS.txt DE_genes.${mom_geno}mom.${F1_geno}.LZ.LRT.parentTableAll.txt DE_genes.${mom_geno}mom.${F1_geno}.RS.LRT.parentTableAll.txt
#Rscript 13_divergence_by_reg_category.r reg_category_output.${mom_geno}mom.${F1_geno}.WT.txt reg_category_output.${mom_geno}mom.${F1_geno}.WT.txt DE_genes.${mom_geno}mom.${F1_geno}.WT.LRT.parentTableAll.txt DE_genes.${mom_geno}mom.${F1_geno}.WT.LRT.parentTableAll.txt

echo "Done!"
