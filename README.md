# cis-trans-pipeline

##Scripts to assign genes to regulatory categories (e.g., cis, trans, cis X trans) based on allele-specific expression data.
Pipeline from Kopania et al.: https://www.biorxiv.org/content/10.1101/2021.08.04.455131v2

###Step 1: Prepare modfiles and pseudogenomes

###Step 2: Prepare RNAseq data
  Raw RNAseq reads were trimmed with Trimmomatic and mapped to pseudogenomes using Tophat2
  Note: Mapping with Tophat2 is required for lapels and suspenders pipeline

###Step 3: Lapels and suspenders pipeline
  These scripts are written to run with slurm; edit SBATCH commands to work with different servers or delete SBATCH commands to run without slurm
  Modify script 02_run_lapels_slurm.sh to read in the proper mapped reads file (bam output from Tophat2) and the proper modfile; runs lapels to get a common coordinate system
  Modify script 03_run_suspenders_slurm.sh to read in the proper bam files to merge (output from lapels); runs suspenders to merge mappings to two different pseudogenomes and assign reads to parents
  Modify script 04_split_F1_bam_by_parent_slurm.sh to read in the proper bam file to split (output from suspenders); extracts only reads that are assigned to one parent or the other by suspenders and saves these in separate bam files
  
###Step 4: Subsample parent data
  If all samples are sequenced at similar coverage, there will be more power to detect differential expression between parents than between alleles within an F1. These scripts downsample parent reads to address this issue.
  Modify 05_run_PEtoSE_slurm.sh to input the proper bam files (output from 04_split_F1_bam_by_parent_slurm.sh); some input data were from paired-end sequencing and some were from single-end, so this script runs 06_PEtoSEbam.R1only.py, which deletes the second read from paired-end files and converts all flags to single-end flags
  Modify 07_subsample_parents_only_slurm.sh to read in the proper bam files (output from 05_run_PEtoSE_slurm.sh); calculates the median number of mapped reads that map to one parent for each F1 sample and downsamples parent samples to match this reacount
  
###Step 5: Quantify number of reads mapping to each gene
  Modify 08_run_HTSeqCount_slurm.sh to read in proper bam files (output from 07_subsample_parents_only_slurm.sh); uses HTSeq count to count number of reads mapping to each gene

###Step 6: Run cis-trans pipeline
  Modify 09_cis_trans_wrapper_slurm.sh to input proper file names. This is a wrapper that will run several scripts:
    10_generate_pca.r Generates a PCA of allele-specific read counts
    11_ASE_DE_edgeR.r Runs EdgeR to identify genes with significant differential expression (DE) or allele-specific expression (ASE) based on a negative-binomial test
    12_fisher_and_logFCplot.r Runs a Fisher's exact test for significant differences between parents and alleles in relative expression levels; plots F1 expression log fold-change (logFC) versus parent logFC
    13_divergence_by_reg_category.r Plots expression divergence for genes separated by regulatory category
