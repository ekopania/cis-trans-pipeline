#PURPOSE: Use edgeR to determine genes DE between parent genotypes and between alleles within a hybrid F1 for ASE analyses

args<-commandArgs(TRUE)
if(length(args) != 7){
        stop("Missing command line arguments.\nArgument 1: minimum RPKM (ex: 1)\nArgument 2: mom genotype (ex: WWLL)\nArgument 3: dad genotype (ex: CCPP)\nArgument 4: F1 genotype (ex: LLPP)\nArgument 5: allele 1 (ex: LEWES)\nArgument 6: allele 2 (ex: PWK)\nArgument 7: cell type (ex: LZ)")
}

min_rpkm<-args[1]
out_pref<-paste0("DE_genes.",args[2],"mom.",args[4],".",args[7],".LRT")
files_mom<-list.files(pattern=paste0(args[2],".*",args[7],".*SE_R1only\\.",args[4],"\\.NEWdownsample\\.htseq_count\\.out$"))
files_dad<-list.files(pattern=paste0(args[3],".*",args[7],".*SE_R1only\\.",args[4],"\\.NEWdownsample\\.htseq_count\\.out$"))
files_f1<-list.files(pattern=paste0(args[4],".*",args[7],".*SE_R1orMateUnmappedOnly.*",args[4],"\\.NEWdownsample\\.htseq_count\\.out$"))
files<-c(files_mom,files_dad,files_f1)
print(min_rpkm)
print(out_pref)
print(files)

print("Loading EdgeR and setting up DGE object...")
library(edgeR)

min_samples<-length(files)
#Assign group based on genotype
mygroups<-c()
for(f in files){
#	geno<-gsub("_SRR.*$","",f) #WT
	geno<-gsub("\\d*M.*$","",f) #LZ and RS
#	if(geno=="WWLL" || geno=="LLWW" || geno=="LLLL"){
#	if(geno=="CCCC"){
	if(geno==args[2]){
		mygroups<-c(mygroups,1)
#	} else if(geno=="CCPP" || geno=="PPCC" || geno=="PPPP"){
#	} else if(geno=="PPPP"){
	} else if(geno==args[3]){
		mygroups<-c(mygroups,2)
#	} else if(geno=="LLPP" || geno=="PPLL"){
#	} else if(geno=="CCPP"){
	} else if(geno==args[4]){
		#Figure out alleles
#		allele<-gsub("allelesONLY.*$","",gsub("LLPP_SRR\\d{7}\\.trimmomatic\\.PwkLewes\\.","",f)) #WT; LLPP
#		allele<-gsub("allelesONLY.*$","",gsub("PPLL_SRR\\d{7}\\.trimmed\\.PwkLewes\\.","",f)) #WT; PPLL
#		allele<-gsub("allelesONLY\\.SE_R1.*$","",gsub("CCPP\\d{3}M.*\\.PwkCzii\\.","",f)) #CCPP
		allele<-gsub("allelesONLY\\.SE_R1.*$","",gsub("[[:alpha:]]{4}\\d*M[[:punct:]]?[[:alpha:]]{2}\\.[[:alpha:]]*\\.","",f)) #hopefully this is generic, definitely works for cell sorted LLPP and PPLL
#		if(allele=="LEWES"){
#		if(allele=="CZII"){
		if(allele==args[5]){
			print(allele)
			mygroups<-c(mygroups,3)
#		} else if(allele=="PWK"){
		} else if(allele==args[6]){
			mygroups<-c(mygroups,4)
		} else{
			stop(paste("Invalid allele:",allele))
		}
	} else{
		stop(paste("Invalid genotype:",geno))
	}
}
print(mygroups)
myDGE<-readDGE(files,columns=c(1,2),group=mygroups,comment.char="#",header=FALSE)

print(dim(myDGE))
#print(head(myDGE$counts))
#print(head(myDGE$samples))

#Filter out things below certain threshold level
#Need to get gene lengths to figure out fpkm to figure out what to filter out
#Using feature counts output from a single sample to get gene lengths
featureCountTable<-read.table("WWLL31MLZ.PwkLewes.LEWESallelesONLY.downsample_.950.SE_R1only.LLPP.NEWdownsample.featurcounts.out",header=TRUE)
featureCountTable_sorted<-featureCountTable[order(featureCountTable$Geneid),]
#print(head(featureCountTable_sorted))
lens<-featureCountTable_sorted$Length
print(head(lens))
print(length(lens))
keep<-rowSums(rpkm(myDGE, gene.length=as.numeric(lens))>min_rpkm) >= min_samples
myDGE<-myDGE[keep, , keep.lib.sizes=FALSE]

#Normalize by library size
myDGE<-calcNormFactors(myDGE)
print(dim(myDGE$counts))
print(head(myDGE$samples))

#Generate design table
spec<-c()
for(i in myDGE$samples$group){
        if(i==1){
                spec<-c(spec,"dom")
        }
        else if(i==2){
                spec<-c(spec,"mus")
        }
        else if(i==3){
                spec<-c(spec,"dom_allele")
        }
        else{
                spec<-c(spec,"mus_allele")
        }
}

design_table<-as.data.frame(cbind(rownames(myDGE$samples), spec))
colnames(design_table)<-c("file","grp")
design_matrix<-model.matrix(~0+design_table$grp, data=myDGE$samples)
colnames(design_matrix)<-levels(design_table$grp)
myDGE<-estimateDisp(myDGE,design_matrix)
print(design_matrix)

#Set up contrasts
print("Setting up contrasts...")
myContrasts<-makeContrasts(domVmus="dom-mus", domAlleleVmusAllele="dom_allele-mus_allele", levels=design_table$grp)
print(myContrasts)

#In my original DE EdgeR analysis (phylogenetic exp evo) I had a whole thing to determine if a gene is exp in a particular species/cell type but I don't think that's necessary here b/c only looking at one cell type at a time and also using a more conservative estimate of above a certain fpkm in all samples

#Only look at autos
#source("http://bioconductor.org/biocLite.R")
#biocLite("ensembldb")
#biocLite("EnsDb.Mmusculus.v75")
print("Getting autosomal genes...")
library(EnsDb.Mmusculus.v75)
edb<-EnsDb.Mmusculus.v75
edb_auto<-addFilter(edb, SeqNameFilter(c(1:19)))
auto_genes<-genes(edb_auto)
write(intersect(rownames(myDGE$counts),auto_genes$gene_id), paste(out_pref,"expGenes.txt",sep="."))

#Test for DE genes using the QLFtest
print("Starting QLFtest...")
#qlfit<-glmQLFit(myDGE,design_matrix)
fit<-glmFit(myDGE,design_matrix) #Calling this qlfit so I don't have to change variables but it's really a regular glmFit instead of glmQLFit
qlf_numDE_table_auto<-c()
#lfc_cutoff_qlf_numDE_table<-c()
qlf_tables<-NULL
for(i in 1:ncol(myContrasts)){
#        result<-glmQLFTest(qlfit,contrast=myContrasts[,i])
	result<-glmLRT(fit,contrast=myContrasts[,i])
        assign(paste("qlf",colnames(myContrasts)[i], sep="."),result)

        compare<-sub("vs", "", colnames(myContrasts)[i])
        #print(compare)

        tt<-topTags(result,n=nrow(result),adjust.method="fdr",p.value=0.05)
        assign(paste("qlf.tt",colnames(myContrasts)[i], sep="."), tt)
#        lfc_cutoff_qlf_numDE_table<-rbind(lfc_cutoff_qlf_numDE_table, cbind(paste("lfc.qlf.numDE",colnames(myContrasts)[i], sep="."), length(intersect(intersect(which(tt$table$FDR<0.05),which(rownames(tt$table) %in% rownames(myDGE$counts))), which(abs(tt$table$logFC)>min_logFC)))))
        qlf_tables[[colnames(myContrasts)[i]]]<-tt$table

        #auto only
        assign(paste("qlf.DE",colnames(myContrasts)[i],"auto",sep="."), tt$table[intersect(which(tt$table$FDR<0.05), which(rownames(tt$table) %in% auto_genes$gene_id)),])
        assign(paste("qlf.numDE",colnames(myContrasts)[i],"auto",sep="."), length(intersect(which(tt$table$FDR<0.05), which(rownames(tt$table) %in% auto_genes$gene_id))))
        qlf_numDE_table_auto<-rbind(qlf_numDE_table_auto,cbind(paste("qlf.numDE",colnames(myContrasts)[i],"auto",sep="."),length(intersect(intersect(which(tt$table$FDR<0.05), which(rownames(tt$table) %in% auto_genes$gene_id)),which(rownames(tt$table) %in% rownames(myDGE$counts))))))
}

#print(lfc_cutoff_qlf_numDE_table)

print(qlf_numDE_table_auto)

print(head(qlf.DE.domVmus.auto))
print(head(qlf.DE.domAlleleVmusAllele.auto))

#Write qlf.DE tables to output files
print("Writing outputs...")
write.table(qlf.domVmus$table,paste(out_pref,"parentTableAll.txt",sep="."),quote=FALSE,sep="\t")
write.table(qlf.domAlleleVmusAllele$table,paste(out_pref,"hybridTableAll.txt",sep="."),quote=FALSE,sep="\t")
write.table(qlf.DE.domVmus.auto, paste(out_pref,"parents.txt",sep="."), quote=FALSE, sep="\t")
write.table(qlf.DE.domAlleleVmusAllele.auto, paste(out_pref,"hybrids.txt",sep="."), quote=FALSE, sep="\t")

print("Done!")
