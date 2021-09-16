#PURPOSE: Determine regulatory categories for each gene and plot logFC

args<-commandArgs(TRUE)
if(length(args) != 8){
        stop("Missing command line arguments.\nArgument 1: mom genotype (ex: WWLL)\nArgument 2: dad genotype (ex: CCPP)\nArgument 3: F1 genotype (ex: LLPP)\nArgument 4: comparison (ex: PwkLewes)\nArgument 5: allele 1 (ex: LEWES)\nArgument 6: allele 2 (ex: PWK)\nArgument 7: cell type (ex: LZ)\nArgument 8: look at all genes or just a subset? (options: \'LZind\',\'RS_ind\',\'none\')")
}

print("Here are the arguments - check to make sure parents and offspring correspond!")
print(args)

#Read in R packages
library(scales)
library(ggplot2)
library(topGO)
library(praise)

mom_files<-list.files(pattern=paste0(args[1],".*",args[7],".*SE_R1only.",args[3],"\\.NEWdownsample\\.htseq_count\\.out$"))
dad_files<-list.files(pattern=paste0(args[2],".*",args[7],".*SE_R1only.",args[3],"\\.NEWdownsample\\.htseq_count\\.out$"))
momAllele_files<-list.files(pattern=paste0(args[3],".*",args[7],"\\.",args[4],"\\.",args[5],"allelesONLY\\.SE_R1orMateUnmappedOnly\\.",args[3],"\\.NEWdownsample\\.htseq_count\\.out$"))
dadAllele_files<-list.files(pattern=paste0(args[3],".*",args[7],"\\.",args[4],"\\.",args[6],"allelesONLY\\.SE_R1orMateUnmappedOnly\\.",args[3],"\\.NEWdownsample\\.htseq_count\\.out$"))

#Prefix for EdgeR output file to read in
out_pref<-paste0("DE_genes.",args[1],"mom.",args[3],".",args[7],".LRT")

print(mom_files)
print(dad_files)
print(momAllele_files)
print(dadAllele_files)
print(out_pref)

#Read in DE genes from script 10 output
parent_DE<-read.table(paste(out_pref,"parents.txt",sep="."),header=TRUE)
hybrid_DE<-read.table(paste(out_pref,"hybrids.txt",sep="."),header=TRUE)
parent_sig_genes<-rownames(parent_DE)
f1_sig_genes<-rownames(hybrid_DE)
#print(dim(parent_DE))
#print(dim(hybrid_DE))
#print(head(parent_DE))
#print(head(hybrid_DE))

#Get genes expressed in all samples
#Maybe outputing this from edgeR after fpkm filtering is best way to do this?
overlap_genes<-scan(paste(out_pref,"expGenes.txt",sep="."),what=character())
##OPTIONAL: only look at induced genes for each cell type
opt<-args[8]
if(opt=="LZind"){
        print("Looking at LZ induced genes only")
        LZind<-scan("gene_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.txt",what=character())
        gene_set<-intersect(overlap_genes,LZind)
} else if(opt=="RSind"){
        print("Looking at RS induced genes only")
        RSind<-scan("gene_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.txt",what=character())
        gene_set<-intersect(overlap_genes,RSind)
} else{
        print("Looking at all expressed genes")
        gene_set<-overlap_genes
}
print(head(gene_set))

#TRYING SEVERAL THINGS HERE:
#As of Apr 2020 I'm doing option 2
#Option 2: Randomly choose 3 samples from the genotype(s) w/ >3 samples; pool read counts from all samples of same geno
#Option 3: Take mean of readcount across all samples
#Option 4: Take median of reacount across all samples

#Implementing option 2: pooling from 3 samples
#MOM
if(length(mom_files)>3){
        mom_subset<-sample(mom_files,2)
	print(mom_subset)
} else{
        mom_subset<-mom_files
}
mom1<-read.table(mom_subset[1],header=FALSE)
mom2<-read.table(mom_subset[2],header=FALSE)
mom3<-read.table(mom_subset[3],header=FALSE)
mom_all<-cbind(as.character(mom1[,1]),mom1[,2],mom2[,2],mom3[,2])
#mom_all<-cbind(as.character(mom1[,1]),mom1[,2],mom2[,2])
keep<-which(mom_all[,1] %in% gene_set)
mom_all_filtered<-mom_all[keep,]
mom_reads<-as.data.frame(cbind(genes=as.character(mom_all_filtered[,1]),readcounts=rowSums(apply(mom_all_filtered[,2:4],c(1,2),as.numeric))))
#mom_reads<-as.data.frame(cbind(genes=as.character(mom_all_filtered[,1]),readcounts=rowSums(apply(mom_all_filtered[,2:3],c(1,2),as.numeric))))
print(dim(mom_reads))
print(head(mom_reads))
#DAD
if(length(dad_files)>3){
	dad_subset<-sample(dad_files,3)
#if(length(dad_files)>2){
#	dad_subset<-sample(dad_files,2)
#	print(dad_subset)
} else{
        dad_subset<-dad_files
}
dad1<-read.table(dad_subset[1],header=FALSE)
dad2<-read.table(dad_subset[2],header=FALSE)
dad3<-read.table(dad_subset[3],header=FALSE)
dad_all<-cbind(as.character(dad1[,1]),dad1[,2],dad2[,2],dad3[,2])
#dad_all<-cbind(as.character(dad1[,1]),dad1[,2],dad2[,2])
keep<-which(dad_all[,1] %in% gene_set)
dad_all_filtered<-dad_all[keep,]
dad_reads<-as.data.frame(cbind(genes=as.character(dad_all_filtered[,1]),readcounts=rowSums(apply(dad_all_filtered[,2:4],c(1,2),as.numeric))))
#dad_reads<-as.data.frame(cbind(genes=as.character(dad_all_filtered[,1]),readcounts=rowSums(apply(dad_all_filtered[,2:3],c(1,2),as.numeric))))
print(dim(dad_reads))
print(head(dad_reads))
#MOM ALLELE
if(length(momAllele_files)>3){
	momAllele_subset<-sample(momAllele_files,3)
} else{
	momAllele_subset<-momAllele_files
}
momAllele1<-read.table(momAllele_subset[1],header=FALSE)
momAllele2<-read.table(momAllele_subset[2],header=FALSE)
momAllele3<-read.table(momAllele_subset[3],header=FALSE)
#momAllele4<-read.table(momAllele_subset[4],header=FALSE)
#momAllele5<-read.table(momAllele_subset[5],header=FALSE)
momAllele_all<-cbind(as.character(momAllele1[,1]),momAllele1[,2],momAllele2[,2],momAllele3[,2])
#momAllele_all<-cbind(as.character(momAllele1[,1]),momAllele1[,2],momAllele2[,2])
#momAllele_all<-cbind(as.character(momAllele1[,1]),momAllele1[,2],momAllele2[,2],momAllele3[,2],momAllele4[,2],momAllele5[,2])
keep<-which(momAllele_all[,1] %in% gene_set)
momAllele_all_filtered<-momAllele_all[keep,]
momAllele_reads<-as.data.frame(cbind(genes=as.character(momAllele_all_filtered[,1]),readcounts=rowSums(apply(momAllele_all_filtered[,2:4],c(1,2),as.numeric))))
#momAllele_reads<-as.data.frame(cbind(genes=as.character(momAllele_all_filtered[,1]),readcounts=rowSums(apply(momAllele_all_filtered[,2:3],c(1,2),as.numeric))))
#momAllele_reads<-as.data.frame(cbind(genes=as.character(momAllele_all_filtered[,1]),readcounts=rowSums(apply(momAllele_all_filtered[,2:6],c(1,2),as.numeric))))
print(dim(momAllele_reads))
print(head(momAllele_reads))
#DAD ALLELE
if(length(dadAllele_files)>3){
        dadAllele_subset<-sample(dadAllele_files,3)
} else{
        dadAllele_subset<-dadAllele_files
}
dadAllele1<-read.table(dadAllele_subset[1],header=FALSE)
dadAllele2<-read.table(dadAllele_subset[2],header=FALSE)
dadAllele3<-read.table(dadAllele_subset[3],header=FALSE)
#dadAllele4<-read.table(dadAllele_subset[4],header=FALSE)
#dadAllele5<-read.table(dadAllele_subset[5],header=FALSE)
dadAllele_all<-cbind(as.character(dadAllele1[,1]),dadAllele1[,2],dadAllele2[,2],dadAllele3[,2])
#dadAllele_all<-cbind(as.character(dadAllele1[,1]),dadAllele1[,2],dadAllele2[,2])
#dadAllele_all<-cbind(as.character(dadAllele1[,1]),dadAllele1[,2],dadAllele2[,2],dadAllele3[,2],dadAllele4[,2],dadAllele5[,2])
keep<-which(dadAllele_all[,1] %in% gene_set)
dadAllele_all_filtered<-dadAllele_all[keep,]
dadAllele_reads<-as.data.frame(cbind(genes=as.character(dadAllele_all_filtered[,1]),readcounts=rowSums(apply(dadAllele_all_filtered[,2:4],c(1,2),as.numeric))))
#dadAllele_reads<-as.data.frame(cbind(genes=as.character(dadAllele_all_filtered[,1]),readcounts=rowSums(apply(dadAllele_all_filtered[,2:3],c(1,2),as.numeric))))
#dadAllele_reads<-as.data.frame(cbind(genes=as.character(dadAllele_all_filtered[,1]),readcounts=rowSums(apply(dadAllele_all_filtered[,2:6],c(1,2),as.numeric))))
print(dim(dadAllele_reads))
print(head(dadAllele_reads))

print("Running Fisher's exact test...")
fisher_genes<-c()
fisherResults<-c()
for(i in 1:nrow(mom_reads)){
        m<-as.numeric(as.character(mom_reads$readcounts[i]))
        d<-as.numeric(as.character(dad_reads$readcounts[i]))
        ma<-as.numeric(as.character(momAllele_reads$readcounts[i]))
        da<-as.numeric(as.character(dadAllele_reads$readcounts[i]))
        if(!(all(sapply(list(as.character(dad_reads$genes[i]),as.character(momAllele_reads$genes[i]),as.character(dadAllele_reads$genes[i])), FUN=identical, as.character(mom_reads$genes[i])))))
        {
                stop(paste("Gene names different across files! Gene:", mom_reads$genes[i]))
        }
        gene<-as.character(mom_reads$genes[i])
        fisher_genes<-c(fisher_genes,gene)
        tab<-cbind(rbind(m,d),rbind(ma,da))
        result<-fisher.test(tab)
        fisherResults<-c(fisherResults,result$p.value)
}
print(length(fisherResults))
print(head(fisherResults))
fisherResults_fdr<-p.adjust(fisherResults,method="fdr")
head(fisher_genes)

#assign reg categories based on DE (output from script 10) and Fisher's Exact Test above
total_expressed<-length(gene_set)
cis<-c()
trans<-c()
compensatory<-c()
cisPlusTrans<-c()
cisPlusTransSame<-c()
cisPlusTransOpp<-c()
cisXtrans<-c()
notDEorComp<-c()
other<-c()
for(i in 1:length(fisher_genes)){
        gene<-as.character(fisher_genes[i])
        fisher_p<-fisherResults_fdr[i]
        m<-as.numeric(as.character(mom_reads$readcounts[i]))
        d<-as.numeric(as.character(dad_reads$readcounts[i]))
        ma<-as.numeric(as.character(momAllele_reads$readcounts[i]))
        da<-as.numeric(as.character(dadAllele_reads$readcounts[i]))
        if(!(all(sapply(list(as.character(mom_reads$genes[i]),as.character(dad_reads$genes[i]),as.character(momAllele_reads$genes[i]),as.character(dadAllele_reads$genes[i])), FUN=identical, gene))))
        {
                stop(paste("Gene names different across files! Gene:", mom_reads$genes[i]))
        }
        if((gene %in% parent_sig_genes) && (gene %in% f1_sig_genes) && (fisher_p>0.05)){
                cis<-c(cis,gene)
        } else if((gene %in% parent_sig_genes) && !(gene %in% f1_sig_genes) && (fisher_p<=0.05)){
                trans<-c(trans,gene)
        } else if(!(gene %in% parent_sig_genes) && (gene %in% f1_sig_genes) && (fisher_p<=0.05)){
                compensatory<-c(compensatory,gene)
        } else if((gene %in% parent_sig_genes) && (gene %in% f1_sig_genes) && (fisher_p<=0.05)){
                if(((m>d) && (ma>da)) || ((m<d) && (ma<da))){
                        cisPlusTrans<-c(cisPlusTrans,gene)
                        if(abs(log2(d/m)) > abs(log2(da/ma))){
                              cisPlusTransSame<-c(cisPlusTransSame,gene)
                        } else if (abs(log(d/m)) < abs(log(da/ma))){
                              cisPlusTransOpp<-c(cisPlusTransOpp,gene)
                        }
                } else{
                        cisXtrans<-c(cisXtrans,gene)
                }
        } else if(!(gene %in% parent_sig_genes)){
		notDEorComp<-c(notDEorComp,gene)
	} else{
                other<-c(other,gene)
        }
}

num_DE_or_ASE<-total_expressed-length(notDEorComp)
num_DE_or_ASE2<-length(cis)+length(trans)+length(cisXtrans)+length(compensatory)+length(cisPlusTrans)+length(other)
if(num_DE_or_ASE!=num_DE_or_ASE2){
	stop("ERROR: disagreement about # of genes that are no conserved (DE or ASE)")
}
print(paste("# Genes expressed in at least one genotype:",total_expressed))
print(paste("# cis:",length(cis)))
print(paste("# trans:",length(trans)))
print(paste("# cisXtrans:",length(cisXtrans)))
print(paste("# compensatory:",length(compensatory)))
print(paste("# cis+trans:",length(cisPlusTrans)))
print(paste("# cis+trans - opposite directions:",length(cisPlusTransOpp)))
print(paste("# cis+trans - same direction:",length(cisPlusTransSame)))
print(paste("# not DE in parents or compensatory:", length(notDEorComp)))
print(paste("# other:",length(other)))
print("Proportion of total # genes expressed:")
print(paste("proportion cis:",length(cis)/total_expressed))
print(paste("proportion trans:",length(trans)/total_expressed))
print(paste("proportion cisXtrans:",length(cisXtrans)/total_expressed))
print(paste("proportion compensatory:",length(compensatory)/total_expressed))
print(paste("proportion cis+trans:",length(cisPlusTrans)/total_expressed))
print(paste("proportion cis+trans - opposite directions:",length(cisPlusTransOpp)/total_expressed))
print(paste("proportion cis+trans - same direction:",length(cisPlusTransSame)/total_expressed))
print(paste("proportion not DE in parents or compensatory:", length(notDEorComp)/total_expressed))
print(paste("proportion other:",length(other)/total_expressed))
print("Proportion of # genes w/ DE or ASE")
print(paste("proportion cis:",length(cis)/num_DE_or_ASE))
print(paste("proportion trans:",length(trans)/num_DE_or_ASE))
print(paste("proportion cisXtrans:",length(cisXtrans)/num_DE_or_ASE))
print(paste("proportion compensatory:",length(compensatory)/num_DE_or_ASE))
print(paste("proportion cis+trans:",length(cisPlusTrans)/num_DE_or_ASE))
print(paste("proportion cis+trans - opposite directions:",length(cisPlusTransOpp)/num_DE_or_ASE))
print(paste("proportion cis+trans - same direction:",length(cisPlusTransSame)/num_DE_or_ASE))
print(paste("proportion other:",length(other)/num_DE_or_ASE))

#Plot logFC, color by gene category
parent_logFC_all<-read.table(paste(out_pref,"parentTableAll.txt",sep="."),header=TRUE)
hybrid_logFC_all<-read.table(paste(out_pref,"hybridTableAll.txt",sep="."),header=TRUE)
parent_logFC_table<-parent_logFC_all[which(rownames(parent_logFC_all) %in% gene_set),]
hybrid_logFC_table<-hybrid_logFC_all[which(rownames(hybrid_logFC_all) %in% gene_set),]
genes_p<-rownames(parent_logFC_table)
genes_f<-rownames(hybrid_logFC_table)
print(head(genes_p))
print(head(genes_f))
gene_categories_p<-c()
for(g in genes_p){
        if(g %in% cis){
                gene_categories_p<-c(gene_categories_p,"cis")
        } else if(g %in% trans){
                gene_categories_p<-c(gene_categories_p,"trans")
        } else if(g %in% compensatory){
                gene_categories_p<-c(gene_categories_p,"comp")
        } else if(g %in% cisPlusTransSame){
                gene_categories_p<-c(gene_categories_p,"cPlusTSame")
        } else if(g %in% cisPlusTransOpp){
                gene_categories_p<-c(gene_categories_p,"cPlusTOpp")
        } else if(g %in% cisXtrans){
                gene_categories_p<-c(gene_categories_p,"cXt")
	} else if(g %in% notDEorComp){
		gene_categories_p<-c(gene_categories_p,"conserved")
        } else if(g %in% other){
                gene_categories_p<-c(gene_categories_p,"other")
        } else{
                stop(paste("ERROR:",g,"not in any regulatory category!"))
        }
}
gene_categories_f<-c()
for(g in genes_f){
        if(g %in% cis){
                gene_categories_f<-c(gene_categories_f,"cis")
        } else if(g %in% trans){
                gene_categories_f<-c(gene_categories_f,"trans")
        } else if(g %in% compensatory){
                gene_categories_f<-c(gene_categories_f,"comp")
        } else if(g %in% cisPlusTransSame){
                gene_categories_f<-c(gene_categories_f,"cPlusTSame")
        } else if(g %in% cisPlusTransOpp){
                gene_categories_f<-c(gene_categories_f,"cPlusTOpp")
        } else if(g %in% cisXtrans){
                gene_categories_f<-c(gene_categories_f,"cXt")
        } else if(g %in% notDEorComp){
                gene_categories_f<-c(gene_categories_f,"conserved") 
	} else if(g %in% other){
                gene_categories_f<-c(gene_categories_f,"other")
        } else{
                stop(paste("ERROR:",g,"not in any regulatory category!"))
        }
}
logfc_p_table<-as.data.frame(cbind(genes_p,logfc_p=parent_logFC_table$logFC,gene_categories_p))
logfc_f_table<-as.data.frame(cbind(genes_f,logfc_f=hybrid_logFC_table$logFC,gene_categories_f))
head(logfc_p_table)
head(logfc_f_table)

print("Plotting fold changes...")
logfc_p_filtered<-logfc_p_table[which(logfc_p_table$genes_p %in% logfc_f_table$genes_f),]
logfc_p_final<-logfc_p_filtered[order(logfc_p_filtered$genes_p),]
logfc_f_filtered<-logfc_f_table[which(logfc_f_table$genes_f %in% logfc_p_table$genes_p),]
logfc_f_final<-logfc_f_filtered[order(logfc_f_filtered$genes_f),]
print(dim(logfc_p_final))
print(dim(logfc_f_final))
print(head(logfc_p_final))
print(head(logfc_f_final))
if(!(all(as.character(logfc_p_final$genes_p)==as.character(logfc_f_final$genes_f)))){
        print("ERROR: gene names between parent and f1 not the same!")
}
if(!(all(as.character(logfc_p_final$gene_categories_p)==as.character(logfc_f_final$gene_categories_f)))){
        print("ERROR: gene categories between parent and f1 not the same!")
}
combo_data<-as.data.frame(cbind(as.character(logfc_p_final$genes_p),as.numeric(as.character(logfc_p_final$logfc_p)),as.numeric(as.character(logfc_f_final$logfc_f)),as.character(logfc_p_final$gene_categories_p)))
colnames(combo_data)<-cbind("gene","logfc_p","logfc_f","gene_category")
write.table(combo_data,file="reg_category_output.txt",sep="\t",quote=FALSE)
#print(head(as.character(combo_data[,4])))
final_data<-combo_data[which(as.character(combo_data[,4])!="conserved"),]
colnames(final_data)<-cbind("gene","logfc_p","logfc_f","gene_category")
#print(dim(combo_data))
#print(dim(final_data))
final_data$gene_category<-factor(final_data$gene_category,levels=c("cis","trans","cXt","comp","cPlusTOpp","cPlusTSame","other"))
print(head(final_data))
p<-ggplot(final_data, aes(x=as.numeric(as.character(logfc_p)), y=as.numeric(as.character(logfc_f)), color=gene_category)) + geom_point()
p<- p + scale_color_manual(values=c("red","brown","blue","orange","green","purple","yellow"))
p<- p + xlab("parent logFC")
p<- p + ylab("F1 logFC")
p<- p + ggtitle("Parent vs F1 expression logFC")
p<- p + theme_minimal() + geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") + theme(text=element_text(size=16))
p<- p + xlim(-4,4) + ylim(-5,5)

pdf(paste("logFC_plot_",args[3],"_",args[7],".pdf",sep=""),width=11,height=8.5)
p
dev.off()

#GO enrichment analysis
#groups<-c("cis","trans","compensatory","cisPlusTrans","cisPlusTransSame","cisPlusTransOpp","cisXtrans","notDEorComp","other")
#for(g in groups){
#        print(g)
#        all<-get(g)
#        inAll<-c()
#        for(gene in gene_set){
#                if(gene %in% all){
#                        inAll<-c(inAll,1)
#                }
#                else{
#                        inAll<-c(inAll,0)
#                }
#        }
#        if(length(all) != length(which(inAll==1))){
#                stop("ERROR: length of 'all' file and 'inAll' vector don't match up!")
#        }
#        names(inAll)<-gene_set
#        geneSelFunc<-function(iO){
#                return(iO==1)
#        }
#        grpVtotal_exp<-new("topGOdata",description=paste0(g,"Vtotal_exp"),ontology="BP",allGenes=inAll,geneSel=geneSelFunc,annot=annFUN.org, mapping = "org.Mm.eg.db", ID="ensembl")
#        result_grpVtotal_exp <- runTest(grpVtotal_exp, algorithm = "classic", statistic = "fisher")
#        resultData_grpVtotal_exp<- GenTable(grpVtotal_exp, classFisher=result_grpVtotal_exp,orderBy="class$Fisher",ranksOf="classicFisher",topNodes=10)
#        print(resultData_grpVtotal_exp)
#        write.table(resultData_grpVtotal_exp,paste0("GO_",g,"Vtotal_exp.",args[5],".txt"),quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")
#}

#Calculate compensatory change using Hunter Fraser's method (Fraser 2019 https://www.sciencedirect.com/science/article/pii/S0168952518301550)
#Step 1: parental difference = logFC(P1_exp/P2_exp)
#Step 2: cis = logFC(A1_exp/A2_exp)
#Step 3: trans = parental - cis <-- calculate cis with DIFFERENT HYBRID INDIVIDUAL than in step 2 to prevent neg corr between cis and trans inherent in using same individual (novel improvement from Fraser 2019)
#Step 4: Compare cis to trans (so this isn't actually in Fraser 2019 but I think intuitively cis=parental means all cis for this gene, trans=parental means all trans, cis>parental and trans<0 means compensatory, cis<parental and trans>0 means reinforcing; I think if the two hybrids give contrasting results such as cis>parental but trans>0 it's ambiguous and we shouldn't assign that gene to a category; also I think we can only do cis vs trans vs compensator vs reinforcing here and can't break down cis and trans together into more nuanced categories like cXt)

#parental difference already in logfc_p_final as long as ok to use combined parental samples - I think so?
#print("Applying Fraser 2019 method...")
#DE_genes<-union(parent_sig_genes,f1_sig_genes)
#DE_genes_table<-logfc_p_final[which(logfc_p_final[,1] %in% DE_genes),]
#print(length(DE_genes))
#print(dim(DE_genes_table))
#fraser2019_propCis<-c()
#fraser2019_propTrans<-c()
#fraser2019_propComp<-c()
#fraser2019_propRein<-c()
#fraser2019_propAmbig<-c()
#for(j in 1:100){
#	reg_cat_fraser2019method<-c()
#	#for(i in 1:nrow(logfc_p_final)){
#	for(i in 1:nrow(DE_genes_table)){
#		#Using logFC between parents based on edgeR - Not sure if this will work b/c it doesn't give the actual values for each parent, which we may need to do a binomial test below
#		#parent_diff<-as.numeric(as.character(gene$logfc_p))
#		#Alternative: randomly sample one mom and one dad; use logFC between the two of them
#		mom_sample<-sample(2:ncol(mom_all_filtered),1)
#		dad_sample<-sample(2:ncol(dad_all_filtered),1)
#	#	mom<-as.numeric(as.character(mom_all_filtered[which(mom_all_filtered[,1]==logfc_p_final[i,1]),mom_sample]))
#		mom<-as.numeric(as.character(mom_all_filtered[which(mom_all_filtered[,1]==DE_genes_table[i,1]),mom_sample]))
#	#	dad<-as.numeric(as.character(dad_all_filtered[which(dad_all_filtered[,1]==logfc_p_final[i,1]),dad_sample]))
#		dad<-as.numeric(as.character(dad_all_filtered[which(dad_all_filtered[,1]==DE_genes_table[i,1]),dad_sample]))
#		parent_diff<-log(mom/dad)
#		#COME UP WITH A WAY TO INCORPORATE VARIANCE/BIO REPS!!!
#		#Fraser 2019 doesn't really address this and instead emphasizes the fact that you don't really need bioreps (just two hybrid F1s)
#		#Maybe iterate through this loop several times, randomly sampling different indivs each time, and see if there is a consensus? ie only report the reg category for a gene if it shows up in the same reg category in over n% of iterations; Could probably figure out the correct n% based on null expectation that a given gene will fall into any of the four categories 25% of the time, so if a gene is associated w/ a specific category significantly more than 25% of the time we can confidently say it belongs in that category
#	
#		#cis technically in logfc_f_final but need to split up by indiv so can't use this...
#		#for cis, try randomly choosing a sample? And then for trans we'll randomly choose a different sample?
#		two_sample<-sample(2:ncol(momAllele_all_filtered),2)
#		cis_sample<-two_sample[1]
#		trans_sample<-two_sample[2]
#	#	cis_mom<-as.numeric(as.character(momAllele_all_filtered[which(momAllele_all_filtered[,1]==logfc_p_final[i,1]),cis_sample]))
#	#	cis_dad<-as.numeric(as.character(dadAllele_all_filtered[which(dadAllele_all_filtered[,1]==logfc_p_final[i,1]),cis_sample]))
#	       cis_mom<-as.numeric(as.character(momAllele_all_filtered[which(momAllele_all_filtered[,1]==DE_genes_table[i,1]),cis_sample]))
#	       cis_dad<-as.numeric(as.character(dadAllele_all_filtered[which(dadAllele_all_filtered[,1]==DE_genes_table[i,1]),cis_sample]))
#		cis<-log(cis_mom/cis_dad)
#	#	trans_mom<-as.numeric(as.character(momAllele_all_filtered[which(momAllele_all_filtered[,1]==logfc_p_final[i,1]),trans_sample]))
#	#	trans_dad<-as.numeric(as.character(dadAllele_all_filtered[which(dadAllele_all_filtered[,1]==logfc_p_final[i,1]),trans_sample]))
#	       trans_mom<-as.numeric(as.character(momAllele_all_filtered[which(momAllele_all_filtered[,1]==DE_genes_table[i,1]),trans_sample]))
#	       trans_dad<-as.numeric(as.character(dadAllele_all_filtered[which(dadAllele_all_filtered[,1]==DE_genes_table[i,1]),trans_sample]))
#		trans<-parent_diff - log(trans_mom/trans_dad)
#	
#		###NEED TO FIGURE OUT A WAY TO TEST SIG DIF - specific method for sig dif logFC? Maybe need to not use EdgeR logFC values for parents and randomly choose a sample from each parent instead? Could do a binomial test for parent dif and allele dif this way, although we still have the issue of way too liberal/not using a good test for overdisp data (unless we narrowed down to just sig DE between parents genes? Although this would probably kick out everything compensatory...)
#		
#		#If parent logFC NOT sig dif from cis and/or trans not sig dif from zero, all cis
#		#If cis close to zero and/or trans NOT sig dif from parents, all trans
#		#If cis > parental and trans neg, compensatory
#		#If cis < parental and trans pos, reinforcing
#		
#		#Comparing parent fold change and ASE fold change with a proportion test for now, although I'm not sure if this is ok for gene expression data given its weird distribution
#	#	print(mom)
#	#	print(cis_mom)
#	#	print(sum(mom,dad))
#	#	print(sum(cis_mom,cis_dad))
#		result<-prop.test(x=c(mom,cis_mom),n=c(mom+dad,cis_mom+cis_dad))
#		parentVScis<-result$p.value
#		result<-prop.test(x=c(mom,trans_mom),n=c(mom+dad,trans_mom+trans_dad))
#		parentVStrans<-result$p.value
#		if(parentVScis > 0.05){
#			reg_cat_fraser2019method<-c(reg_cat_fraser2019method, "cis")
#		} else if(parentVStrans > 0.05){
#			reg_cat_fraser2019method<-c(reg_cat_fraser2019method, "trans")
#		} else if((cis>parent_diff) && (trans<0)){
#			reg_cat_fraser2019method<-c(reg_cat_fraser2019method, "compensatory")
#		} else if((cis<parent_diff) && (trans>0)){
#			reg_cat_fraser2019method<-c(reg_cat_fraser2019method, "reinforcing")
#		} else{
#			reg_cat_fraser2019method<-c(reg_cat_fraser2019method, "ambiguous")
#		}
#	}
#	fraser2019_propCis<-c(fraser2019_propCis,length(which(reg_cat_fraser2019method=="cis"))/length(reg_cat_fraser2019method))
#	fraser2019_propTrans<-c(fraser2019_propTrans,length(which(reg_cat_fraser2019method=="trans"))/length(reg_cat_fraser2019method))
#	fraser2019_propComp<-c(fraser2019_propComp,length(which(reg_cat_fraser2019method=="compensatory"))/length(reg_cat_fraser2019method))
#	fraser2019_propRein<-c(fraser2019_propRein,length(which(reg_cat_fraser2019method=="reinforcing"))/length(reg_cat_fraser2019method))
#	fraser2019_propAmbig<-c(fraser2019_propAmbig,length(which(reg_cat_fraser2019method=="ambiguous"))/length(reg_cat_fraser2019method))
##	print("Regulatory category results from Fraser 2019 method:")
##	print(paste("Proportion cis:",length(which(reg_cat_fraser2019method=="cis"))/length(reg_cat_fraser2019method)))
##	print(paste("Proportion trans:",length(which(reg_cat_fraser2019method=="trans"))/length(reg_cat_fraser2019method)))
##	print(paste("Proportion compensatory:",length(which(reg_cat_fraser2019method=="compensatory"))/length(reg_cat_fraser2019method)))
##	print(paste("Proportion reinforcing:",length(which(reg_cat_fraser2019method=="reinforcing"))/length(reg_cat_fraser2019method)))
##	print(paste("Proportion ambiguous:",length(which(reg_cat_fraser2019method=="ambiguous"))/length(reg_cat_fraser2019method)))
#}
#print("Fraser 2019 cis:")
#print(paste(mean(fraser2019_propCis),median(fraser2019_propCis),sd(fraser2019_propCis)))
#print(paste(mean(fraser2019_propTrans),median(fraser2019_propTrans),sd(fraser2019_propTrans)))
#print(paste(mean(fraser2019_propComp),median(fraser2019_propComp),sd(fraser2019_propComp)))
#print(paste(mean(fraser2019_propRein),median(fraser2019_propRein),sd(fraser2019_propRein)))
#print(paste(mean(fraser2019_propAmbig),median(fraser2019_propAmbig),sd(fraser2019_propAmbig)))


print(paste("Done with 12_fisher_and_logFCplot!",praise()))
