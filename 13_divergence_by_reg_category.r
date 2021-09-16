#PURPOSE: Test for differences in the median divergence level across different expression categories (ex: ask if cis-regulatory changes are associated with higer divergence at the phylogenetic level than trans-regulatory changes, do this across all categories)

#Read in command line args
args<-commandArgs(TRUE)
if(length(args) != 4){
        stop("Missing command line arguments. Enter two files of regulatory categories, one for LZ and one for RS (output from script 12), and enter two files of logFC (DE_genes*parentTableAll.txt output from script 11).")
}
print("Here are the input regulatory category files:")
print(args)

geno<-gsub(".LZ.txt","",gsub("reg_category_output.[[:alpha:]]{4}mom.","",args[1]))
print(geno)

#R Libraries
library(ggplot2)

#Read in divergence values from EVE (beta)
print("Reading in divergence values...")
mydata<-read.table("../mus_expression_analysis/MULTI_MAP/results/indivBetaMLparams_ensemblOrthos_BOTHinduced.res")
sharedData<-read.table("../mus_expression_analysis/MULTI_MAP/results/sharedBetaMLparams_ensemblOrthos_BOTHinduced.res")
print(head(mydata))

betai<-cbind(betai=0-log(mydata[,4]),group=c(rep("LZ",3375),rep("RS",2769)))
betaShared<-0-log(sharedData[1,4])
print(dim(betai))
print(head(betai))
print(betaShared)

#TODO: include eve div values for genes expressed (not just induced) in each category (from ../mus_expression_analysis/MULTI_MAP/results/indivBetaMLparams_ensemblOrthos_LZ.res and ../mus_expression_analysis/MULTI_MAP/results/indivBetaMLparams_ensemblOrthos_RS.res)
#TODO: include logFC values (plot this as another form of "expression divergence") (from DE_genes*.LRT.parentTableAll.txt)
#Read in logFC tables
print("Reading in logFC table...")
logFCdata_LZ<-read.table(args[3], header=TRUE)
logFCdata_RS<-read.table(args[4], header=TRUE)
#TODO: include dN/dS values (plot protein coding divergence by reg category) (from ../PAML_runs/WHOLE_GENOMES/omega_list.ensemblOrthos.txt) 
#Read in dN/dS table
print("Reading in dN/dS values...")
dndsData<-read.table("../PAML_runs/WHOLE_GENOMES/omega_list.ensemblOrthos.txt")

#Read in regulatory categories (output from script 11)
print("Reading in regulatory categories...")
reg_cats_LZ<-read.table(args[1],header=TRUE) 
print(dim(reg_cats_LZ))
print(head(reg_cats_LZ))
reg_cats_RS<-read.table(args[2],header=TRUE)
print(dim(reg_cats_RS))
print(head(reg_cats_RS))

#Make a giant table that associates each gene with its divergence value and reg category
LZgenes<-scan("gene_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.txt",what=character())
RSgenes<-scan("gene_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.txt",what=character())
gene_names<-c(LZgenes,RSgenes)
betaFinal<-as.data.frame(cbind(gene_names, betai),stringsAsFactors=FALSE)
print(dim(betaFinal))
print(length(intersect(betaFinal$gene_names,c(reg_cats_LZ$gene,reg_cats_RS$gene))))
comboData<-c()
#for(i in 1:100){ #For testing
for(i in 1:nrow(betaFinal)){
	this_gene<-betaFinal$gene_name[i]
	if( (this_gene %in% reg_cats_LZ$gene) && (betaFinal$group[i]=="LZ") ){
		cat<-as.character(reg_cats_LZ$gene_category[which(reg_cats_LZ$gene==this_gene)])
		temp_line<-c(unlist(as.vector(betaFinal[i,])),cat,paste(betaFinal$group[i],cat,sep="."))
		logFC<-as.character(logFCdata_LZ$logFC[which(rownames(logFCdata_LZ)==this_gene)])
		logFC_abs<-abs(as.numeric(logFC))
		temp_line<-c(temp_line,logFC,logFC_abs)
	} else if( (this_gene %in% reg_cats_RS$gene) && (betaFinal$group[i]=="RS") ){
		cat<-as.character(reg_cats_RS$gene_category[which(reg_cats_RS$gene==this_gene)])
                temp_line<-c(unlist(as.vector(betaFinal[i,])),cat,paste(betaFinal$group[i],cat,sep="."))
		logFC<-as.character(logFCdata_RS$logFC[which(rownames(logFCdata_RS)==this_gene)])
		logFC_abs<-abs(as.numeric(logFC))
                temp_line<-c(temp_line,logFC,logFC_abs)
	} else{
		temp_line<-c(unlist(as.vector(betaFinal[i,])),"none",paste(betaFinal$group[i],"none",sep="."),NA,NA)
	}
	if(this_gene %in% dndsData$V1){
		dnds<-dndsData$V4[which(dndsData$V1==this_gene)]
		temp_line<-c(temp_line,dnds)
	} else{
		temp_line<-c(temp_line,NA)
	}
#	print(temp_line)
	comboData<-rbind(comboData,temp_line)
}
dataFinal<-as.data.frame(comboData)
colnames(dataFinal)<-c("gene","betai","group","categ","group_cat","logFC","abs_logFC","dNdS")
print(dim(dataFinal))
print(head(dataFinal))
dataFinal_LZ<-dataFinal[which(dataFinal$group=="LZ"),]
dataFinal_RS<-dataFinal[which(dataFinal$group=="RS"),]
print(dim(dataFinal_LZ))
print(dim(dataFinal_RS))

####ABSOLUTE VALUE OF LOGFC!!!

#These are necessary for the pairwise wilcoxon test, otherwise only "NA" is associated with the "none" category and it causes an error in the test
dataFinal_LZ_noNone<-dataFinal_LZ[which(dataFinal_LZ$categ!="none"),]
dataFinal_RS_noNone<-dataFinal_RS[which(dataFinal_RS$categ!="none"),]
#Get rid of ridiculous dN/dS values
dataFinal_LZ_filterDNDS<-dataFinal_LZ[which(as.numeric(as.character(dataFinal_LZ$dNdS))<1.5),]
dataFinal_RS_filterDNDS<-dataFinal_RS[which(as.numeric(as.character(dataFinal_RS$dNdS))<1.5),]

#Perform pairwise wilcoxon test
print("Performing pairwise wilcoxon test...")
result_div_LZ<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_LZ$betai)), dataFinal_LZ$group_cat, p.adjust.method = "fdr")
result_div_RS<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_RS$betai)), dataFinal_RS$group_cat, p.adjust.method = "fdr")
print("logFCLZ")
result_logFC_LZ<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_LZ_noNone$logFC)), dataFinal_LZ_noNone$group_cat, p.adjust.method = "fdr")
print("logFCRS")
result_logFC_RS<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_RS_noNone$logFC)), dataFinal_RS_noNone$group_cat, p.adjust.method = "fdr")
print("abslogFCLZ")
result_logFCabs_LZ<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_LZ_noNone$abs_logFC)), dataFinal_LZ_noNone$group_cat, p.adjust.method = "fdr")
print("abslogFCRS")
result_logFCabs_RS<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_RS_noNone$abs_logFC)), dataFinal_RS_noNone$group_cat, p.adjust.method = "fdr")
result_dNdS_LZ<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_LZ_filterDNDS$dNdS)), dataFinal_LZ_filterDNDS$group_cat, p.adjust.method = "fdr")
result_dNdS_RS<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_RS_filterDNDS$dNdS)), dataFinal_RS_filterDNDS$group_cat, p.adjust.method = "fdr")
print("LZ EVE divergence:")
print(result_div_LZ)
print("RS EVE divergence:")
print(result_div_RS)
print("LZ logFC:")
print(result_logFC_LZ)
print("RS logFC:")
print(result_logFC_RS)
print("LZ absolute value of logFC:")
print(result_logFCabs_LZ)
print("RS absolute value of logFC:")
print(result_logFCabs_RS)
print("LZ dN/dS:")
print(result_dNdS_LZ)
print("RS dN/dS:")
print(result_dNdS_RS)

#Plot - maybe violin plot for each reg category and each cell type
print("Plotting LZ...")
#get rid of outliers for nicer plots
dataFinal_LZ<-dataFinal_LZ[which(as.numeric(as.character(dataFinal_LZ$betai)) > -5),]
p_div<-ggplot(dataFinal_LZ, aes(x=categ, y=as.numeric(as.character(betai)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Gene expression divergence (EVEmodel) for different regulatory categories - early", x="Regulatory Category", y="Expression Divergence")
p_div<-p_div + theme(axis.text.y = element_text(size=20))
p_div<-p_div + theme_minimal()
p_div<-p_div + scale_fill_manual(values=c("chocolate1"))
p_logFC<-ggplot(dataFinal_LZ_noNone, aes(x=categ, y=as.numeric(as.character(logFC)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Expression logFC for different regulatory categories - early", x="Regulatory Category", y="logFC")
p_logFC<-p_logFC + theme(axis.text.y = element_text(size=20))
p_logFC<-p_logFC + theme_minimal()
p_logFC<-p_logFC + scale_fill_manual(values=c("chocolate1"))
p_logFC_abs<-ggplot(dataFinal_LZ_noNone, aes(x=categ, y=as.numeric(as.character(abs_logFC)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Absolute value of logFC for different regulatory categories - early", x="Regulatory Category", y="absolute value of logFC")
p_logFC_abs<-p_logFC_abs + theme(axis.text.y = element_text(size=20))
p_logFC_abs<-p_logFC_abs + theme_minimal()
p_logFC_abs<-p_logFC_abs + scale_fill_manual(values=c("chocolate1"))
p_dnds<-ggplot(dataFinal_LZ_filterDNDS, aes(x=categ, y=as.numeric(as.character(dNdS)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="dN/dS for different regulatory categories - early", x="Regulatory Category", y="dN/dS")
p_dnds<-p_dnds + theme(axis.text.y = element_text(size=20))
p_dnds<-p_dnds + theme_minimal()
p_dnds<-p_dnds + scale_fill_manual(values=c("chocolate1"))
pdf(paste0("divergence_violins.byRegCategory.",geno,".LZ.pdf"),onefile=TRUE,width=11,height=8.5)
p_div
p_logFC
p_logFC_abs
p_dnds
dev.off()
#ggsave("EVEmodel_divergence_violins.byRegCategory.LZ.pdf",p,width=11,height=8.5,units="in")

print("Plotting RS...")
#get rid of outliers for nicer plots
dataFinal_RS<-dataFinal_RS[which(as.numeric(as.character(dataFinal_RS$betai)) > -5),]
p_div<-ggplot(dataFinal_RS, aes(x=categ, y=as.numeric(as.character(betai)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Gene expression divergence (EVEmodel) for different regulatory categories - late", x="Regulatory Category", y="Expression Divergence")
p_div<-p_div + theme(axis.text.y = element_text(size=20))
p_div<-p_div + theme_minimal()
p_div<-p_div + scale_fill_manual(values=c("lightsteelblue1"))
p_logFC<-ggplot(dataFinal_RS_noNone, aes(x=categ, y=as.numeric(as.character(logFC)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Expression logFC for different regulatory categories - late", x="Regulatory Category", y="logFC")
p_logFC<-p_logFC + theme(axis.text.y = element_text(size=20))
p_logFC<-p_logFC + theme_minimal()
p_logFC<-p_logFC + scale_fill_manual(values=c("lightsteelblue1"))
p_logFC_abs<-ggplot(dataFinal_RS_noNone, aes(x=categ, y=as.numeric(as.character(abs_logFC)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Absolute value of logFC for different regulatory categories - late", x="Regulatory Category", y="absolute value of logFC")
p_logFC_abs<-p_logFC_abs + theme(axis.text.y = element_text(size=20))
p_logFC_abs<-p_logFC_abs + theme_minimal()
p_logFC_abs<-p_logFC_abs + scale_fill_manual(values=c("lightsteelblue1"))
p_dnds<-ggplot(dataFinal_RS_filterDNDS, aes(x=categ, y=as.numeric(as.character(dNdS)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="dN/dS for different regulatory categories - late", x="Regulatory Category", y="dN/dS")
p_dnds<-p_dnds + theme(axis.text.y = element_text(size=20))
p_dnds<-p_dnds + theme_minimal()
p_dnds<-p_dnds + scale_fill_manual(values=c("lightsteelblue1"))
pdf(paste0("divergence_violins.byRegCategory.",geno,".RS.pdf"),onefile=TRUE,width=11,height=8.5)
p_div
p_logFC
p_logFC_abs
p_dnds
dev.off()

#Plot and perform pairwise wilcoxon test for more "cleaned up" versions - no other, none, or cXt; just abs(logFC) and eve exp div
#keep<-c(which(dataFinal_LZ$categ=="cis"),which(dataFinal_LZ$categ=="trans"),which(dataFinal_LZ$categ=="conserved"),which(dataFinal_LZ$categ=="comp"),which(dataFinal_LZ$categ=="cPlusTOpp"),which(dataFinal_LZ$categ=="cPlusTSame"))
#dataFinal_LZ_cleaned<-dataFinal_LZ[keep,]
#keep<-c(which(dataFinal_RS$categ=="cis"),which(dataFinal_RS$categ=="trans"),which(dataFinal_RS$categ=="conserved"),which(dataFinal_RS$categ=="comp"),which(dataFinal_RS$categ=="cPlusTOpp"),which(dataFinal_RS$categ=="cPlusTSame"))
#dataFinal_RS_cleaned<-dataFinal_RS[keep,]
#x_order<-c("conserved","cis","trans","compensatory","cPlusTOpp","cPlusTSame")
#To plot cis and trans only
#keep<-c(which(dataFinal_LZ$categ=="cis"),which(dataFinal_LZ$categ=="trans"))
#dataFinal_LZ_cleaned<-dataFinal_LZ[keep,]
#keep<-c(which(dataFinal_RS$categ=="cis"),which(dataFinal_RS$categ=="trans"))
#dataFinal_RS_cleaned<-dataFinal_RS[keep,]
#x_order<-c("cis","trans")
#Cis, trans, cis+trans opp, cis+trans same
keep<-c(which(dataFinal_LZ$categ=="cis"),which(dataFinal_LZ$categ=="trans"),which(dataFinal_LZ$categ=="cPlusTOpp"),which(dataFinal_LZ$categ=="cPlusTSame"))
dataFinal_LZ_cleaned<-dataFinal_LZ[keep,]
keep<-c(which(dataFinal_RS$categ=="cis"),which(dataFinal_RS$categ=="trans"),which(dataFinal_RS$categ=="cPlusTOpp"),which(dataFinal_RS$categ=="cPlusTSame"))
dataFinal_RS_cleaned<-dataFinal_RS[keep,]
x_order<-c("cis","trans","cPlusTOpp","cPlusTSame")

print(dim(dataFinal_LZ_cleaned))
print(dim(dataFinal_RS_cleaned))
print(head(dataFinal_LZ_cleaned))
print(head(dataFinal_RS_cleaned))

#Plot LZ
p_div<-ggplot(dataFinal_LZ_cleaned, aes(x=categ, y=as.numeric(as.character(betai)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Gene expression divergence (EVEmodel) for different regulatory categories - early", x="Regulatory Category", y="Expression Divergence")
p_div<-p_div + ylim(-6, 5) + scale_x_discrete(limits=x_order) #ylim(-30,0)
p_div<-p_div + theme(axis.text.y = element_text(size=20))
p_div<-p_div + theme_minimal()
p_div<-p_div + scale_fill_manual(values=c("chocolate1"))
p_logFC_abs<-ggplot(dataFinal_LZ_cleaned, aes(x=categ, y=as.numeric(as.character(abs_logFC)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Absolute value of logFC for different regulatory categories - early", x="Regulatory Category", y="absolute value of logFC")
p_logFC_abs<-p_logFC_abs + ylim(0,4) + scale_x_discrete(limits=x_order)
p_logFC_abs<-p_logFC_abs + theme(axis.text.y = element_text(size=20))
p_logFC_abs<-p_logFC_abs + theme_minimal()
p_logFC_abs<-p_logFC_abs + scale_fill_manual(values=c("chocolate1"))
pdf(paste0("divergence_violins.byRegCategory.cleaned.noOutliers.",geno,".LZ.pdf"),onefile=TRUE,width=11,height=8.5)
p_div
p_logFC_abs
dev.off()

#Plot RS
p_div<-ggplot(dataFinal_RS_cleaned, aes(x=categ, y=as.numeric(as.character(betai)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Gene expression divergence (EVEmodel) for different regulatory categories - late", x="Regulatory Category", y="Expression Divergence")
p_div<-p_div + ylim(-6, 5) + scale_x_discrete(limits=x_order) #ylim(-30,0)
p_div<-p_div + theme(axis.text.y = element_text(size=20))
p_div<-p_div + theme_minimal()
p_div<-p_div + scale_fill_manual(values=c("lightsteelblue1"))
p_logFC_abs<-ggplot(dataFinal_RS_cleaned, aes(x=categ, y=as.numeric(as.character(abs_logFC)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Absolute value of logFC for different regulatory categories - late", x="Regulatory Category", y="absolute value of logFC")
p_logFC_abs<-p_logFC_abs + ylim(0,4) + scale_x_discrete(limits=x_order) 
p_logFC_abs<-p_logFC_abs + theme(axis.text.y = element_text(size=20))
p_logFC_abs<-p_logFC_abs + theme_minimal()
p_logFC_abs<-p_logFC_abs + scale_fill_manual(values=c("lightsteelblue1"))
pdf(paste0("divergence_violins.byRegCategory.cleaned.noOutliers.",geno,".RS.pdf"),onefile=TRUE,width=11,height=8.5)
p_div
p_logFC_abs
dev.off()

print("Performing pairwise wilcoxon test on subset of reg categories (cis, trans, comp, cons, c+tOpp, c+tSame)...")
result_div_LZ<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_LZ_cleaned$betai)), dataFinal_LZ_cleaned$group_cat, p.adjust.method = "fdr")
result_div_RS<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_RS_cleaned$betai)), dataFinal_RS_cleaned$group_cat, p.adjust.method = "fdr")
result_logFCabs_LZ<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_LZ_cleaned$abs_logFC)), dataFinal_LZ_cleaned$group_cat, p.adjust.method = "fdr")
result_logFCabs_RS<-pairwise.wilcox.test(as.numeric(as.character(dataFinal_RS_cleaned$abs_logFC)), dataFinal_RS_cleaned$group_cat, p.adjust.method = "fdr")
print("LZ EVE divergence:")
print(result_div_LZ)
print("RS EVE divergence:")
print(result_div_RS)
print("LZ abs(logFC):")
print(result_logFCabs_LZ)
print("RS abs(logFC):")
print(result_logFCabs_RS)


print("Done with 14_divergence_by_reg_category")
