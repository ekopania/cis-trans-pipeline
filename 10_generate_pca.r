#PURPOSE: generate PCAs of readcounts

args<-commandArgs(TRUE)
if(length(args) != 7){
        stop("Missing command line arguments.\nArgument 1: mom genotype (ex: WWLL)\nArgument 2: dad genotype (ex: CCPP)\nArgument 3: F1 genotype (ex: LLPP)\nArgument 4: comparison (ex: PwkLewes)\nArgument 5: allele 1 (ex: LEWES)\nArgument 6: allele 2 (ex: PWK)\nArgument 7: cell type (ex: LZ)")
}

LLLL<-list.files(pattern=paste0(args[1],"\\d*M",args[7],"\\.",args[4],"\\.",args[5],"allelesONLY[[:print:]]*\\.SE_R1only\\.",args[3],"\\.NEWdownsample\\.htseq_count\\.out$"))
PPPP<-list.files(pattern=paste0(args[2],"\\d*M",args[7],"\\.",args[4],"\\.",args[6],"allelesONLY[[:print:]]*\\.SE_R1only\\.",args[3],"\\.NEWdownsample\\.htseq_count\\.out$"))
F1s<-list.files(pattern=paste0(args[3],"\\d*M",args[7],".*allelesONLY.*\\.SE_R1orMateUnmappedOnly\\.",args[3],"\\.NEWdownsample\\.htseq_count\\.out$"))

pdf(paste("ASE_readcounts_scaled_pca_plots",args[3],args[7],"pdf",sep="."), height=8.5,width=11)

print(LLLL)
print(PPPP)
print(F1s)

exp_data<-c()
samples<-c()
genos<-c()
celltypes<-c()
for(i in 1:length(LLLL)){
	temp.tab<-read.table(LLLL[i])
	#samp=gsub("_toLEWES(.*)","",LLLL[i])
#	samp=gsub("[[:punct:]]katyaData(.*)","",LLLL[i])
#	samp=gsub("trimmomatic(.*)","",LLLL[i])
#	samp=gsub("PwkLewes(.*)","",LLLL[i])
#	samp=gsub("PwkCzii(.*)","",LLLL[i])
	samp=gsub("downsample(.*)","",LLLL[i])
	print(samp)
	celltype=gsub("\\..*","",gsub(".*\\d*M","",samp))
#	celltype=sub(".*[_]([[:alpha:]]{2})[_].*", "\\1", samp)
	print(celltype)
	exp_data<-cbind(exp_data,temp.tab[,2])
	genos<-c(genos,"LLLL")
	celltypes<-c(celltypes,celltype)
	if(i==1){
		rownames(exp_data)<-temp.tab[,1]
	}
	else if(!(all(rownames(exp_data)==temp.tab[,1]))){
		stop("Gene names don't match up!")
	}
	samples<-c(samples,samp)
}
for(i in 1:length(PPPP)){
        temp.tab<-read.table(PPPP[i])
        #samp=gsub("_toPWK(.*)","",PPPP[i])
#	samp=gsub("[[:punct:]]katyaData(.*)","",PPPP[i])
#	samp=gsub("trimmomatic(.*)","",PPPP[i])
#	samp=gsub("PwkLewes(.*)","",PPPP[i])
#	samp=gsub("PwkCzii(.*)","",PPPP[i])
	samp=gsub("downsample(.*)","",PPPP[i])
        print(samp)
#	celltype=sub(".*[_]([[:alpha:]]{2})[_].*", "\\1", samp)
	celltype=gsub("\\..*","",gsub(".*\\d*M","",samp))
        print(celltype)
	exp_data<-cbind(exp_data,temp.tab[,2])
	genos<-c(genos,"PPPP")
	celltypes<-c(celltypes,celltype)
        if(!(all(rownames(exp_data)==temp.tab[,1]))){
                stop("Gene names don't match up!")
        }
	samples<-c(samples,samp)
}
for(i in 1:length(F1s)){
        temp.tab<-read.table(F1s[i])
#	samp=gsub("[[:punct:]]katyaData(.*)","",F1s[i])
	samp=gsub("\\.SE_R1(.*)","",F1s[i])
        print(samp)
	celltype=sub("[[:alnum:]]*[M_]([[:alpha:]]{2}).*", "\\1", samp)
	geno=substr(samp,1,4)
        print(celltype)
        exp_data<-cbind(exp_data,temp.tab[,2])
	genos<-c(genos,geno)
        celltypes<-c(celltypes,celltype)
	if(!(all(rownames(exp_data)==temp.tab[,1]))){
                stop("Gene names don't match up!")
        }        
        samples<-c(samples,samp)
}

colnames(exp_data)<-samples
print(head(exp_data))
exp_data_tab<-t(exp_data)
print(dim(exp_data_tab))
#Generate shortened labels for each sample instead of printing whole htseq output filename to PCA
labs<-unlist(lapply(row.names(exp_data_tab), function(x) gsub(".htseq(.*)","",x)))
print(labs)
#remove genes with zero variance, otherwise can't calculate scaled PCA
new_exp_data<-exp_data_tab[ , apply(exp_data_tab, 2, var) != 0]
print(dim(new_exp_data))
print(new_exp_data[1:5,1:5])
exp_pca<-prcomp(new_exp_data,center=TRUE,scale.=TRUE)
print(exp_pca$x)
print(summary(exp_pca))
library(ggplot2)
###To label with sample names
ggplot(as.data.frame(exp_pca$x),aes(x=PC1,y=PC2,label=labs))+geom_text(size=3)+labs(title="Gene Expression Read Count PCA")
###To color by genotype and shape by celltype
#ggplot(as.data.frame(exp_pca$x),aes(x=PC1,y=PC2,shape=celltypes,color=genos))+geom_point(size=3)+labs(title="Gene Expression Read Count PCA")
dev.off()
