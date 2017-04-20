setwd("C:/Users/Devin/Documents/GradSchool/IowaState/EEOB_546X/BCB546X-Spring2017_Valerin_final_project")
getwd()
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("edgeR")
library(edgeR)
counts_data<-read.delim("Asparagus.RSEM_genecounts.txt", check.names=FALSE, stringsAsFactors = FALSE) #read in gene counts
head(counts_data)
#group<-c(Line_8A,Line_8A,Line_8A,Line_8B,Line_8B,Line_8B,Line_10, Line_10, Line_10, Line_9, Line_9)
cd<-DGEList(counts=counts_data) #convert gene count data to DGEList
keep <- rowSums(cpm(cd)>1) >= 3 #filters out lowly expressed genes
cd <- cd[keep, , keep.lib.sizes=FALSE]
cd<-calcNormFactors(cd) #normalizes
#line<-factor(c(88,88,88,89,89,89,103,103,103,9,9))
colors <- c("black","black","black","darkgreen","darkgreen","darkgreen","blue","blue","blue", "red", "red")
names<-c("female", "male","supermale","female", "male","supermale","female", "male")
plotMDS(cd, col=colors, labels=names, top = 250, gene.selection = "pairwise")
#legend('center',legend=c("Line_8A", "Line_8B","Line_10", "Line_9"), col= c("black","darkgreen","blue", "red"), ncol=2, title="Key" )


#design matrix of sex and line
rm(design)
line<-factor(c(88,88,88,89,89,89,103,103,103,9,9))
sex<-factor(c("XX","XY","YY","XX","XY","YY","XX","XY","YY","XX","XY"))
data.frame(Sample=colnames(cd),line,sex)
design <- model.matrix(~line+sex)
rownames(design) <- colnames(cd)
design

#design matrix by by sex
#rm(design)
#sex<-factor(c("XX","XY","YY","XX","XY","YY","XX","XY","YY","XX","XY"))
#data.frame(Sample=colnames(cd),sex)
#design<-model.matrix(~0+sex)
#rownames(design) <- colnames(cd)
#design

#estimates dispersion
cd <- estimateDisp(cd, design)
#gits using glm QL fit (recommended glm)
fit <- glmQLFit(cd, design)
#plotQLDisp(fit)

lrt_FSupM <- glmLRT(fit, coef=6)  #Female VS Supermale
summary(de_FemvsSupMal <- decideTestsDGE(lrt_FSupM, adjust.method="fdr"))
FDR_FSupM <- p.adjust(lrt_FSupM$table$PValue, method="fdr") 
sum(FDR_FSupM < 0.05)

detags <- rownames(cd)[as.logical(de_FemvsSupMal)]
plotSmear(lrt_FSupM, de_FemvsSupMal.tags=detags)
abline(h=c(-1, 1))

lrt_FM <- glmLRT(fit, coef=5) #female VS male
summary(de_FemvsMal <- decideTestsDGE(lrt_FM, adjust.method="fdr"))
FDR_FMale<-p.adjust(lrt_FM$table$PValue,method="fdr")
sum(FDR_FMale<0.05)
plotMD(lrt_FM)

lrt_MSupM <- glmLRT(fit, contrast=c(0,0,0,0,-1,1))
summary(de_MalvsSupMal <- decideTestsDGE(lrt_MSupM, adjust.method="fdr"))
FDR_MSupMale<-p.adjust(lrt_MSupM$table$PValue,method="fdr")
sum(FDR_MSupMale<0.05)

#qlf<-glmQLFTest(fit,coef=5:6)
#FDR<-p.adjust(qlf$table$PValue, method="BH")
#sum(FDR<0.05)
#summary(decideTests(qlf))
#qlf<-glmQLFTest(fit, contrast=c(0,0,0,-1,1)) #male VS supermale
#topTags(qlf)
#plotMD(qlf)
#design
#qlf2<-glmQLFTest(fit, coef=3)#female vs male
#FDR2<-p.adjust(qlf2$table$PValue, method="fdr")
#summary(de_FemaleVSMale <- decideTestsDGE(qlf2, adjust.method="fdr"))
#FDR_FSupM <- p.adjust(qlf2$table$PValue, method="fdr") 
#sum(FDR2<0.05)
#topTags(qlf2)
#plotMD(qlf2)
#qlf3<-glmQLFTest(fit, coef=6)
#FDR3<-p.adjust(qlf3$table$PValue, method="BH")
#sum(FDR<0.05)
#topTags(qlf3)
#plotMD(qlf3)

