##EEOB_546X
##Final Project
##Members:Devin,Valeria,Katerina
#The following is the script for expression data analysis

#set directory
setwd("C:/Users/Devin/Documents/GradSchool/IowaState/EEOB_546X/BCB546X-Spring2017_Valerin_final_project")
getwd()
#import bioconductor and required libraries
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("edgeR")
library(edgeR)

#read in the gene counts from paper
counts_data<-read.delim("Asparagus.RSEM_genecounts.txt", check.names=FALSE, stringsAsFactors = FALSE, row.names = TRUE) #read in gene counts
head(counts_data) #take a look at counts_data

cd<-DGEList(counts=counts_data) #convert gene count data to DGEList
keep <- rowSums(cpm(cd)>1) >= 3 #filters out lowly expressed genes
cd <- cd[keep, , keep.lib.sizes=FALSE]
cd<-calcNormFactors(cd) #normalizes
nrow(cd) #verify that we have 19980 genes
#generates first figure (beautified by Katerina because she is awesome)
colors <- c("black","black","black","darkgreen","darkgreen","darkgreen","blue","blue","blue", "red", "red")
names<-c("female", "male","supermale","female", "male","supermale","female", "male")
plotMDS(cd, col=colors, labels=names, top = 250, gene.selection = "pairwise")
#legend('center',legend=c("Line_8A", "Line_8B","Line_10", "Line_9"), col= c("black","darkgreen","blue", "red"), ncol=2, title="Key" )


#design matrix of sex and line
rm(design)
line<-factor(c(88,88,88,89,89,89,103,103,103,9,9))
sex<-factor(c("XX","XY","YY","XX","XY","YY","XX","XY","YY","XX","XY"))
#data.frame(Sample=colnames(cd),line,sex)
design <- model.matrix(~line+sex)
rownames(design) <- colnames(cd)
design

#estimates dispersion
#cd <- estimateDisp(cd, design)
cd <- estimateGLMCommonDisp(cd, design)
cd <- estimateGLMTrendedDisp(cd, design)
cd <- estimateGLMTagwiseDisp(cd, design)
#gits using glm QL fit (recommended glm)
#fit <- glmQLFit(cd, design) #preferred but lenient
fit <- glmFit(cd, design) #more restrictive glm
#plotQLDisp(fit)

###COMPARISONS####
lrt_FSupM <- glmLRT(fit, coef=6)  #Female VS Supermale
summary(de_FemvsSupMal <- decideTestsDGE(lrt_FSupM, adjust.method="fdr"))
FDR_FSupM <- p.adjust(lrt_FSupM$table$PValue, method="fdr") 
sum(FDR_FSupM < 0.05)
plotMD(lrt_FSupM)

lrt_FM <- glmLRT(fit, coef=5) #female VS male
summary(de_FemvsMal <- decideTestsDGE(lrt_FM, adjust.method="fdr"))
FDR_FMale<-p.adjust(lrt_FM$table$PValue,method="fdr")
sum(FDR_FMale<0.05)
plotMD(lrt_FM)

lrt_MSupM <- glmLRT(fit, contrast=c(0,0,0,0,-1,1)) #Male VS Supermale
summary(de_MalvsSupMal <- decideTestsDGE(lrt_MSupM, adjust.method="fdr"))
FDR_MSupMale<-p.adjust(lrt_MSupM$table$PValue,method="fdr")
sum(FDR_MSupMale<0.05)
plotMD(lrt_MSupM)


##Get the names of all of the genes that have a Pvalue<0.05 for FDR_FSupM
total_FDR_FSupM<-length(FDR_FSupM)
total_FDR_FSupM
names_of_FSupM_hit<-c()
counter<-0
for (i in 1:total_FDR_FSupM){
  if (FDR_FSupM[i]<0.05 ){
    counter<-counter+1
    names_of_FSupM_hit[counter]=rownames(lrt_FSupM)[i]
  }
}
names_of_FSupM_hit
length(names_of_FSupM_hit)

##Get the names of all of the genes that have a Pvalue<0.05 for FDR_FMale
total_FDR_FMale<-length(FDR_FMale)
total_FDR_FMale
names_of_FMale_hit<-c()
counter<-0
for (i in 1:total_FDR_FMale){
  if (FDR_FMale[i]<0.05 ){
    counter<-counter+1
    names_of_FMale_hit[counter]=rownames(lrt_FM)[i]
  }
}
names_of_FMale_hit
length(names_of_FMale_hit)


##Get the names of all of the genes that have a Pvalue<0.05 for FDR_FMale
total_FDR_MSupMale<-length(FDR_MSupMale)
total_FDR_MSupMale
names_of_MSupMale_hit<-c()
counter<-0
for (i in 1:total_FDR_MSupMale){
  if (FDR_MSupMale[i]<0.05 ){
    counter<-counter+1
    names_of_MSupMale_hit[counter]=rownames(lrt_MSupM)[i]
  }
}
names_of_MSupMale_hit
length(names_of_MSupMale_hit)



###Venn Diagram###
library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 = 22, area2 = 20, area3 = 13, n12 = 11, n23 = 4, n13 = 5, 
                 n123 = 1, category = c("Dog People", "Cat People", "Lizard People"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))