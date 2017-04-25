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
counts_data<-read.delim("Asparagus.RSEM_genecounts.txt", check.names=FALSE, stringsAsFactors = FALSE) #read in gene counts
head(counts_data) #take a look at counts_data

cd<-DGEList(counts=counts_data) #convert gene count data to DGEList
keep <- rowSums(cpm(cd)>1) >= 3 #filters out lowly expressed genes
cd <- cd[keep, , keep.lib.sizes=FALSE]
#cd<-calcNormFactors(cd) #normalizes
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
sum(FDR_FSupM < 0.05) #462 differentially expressed genes
plotMD(lrt_FSupM)

lrt_FM <- glmLRT(fit, coef=5) #female VS male
summary(de_FemvsMal <- decideTestsDGE(lrt_FM, adjust.method="fdr"))
FDR_FMale<-p.adjust(lrt_FM$table$PValue,method="fdr")
sum(FDR_FMale<0.05) #263
plotMD(lrt_FM)

lrt_MSupM <- glmLRT(fit, contrast=c(0,0,0,0,-1,1)) #Male VS Supermale
summary(de_MalvsSupMal <- decideTestsDGE(lrt_MSupM, adjust.method="fdr"))
FDR_MSupMale<-p.adjust(lrt_MSupM$table$PValue,method="fdr")
sum(FDR_MSupMale<0.05) #29
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
num_FSupM_hit<-length(names_of_FSupM_hit)
num_FSupM_hit #462

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
num_FMale_hit<-length(names_of_FMale_hit)
num_FMale_hit #263

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
num_MSupMale_hit<-length(names_of_MSupMale_hit)
num_MSupMale_hit #29

total_number_of_DF_genes<-sum(c(num_FMale_hit,num_FSupM_hit,num_MSupMale_hit), na.is=FALSE)
total_number_of_DF_genes


#FMale_FSupM<-table(names_of_FMale_hit[names_of_FMale_hit %in% names_of_MSupMale_hit])
#FMale_FSupM
FMale_FSupM<-names_of_FMale_hit %in% names_of_FSupM_hit
sum(FMale_FSupM)
num_FMale_FSupM<-sum(FMale_FSupM) #166
#length(FMale_FSupM[FMale_FSupM==TRUE])
length(FMale_FSupM) #263

FMale_MSupMale<-names_of_FMale_hit %in% names_of_MSupMale_hit
sum(FMale_MSupMale)
num_FMale_MSupMale<-sum(FMale_MSupMale) #2
length(FMale_MSupMale) #263 #check length of total to make sure it matches
which_FMale_MSupMale<-which(names_of_FMale_hit %in% names_of_MSupMale_hit)
which_FMale_MSupMale # 177 251, indexes of the found names
#grep_Fmale_MSupMale<-grep(TRUE,names_of_FMale_hit %in% names_of_MSupMale_hit)
#grep_Fmale_MSupMale
names_of_Fmale_MSupMale<-c(names_of_FMale_hit[177],names_of_FMale_hit[251])
names_of_Fmale_MSupMale #names of the matching genes

FSupM_MSupM<-names_of_FSupM_hit %in% names_of_MSupMale_hit
sum(FSupM_MSupM) 
num_FSupM_MSupM<-sum(FSupM_MSupM) #24

FSupM_MSupM_FMale<-names_of_FMale_hit %in% names_of_FSupM_hit %in% names_of_MSupMale_hit
sum(FSupM_MSupM_FMale)
num_FSupM_MSupM_FMale<-sum(FSupM_MSupM_FMale) #0

compare_all<-names_of_Fmale_MSupMale %in% names_of_FSupM_hit
sum(compare_all) #0
###Checking for unique names###
num_FMale_hit+num_FSupM_hit+num_MSupMale_hit #754 genes total with duplicates
all_DF_genes<-c(names_of_FMale_hit,names_of_FSupM_hit, names_of_MSupMale_hit)
unique_genes<-unique(all_DF_genes)
unique_genes
length(unique_genes) #481 unique differentially expressed genes

###Compare to genes in paper### #script from Kat
paper_fsupm <- read.delim("paper_fsupm.txt") #extracted gene_id column from supplementary table 1
paper_fm <- read.delim("paper_fm.txt") # extracted gene_id column from supplementary table 1
paper_msupm <- read.delim("paper_msupm.txt") #extracted gene_id column from supplementary table 1

paper_detags <- rbind(paper_fsupm, paper_fm, paper_msupm) #combine all names from DF data of paper
typeof(paper_detags) #type list, we want a vector
uniq_paper_detags <- unique(paper_detags) #grab just the unique names
uniq_paper_detags <-as.vector(uniq_paper_detags) #convert the list to a vector
uniq_paper_detags
typeof(uniq_paper_detags) #vector now!
length(uniq_paper_detags) #570 DF genes from paper

###Compare our unique genes to the paper's###
compare_to_paper<-uniq_paper_detags %in% unique_genes 
sum(compare_to_paper) #516 genes in the 570 unique from the paper
#when normalized, 476 of our 481 differential genes match the paper's 570 differential genes

###Venn Diagram###
install.packages('VennDiagram')
library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 = num_FMale_hit, area2 = num_FSupM_hit, area3 = num_MSupMale_hit, n12 = num_FMale_FSupM, n23 = num_FSupM_MSupM, n13 = num_FMale_MSupMale, 
                 n123 = num_FSupM_MSupM_FMale, category = c("Male vs Female", "Supermale VS female", "Supermale VS male"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))

