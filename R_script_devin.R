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
cd<-calcNormFactors(cd) #normalizes
<<<<<<< HEAD
<<<<<<< HEAD

=======
>>>>>>> 6525f3c85b3ff9749df8bf943f406d46b74f65f8
>>>>>>> 7855947ecb9001ee36647d07fdfecb2d38cad1e5
=======

>>>>>>> b6a7f296419edf6307028d7363142d2b54372506
nrow(cd) #verify that we have 19980 genes
#generates first figure (legend perfected by Kat)
colors <- c("black","black","black","darkgreen","darkgreen","darkgreen","blue","blue","blue", "red", "red")
plotMDS(cd, col=colors, top = 250, gene.selection = "pairwise", main="MDS plot of various lines and sex of Asparagus officinalis")  #MDS plot without legend but the names indicate line and sex
#names<-c("female", "male","supermale","female", "male","supermale","female", "male")
#plotMDS(cd, col=colors, labels=names, top = 250, gene.selection = "pairwise")
#my initial legend below
#legend("center",legend=c("Line_8A", "Line_8B","Line_10", "Line_9"), col= c("black","darkgreen","blue", "red"), ncol=2, title="Key" )


#design matrix of sex and line
rm(design)
line<-factor(c(88,88,88,89,89,89,103,103,103,9,9))
#line<-factor(c("Line_8A","Line_8A","Line_8A", "Line_8B","Line_8B","Line_8B","Line_10","Line_10","Line_10","Line_9","Line_9"))
sex<-factor(c("XX","XY","YY","XX","XY","YY","XX","XY","YY","XX","XY"))
df<-data.frame(Sample=colnames(cd),line,sex)
df
is.factor(df$line)
is.factor(df$sex)
#design <- model.matrix(~sex+sex:line, data=df) #nested data frame
#design <- model.matrix(~line+line:sex, data=df) #nested data frame
#design <- model.matrix(~sex*line, data=df)
design <- model.matrix(~line+sex,data=df)
summary(design)
rownames(design) <- colnames(cd)
design

colnames(design)




#estimates dispersion
#cd <- estimateDisp(cd, design) #if we use this, even though it is supposed to be equivalent to the three functions below, we will get far lower numbers
cd <- estimateGLMCommonDisp(cd, design)
cd <- estimateGLMTrendedDisp(cd, design)
cd <- estimateGLMTagwiseDisp(cd, design)
#gits using glm QL fit (recommended glm)
#fit <- glmQLFit(cd, design) #preferred but lenient
fit <- glmFit(cd, design) #more restrictive glm


###COMPARISONS####
#lrt_FSupM <- glmQLFTest(fit, coef=6) #gives us just one differentially expressed gene
lrt_FSupM <- glmLRT(fit, coef=6)  #Female VS Supermale
summary(de_FemvsSupMal <- decideTestsDGE(lrt_FSupM, adjust.method="fdr"))
FDR_FSupM <- p.adjust(lrt_FSupM$table$PValue, method="fdr") 
sum(FDR_FSupM < 0.05) #408  differentially expressed genes (462 if not normalized)
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> b6a7f296419edf6307028d7363142d2b54372506
status<-rep("not significant",length(FDR_FSupM))
status[FDR_FSupM<0.05 & lrt_FSupM$table$logFC>0]<-"up"
status[FDR_FSupM<0.05 & lrt_FSupM$table$logFC<0]<-"down"
plotMD(lrt_FSupM, main="SuperMale VS Female",hl.col=c("red","blue"),status = status, bg.col="grey")
#abline(h=c(-1,1), col="blue")
sum(FDR_FSupM<0.05 & lrt_FSupM$table$logFC>0)
sum(FDR_FSupM<0.05 & lrt_FSupM$table$logFC<0)
sum(FDR_FSupM<0.05 & lrt_FSupM$table$logFC>0)+sum(FDR_FSupM<0.05 & lrt_FSupM$table$logFC<0)
<<<<<<< HEAD

=======
<<<<<<< HEAD

=======
plotMD(lrt_FSupM, main="Female VS SuperMale",hl.col=c("blue","red"), bg.col="grey")
abline(h=c(-1,1), col="blue")
>>>>>>> 6525f3c85b3ff9749df8bf943f406d46b74f65f8
>>>>>>> 7855947ecb9001ee36647d07fdfecb2d38cad1e5
=======

>>>>>>> b6a7f296419edf6307028d7363142d2b54372506


lrt_FM <- glmLRT(fit, coef=5) #female VS male
#lrt_FM <- glmLRT(fit, contrast=c(-1,0,0,0,1,0))
summary(de_FemvsMal <- decideTestsDGE(lrt_FM, adjust.method="fdr"))
FDR_FMale<-p.adjust(lrt_FM$table$PValue,method="fdr") #adjust Pvalue
sum(FDR_FMale<0.05) #221 significant differentially expressed genes
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> b6a7f296419edf6307028d7363142d2b54372506
status<-rep("not significant",length(FDR_FMale))
status[FDR_FMale<0.05 & lrt_FM$table$logFC>0]<-"up"
status[FDR_FMale<0.05 & lrt_FM$table$logFC<0]<-"down"
plotMD(lrt_FM, main="Male VS Female",hl.col=c("red","blue"),status = status, bg.col="grey")
#abline(h=c(-1,1), col="blue")
<<<<<<< HEAD

=======
plotMD(lrt_FM, main="Female VS Male")
<<<<<<< HEAD

=======
>>>>>>> 6525f3c85b3ff9749df8bf943f406d46b74f65f8
>>>>>>> 7855947ecb9001ee36647d07fdfecb2d38cad1e5
=======

>>>>>>> b6a7f296419edf6307028d7363142d2b54372506

lrt_MSupM <- glmLRT(fit, contrast=c(0,0,0,0,-1,1)) #Male VS Supermale
summary(de_MalvsSupMal <- decideTestsDGE(lrt_MSupM, adjust.method="fdr"))
FDR_MSupMale<-p.adjust(lrt_MSupM$table$PValue,method="fdr")
sum(FDR_MSupMale<0.05) #27
status<-rep("not significant",length(FDR_MSupMale))
status[FDR_MSupMale<0.05 & lrt_MSupM$table$logFC>0]<-"up"
status[FDR_MSupMale<0.05 & lrt_MSupM$table$logFC<0]<-"down"
plotMD(lrt_MSupM, main="Supermale VS Male",hl.col=c("red","blue"),status = status, bg.col="grey")
#abline(h=c(-1,1), col="blue")

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
num_FSupM_hit #408 

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
num_FMale_hit #221 DF genes in Female/Male pairwise comparison

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
num_MSupMale_hit #27

total_number_of_DF_genes<-sum(c(num_FMale_hit,num_FSupM_hit,num_MSupMale_hit), na.is=FALSE)
total_number_of_DF_genes #656 differentially expressed genes (including duplicates)

###Checks the number of hits between FSupM and MSupMale
FSupM_MSupM<-names_of_FSupM_hit %in% names_of_MSupMale_hit
sum(FSupM_MSupM) 
num_FSupM_MSupM<-sum(FSupM_MSupM) #23 common DF genes shared


FMale_FSupM<-names_of_FMale_hit %in% names_of_FSupM_hit
sum(FMale_FSupM)
num_FMale_FSupM<-sum(FMale_FSupM) #150 differential genes shared between Female/Male and Female/Supermale
#length(FMale_FSupM[FMale_FSupM==TRUE])
length(FMale_FSupM) #221 genes in vector, checking the length of the logic vector is equivalent to the number of differentially expressed genes in FMale hits

#FMale_MSupMale mini vector for all_compare
FMale_MSupMale<-names_of_FMale_hit %in% names_of_MSupMale_hit
sum(FMale_MSupMale)
num_FMale_MSupMale<-sum(FMale_MSupMale) #2


length(FMale_MSupMale) #221 #check length of total to make sure it matches
which_FMale_MSupMale<-which(names_of_FMale_hit %in% names_of_MSupMale_hit)
which_FMale_MSupMale # returns indexes of the found names
names_of_Fmale_MSupMale<-c(names_of_FMale_hit[which_FMale_MSupMale[1]],names_of_FMale_hit[which_FMale_MSupMale[2]])
names_of_Fmale_MSupMale #names of the matching genes between Fmale and MSupMale groups


###compare all groups###
FSupM_MSupM_FMale<-names_of_FMale_hit %in% names_of_FSupM_hit %in% names_of_MSupMale_hit
sum(FSupM_MSupM_FMale)
num_FSupM_MSupM_FMale<-sum(FSupM_MSupM_FMale) #0 DF genes shared between all three groups
###verify the comparison of all groups###
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
uniq_paper_detags<-unlist(uniq_paper_detags)
uniq_paper_detags <-as.vector(uniq_paper_detags) #convert the list to a vector
uniq_paper_detags
typeof(uniq_paper_detags) #vector now!
length(uniq_paper_detags) #570 DF genes from paper

###Compare our unique genes to the paper's###
compare_to_paper<-uniq_paper_detags %in% unique_genes 
sum(compare_to_paper) #516 genes (when not normalized)in the 570 unique from the paper
#when normalized, 476 of our 481 differential genes match the paper's 570 differential genes

###Venn Diagram###
install.packages('VennDiagram')
library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 = num_FMale_hit, area2 = num_FSupM_hit, area3 = num_MSupMale_hit, n12 = num_FMale_FSupM, n23 = num_FSupM_MSupM, n13 = num_FMale_MSupMale, 
                 n123 = num_FSupM_MSupM_FMale, category = c("Male vs Female", "Supermale VS female", "Supermale VS male"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))

