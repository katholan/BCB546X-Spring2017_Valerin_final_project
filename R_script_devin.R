setwd("C:/Users/Devin/Documents/GradSchool/IowaState/EEOB_546X/BCB546X-Spring2017_Valerin_final_project")
getwd()
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("edgeR")
library(edgeR)
counts_data<-read.delim("Asparagus.RSEM_genecounts.txt", check.names=FALSE, stringsAsFactors = FALSE) #read in gene counts
head(counts_data)
cd<-DGEList(counts=counts_data) #convert gene count data to DGEList
keep <- rowSums(cpm(cd)>1) >= 2 #filters out lowly expressed genes
cd <- cd[keep, , keep.lib.sizes=FALSE]
cd<-calcNormFactors(cd) #normalizes
plotMDS(cd)

#design matrix of sex and line
rm(design)
line<-factor(c(88,88,88,89,89,89,103,103,103,9,9))
sex<-factor(c("XX","XY","YY","XX","XY","YY","XX","XY","YY","XX","XY"))
data.frame(Sample=colnames(cd),line,sex)
design <- model.matrix(~line+sex)
rownames(design) <- colnames(cd)
design

#design matrix by by sex
rm(design)
sex<-factor(c("XX","XY","YY","XX","XY","YY","XX","XY","YY","XX","XY"))
data.frame(Sample=colnames(cd),sex)
design<-model.matrix(~sex)
rownames(design) <- colnames(cd)
design

#estimates dispersion
cd <- estimateDisp(cd, design)
#gits using glm QL fit (recommended glm)
fit <- glmQLFit(cd, design)
plotQLDisp(fit)

qlf<-glmQLFTest(fit,coef=2:3)
plotMD(qlf)

tr<-glmTreat(fit) 
plotMDS()
logcpm <- cpm(y, prior.count=2, log=TRUE) #for heat map

cd
et
