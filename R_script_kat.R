###EEOB 546X Final Project
###Group: Valerin
###Members: Valeria, Katerina, and Devin
###Paper: Sex-based gene expression in dioecious garden asparagus (Asparagus officinalis), Harkess et. al, New Phytologist, 2015

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
library("dplyr")
###----------------------------------------------------------------------------------------------------------------
###Import csv of total gene counts

genecounts <- read.csv("genecounts.csv", stringsAsFactors = F, row.names = "gene_ID")
#dim(genecounts)
###120061 rows (aka genes/transcripts)
###12 columns (first column gene id, other columns line and sex libraries)
#class(genecounts$X88F)

###----------------------------------------------------------------------------------------------------------------
###From Devin
counts_data<-read.delim("Asparagus.RSEM_genecounts.txt", check.names=FALSE, stringsAsFactors = FALSE)

cd<-DGEList(counts=counts_data) 
keep <- rowSums(cpm(cd)>1) >= 3 
cd <- cd[keep, , keep.lib.sizes=FALSE]
cd<-calcNormFactors(cd) 
colors <- c("black","black","black","darkgreen","darkgreen","darkgreen","blue","blue","blue", "red", "red")
plotMDS(cd, col=colors, top = 250, gene.selection = "pairwise", main="MDS plot ")
names<-c("female", "male","supermale","female", "male","supermale","female", "male")
plotMDS(cd, col=colors, labels=names, top = 250, gene.selection = "pairwise")

#design matrix of sex and line
line<-factor(c(88,88,88,89,89,89,103,103,103,9,9))
sex<-factor(c("XX","XY","YY","XX","XY","YY","XX","XY","YY","XX","XY"))
design <- model.matrix(~line+sex)
rownames(design) <- colnames(cd)

#estimates dispersion
cd <- estimateGLMCommonDisp(cd, design)
cd <- estimateGLMTrendedDisp(cd, design)
cd <- estimateGLMTagwiseDisp(cd, design)
fit <- glmFit(cd, design) 

###Comparisons

#Female vs supermale
lrt_FSupM <- glmLRT(fit, coef=6)
summary(de_FemvsSupMal <- decideTestsDGE(lrt_FSupM, adjust.method="fdr"))
FDR_FSupM <- p.adjust(lrt_FSupM$table$PValue, method="fdr") 
sum(FDR_FSupM < 0.05) 

#Female vs male
lrt_FM <- glmLRT(fit, coef=5) #female VS male
summary(de_FemvsMal <- decideTestsDGE(lrt_FM, adjust.method="fdr"))
FDR_FMale<-p.adjust(lrt_FM$table$PValue,method="fdr") #adjust Pvalue
sum(FDR_FMale<0.05)

#Male vs supermale
lrt_MSupM <- glmLRT(fit, contrast=c(0,0,0,0,-1,1)) #Male VS Supermale
summary(de_MalvsSupMal <- decideTestsDGE(lrt_MSupM, adjust.method="fdr"))
FDR_MSupMale<-p.adjust(lrt_MSupM$table$PValue,method="fdr")
sum(FDR_MSupMale<0.05)
###----------------------------------------------------------------------------------------------------------------



###----------------------------------------------------------------------------------------------------------------
###Heatmap
###workflow from https://github.com/alexharkess/PBIO_6550/wiki/Lab-8-:-Data-munging-and-visualization

library(ggplot2)
#install.packages(c("reshape", "gplots"))
library(reshape)
library(gplots)

FPKM_matrix <- read.table("TMM_normalized_FPKM_matrix.txt", header = T, row.names = 1)
FPKM_matrix_melt <- melt(FPKM_matrix)
plot <- ggplot(FPKM_matrix_melt, aes(x = variable, y = value)) + geom_boxplot() +
  geom_jitter()
plot
plot + ylim(0,100)

FPKM_matrix <- as.matrix(FPKM_matrix)

my_palette <- colorRampPalette(c("turquoise4","white","maroon2"))(n = 1000)
heatmap.2(FPKM_matrix, dendrogram="both", trace="none", scale="row", density.info="none", col=my_palette)
###tutorial heatmap of 570 genes, only 9 lines (doesn't include 8A male or 10 supermale)
###---------------------------------------------------------------------------------------------------


###---------------------------------------------------------------------------------------------------
###female vs supermale
detags_FSupM <- rownames(cd)[as.logical(de_FemvsSupMal)]

plotMD(lrt_FSupM, main="Female VS SuperMale", bg.col = "gray", hl.col = c("blue", "red"), status = (row.names(lrt_FSupM$table) %in% detags_FSupM), legend = "topright")
abline(h=c(-1, 1), col = "red")

###female vs male
detags_FM <- rownames(cd)[as.logical(de_FemvsMal)]

plotMD(lrt_FM, main="Female VS Male", bg.col = "gray", hl.col = c("blue", "red"), status = (row.names(lrt_FSupM$table) %in% detags_FM), legend = "topright")
abline(h=c(-1, 1), col = "red")

###male vs supermale
detags_MSupM <- rownames(cd)[as.logical(de_MalvsSupMal)]

plotMD(lrt_MSupM, main="Male VS Supermale", bg.col = "gray", hl.col = c("blue", "red"), status = (row.names(lrt_FSupM$table) %in% detags_MSupM), legend = "topright")
abline(h=c(-1, 1), col = "red")

###---------------------------------------------------------------------------------------------------


###---------------------------------------------------------------------------------------------------
###total number of unique DEGs in our analysis
detags <- c(detags_FM, detags_FSupM, detags_MSupM)
detags_df <- as.data.frame(detags)
uniq_detags <- unique(detags_df)
colnames(uniq_detags) <- "gene_id"
###---------------------------------------------------------------------------------------------------


###---------------------------------------------------------------------------------------------------
###heatmap
###pulling out only genes we found to be DE from the counts data
cd_matrix <- as.matrix(cd)
cd_df <- as.data.frame(cd_matrix)
genenames <- uniq_detags$gene_id

library(data.table)
setDT(cd_df, keep.rownames = TRUE)[]

DEG_counts <- cd_df[cd_df$rn %in% genenames,]
write.csv(DEG_counts, file = "DEG_counts.csv")

heatmap_counts <- read.csv("DEG_counts.csv", header = T, row.names = "rn")

###introduces column with row numbers from original counts data frame
heatmap_counts <- heatmap_counts[,2:12]
heatmap_counts_matrix <- as.matrix(heatmap_counts)

my_palette <- colorRampPalette(c("turquoise4","white","maroon2"))(n = 1000)

heatmap.2(heatmap_counts_matrix, dendrogram="both", trace="none", scale="row", density.info="none", col=my_palette)
###---------------------------------------------------------------------------------------------------


###---------------------------------------------------------------------------------------------------
###DEGs from paper
paper_fsupm <- read.delim("paper_fsupm.txt")
paper_fm <- read.delim("paper_fm.txt")
paper_msupm <- read.delim("paper_msupm.txt")

paper_detags <- rbind(paper_fsupm, paper_fm, paper_msupm)
uniq_paper_detags <- unique(paper_detags)

in_both <- sum(sort(uniq_detags[,]) %in% sort(uniq_paper_detags[,]))
###475 common between our analysis and theirs
###---------------------------------------------------------------------------------------------------


###---------------------------------------------------------------------------------------------------
###Venn Diagrams

###Venn Diagram of our genes
###from devin (edited by Kat)
FM <- length(detags_FM)
FSupM <- length(detags_FSupM)
MSupM <- length(detags_MSupM)

FM_FSupM <- sum(detags_FM %in% detags_FSupM)
FM_MSupM <- sum(detags_FM %in% detags_MSupM)
FSupM_MSupM <- sum(detags_FSupM %in% detags_MSupM)

FM_FSupM_MSupM <- sum(detags_FM %in% detags_FSupM %in% detags_MSupM)

###Venn Diagram###
###From devin
#install.packages('VennDiagram')
library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 = FM, area2 = FSupM, area3 = MSupM, 
                 n12 = FM_FSupM, n23 = FSupM_MSupM, n13 = FM_MSupM, 
                 n123 = FM_FSupM_MSupM, 
                 category = c("Male vs Female", "Supermale VS female", "Supermale VS male"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"),
                 cex = 2, cat.cex = 2)


###Venn Diagram of genes that matched between our analysis and the paper's
###Kat
paper_fsupm <- read.delim("paper_fsupm.txt")
paper_fm <- read.delim("paper_fm.txt")
paper_msupm <- read.delim("paper_msupm.txt")
paper_fsupm <- as.vector(paper_fsupm)
paper_fm <- as.vector(paper_fm)
paper_msupm <- as.vector(paper_msupm)

#p_fm_fsupm <- paper_fm %in% paper_fsupm
#p_fm_msupm <- paper_fm %in% paper_msupm
#p_fsupm_msupm <- paper_fsupm %in% paper_msupm
#p_all <- paper_fm %in% paper_fsupm %in% paper_msupm

compare_FM <- sum(detags_FM %in% paper_fm)
compare_FSupM <- sum(detags_FSupM %in% paper_fsupm)
compare_MSupM <- sum(detags_MSupM %in% paper_msupm)

#FM_FSupM <- sum(detags_FM %in% detags_FSupM %in% p_fm_fsupm)
#FM_MSupM <- sum(detags_FM %in% detags_MSupM %in% p_fm_msupm)
#FSupM_MSupM <- sum(detags_FSupM %in% detags_MSupM %in% p_fsupm_msupm)

compare_all <- sum(detags_FM %in% detags_FSupM %in% detags_MSupM %in% p_all)


grid.newpage()
draw.pairwise.venn(area1 = FM, area2 = paper_fm, cross.area = compare_FM,
                 category = c("Ours", "Paper"), lty = "blank", 
                 fill = c("skyblue", "pink1"),
                 cex = 2, cat.cex = 2)


###---------------------------------------------------------------------------------------------------
