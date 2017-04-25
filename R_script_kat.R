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



###----------------------------------------------------------------------------------------------------------------
###From Devin/Valeria
counts_data <- read.delim("Asparagus.RSEM_genecounts.txt", check.names = F, stringsAsFactors = F)
cd<-DGEList(counts=counts_data)
keep <- rowSums(cpm(cd)>1) >= 3
cd <- cd[keep, , keep.lib.sizes=FALSE]
group_line<-factor(c(88,88,88,89,89,89,103,103,103,9,9))
group_sex <- factor(c("Female", "Male",	"Super Male",	"Female",	"Male",	"Super Male",	"Female",	"Male",	"Super Male",	"Female",	"Male"))
design <- model.matrix(~group_line + group_sex)
rownames(design) <- colnames(cd)
cd <- estimateGLMCommonDisp(cd, design)
cd <- estimateGLMTrendedDisp(cd, design)
cd <- estimateGLMTagwiseDisp(cd, design)
fit <- glmFit(cd, design)

###19980 genes.  almost the same
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

my_palette <- colorRampPalette(c("green","black","red"))(n = 1000)
heatmap.2(FPKM_matrix, dendrogram="both", trace="none", scale="row", density.info="none", col=my_palette)
###tutorial heatmap of 570 genes, only 9 lines (doesn't include 8A male or 10 supermale)
###---------------------------------------------------------------------------------------------------


###---------------------------------------------------------------------------------------------------
###female vs supermale
lrt_FSupM <- glmLRT(fit, coef=6)
summary(de_FemvsSupMal <- decideTestsDGE(lrt_FSupM, adjust.method="fdr"))
FDR_FSupM <- p.adjust(lrt_FSupM$table$PValue, method="fdr")

detags_FSupM <- rownames(cd)[as.logical(de_FemvsSupMal)]
plotSmear(lrt_FSupM, de_FemvsSupMal.tags=detags)
abline(h=c(-1, 1))

sum(FDR_FSupM < 0.05)


###female vs male
lrt_FM <- glmLRT(fit, coef=5)
summary(de_FemvsMal <- decideTestsDGE(lrt_FM, adjust.method="fdr"))
FDR_FM <- p.adjust(lrt_FM$table$PValue, method="fdr")

detags_FM <- rownames(cd)[as.logical(de_FemvsMal)]
plotSmear(lrt_FM, de_FemvsMal.tags=detags_FM)
abline(h=c(-1, 1))

sum(FDR_FM < 0.05)


###male vs supermale
lrt_MSupM <- glmLRT(fit, contrast=c(0,0,0,0,-1,1))
summary(de_MalvsSupMal <- decideTestsDGE(lrt_MSupM, adjust.method="fdr"))
FDR_MSupM <- p.adjust(lrt_MSupM$table$PValue, method="fdr")
sum(FDR_MSupM < 0.05)

detags_MSupM <- rownames(cd)[as.logical(de_MalvsSupMal)]
plotSmear(lrt_FM, de_MalvsSupMal.tags=detags_MSupM)
abline(h=c(-1, 1))

sum(FDR_MSupM < 0.05)
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
counts2 <- counts_data
counts2 <- read.csv("genecounts.csv", stringsAsFactors = F)
genenames <- uniq_detags$gene_id
DEG_counts <- counts2[counts2$gene_ID %in% genenames,]
write.csv(DEG_counts, file = "DEG_counts.csv")

heatmap_counts <- read.csv("DEG_counts.csv", header = T, row.names = "gene_ID")

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

sort(uniq_detags[,1])
sort(uniq_paper_detags[,1])
in_both <- sum(sort(uniq_detags[,]) %in% sort(uniq_paper_detags[,]))
###516 common between our analysis and theirs
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



###---------------------------------------------------------------------------------------------------
###from devin
library(edgeR)
counts_data<-read.delim("Asparagus.RSEM_genecounts.txt", check.names=FALSE, stringsAsFactors = FALSE) #read in gene counts
cd<-DGEList(counts=counts_data) #convert gene count data to DGEList
keep <- rowSums(cpm(cd)>1) >= 3 #filters out lowly expressed genes
cd <- cd[keep, , keep.lib.sizes=FALSE]
cd<-calcNormFactors(cd) #normalizes

par(xpd = T)
colors <- c("black","black","black","darkgreen","darkgreen","darkgreen","blue","blue","blue", "red", "red")
names<-c("female", "male","supermale","female", "male","supermale","female", "male")
plotMDS(cd, col=colors, labels=names, top = 250, gene.selection = "pairwise", main = "Multidimensional Scaling Plot")
legend(3.5, 2.57, legend=c("Line 8A", "Line 8B","Line 10", "Line 9"), col= c("black","darkgreen","blue", "red"), title="Key", ncol = 1, pch = 19 )


###make plots pretty/colored/etc
###compare gene lists from comparisons to paper gene lists (detags), extract gene names that are different
###double checking total number of unique DEGs from paper


