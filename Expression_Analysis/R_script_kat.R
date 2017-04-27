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

#design matrix of sex and line
line<-factor(c(88,88,88,89,89,89,103,103,103,9,9))
sex<-factor(c("XX","XY","YY","XX","XY","YY","XX","XY","YY","XX","XY"))

design <- model.matrix(~line+sex, data = df)
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
heatmap.2(FPKM_matrix, dendrogram="both", trace="none", scale="row", density.info="none", col=my_palette,
          srtCol=45, cexCol = 1)
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
library(ggplot2)
#install.packages(c("reshape", "gplots"))
library(reshape)
library(gplots)

cd_matrix <- as.matrix(cd)
cd_df <- as.data.frame(cd_matrix)
genenames <- uniq_detags$gene_id

library(data.table)
setDT(cd_df, keep.rownames = TRUE)[]

#DEG_counts <- cd_df[cd_df$rn %in% genenames,]
#write.csv(DEG_counts, file = "DEG_counts.csv")

heatmap_counts <- read.csv("DEG_counts.csv", header = T, row.names = "rn")

###introduces column with row numbers from original counts data frame
heatmap_counts <- heatmap_counts[,2:12]
heatmap_counts_matrix <- as.matrix(heatmap_counts)

my_palette <- colorRampPalette(c("turquoise4","white","maroon2"))(n = 1000)

heatmap.2(heatmap_counts_matrix, dendrogram="both", trace="none", scale="row", density.info="none", col=my_palette, srtCol=45, cexCol = 1)
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
###From Devin (edited by Kat)
#install.packages('VennDiagram')
library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 = FM, area2 = FSupM, area3 = MSupM, 
                 n12 = FM_FSupM, n23 = FSupM_MSupM, n13 = FM_MSupM, 
                 n123 = FM_FSupM_MSupM, 
                 category = c("Male vs Female", "Supermale vs female", "Supermale vs male"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"),
                 cex = 2, cat.cex = 2,
                 cat.just=list(c(-0.2,4) , c(1,4) , c(0.5,-4)))


###Breakdown of genes that matched between our analysis and the paper's
compare_FM <- sum(detags_FM %in% paper_fm$gene_id)
compare_FSupM <- sum(detags_FSupM %in% paper_fsupm$gene_id)
compare_MSupM <- sum(detags_MSupM %in% paper_msupm$gene_id)
###---------------------------------------------------------------------------------------------------




###---------------------------------------------------------------------------------------------------
###Creating heatmap from samples that went through transcriptome pipeline

stringsAsFactors = F
#read in files with FPKM
samp1 <- read.delim("SRR1639681_counting.genes.results", stringsAsFactors = FALSE)
samp2 <- read.delim("SRR1642915_counting.genes.results", stringsAsFactors = FALSE)
 
#join fpkm columns together
#rename fpkm columns with sample name
colnames(samp1)[6] <- "FPKM_1"
colnames(samp2)[6] <- "FPKM_2"

joined_fpkm <- cbind(samp1$gene_id, samp1$FPKM_1, samp2$FPKM_2)
joined_fpkm <- as.data.frame(joined_fpkm)
fpkm_rownames <- data.frame(joined_fpkm[,-1], row.names=joined_fpkm[,1])
fpkm_colnames <- fpkm_rownames
colnames(fpkm_colnames)[1] <- "FPKM_1"
colnames(fpkm_colnames)[2] <- "FPKM_2"

###gene id names are different format in uniq_detags and fpkm data frame 
###change uniq_detags to match those in fpkm
uniq_detags_vector <- as.vector(uniq_detags$gene_id)
uniq_detags_aspof <- gsub("comp", "Aspof_comp", uniq_detags_vector)

###subset fpkm values by only those determined to be DEG
deg_fpkm <- subset(fpkm_colnames, rownames(fpkm_colnames) %in% uniq_detags_aspof)
deg_fpkm$FPKM_1 <- as.numeric(as.character(deg_fpkm$FPKM_1))
deg_fpkm$FPKM_2 <- as.numeric(as.character(deg_fpkm$FPKM_2))

###make into matrix and generate heatmap
#since FPKM values are strongly right-tailed with few very large numbers, a log2 transformation is done to better visualize the data.  1 is added the each value to prevent the 0s from being changed to infinity values
deg_fpkm_log <- log2(deg_fpkm[,] + 1)
#use summary to decide where you want to have your breaks.
summary(deg_fpkm_log)
#breaks sets how many colors in your gradient and what values get assigned what colors.  Max is 10.567
breaks=seq(0, 11, by=0.1)

deg_fpkm_matrix <- as.matrix(deg_fpkm_log)

my_col <- colorpanel(n=length(breaks)-1,low="gray75",high="black")

heatmap.2(deg_fpkm_matrix, breaks = breaks, dendrogram = "none", trace="none", density.info="none", col=my_col, srtCol=45, cexCol = 1)





