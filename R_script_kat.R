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
###filtered genecounts to only include transcripts >1 CPM (counts per million) in at least 3 libraries
###paper gets 19,986 genes with these parameters

###round integers to whole numbers
genecounts_round <- round(genecounts)

#genecounts_round$total <- rowSums(genecounts_round[,])
#filtered_genecounts <- genecounts_round[genecounts_round$total >= 3,]

###results in 97630 total genes, but this includes genes that might have 2 or 3 CPM in a single library
filtered_genecounts <- genecounts_round %>% 
  filter((rowSums(genecounts_round[,,drop = F] > 1)) >= 3) 

###results in 56482, which is more than 2x the number in the paper.
###this code keeps dropping the rownames, even when drop = F


###alternative code to filter the genecounts dataframe.  does the same as the above but doesn't require the dplyr package
###doesn't drop row names
filtered_genecounts <- genecounts_round[rowSums(genecounts_round[,,drop = F] > 1) >= 3,]

###possible explanation for discrepency of gene numbers could be because I don't know what the units are in the file we received
###the paper dropped any less than 1 CPM, but if the RSEM_genecounts file is NOT in CPM, that could be where the problem lies

write.csv(filtered_genecounts, file = "filtered_genecounts")

###----------------------------------------------------------------------------------------------------------------
###From Devin
counts <- read.delim("Asparagus.RSEM_genecounts.txt", check.names = F, stringsAsFactors = F)
cd <- DGEList(counts = counts)
filtered <- rowSums(cpm(cd)>1) >=3
filtered_cd <- cd[filtered, , keep.lib.sizes = F]
nrow(filtered_cd)
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

breaks=seq(-3, 3, by=0.5)
my_palette <- colorpanel(n=length(breaks)-1,low="turquoise4",mid="white",high="maroon2")
heatmap.2(FPKM_matrix, dendrogram="both", trace="none", scale="row", density.info="none", col=my_palette)
###tutorial heatmap of 570 genes, only 9 lines (doesn't include 8A male or 10 supermale)


###adapted workflow for 570 DEGs
heatmap.2(plot.data,density.info="none",col=mycol, key = F, trace='none', breaks=breaks, margins = c(1,10), dendrogram="none",
          Rowv=NA, Colv=NA, cexRow=1, cexCol=1,
          lmat=rbind(c(0,3), c(2,1), c(0,4)), lwid=c(0.5,4), lhei=c(0.5, 4,1))



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
###compare gene lists from comparisons to paper gene lists (detags)
###venn diagram of de genes
###heatmap of total DEGs, as well as smaller heatmaps for comparisons??
lrt_FSupM <- glmLRT(fit, coef=4)
summary(de_FemvsSupMal <- decideTestsDGE(lrt_FSupM, adjust.method="fdr"))
FDR_FSupM <- p.adjust(lrt_FSupM$table$PValue, method="fdr")
sum(FDR_FSupM < 0.05)
detags_FSupM <- rownames(cd)[as.logical(de_FemvsSupMal)]
plotSmear(lrt_FSupM, de_FemvsSupMal.tags=detags)
abline(h=c(-1, 1))
sum(FDR_FSupM < 0.05)

lrt_FM <- glmLRT(fit, coef=3)
summary(de_FemvsMal <- decideTestsDGE(lrt_FM, adjust.method="fdr"))
detags_FM <- rownames(cd)[as.logical(de_FemvsMal)]
FDR_FM <- p.adjust(lrt_FM$table$PValue, method="fdr")
sum(FDR_FM < 0.05)

lrt_MSupM <- glmLRT(fit, contrast=c(0,0,-1,1))
summary(de_MalvsSupMal <- decideTestsDGE(lrt_MSupM, adjust.method="fdr"))
detags_MSupM <- rownames(cd)[as.logical(de_MalvsSupMal)]
FDR_MSupM <- p.adjust(lrt_MSupM$table$PValue, method="fdr")
sum(FDR_MSupM < 0.05)

###total number of unique DEGs in our analysis
detags <- c(detags_FM, detags_FsupM, detags_MsupM)
detags_df <- as.data.frame(detags)
uniq_detags <- unique(detags_df)
colnames(uniq_detags) <- "gene_id"

###double checking total number of unique DEGs from paper
paper_fsupm <- read.delim("paper_fsupm.txt")
paper_fm <- read.delim("paper_fm.txt")
paper_msupm <- read.delim("paper_msupm.txt")

paper_detags <- rbind(paper_fsupm, paper_fm, paper_msupm)
uniq_paper_detags <- unique(paper_detags)


###comparing paper DEGs to our DEGs
total_detags <- rbind(uniq_detags, uniq_paper_detags) #931
same_detags <- total_detags != unique(total_detags) #615
