# This directs to DESeq2 and Rsubread manual
# browseVignettes("DESeq2")
# RsubreadUsersGuide()

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsubread")
# BiocManager::install("DESeq2")
# BiocManager::install("ggplot2")
# BiocManager::install("GGally")
# BiocManager::install("pheatmap")

library(pheatmap)
library(Rsubread)
library(DESeq2)
library(ggplot2)
library(GGally)
library(psych)
#!!!input path to GTF annotation file
annotFile = "/Volumes/Mitch_Home_Backup/Mitch_Ellison_Data_All_UPitt_2015_to_2020/ChIP_seq_Analysis_Docker_started_on_11012020/DEseq_Analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.gff"

dataset = c("IgG")

features = c("non_overlapping_gene")



for(feature in features){
  
dir.create(paste0("/Volumes/Mitch_Home_Backup/Arndt_Lab_Consulting/ChIPseq_Redo_Spt6_AID_112022/DESeq_Analysis_ChIPseq_Redo_11192022/output/output_data_", feature, "_", dataset))
  
setwd(paste0("/Volumes/Mitch_Home_Backup/Arndt_Lab_Consulting/ChIPseq_Redo_Spt6_AID_112022/DESeq_Analysis_ChIPseq_Redo_11192022/output/output_data_", feature, "_", dataset))


# #!!!Ensure the you make your featurecount files is consistent with the order of the table you make








fCount1 <- featureCounts(files="/Volumes/Mitch_Home_Backup/Arndt_Lab_Consulting/ChIPseq_Redo_Spt6_AID_112022/DESeq_Analysis_ChIPseq_Redo_11192022/bams/SN126_S13_nodups.bam",
                         annot.ext=annotFile, isGTFAnnotationFile=TRUE, GTF.featureType=feature, GTF.attrType="gene_id")
fCount2 <- featureCounts(files="/Volumes/Mitch_Home_Backup/Arndt_Lab_Consulting/ChIPseq_Redo_Spt6_AID_112022/DESeq_Analysis_ChIPseq_Redo_11192022/bams/SN127_S14_nodups.bam",
                         annot.ext=annotFile, isGTFAnnotationFile=TRUE, GTF.featureType=feature, GTF.attrType="gene_id")
fCount3 <- featureCounts(files="/Volumes/Mitch_Home_Backup/Arndt_Lab_Consulting/ChIPseq_Redo_Spt6_AID_112022/DESeq_Analysis_ChIPseq_Redo_11192022/bams/SN128_S15_nodups.bam",
                         annot.ext=annotFile, isGTFAnnotationFile=TRUE, GTF.featureType=feature, GTF.attrType="gene_id")
fCount4 <- featureCounts(files="/Volumes/Mitch_Home_Backup/Arndt_Lab_Consulting/ChIPseq_Redo_Spt6_AID_112022/DESeq_Analysis_ChIPseq_Redo_11192022/bams/SN129_S16_nodups.bam",
                         annot.ext=annotFile, isGTFAnnotationFile=TRUE, GTF.featureType=feature, GTF.attrType="gene_id")


#Everything below is for formating featurecount ouput for Deseq2 input and further processing

#"fCount1$annotation[,c("GeneID")]" makes first row GeneID, the others are count inputs for each featureCounts file
#!!!add or remove fCountX$counts repeats as necessary
write.table(x=data.frame(fCount1$annotation[,c("GeneID")],fCount1$counts,fCount2$counts,fCount3$counts,fCount4$counts,
                         stringsAsFactors=FALSE),file="counts.txt",
                         quote=FALSE,sep="\t",row.names=FALSE)


conditions <- c("IgG_DMSO", "IgG_DMSO",
                "IgG_IAA", "IgG_IAA")

sampleNum = length(conditions)
#if your data is not single read indicate so here
coldata <- matrix(c(conditions,rep("paired_end", sampleNum)), nrow = sampleNum, ncol = 2)



# }
head(coldata)

#!!!assigns sample ID to each condition
rownames(coldata) <- (c("IgG_DMSO", "IgG_DMSO",
                        "IgG_IAA", "IgG_IAA"))
#adds header to coldata file
colnames(coldata) <- c("condition","type")
head(coldata)

#sets variable cts to the counts table generated above, automatically sets column names to row names 
#!!!MAKE SURE THE SAMPLE ORDER IN THE COUNT DATA FILE IS THE SAME AS THE COLUMN DATA FILE!!!
cts <- as.matrix(read.csv("counts.txt", header = TRUE, 
  sep="\t", colClasses = c("NULL", rep(NA,sampleNum)), col.names = rownames(coldata), row.names= fCount1$annotation[,c("GeneID")]))

head(cts)

#generate realtive read depth
#estimateSizeFactorsForMatrix(cts)

# dim(cts)
# keep <- rowSums(cts) >= 100
# cts <- cts[keep,]
# dim(cts)

# take only genes with a raw read count of greater than 20 meaning that at least 10 fragments were mapped to that gene
cts = cts[apply(cts[,1:4], 1, function(x) all(x>=10)),]


#Differential Expression Analysis
# library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

#dds <- DESeq(dds)
resultsNames(dds)
sizeFactors(dds) <- c(1/0.129471, 1/0.0896311, 1/0.126794, 1/0.0849086)

dds = estimateDispersions(dds) 
dds = nbinomWaldTest(dds)

#!!!input contrasting conditions here 
contrast=c("condition","IgG_IAA","IgG_DMSO")
IgG_IAA_vs_IgG_DMSO <- results(dds, contrast, pAdjustMethod = "fdr")

write.table(subset(IgG_IAA_vs_IgG_DMSO, padj < 0.05), file="IgG_IAA_vs_IgG_DMSO_significant.tsv", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


write.table(subset(IgG_IAA_vs_IgG_DMSO, log2FoldChange < -0.58), file="IgG_IAA_vs_IgG_DMSO_l2fc_less_than_-0.58.tsv", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)



#output is a results file comparing two conditions
png("IgG_IAA_vs_IgG_DMSO_MA_plot.png", height = 600, width = 600)
plotMA(IgG_IAA_vs_IgG_DMSO, ylim=c(-5,5), alpha = 0.05, xlab = "NCIS Normalized Read Counts")
dev.off()



#re-defines counts as a dataframe compatible with pairs and ggpairs functions
ddsdf <- as.data.frame(counts(dds, normalized=TRUE))

library(corrplot)
M = cor(ddsdf)

png("corrplot_numbers.png")
corrplot(M, method = "number")
dev.off()

png("corrplot_color.png")
corrplot(M, method = "color")
dev.off()

}

library(psych)

png("pairs.png", height = 800, width = 800)
pairs.panels(ddsdf, pch = 19)
dev.off()

png("log2_pairs.png", height = 800, width = 800)
pairs.panels(log2(ddsdf))
dev.off()

png("log2_pairs_custom.png", height = 800, width = 800)
pairs.panels(log2(ddsdf), xlim = c(0,15), ylim = c(0,15),
             smooth = FALSE, scale = FALSE, density = FALSE,
             ellipses = FALSE, method = "pearson", hist.col = "blue",
             breaks = 50)
dev.off()
