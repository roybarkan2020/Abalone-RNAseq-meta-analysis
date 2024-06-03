# The transcriptomic response to heat-stress in abalone - A meta-analysis 

You can follow the steps we used to perform the meta-analysis here. This study aimed to test the "core" response to heat stress across multiple abalone species. 
The general steps are:
1. A collection of appropriate RNA-seq datasets from multiple studies was downloaded.
2. Each dataset was mapped to the same reference transcriptome (red abalone GCF_023055435 accession), and gene counts were used for differential gene expression analysis. 
3. Normalized gene counts were then merged and used for weighted gene co-expression network analysis (WGCNA).

## Step 1

### Import meta_table of all subsets:

This can be downloaded from NCBI BioProject(s) of interest. 

```
metable_sra_all <- read.csv("metadata_table_file.csv")
```

### Downloading reads from different studies (wget) and converting to fastq files (fastq-dump) :

This code downloads the SRA data:

```
wget --input [SRA_list_from_metable]
```

This code converts SRA to fastq:

```
fastq-dump --split-files [downloaded_SRA_files]
```

### QC - fastqc and multiqc (do this for each dataset) 

```
fastqc "fastq_file.fastq"
```

```
multiqc "path/to/the/fastqc/folder.zip" 
```

* If the data passes the QC you can move forward. If not, consider trimming and processing the "problematic data" with Trimmomatic or other available tools. 


### Build your reference for mapping (in this case, the red abalone transcriptome). Generating 'transcript-to-gene-map.txt' and 'mRNAonly.fna' files from the red abalone genomic files:

```
cat GCF_023055435.1_xgHalRufe1.0.p_genomic.gff | awk -F '\t' '$3=="mRNA"' | sed -E 's/.*Parent=(gene-LOC[^;]+).*transcript_id=([^;]+).*/\1 \2/' > gene_to_transcript_map.txt
```
```
cat gene_to_transcript_map.txt | awk '{print $2}' | xargs samtools faidx GCF_023055435.1_xgHalRufe1.0.p_rna.fna > GCF_023055435.1_xgHalRufe1.0.p_mRNAonly.fna
```

### Reference index (transcripts) for rsem-caclulate-expression:

```
singularity run ../rsem.sif rsem-prepare-reference GCF_023055435.1_xgHalRufe1.0.p_mRNAonly.fna --transcript-to-gene-map gene_to_transcript_map.txt --bowtie2 Rufref
```

### Map and align one set of paired-end reads using rsem-calculate-expression:

```
singularity run ../rsem.sif rsem-calculate-expression -p 8 --bowtie2 --no-bam-output --paired-end "/path/for/forward/read.fastq" "/path/for/reverse/read.fastq" /home/jc748673/test/Rufref  <output file name> 
```

## Step 2 (Moving to R) 

This analysis was performed on study PRJNA597237 which includes two control samples (SRR_830 and SRR_833) and two heat-stressed samples (SRR_829 and SRR_832)


### Loading packages/tools

```{r}

library(DESeq2)
library(tximport)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

```


### Import rsem count files

```
list_files <- c("829_ruf.genes.results", "830_ruf.genes.results", "832_ruf.genes.results", "833_ruf.genes.results")
names(list_files) <- c("829_H", "830_C", "832_H", "833_C")
txData <- tximport(list_files, 'rsem')
txData$length[txData$length == 0] <- 1
colums <- read.table("coldata_237.txt", header = T, sep = '\t')
```

### Build a data structure using DESeq2 that defines how we are going to compare values

```
dds_237 <- DESeqDataSetFromTximport(txData, columns, ~treatment)
```

### Remove lowly expressed genes

```
keep_237 <- rowSums(counts(dds_237)) >= 10
dds_237 <- dds_237[keep_237,]
```

### DESeq2 (estimating: size factor, dispersions, model testing)

```
ddsDE_237 <- DESeq(dds_237)
```
### Export normalized counts

```
normCounts_237 <- counts(ddsDE_237, normalized = T)
write.csv(normCounts_237, "normCounts_237.csv")
```

### DESeq2 results

```
DESeq_res_237 <- results(ddsDE_237, alpha = 0.05)
DESeq_resOrd_237 <- DESeq_res_237[order(DESeq_res_237$padj),]
write.csv(DESeq_resOrd_237, "DEseqResults_ordered.csv")
```
### Look at the percentage of DE genes, outliers, etc

```
summary(DESeq_res_237)
```

### Further data manipulation and plotting

```
plotMA(ddsDE_237)
normCountMatrix_237 <- read.csv("normCounts_237.csv", row.names = 1)
deseq_res_matrix_237 <- read.csv("DEseqResults_ordered.csv", row.names = 1)
deseq_res_matrix_237$sig <- ifelse(deseq_res_matrix_237$padj <= 0.05, "yes", "no")
deseq_res_matrix_237 <- na.omit(deseq_res_matrix_237)
ggplot(deseq_res_matrix_237, aes(x = log10(baseMean), y= log2FoldChange, color = sig)) + geom_point()
```
### Volcano plot

```
ggplot(deseq_res_matrix_237, aes(x = log2FoldChange, y= -log10(padj), color = sig)) + geom_point()
```

### Filter all non-sig values

```
sig_only_237 <- subset(deseq_res_matrix_237, padj <= 0.05)
```

### Merge both data frames (sig_only) 

```
allSig_237 <- merge(normCountMatrix_237, sig_only_237, by = 0)
allSig_237_filter <- allSig_237[,2:5]
row.names(allSig_237_filter) <- allSig_237$Row.names
```
### Heatmap

```
pheatmap(log2(allSig_237_filter + 1))
pheatmap(log2(allSig_237_filter + 1), scale = 'row')
pheatmap(log2(allSig_237_filter + 1), scale = 'row', show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T)
pheatmap(log2(allSig_237_filter + 1), scale = 'row', show_colnames = F, show_rownames = F, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
pheatmap(log2(allSig_237_filter + 1), scale = 'row', show_colnames = T, show_rownames = F, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
```

### More heatmaps (visualizing the most DEG)

```
top_ten_deg_237 <- allSig_237_filter[1:10,]
top_fifty_deg_237 <- allSig_237_filter[1:50,]
top_hundred_deg_237 <- allSig_237_filter[1:100,]
pheatmap(log2(top_ten_deg_237 + 1))
pheatmap(log2(top_ten_deg_237 + 1), scale = 'row', show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
pheatmap(log2(top_fifty_deg_237 + 1), scale = 'row', show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
pheatmap(log2(top_hundred_deg_237 + 1), scale = 'row', show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
```

### Looking/compating top expressed genes

```
top_ten_melt_237 <- melt(as.matrix(top_ten_deg_237))
names(top_ten_melt_237) <- c("gene","treatment","val")
top_ten_melt_237$treatment <- ifelse(grepl("heat", top_ten_melt_237$treatment), "heat", "control")
geneExp_top_ten <- ggplot(top_ten_melt_237, aes(x = treatment, y = log2(val + 1), color = treatment)) + geom_point() + facet_grid(~gene)
top_fifty_melt_237 <- melt(as.matrix(top_fifty_deg_237))
names(top_fifty_melt_237) <- c("gene","treatment","val")
top_fifty_melt_237$treatment <- ifelse(grepl("heat", top_fifty_melt_237$treatment), "heat", "control")
geneExp_top_fifty <- ggplot(top_fifty_melt_237, aes(x = treatment, y = log2(val + 1), color = treatment)) + geom_point() + facet_grid(~gene)
```

### Plot counts - can be done per gene by gene name/row name or for the highest/lowest expressed gene

```
plotCounts(dds_237, gene = which.max(DESeq_res_237$padj), intgroup = "treatment")
plotCounts(dds_237, gene = which.min(DESeq_res_237$padj), intgroup = "treatment")
plotCounts(dds_237, gene = "gene-LOC124120723", intgroup = "treatment")
#normalization for more ploting
vsd_237 <- vst(dds_237, blind = F)
```

### PCA

```
plotPCA(vsd_237, intgroup = c("run" , "treatment"))
plotPCA(vsd_237, intgroup = c("treatment"))
```

### Sample distance

```
sampleDist_237 <- dist(t(assay(vsd_237)))
sampleDist_237_matrix <- as.matrix(sampleDist_237)
rownames(sampleDist_237_matrix) <- vsd_237$treatment
colnames(sampleDist_237_matrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDist_237_matrix, clustering_distance_rows = sampleDist_237, clustering_distance_cols = sampleDist_237, col = colors)
```


```

