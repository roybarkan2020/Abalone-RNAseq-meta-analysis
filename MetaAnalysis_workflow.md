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
metable_sra_all <- read.csv(<metadata_table_file.csv>)
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

### QC - fastqc and multiqc (do this for all reads of each dataset) 

```
fastqc <fastq_file.fastq>
```

```
multiqc <path/to/the/fastqc/folder.zip> 
```

* If the data passes the QC you can move forward. If not, consider trimming and processing the "problematic data" with Trimmomatic or other available tools. 


### Build your reference for mapping (in this case, the red abalone transcriptome). Generating 'transcript-to-gene-map.txt' and 'mRNAonly.fna' files from the red abalone genomic files:

```
cat <gff_file.gff> | awk -F '\t' '$3=="mRNA"' | sed -E 's/.*Parent=(gene-LOC[^;]+).*transcript_id=([^;]+).*/\1 \2/' > <gene_to_transcript_map.txt>
```
```
cat gene_to_transcript_map.txt | awk '{print $2}' | xargs samtools faidx <rna_file.fna> > <mRNAonly_file.fna>
```

### Reference index (transcripts) for rsem-caclulate-expression:

```
rsem rsem-prepare-reference <mRNAonly_file.fna> --transcript-to-gene-map <gene_to_transcript_map.txt> --bowtie2 <perfix_output_reference_file>
```

### Map, align and get gene counts for each set of paired-end reads using rsem-calculate-expression:

```
rsem rsem-calculate-expression -p 8 --bowtie2 --no-bam-output --paired-end </path/for/forward/read.fastq> </path/for/reverse/read.fastq> </path/for/output_reference_file>  <perfix_rsem_output_file_name>
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

```{r}
list_files <- c("829_ruf.genes.results", "830_ruf.genes.results", "832_ruf.genes.results", "833_ruf.genes.results")
names(list_files) <- c("829_H", "830_C", "832_H", "833_C")
txData <- tximport(list_files, 'rsem')
txData$length[txData$length == 0] <- 1
colums <- read.table("coldata_237.txt", header = T, sep = '\t')
```

### Build a data structure using DESeq2 that defines how we are going to compare values

```{r}
dds_237 <- DESeqDataSetFromTximport(txData, columns, ~treatment)
```

### Remove lowly expressed genes

```{r}
keep_237 <- rowSums(counts(dds_237)) >= 10
dds_237 <- dds_237[keep_237,]
```

### DESeq2 (estimating: size factor, dispersions, model testing)

```{r}
ddsDE_237 <- DESeq(dds_237)
```

### Export normalized counts

```{r}
normCounts_237 <- counts(ddsDE_237, normalized = T)
write.csv(normCounts_237, "normCounts_237.csv")
```

### DESeq2 results

```{r}
DESeq_res_237 <- results(ddsDE_237, alpha = 0.05)
DESeq_resOrd_237 <- DESeq_res_237[order(DESeq_res_237$padj),]
write.csv(DESeq_resOrd_237, "DEseqResults_ordered.csv")
```

### Look at the percentage of DE genes, outliers, etc

```{r}
summary(DESeq_res_237)
```

### Further data manipulation and plotting

```{r}
plotMA(ddsDE_237)
normCountMatrix_237 <- read.csv("normCounts_237.csv", row.names = 1)
deseq_res_matrix_237 <- read.csv("DEseqResults_ordered.csv", row.names = 1)
deseq_res_matrix_237$sig <- ifelse(deseq_res_matrix_237$padj <= 0.05, "yes", "no")
deseq_res_matrix_237 <- na.omit(deseq_res_matrix_237)
ggplot(deseq_res_matrix_237, aes(x = log10(baseMean), y= log2FoldChange, color = sig)) + geom_point()
```

### Volcano plot

```{r}
ggplot(deseq_res_matrix_237, aes(x = log2FoldChange, y= -log10(padj), color = sig)) + geom_point()
```

### Filter all non-sig values

```{r}
sig_only_237 <- subset(deseq_res_matrix_237, padj <= 0.05)
```

### Merge both data frames (sig_only) 

```{r}
allSig_237 <- merge(normCountMatrix_237, sig_only_237, by = 0)
allSig_237_filter <- allSig_237[,2:5]
row.names(allSig_237_filter) <- allSig_237$Row.names
```

### Heatmap

```{r}
pheatmap(log2(allSig_237_filter + 1))
pheatmap(log2(allSig_237_filter + 1), scale = 'row')
pheatmap(log2(allSig_237_filter + 1), scale = 'row', show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T)
pheatmap(log2(allSig_237_filter + 1), scale = 'row', show_colnames = F, show_rownames = F, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
pheatmap(log2(allSig_237_filter + 1), scale = 'row', show_colnames = T, show_rownames = F, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
```

### More heatmaps (visualizing the most DEG)

```{r}
top_ten_deg_237 <- allSig_237_filter[1:10,]
top_fifty_deg_237 <- allSig_237_filter[1:50,]
top_hundred_deg_237 <- allSig_237_filter[1:100,]
pheatmap(log2(top_ten_deg_237 + 1))
pheatmap(log2(top_ten_deg_237 + 1), scale = 'row', show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
pheatmap(log2(top_fifty_deg_237 + 1), scale = 'row', show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
pheatmap(log2(top_hundred_deg_237 + 1), scale = 'row', show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, treeheight_row = 0, treeheight_col = 0)
```

### Looking at the top most expressed genes

```{r}
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

```{r}
plotCounts(dds_237, gene = which.max(DESeq_res_237$padj), intgroup = "treatment")
plotCounts(dds_237, gene = which.min(DESeq_res_237$padj), intgroup = "treatment")
plotCounts(dds_237, gene = "gene-LOC124120723", intgroup = "treatment")
#normalization for more ploting
vsd_237 <- vst(dds_237, blind = F)
```

### PCA

```{r}
plotPCA(vsd_237, intgroup = c("run" , "treatment"))
plotPCA(vsd_237, intgroup = c("treatment"))
```

### Sample distance

```{r}
sampleDist_237 <- dist(t(assay(vsd_237)))
sampleDist_237_matrix <- as.matrix(sampleDist_237)
rownames(sampleDist_237_matrix) <- vsd_237$treatment
colnames(sampleDist_237_matrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDist_237_matrix, clustering_distance_rows = sampleDist_237, clustering_distance_cols = sampleDist_237, col = colors)
```

## Step 3 - WGCNA (adopted from https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html)

### Load tools 

```{r}
library(DESeq2)
library(magrittr)
library(WGCNA)
library(ggplot2)
```

### Prep files for WGCNA 

Merge (by gene) all normalized gene counts from each study into one dataframe: 

```{r}
merged <- left_join(norm_counts1, norm_counts2, by = join_by(Gene == Gene))
```

Store the gene IDs as row names: 

```{r}
df <- merged %>% tibble::column_to_rownames("Gene")
```

Extract relevant data for the metadata:

```{r}
meta_df <- data.frame( Sample = names(df)) %>% mutate(Type = gsub("_.*","", Sample) %>% gsub("[.].*","", .))
```

Test that the metadata fits the counts table:

```{r}
all.equal(colnames(df), meta_df$Sample)
```

Only keep rows that have total counts above the cutoff (50 in the example below)

```{r}
df <- round(df) %>% as.data.frame() %>% dplyr::filter(rowSums(.) >= 50)
```

Group samples by treatments (in this case - CTRL for "control samples" and HEAT for "heat-stressed samples")

```{r}
metadata <- meta_df %>% dplyr::mutate(Type = dplyr::case_when(stringr::str_detect(Sample, "CTRL_") ~ "CTRL",stringr::str_detect(Sample, "HEAT_") ~ "HEAT"),Type = as.factor(Type))
```

Test it:

```{r}
levels(metadata$Type)
```

Create a `DESeqDataSet` object (now done for all datasets combined, not as in step 2, which was done for each dataset separately) 

```{r}
dds <- DESeqDataSetFromMatrix(countData = df,colData = metadata, design = ~1)
```

Normalize data:

```{r}
dds_norm <- vst(dds)
normalized_counts <- assay(dds_norm) %>% t()
```

Determine parameters for WGCNA:

The pickSoftThreshold() function to help identify good threshold choice for the downstream analysis:

```{r}
sft <- pickSoftThreshold(normalized_counts, dataIsExpr = TRUE, corFnc = cor, networkType = "signed")
```

```{r}
sft_df <- data.frame(sft$fitIndices) %>% dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)
```

Create a plot of the model fitting by the "power threshold" to get a better view of choosing an appropriate threshold for power (normally, you will aim for values above y=0.8):

```{r}
ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +  geom_point() + geom_text(nudge_y = 0.1) + geom_hline(yintercept = 0.80, col = "red") + ylim(c(min(sft_df$model_fit), 1.05)) + xlab("Soft Threshold (power)") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale independence") + theme_classic()
```

Find gene co-expression modules, using the selected threshold (from the plot above, in this case power=6): 

```{r}
bwnet <- blockwiseModules(normalized_counts,
  maxBlockSize = 5000, 
  TOMType = "signed", 
  power =6, 
  numericLabels = TRUE, 
  randomSeed = 1234,
)
```

Create main results to file:

```{r}
readr::write_rds(bwnet,"wgcna_results.RDS")
```

Extract the most relevant parts from bwnet (for plotting):

```{r}
module_eigengenes <- bwnet$MEs
```

Check that samples are in the right order:

```{r}
all.equal(metadata$Sample, rownames(module_eigengenes))
```

Create the designMatrix from the `Type` (i.e. HEAT/CTRL) variable

```{r}
des_mat <- model.matrix(~ metadata$Type)
```

Run linear model (limma):

```{r}
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
fit <- limma::eBayes(fit)
```

Apply multiple testing corrections, obtain stats and see which modules are significant (look at the adjusted P value):

```{r}
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>% tibble::rownames_to_column("module")
head(stats_df)
```

Set up the module eigengene for a chosen module (in this case, module 5) with the sample metadata labels:

```{r}
module_5_df <- module_eigengenes %>% tibble::rownames_to_column("Sample") %>% dplyr::inner_join(metadata %>% dplyr::select(Sample, Type), by = c("Sample" = "Sample"))
```

Plot:

```{r}
ggplot(module_5_df,aes(x = Type, y = ME5, color = Type)) + geom_boxplot(width = 0.2, outlier.shape = NA) +  ggforce::geom_sina(maxwidth = 0.3) + theme_classic()
```

Create a file with the genes associated with the chosen module(s):

```{r}
gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>% dplyr::mutate(module = paste0("ME", module))
gene_module_key %>% dplyr::filter(module == "ME5")
```

Create a function to produce heatmaps of modules of choice: 

```{r}
make_module_heatmap <- function(module_name, expression_mat = normalized_counts, metadata_df = metadata, gene_module_key_df = gene_module_key, module_eigengenes_df = module_eigengenes) {
module_eigengene <- module_eigengenes_df %>% dplyr::select(all_of(module_name)) %>% tibble::rownames_to_column("Sample")
+     col_annot_df <- metadata_df %>% dplyr::select(Sample, Type, Sample) %>% dplyr::inner_join(module_eigengene, by = "Sample") %>%  dplyr::arrange(Type, Sample) %>% tibble::column_to_rownames("Sample")
+     
col_annot <- ComplexHeatmap::HeatmapAnnotation(Treatment = col_annot_df$Type, module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)), col = list(Treatment = c("CTRL" = "#f1a340", "HEAT" = "#998ec3"))
)
+     module_genes <- gene_module_key_df %>% dplyr::filter(module == module_name) %>% dplyr::pull(gene)
+     mod_mat <- expression_mat %>% t() %>% as.data.frame() %>% dplyr::filter(rownames(.) %in% module_genes) %>% dplyr::select(rownames(col_annot_df)) %>% as.matrix()
+     mod_mat <- mod_mat %>% t() %>% scale() %>% t()
+     color_func <- circlize::colorRamp2( c(-2, 0, 2), c("#67a9cf", "#f7f7f7", "#ef8a62"))
+      heatmap <- ComplexHeatmap::Heatmap(mod_mat, name = module_name,col = color_func, bottom_annotation = col_annot, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE
)
return(heatmap)
}
```

Create a heatmap of your chosen module(s) - module 5, in this case, was the most significant:

```{r}
mod_5_heatmap <- make_module_heatmap(module_name = "ME5")
```

Look at the heatmap

```{r}
mod_5_heatmap
```


