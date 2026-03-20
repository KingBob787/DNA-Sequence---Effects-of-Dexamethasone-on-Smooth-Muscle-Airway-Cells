# ============================================
# RNA-seq Analysis of Dexamethasone Treatment
# Dataset: Airway smooth muscle cells satvik desaiauauuuua
# ============================================

# --- 1. Load Libraries ---
# These are all the packages we need for the analysis
library(DESeq2)          # Main differential expression package
library(airway)          # Example dataset
library(ggplot2)         # General plotting
library(pheatmap)        # Heatmap visualization
library(EnhancedVolcano) # Volcano plot visualization
library(org.Hs.eg.db)   # Human gene annotation database

# --- 2. Load Data ---
# Load the airway dataset - an RNA-seq experiment testing
# the effect of dexamethasone on airway smooth muscle cells
data("airway")

# --- 3. Set Up DESeq2 Object ---
# Create a DESeq2 dataset
# design = ~ cell + dex means:
# - compare treated vs untreated (dex)
# - while accounting for different cell lines (cell)
dds <- DESeqDataSet(airway, design = ~ cell + dex)

# --- 4. Filter Low Count Genes ---
# Remove genes with 10 or fewer total counts across all samples
# These genes are too lowly expressed to give reliable results
# Reduces dataset from 63,677 genes to ~22,008 genes
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]

# --- 5. Run DESeq2 Analysis ---
# This does three things automatically:
# 1. Normalizes counts so samples are fairly comparable
# 2. Estimates variability for each gene
# 3. Tests each gene for significant changes
dds <- DESeq(dds)

# --- 6. Extract Results ---
# Get results comparing treated vs untreated
# contrast = c("dex", "trt", "untrt") means:
# look at the dex column, compare trt against untrt
results <- results(dds, contrast = c("dex", "trt", "untrt"))

# Quick summary of results
summary(results)

# How many genes are significant at padj < 0.05?
sum(results$padj < 0.05, na.rm = TRUE)

# --- 7. Sort Results by Significance ---
# Order genes from most to least significant
resultsOrdered <- results[order(results$padj),]

# --- 8. Add Human Readable Gene Names ---
# Convert ENSG codes to gene symbols using human genome database
resultsOrdered$symbol <- mapIds(org.Hs.eg.db,
                                keys = rownames(resultsOrdered),
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "first")

# View top 10 most significant genes
head(resultsOrdered, 10)

# --- 9. Volcano Plot ---
# Shows all genes plotted by fold change (x) vs significance (y)
# Red dots = significant AND large fold change (our hits)
# Top corners = most interesting genes
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Dexamethasone vs Untreated',
                pCutoff = 0.05,
                FCcutoff = 1.0)

# --- 10. Heatmap ---
# Shows expression patterns of top 20 most significant genes
# Red = higher than average expression
# Blue = lower than average expression
# Samples should cluster into treated vs untreated groups
top20 <- head(order(results$padj, na.last = NA), 20)
mat <- counts(dds, normalized = TRUE)[top20,]
pheatmap(mat,
         scale = "row",
         color = colorRampPalette(c("darkblue", "white", "darkred"))(100))

#11. Using datacamp PCA tutorial.
#install corrr package, used for correlation analysis
install.packages("corrr")
library(corrr)

install.packages("ggcorrplot")
library(ggcorrplot)

install.packages("FactoMineR")
library(FactoMineR)

install.packages("factoextra")
library(factoextra)

colSums(is.na(resultsOrdered))

resultsOrdered

resultsTop10 <- head(resultsOrdered, 10)

colSums(is.na(resultsTop10))

View(as.data.frame(resultsOrdered)) # visualize compressed data

counts_table <- as.data.frame(counts(dds, normalized = TRUE)) #visualize raw data of counts
View(counts_table)

#see the column data for cells and whether they were treated or not
colData(dds)[, c("cell", "dex")]

#select top 500 variable genes. 
top_var_genes <- head(order(rowVars(counts(dds, normalized = TRUE)), decreasing = TRUE), 500)

mat_pca <- t(counts(dds, normalized = TRUE)[top_var_genes,])

mat_scaled <- scale(mat_pca)

data.pca <- prcomp(mat_scaled)
summary(data.pca)

library(factoextra)
fviz_pca_ind(data.pca,
             col.ind = colData(dds)$dex,
             palette = c("blue", "red"),
             addEllipses = TRUE,
             legend.title = "Dexamethasone Treatment")

fviz_pca_ind(data.pca,
             col.ind = colData(dds)$cell,
             palette = c("orange", "purple", "green", "red"),
             addEllipses = FALSE,
             legend.title = "Cell Line")

colData(dds)[, c("cell", "dex")]


#12. Enrichment analysis
BiocManager::install("clusterProfiler")
library("clusterProfiler")

sig_genes <- subset(resultsOrdered, padj < 0.05) # sorts all the statistically significant genes

sig_genes$entrez <- mapIds(org.Hs.eg.db,
                           keys = rownames(sig_genes),
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first")
#^ extract the ENTREZ IDs of each significant gene

sig_genes

go_results <- enrichGO(gene = sig_genes$entrez,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)
View(go_results)

dotplot(go_results, showCategory = 20)

resultsOrdered[!is.na(resultsOrdered$symbol) & resultsOrdered$symbol == "CRISPLD2", ]


which(!is.na(resultsOrdered$symbol) & resultsOrdered$symbol == "CRISPLD2")

#13. upregulated vs downregulated gene comparison.


upregulated <- subset(resultsOrdered, padj < 0.05 & log2FoldChange > 0)
downregulated <- subset(resultsOrdered, padj < 0.05 & log2FoldChange < 0)

nrow(upregulated)
nrow(downregulated)

install.packages("AnnotationDbi")
library("AnnotationDbi")

up_entrez <- mapIds(org.Hs.eg.db,
                    keys = rownames(upregulated),
                    column = "ENTREZID",
                    keytype = "ENSEMBL",
                    multiVals = "first")

go_up <- enrichGO(gene = up_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

down_entrez <- mapIds(org.Hs.eg.db,
                      keys = rownames(downregulated),
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")

go_down <- enrichGO(gene = down_entrez,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)

dotplot(go_up, showCategory = 20, title = "Upregulated Genes - GO Biological Processes")

dotplot(go_down, showCategory = 20, title = "Downregulated Genes - GO Biological Processes")




