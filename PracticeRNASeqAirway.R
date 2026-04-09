  # RNA-seq Analysis of Dexamethasone Treatment
  # Dataset: Airway smooth muscle cells
  
  #1. Load libraries
  # These are all the packages we need for the analysis
  library(DESeq2)          # Main differential expression package
  library(airway)          # Example dataset
  library(ggplot2)         # General plotting
  library(pheatmap)        # Heatmap visualization
  library(EnhancedVolcano) # Volcano plot visualization
  library(org.Hs.eg.db)   # Human gene annotation database, contains ENTREZ IDs
  
  #2. Load data
  data("airway")
  
  # 3. Set Up DESeq2 Object
  # Create a DESeq2 dataset
  dds <- DESeqDataSet(airway, design = ~ cell + dex)
  
  #4. Filter Low Count Genes
  # Remove genes with 10 or fewer total counts across all samples
  # These genes are too lowly expressed to give reliable results
  # Reduces dataset from 63,677 genes to 22,008 genes
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep,]
  
  # 5. Run DESeq2 Analysis
  # Does these 3 things:
  # 1. Normalizes counts so samples are fairly comparable
  # 2. Estimates variability for each gene
  # 3. Tests each gene for significant changes
  dds <- DESeq(dds)
  
  # 6. Extract Results
  # Get results comparing treated vs untreated
  results <- results(dds, contrast = c("dex", "trt", "untrt"))
  summary(results)
  
  sum(results$padj < 0.05, na.rm = TRUE)
  
  # 7. Sort Results by Significance
  resultsOrdered <- results[order(results$padj),]
  
  # 8. Add Human Readable Gene Names (ENSEMBL IDs)
  # Convert ENSG codes to gene symbols using human genome database
  resultsOrdered$symbol <- mapIds(org.Hs.eg.db,
                                  keys = rownames(resultsOrdered),
                                  column = "SYMBOL",
                                  keytype = "ENSEMBL",
                                  multiVals = "first")
  
  # View top 10 most significant genes
  head(resultsOrdered, 10)

  
  # 9. Volcano Plot 
  # Shows all genes plotted by fold change (x) vs significance (y)
  # Red dots = significant AND large fold change
  EnhancedVolcano(results,
                  lab = rownames(results),
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = 'Dexamethasone vs Untreated',
                  pCutoff = 0.05,
                  FCcutoff = 1.0)
  
  # 10. Heatmap
  # Shows expression patterns of top 20 most significant genes, organized by sample
  # Red = higher than average expression
  # Blue = lower than average expression
  top20 <- head(order(results$padj, na.last = NA), 20)
  mat <- counts(dds, normalized = TRUE)[top20,]
  pheatmap(mat,
           scale = "row",
           color = colorRampPalette(c("darkblue", "white", "darkred"))(100))
  
  #11. Using datacamp PCA tutorial, conducts Principal Component Analysis
  #determine main sources of variability in samples, reduces confounding variables
  # confirms that treatment (drug) was main source of variability
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
  
  #select top 500 variable genes
  top_var_genes <- head(order(rowVars(counts(dds, normalized = TRUE)), decreasing = TRUE), 500)
  
  #conduct PCA analysis
  mat_pca <- t(counts(dds, normalized = TRUE)[top_var_genes,])
  
  mat_scaled <- scale(mat_pca)
  
  data.pca <- prcomp(mat_scaled)
  summary(data.pca)
  
  #visualize PCA analysis results
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
  library(clusterProfiler)
  install.packages("igraph", type = "binary")
  library(igraph)
  
  sig_genes <- subset(resultsOrdered, padj < 0.05 & abs(log2FoldChange > 1)) # sorts all the statistically significant genes
  
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
  library(AnnotationDbi)
  
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
  
  
  #14. transcription factor analysis to determine which transcription factors are being regulated by dexamethasone
  BiocManager::install("dorothea")
  BiocManager:: install("viper")
  BiocManager:: install("dplyr")
  
  library(dorothea)
  library(viper)
  library(dplyr)
  
  #normalized counts matrix
  counts_matrix <- counts(dds, normalized = TRUE)
  
  #log transform the counts
  log_counts <- log2(counts_matrix + 1)
  
  #load human dorothea regulons
  data(dorothea_hs, package = "dorothea")
  regulons <- dorothea_hs %>%
    filter(confidence %in% c("A", "B"))
  #^^Gives us a list of transcription factors relating to genes with high confidence levels
  
  #correct format to be able to run viper
  viper_regulons <- df2regulon(regulons)
  viper_regulons <- split(regulons, regulons$tf)
  viper_regulons <- lapply(viper_regulons, function(x){
    tfmode<- setNames(as.numeric(x$mor), x$target)
    list(tfmode = tfmode, likelihood = as.numeric(rep(1, length(tfmode))))
  })
  
  #how active is each transcription factor in our dataset?
  tf_activities <- viper(eset = log_counts,
                         regulon = viper_regulons,
                         nes = TRUE,
                         method = "none",
                         minsize = 4,
                         eset.filter = FALSE)
  
  class(regulons$mor)
  head(regulons$mor)
  sum(is.na(regulons$mor))
  #Viper/Dorothea incompatible, switching to decoupler
  
  BiocManager::install("decoupleR")
  library(decoupleR)
  
  #convert to gene symbols 
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = rownames(log_counts),
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  #Remove NAs
  keep_genes <- !is.na(gene_symbols)
  log_counts_symbols <- log_counts[keep_genes,]
  rownames(log_counts_symbols) <- gene_symbols[keep_genes]
  
  head(rownames(log_counts_symbols))
  
  regulons_df <- dorothea_hs %>%
    filter(confidence %in% c("A", "B"))
  
  #finally, do transcription factor analysis using our newly formatted datasets
  tf_activities <- run_ulm(mat = log_counts_symbols,
                           net = regulons_df,
                           .source = "tf",
                           .target = "target",
                           .mor = "mor",
                           minsize = 4)
  head(tf_activities)
  
  #now we filter for statistically significant TF activities. 
  tf_sig <- tf_activities %>%
    filter(p_value < 0.05)
  
  nrow(tf_sig)
  
  treated_samples <- colnames(dds)[dds$dex == "trt"]
  untreated_samples <- colnames(dds)[dds$dex == "untrt"]
  
  BiocManager::install("tidyr")
  library(tidyr)
  #calculate mean activity per TF per treatment group
  tf_summary <- tf_sig %>%
    mutate(treatment = ifelse(condition %in% treated_samples,
                              "treated", "untreated")) %>%
    group_by(source, treatment) %>%
    summarise(mean_score = mean(score)) %>%
    pivot_wider(names_from = treatment,
                        values_from = mean_score) %>%
    mutate(diff = treated-untreated) %>%
    arrange(desc(abs(diff)))
  
  head(tf_summary, 20)

  
  ggplot(head(tf_summary, 20), 
         aes(x = reorder(source, diff), 
             y = diff,
             fill = diff > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("blue", "red"),
                      labels = c("Less active", "More active"),
                      name = "Direction") +
    coord_flip() +
    labs(title = "Top 20 Most Changed Transcription Factors",
         subtitle = "Dexamethasone treated vs untreated",
         x = "Transcription Factor",
         y = "Activity Difference (treated - untreated)") +
    theme_minimal()
  
  
  #15. TF Heatmap Analysis
  
  library(tidyr)
  library(dplyr)
  install.packages("tibble")
  library(tibble)
  #extract names of top 20 affected transcription factors
  top20_tf_names <- head(tf_summary$source, 20)
  
  #modify data to be heatmap compatible
  tf_heatmap_data <- tf_activities %>%
    filter(source %in% top20_tf_names) %>%
    pivot_wider(names_from = condition,
                values_from = score,
                id_cols = source) %>%
    column_to_rownames("source")
  
  #transform data into matrix format
  tf_heatmap_matrix <- as.matrix(tf_heatmap_data)  
  
  #categorize and label each sample as either treated or untreated
  #maybe should add this to the earlier gene expression heatmap as well
  sample_annotations <- data.frame(
    Treatment = ifelse(colnames(tf_heatmap_matrix) %in% treated_samples,
                       "Treated", "Untreated"),
    row.names = colnames(tf_heatmap_matrix)
  )
  
  #create actual heatmap  
  pheatmap(tf_heatmap_matrix,
           annotation_col = sample_annotations,
           scale = "row",
           color = colorRampPalette(c("darkblue", "white", "darkred"))(100),
           main = "Top 20 Transcription Factor Activities")
    
  
  
    
