###### Scripts de Fer


# library(biomaRt)
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# head(listAttributes(ensembl))
# head(listFilters(ensembl))
# ```
# ```{r echo=FALSE}
# library(tximport)
# library(tidyverse)
# library(DESeq2)
# library(ggplot2)
# library(ggrepel) # libreria que evita el overlap de texto en labels
# library(rhdf5)
# ```
#
# ```{r}
# etID_egID <- getBM(attributes=c('ensembl_transcript_id','ensembl_gene_id'),  mart = ensembl) %>% relocate(ensembl_transcript_id)
#
# etID_hgnc <- getBM(attributes=c('ensembl_transcript_id', 'hgnc_symbol'), mart = ensembl)
# ```
#
# ```{r}
# # generar tabla de metadatos
# metadata.tsv <- data.frame("SRA" =c("SRR20114180","SRR20114179", "SRR20114178", "SRR20114177", "SRR20114176",
#                                           "SRR20114175", "SRR20114174", "SRR20114173", "SRR20114172", "SRR20114171"), "sample" = c("control_1","control_2","control_3", "patient_1", "patient_2", "patient_3", "patient_4", "patient_5", "patient_6", "patient_7") , "dex" = c(rep("control",3), rep("patient",7)), "species" = "Homo_sapiens")
# ```
#
# ```{r}
# # Anotacion articulo
# # Cargar los archivos por nombre del archivo (ubicacion empleando la anotacion de los 224 transcriptomas)
# samples <- metadata.tsv
# files   <- file.path("./kallisto_quant", samples$SRA,"abundance.tsv")
# names(files) <- samples$SRA
# ```
#
# ```{r}
# # Load table with trx id and gene id corrspondence
# tx2gene <- etID_egID
# # Load table with trx id and gene name corrspondence
# tx2genename <- etID_hgnc
#
# # Run tximport
# # tx2gene links transcript IDs to gene IDs for summarization
# txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar=TRUE, ignoreTxVersion=TRUE)
# txi.kallisto.name <- tximport(files, type = "kallisto", tx2gene = tx2genename, ignoreAfterBar=TRUE, ignoreTxVersion=TRUE)
# ```
#
#
# ```{r}
# names(txi.kallisto)
# ```
#
# ```{r}
# head(txi.kallisto$counts)
# ```
#
# ```{r}
# # nombre de los transcriptomas
# rownames(samples) <- samples$sample
# colnames(txi.kallisto$counts) <-rownames(samples)
# ```
#
# ```{r echo=FALSE}
# # Importacion de los datos convirtiendolos en un objeto que puede leer Deseq.
# ddsTxi_all <- DESeqDataSetFromTximport(txi.kallisto, samples, design = ~ dex) # Create a DESeq object from the tximport data
# ```
#
# ```{r echo=FALSE}
# # Prefiltrado, eliminacion de genes con bajas cuentas
# keep  <- rowSums(counts(ddsTxi_all)) >= 10
# ddsTxi_all  <- ddsTxi_all[keep,]
# dds_all <- DESeq(ddsTxi_all) # run Differential expression analysis
# ```
#
# ```{r}
# res <- results(dds_all)
# res
# ```
#
# ```{r}
# summary(res)
# ```

```{r}
# resSig <- subset(res, padj < 0.1) # subset
# head(resSig[ order(resSig$log2FoldChange), ]) # more downregulated genes
# head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ]) # genes more upregulated
```




```{r}
# # extraer UP
# all_de_gene_matrix_UP  <- subset(res, padj < 0.05 & log2FoldChange >= 2)
# write.table(all_de_gene_matrix_UP,file ="./all_DEG_kallisto_Patients.tsv", quote=FALSE, sep="\t")
# # Extraer nombres
# all_de_gene_names_UP <- rownames(all_de_gene_matrix_UP)
# 
# # extraer down genes
# all_de_gene_matrix_DOWN  <- subset(res, padj < 0.05 & log2FoldChange < -2)
# write.table(all_de_gene_matrix_DOWN,file ="./all_DEG_kallisto_Control.tsv", quote=FALSE, sep="\t")
# # Extraer nombres
# all_de_gene_names_DOWN <- rownames(all_de_gene_matrix_DOWN)
# 
# # Numero de genes expresados
# length(all_de_gene_names_UP)
# length(all_de_gene_names_DOWN)
# ```
# 
# # Visualización gráfica
# 
# **Cuentas normalizadas para graficas (rlog)**
# ```{r}
# all_normalized <- rlog(dds_all, blind=FALSE) # result rld, vst
# all_normalized_db <- as.data.frame(assay(all_normalized))
# head(all_normalized_db)
# ```
# 
# ## PCA
# ```{r}
# vst = vst(dds, blind = FALSE)
# plotPCA(vst, intgroup = "disease_state")
# ```
# 
# ## Volcano plot
# ```{r}
# de <- as.data.frame(res_all)
# # add a column of NAs
# de$diffexpressed <- "NO"
# 
# # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
# de$diffexpressed[de$log2FoldChange > 2 & de$pvalue < 0.05] <- "UP"
# # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
# de$diffexpressed[de$log2FoldChange < -2 & de$pvalue < 0.05] <- "DOWN"
# # Create a new column "names" to de, that will contain the name of a subset if genes differentially expressed (NA in case they are not)
# de$names <- NA
# # filter for a subset of interesting genes
# filter <- which(de$diffexpressed != "NO" & de$padj < 0.05 & (de$log2FoldChange >= 5  | de$log2FoldChange <= -5))
# de$names[filter] <- rownames(de)[filter]
# 
# # grafica
# png(file = "volcano-COVID.png",
#     width = 800, height = 800) # guardar el plot en formato png
# ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=names)) +
#     geom_point() +
#     scale_color_manual(values=c("blue", "black", "red")) + # cambiar colores de puntos
#     theme_minimal() +
#     geom_text_repel() +
#     xlim(-15,15)
# dev.off()
## Heatmap
# library("pheatmap")
# los primeros 20 genes
select <- order(rowMeans(counts(dds ,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("disease_state")])

#heatmap
pheatmap(assay(vst)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE)


## heatmap attempt

