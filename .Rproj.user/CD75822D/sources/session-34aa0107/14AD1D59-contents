# Differential expression analysis of dexamethasone-treated
# human airway smooth muscle cells
#
# Data: airway Bioconductor package (Himes et al., 2014; GEO: GSE52778)
# Run from the project root (open the .Rproj file first).

# EnhancedVolcano is loaded first to avoid a ggplot2 version conflict
library(EnhancedVolcano)
library(DESeq2)
library(airway)
library(org.Hs.eg.db)
library(pheatmap)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# ---- 1. Load data ----------------------------------------------------------

data(airway)

# 63,677 genes x 8 samples (4 donors, each treated and untreated)
dim(airway)
colData(airway)[, c("cell", "dex")]

# Library sizes range from ~15M to ~31M reads, so raw counts
# are not directly comparable between samples
colSums(assay(airway))

# Use untreated samples as the reference level
airway$dex <- relevel(airway$dex, ref = "untrt")


# ---- 2. Differential expression --------------------------------------------

# Donor is included in the design so that the treatment effect
# is estimated within each donor (paired design)
dds <- DESeqDataSet(airway, design = ~ cell + dex)

# Drop genes with fewer than 10 reads in total
dds <- dds[rowSums(counts(dds)) >= 10, ]
nrow(dds)

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
summary(res)

# Add gene symbols
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

res_ordered <- res[order(res$padj), ]
head(res_ordered[, c("symbol", "log2FoldChange", "padj")], 10)

write.csv(as.data.frame(res_ordered), "results/deseq2_results.csv")

# Genes with padj < 0.05 and at least a 2-fold change
sig <- res_ordered[!is.na(res_ordered$padj) & res_ordered$padj < 0.05, ]
n_up   <- sum(sig$log2FoldChange > 1)
n_down <- sum(sig$log2FoldChange < -1)
cat("Up:", n_up, " Down:", n_down, "\n")

# CRISPLD2 was highlighted in the original paper
res_ordered[which(res_ordered$symbol == "CRISPLD2"),
            c("symbol", "log2FoldChange", "padj")]


# ---- 3. PCA ----------------------------------------------------------------

vsd <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(vsd, intgroup = "dex")
pca_plot
ggplot2::ggsave("figures/pca_plot.png", pca_plot, width = 6, height = 5)

# Same PCA coloured by donor
plotPCA(vsd, intgroup = "cell")


# ---- 4. Volcano plot -------------------------------------------------------

res_plot <- res[!is.na(res$padj), ]

# Label the 8 strongest up- and down-regulated genes, plus CRISPLD2
sig_named  <- sig[!is.na(sig$symbol), ]
up_genes   <- head(sig_named$symbol[sig_named$log2FoldChange > 1], 8)
down_genes <- head(sig_named$symbol[sig_named$log2FoldChange < -1], 8)
label_genes <- c(up_genes, down_genes, "CRISPLD2")

volcano <- EnhancedVolcano(res_plot,
                           lab = res_plot$symbol,
                           x = "log2FoldChange",
                           y = "padj",
                           selectLab = label_genes,
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           pointSize = 1.5,
                           labSize = 3.5,
                           drawConnectors = TRUE,
                           widthConnectors = 0.4,
                           max.overlaps = Inf,
                           legendPosition = "right",
                           ylab = bquote(~-Log[10] ~ "adjusted" ~ italic(P)),
                           title = "Dexamethasone vs Control",
                           subtitle = "Human airway smooth muscle cells")
volcano
ggplot2::ggsave("figures/volcano_plot.png", volcano, width = 10, height = 8)


# ---- 5. Heatmap of the top 30 genes ----------------------------------------

top30 <- head(rownames(res_ordered), 30)
mat <- assay(vsd)[top30, ]
mat <- mat - rowMeans(mat)

gene_names <- res_ordered[top30, "symbol"]
rownames(mat) <- ifelse(is.na(gene_names), top30, gene_names)

anno <- as.data.frame(colData(vsd)[, c("dex", "cell")])

pheatmap(mat, annotation_col = anno)
pheatmap(mat, annotation_col = anno,
         filename = "figures/heatmap_top30.png", width = 7, height = 9)


# ---- 6. Session info -------------------------------------------------------

writeLines(capture.output(sessionInfo()), "results/session_info.txt")
