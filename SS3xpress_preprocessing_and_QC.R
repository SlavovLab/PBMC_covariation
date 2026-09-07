########################################################################################
## Processing raw data and quality control for transcriptomic dataset: ss3xpress
########################################################################################


# ── 1. Load zUMIs output and map Ensembl IDs to gene symbols ─────────────────────────

pbmc.rna <- readRDS("/projects/slavov/LK/ImmPTR/data/raw_data/ss3Xpress/zUMIs_output/Xpress.dgecounts.rds")


## ---- 1.1 Load the Ensembl annotation snapshot ----
metadata_dir     <- "/projects/slavov/LK/ImmPTR/data/meta_data"
ensembl_map_csv  <- file.path(metadata_dir, "ensembl_biotype_mapping.csv")
stopifnot(file.exists(ensembl_map_csv))

mapping_all <- read.csv(ensembl_map_csv, stringsAsFactors = FALSE, na.strings = "")
mapping_all$hgnc_symbol[is.na(mapping_all$hgnc_symbol)] <- ""

cat("\n=== 1.1 Frozen Ensembl annotation map ===\n")
cat(sprintf("  Rows: %d | unique Ensembl IDs: %d | protein-coding IDs: %d\n",
            nrow(mapping_all),
            dplyr::n_distinct(mapping_all$ensembl_gene_id),
            sum(mapping_all$gene_biotype == "protein_coding")))


## ---- 1.2 Extract the exon UMI count matrix — primary analysis matrix ----

rna.counts <- as.matrix(pbmc.rna[["umicount"]][["exon"]][["all"]])
rownames_no_version <- sub("\\.\\d+$", "", rownames(rna.counts))


## ---- 1.3 Restrict the annotation map to the exon matrix ----

mapping <- mapping_all[mapping_all$ensembl_gene_id %in% rownames_no_version, ]

cat("\n=== 1.3 Annotation scoped to exon matrix ===\n")
cat(sprintf("  Exon matrix: %d genes x %d cells\n",
            nrow(rna.counts), ncol(rna.counts)))
cat(sprintf("  IDs resolved in frozen map: %d / %d (%d absent — retired IDs)\n",
            length(intersect(rownames_no_version, mapping_all$ensembl_gene_id)),
            length(rownames_no_version),
            length(setdiff(rownames_no_version, mapping_all$ensembl_gene_id))))
cat(sprintf("  `mapping` rows: %d | protein-coding IDs: %d\n",
            nrow(mapping), sum(mapping$gene_biotype == "protein_coding")))


## ---- 1.4 Apply gene name mapping; keep Ensembl ID where mapping fails ----

rownames(rna.counts) <- rownames_no_version
gene_name_map <- mapping[match(rownames(rna.counts), mapping$ensembl_gene_id), "hgnc_symbol"]
rownames(rna.counts) <- ifelse(is.na(gene_name_map), rownames(rna.counts), gene_name_map)

# Report and remove IDs that remain unmapped (no HGNC symbol found)
unmapped_ids <- rownames(rna.counts)[grepl("^ENSG", rownames(rna.counts))]
cat("\nUnmapped Ensembl IDs remaining:", length(unmapped_ids), "\n")
rna.counts <- rna.counts[!(rownames(rna.counts) %in% unmapped_ids), ]


## ---- 1.5 Biotype composition summary (reference plot) ----

sum_RNA_types <- mapping %>%
  dplyr::group_by(gene_biotype) %>%
  dplyr::summarise(Number_RNA_Species = dplyr::n_distinct(hgnc_symbol), .groups = "drop")

sum_RNA_types$bar_color <- ifelse(sum_RNA_types$gene_biotype == "protein_coding", "red", "black")

ggplot(sum_RNA_types,
       aes(x = reorder(gene_biotype, Number_RNA_Species),
           y = log10(Number_RNA_Species),
           fill = bar_color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red" = "red", "black" = "grey50")) +
  labs(
    title = paste0("PBMC RNA Species Quantified by SS3Xpress, n = ", nrow(mapping)),
    y     = expression("Count, Log"["10"]),
    x     = "RNA Biotype"
  ) +
  theme_light() +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title  = element_text(size = 20),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  )


## ---- 1.6 Restrict to protein-coding genes and aggregate duplicate symbols ----

selected_genes <- mapping[mapping$gene_biotype == "protein_coding", "hgnc_symbol"]
selected_genes <- selected_genes[selected_genes != ""]
rna.counts     <- rna.counts[rownames(rna.counts) %in% selected_genes, ]

dense_counts      <- as.matrix(rna.counts)
aggregated_counts <- rowsum(dense_counts, group = rownames(rna.counts))
rna.counts        <- as(aggregated_counts, "dgCMatrix")

cat("\n=== 1.6 Gene set after protein-coding filter ===\n")
cat(sprintf("  Genes retained: %d | cells: %d\n", nrow(rna.counts), ncol(rna.counts)))
rm(dense_counts, aggregated_counts)



# ── 2. Load donor barcode metadata ───────────────────────────────────────────────────

barcode_donor_meta_labels_full <- read.csv(
  "/projects/slavov/LK/ImmPTR/data/meta_data/ss3Xpress/barcode_annotation_slavov_xpressJMIndex.csv"
)



# ── 3. Seurat object construction and cell-level QC ───────────────────────────────────

pbmc.so <- CreateSeuratObject(counts = rna.counts, min.cells = 3, min.features = 200)

meta_idx <- match(colnames(pbmc.so), barcode_donor_meta_labels_full$XC)

stopifnot(
  !anyNA(meta_idx),
  length(meta_idx) == ncol(pbmc.so)
)

pbmc.so <- AddMetaData(pbmc.so,
                       metadata = barcode_donor_meta_labels_full$donor[meta_idx],         col.name = "donor")
pbmc.so <- AddMetaData(pbmc.so,
                       metadata = barcode_donor_meta_labels_full$PlateName[meta_idx],     col.name = "PlateName")
pbmc.so <- AddMetaData(pbmc.so,
                       metadata = barcode_donor_meta_labels_full$lonza_donorid[meta_idx], col.name = "lonza_donorid")
pbmc.so <- AddMetaData(pbmc.so,
                       metadata = barcode_donor_meta_labels_full$XC_NovaSeq[meta_idx],    col.name = "XC_NovaSeq")

pbmc.so$PlateName <- gsub("Slavov_", "", pbmc.so$PlateName)

pbmc.so[["percent.mt"]]   <- PercentageFeatureSet(pbmc.so, pattern = "^MT-")
pbmc.so[["percent.ribo"]] <- PercentageFeatureSet(pbmc.so, pattern = "^RPL|^RPS")

# QC diagnostic plots
VlnPlot(pbmc.so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)

feature.scat1 <- FeatureScatter(pbmc.so, feature1 = "nCount_RNA", feature2 = "percent.mt")
feature.scat2 <- FeatureScatter(pbmc.so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
feature.scat1 + feature.scat2

# Apply cell-level QC filters
pbmc.so <- subset(pbmc.so, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

# Log-normalize
pbmc.so <- NormalizeData(pbmc.so, normalization.method = "LogNormalize", scale.factor = 10000)


# ── 4. Gene filtering: detection rate + mean expression ──────────────────────────────

counts_filtered  <- GetAssayData(pbmc.so, layer = "counts")
data_normalized  <- GetAssayData(pbmc.so, layer = "data")
n_cells          <- ncol(counts_filtered)

cat("Cells after QC filtering:", n_cells, "\n")
cat("Genes before gene filtering:", nrow(counts_filtered), "\n")

detection_threshold  <- 0.03   # gene must be detected in ≥3% of cells
mean_expr_threshold  <- 0.05   # mean log-normalized expression
min_cells_expressing <- detection_threshold * n_cells

cells_expressing <- Matrix::rowSums(counts_filtered > 0)
mean_norm_expr   <- Matrix::rowMeans(data_normalized)

detection_pass <- cells_expressing >= min_cells_expressing
expression_pass <- mean_norm_expr   >= mean_expr_threshold
genes_to_keep   <- detection_pass & expression_pass

cat("\n=== Gene Filtering Summary ===\n")
cat("Detection rate filter (>=", round(min_cells_expressing), "cells):", sum(detection_pass), "genes pass\n")
cat("Mean expression filter (>=", mean_expr_threshold, "):", sum(expression_pass), "genes pass\n")
cat("Both filters combined:", sum(genes_to_keep), "genes pass\n")
cat("Genes removed:", sum(!genes_to_keep), "\n")
cat("Percentage retained:", round(100 * sum(genes_to_keep) / length(genes_to_keep), 1), "%\n")

# Diagnostic data frame
gene_stats <- data.frame(
  gene               = rownames(counts_filtered),
  n_cells_expressing = cells_expressing,
  detection_rate     = cells_expressing / n_cells,
  mean_log_expr      = mean_norm_expr,
  passes_filter      = genes_to_keep
)

p_detect <- ggplot(gene_stats, aes(x = detection_rate, fill = passes_filter)) +
  geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = detection_threshold, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "#00BA38", "FALSE" = "#F8766D"),
                    labels = c("TRUE" = "Retained", "FALSE" = "Removed")) +
  labs(title = "Gene Detection Rate Distribution",
       x = "Detection Rate (Fraction of Cells)", y = "Number of Genes", fill = "Filter Status") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.position = "top")

p_expr <- ggplot(gene_stats, aes(x = mean_log_expr, fill = passes_filter)) +
  geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = mean_expr_threshold, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "#00BA38", "FALSE" = "#F8766D"),
                    labels = c("TRUE" = "Retained", "FALSE" = "Removed")) +
  labs(title = "Gene Mean Expression Distribution",
       x = "Mean Log-Normalized Expression", y = "Number of Genes", fill = "Filter Status") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.position = "top")

p_scatter <- ggplot(gene_stats, aes(x = detection_rate, y = mean_log_expr, color = passes_filter)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = mean_expr_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = detection_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c("TRUE" = "#00BA38", "FALSE" = "#F8766D"),
                     labels = c("TRUE" = "Retained", "FALSE" = "Removed")) +
  labs(title = "Detection Rate vs Mean Expression",
       x = "Detection Rate (Fraction of Cells)", y = "Mean Log-Normalized Expression",
       color = "Filter Status") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.position = "top")

(p_detect / p_expr) | p_scatter

# Apply gene filter and update Seurat object
pbmc.so <- subset(pbmc.so, features = rownames(counts_filtered)[genes_to_keep])
cat("\nFinal gene count:", nrow(pbmc.so), "\n")
cat("Final cell count:", ncol(pbmc.so), "\n")


# ── 5. Dimensionality reduction and clustering ────────────────────────────────────────

pbmc.so <- FindVariableFeatures(pbmc.so, selection.method = "vst", nfeatures = 3000)
pbmc.so <- ScaleData(pbmc.so, vars.to.regress = "donor")
pbmc.so <- RunPCA(pbmc.so, features = VariableFeatures(pbmc.so))

ElbowPlot(pbmc.so, ndims = 35)
DimHeatmap(pbmc.so, dims = 1:15, cells = 500, balanced = TRUE)

pbmc.so <- FindNeighbors(pbmc.so, dims = 1:50)
pbmc.so <- FindClusters(pbmc.so, resolution = 0.9)
pbmc.so <- RunUMAP(pbmc.so, dims = 1:50)

DimPlot(pbmc.so, reduction = "umap", label = TRUE, pt.size = 0.7)
DimPlot(pbmc.so, reduction = "pca", label = TRUE, pt.size = 0.7)
DimPlot(pbmc.so, reduction = "umap", group.by = "donor", label = TRUE, pt.size = 0.3, alpha = 0.3)
DimPlot(pbmc.so, reduction = "pca",  group.by = "donor", label = TRUE, pt.size = 0.3, alpha = 0.3)


# ── 6. DE marker discovery for cluster annotation ────────────────────────────────────

pbmc.markers <- FindAllMarkers(
  pbmc.so,
  only.pos        = TRUE,
  min.pct         = 0.10,
  logfc.threshold = 0.10
)
significant_markers <- pbmc.markers %>% dplyr::filter(p_val_adj < 0.05)
top6                <- significant_markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)

GeneMarkersUsed_numbered <- significant_markers %>%
  dplyr::filter(grepl(
    "PTPRC|GZMA|GZMK|CD4|CD3E|CD3D|IL7R|CCR7|SELL|GATA3|CD8A|CD8B|GNLY|GZMB|S1PR5|NKG7|MS4A1|CD79A|FCGR3A|MS4A7|ITGAX|CD14|S100A8|S100A9|IL1B|CCR1|FCER1A|CD1C|SPIB|JCHAIN",
    gene))

# Inspect DotPlot before committing to annotation
DotPlot(
  pbmc.so,
  features       = unique(GeneMarkersUsed_numbered$gene),
  cluster.idents = TRUE
) +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "#CE2029", midpoint = 0) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(color = "Average Expression")


# ── 7. First annotation pass — all populations including Dendritic cells ──────────────

new.cluster.ids.all.pops <- c("CD4 T cells", "Monocytes", "B cells", "CD8 T cells",
                              "CD8 T cells", "NK cells", "Dendritic cells",
                              "Monocytes", "Dendritic cells")

names(new.cluster.ids.all.pops) <- levels(pbmc.so)

pbmc.so.all.pops <- pbmc.so
pbmc.so.all.pops <- RenameIdents(pbmc.so.all.pops, new.cluster.ids.all.pops)

celltype_colors_all_pops <- c(
  "CD4 T cells"    = "#FFB347",
  "CD8 T cells"    = "#5A8DEE",
  "B cells"        = "#9B3D6A",
  "NK cells"       = "#008B72",
  "Monocytes"      = "#D2BD92",
  "Dendritic cells" = "#C43E00"
)

umap_df_allpops    <- as.data.frame(Embeddings(pbmc.so.all.pops, "umap"))
umap_df_allpops$celltype <- pbmc.so.all.pops@active.ident
centroids_allpops  <- umap_df_allpops %>%
  group_by(celltype) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

sct.umap.allpops <- ggplot(umap_df_allpops, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 1, alpha = 0.4) +
  scale_color_manual(values = celltype_colors_all_pops) +
  geom_text_repel(
    data         = centroids_allpops,
    aes(label    = celltype),
    size         = 4,
    color        = "black",
    box.padding  = 1.5,
    max.overlaps = Inf
  ) +
  labs(title = "PBMCs", x = "UMAP 1", y = "UMAP 2") +
  theme_light(base_size = 16) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "none"
  )

sct.umap.allpops

# Factor idents to enforce consistent cell type ordering before DotPlot

Idents(pbmc.so.all.pops) <- factor(
  Idents(pbmc.so.all.pops),
  levels = c("Monocytes", "NK cells", "CD8 T cells", "CD4 T cells", "B cells", "Dendritic cells")
)

DotPlot(
  pbmc.so.all.pops,
  features       = unique(GeneMarkersUsed_numbered$gene),
  cluster.idents = FALSE
) +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "#CE2029", midpoint = 0) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(color = "Average Expression")


# ── 8. Final annotation — rename pbmc.so, recompute markers ──────────────────────────

names(new.cluster.ids.all.pops) <- levels(pbmc.so)
pbmc.so <- RenameIdents(pbmc.so, new.cluster.ids.all.pops)

# Recompute markers on annotated cell types
pbmc.markers <- FindAllMarkers(
  pbmc.so,
  only.pos        = TRUE,
  min.pct         = 0.10,
  logfc.threshold = 0.10
)
significant_markers <- pbmc.markers %>% dplyr::filter(p_val_adj < 0.05)
top6                <- significant_markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)

GeneMarkersUsed <- significant_markers %>%
  dplyr::filter(grepl(
    "CD3E|CD8A|GZMA|MS4A1|CD74|CD37|CD14|GNLY|ITGAM|ZAP70|FCER1G|CD247|CD1C",
    gene))

DotPlot(
  pbmc.so,
  features       = unique(GeneMarkersUsed$gene),
  cluster.idents = TRUE
) +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "#CE2029", midpoint = 0) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(color = "Average Expression", title = "RNA Markers") +
  theme(plot.title = element_text(hjust = 0.5))

# Annotated UMAP — main populations (including Dendritic cells)
celltype_colors <- c(
  "CD4 T cells"    = "#FFB347",
  "CD8 T cells"    = "#5A8DEE",
  "B cells"        = "#9B3D6A",
  "NK cells"       = "#008B72",
  "Monocytes"      = "#D2BD92",
  "Dendritic cells" = "#C43E00"
)

umap_df          <- as.data.frame(Embeddings(pbmc.so, "umap"))
umap_df$celltype <- pbmc.so@active.ident
centroids        <- umap_df %>%
  group_by(celltype) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

sct.umap.mainpops <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(values = celltype_colors) +
  geom_text_repel(
    data         = centroids,
    aes(label    = celltype),
    size         = 6,
    color        = "black",
    box.padding  = 0.8,
    max.overlaps = Inf
  ) +
  labs(title = "mRNA", x = "UMAP 1", y = "UMAP 2") +
  theme_light(base_size = 16) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "none"
  )

sct.umap.mainpops

# Cell-type order for all subsequent per-cell-type plots
celltype_order_rna <- c("Monocytes", "NK cells", "CD8 T cells", "CD4 T cells", "B cells")


# ── 9. Detection sensitivity analysis ────────────────────────────────────────────────

## 9a. Define the final QC-passed cells (from pbmc.so, excluding Dendritic cells)
##     These are the cells used in all downstream analyses
final_cells <- colnames(pbmc.so)[pbmc.so@active.ident != "Dendritic cells"]

## 9b. Helper: from a zUMIs sparse matrix, count unique genes detected per cell restricted to the final cell set
count_genes_per_cell <- function(mat, cell_ids) {
  # Subset to final cells that exist in this matrix
  cells_in_mat <- intersect(cell_ids, colnames(mat))
  mat_sub      <- mat[, cells_in_mat, drop = FALSE]
  Matrix::colSums(mat_sub > 0)
}

## 9c. Extract all four raw zUMIs count matrices and restrict to final cells

# readcount / exon
rc_exon_raw  <- pbmc.rna[["readcount"]][["exon"]][["all"]]
rc_exon_genes <- count_genes_per_cell(rc_exon_raw, final_cells)

# readcount / inex (intron + exon)
rc_inex_raw  <- pbmc.rna[["readcount"]][["inex"]][["all"]]
rc_inex_genes <- count_genes_per_cell(rc_inex_raw, final_cells)

# umicount / exon
umi_exon_raw  <- pbmc.rna[["umicount"]][["exon"]][["all"]]
umi_exon_genes <- count_genes_per_cell(umi_exon_raw, final_cells)

# umicount / inex (intron + exon)
umi_inex_raw  <- pbmc.rna[["umicount"]][["inex"]][["all"]]
umi_inex_genes <- count_genes_per_cell(umi_inex_raw, final_cells)

## 9d. UMI exon with all gene filters applied — use the processed Seurat data layer
##     (protein-coding, detection-rate filtered, mean-expression filtered)
expr_mat_filtered   <- GetAssayData(pbmc.so, layer = "data")
cells_no_dc         <- intersect(final_cells, colnames(expr_mat_filtered))
umi_exon_filtered_genes <- Matrix::colSums(
  expr_mat_filtered[, cells_no_dc, drop = FALSE] > 0
)

## 9e. Assemble long-format data frame for plotting
##     Label order reflects increasing stringency of quantification
matrix_labels <- c(
  "Reads\nIntron+Exon",
  "Reads\nExon",
  "UMIs\nIntron+Exon",
  "UMIs\nExon",
  "UMIs\nExon\nProtein-coding,\nfiltered"
  )

sensitivity_df <- bind_rows(
  data.frame(Matrix = matrix_labels[1], Genes_Detected = as.numeric(rc_inex_genes)),
  data.frame(Matrix = matrix_labels[2], Genes_Detected = as.numeric(rc_exon_genes)),
  data.frame(Matrix = matrix_labels[3], Genes_Detected = as.numeric(umi_inex_genes)),
  data.frame(Matrix = matrix_labels[4], Genes_Detected = as.numeric(umi_exon_genes)),
  data.frame(Matrix = matrix_labels[5], Genes_Detected = as.numeric(umi_exon_filtered_genes))
)
sensitivity_df$Matrix <- factor(sensitivity_df$Matrix, levels = matrix_labels)

# Color the gene-filtered UMI exon panel to distinguish it as the analytical matrix
sensitivity_fills <- c(
  "Reads\nIntron+Exon"          = "grey40",
  "Reads\nExon"                 = "grey40",
  "UMIs\nIntron+Exon"           = "grey40",
  "UMIs\nExon"                  = "grey40",
  "UMIs\nExon\nProtein-coding,\nfiltered" = "#0072B2"
)

SCT_Detection_Sensitivity <- ggplot(
  sensitivity_df,
  aes(x = Matrix, y = Genes_Detected, fill = Matrix)
) +
  geom_violin(trim = FALSE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "black",
               fill = "white", alpha = 0.8, linewidth = 0.3) +
  scale_fill_manual(values = sensitivity_fills) +
  labs(
    x = NULL,
    y = "# Genes / Cell"
  ) +
  theme_light(base_size = 8) +
  theme(
    axis.title      = element_text(size = 8),
    axis.title.y    = element_text(size = 8),
    axis.text.y     = element_text(size = 7),
    axis.text.x = element_text(size = 5, lineheight = 0.75),
    legend.position = "none",
    plot.margin     = margin(2, 2, 2, 2)
  )

SCT_Detection_Sensitivity

## 9f. Genes detected per cell by cell type (UMI exon, gene-filtered — analysis matrix)
celltype_df <- data.frame(
  Cell            = names(umi_exon_filtered_genes),
  genes_quantified = as.numeric(umi_exon_filtered_genes)
)
celltype_df$CellType <- pbmc.so@active.ident[celltype_df$Cell]
celltype_df          <- celltype_df %>% filter(CellType != "Dendritic cells")
celltype_df$CellType <- factor(celltype_df$CellType, levels = celltype_order_rna)

# SCT_IDs_Across_CellTypes

## 9g. Genes detected per cell by cell type (UMI exon, gene-filtered)

SCT_IDs_Across_CellTypes <- ggplot(
  celltype_df,
  aes(x = CellType, y = genes_quantified, fill = CellType)
) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", alpha = 0.8) +
  scale_fill_manual(values = celltype_colors) +
  ylim(0, 3200) +
  labs(x = "", y = "# Unique Genes Detected") +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 12),
    axis.title.y    = element_text(size = 12)
  )

SCT_IDs_Across_CellTypes


# ── 10. Data completeness — fraction of cells each gene is detected in ────────────────

expr_mat      <- GetAssayData(pbmc.so, layer = "data")
total_cells   <- ncol(expr_mat)
gene_pct_det  <- (Matrix::rowSums(expr_mat > 0) / total_cells) * 100

completeness_df <- data.frame(
  Gene                 = rownames(expr_mat),
  PercentCellsDetected = gene_pct_det
)

Data_Completeness_RNA_Dataset <- ggplot(completeness_df, aes(x = PercentCellsDetected)) +
  geom_histogram(bins = 40, alpha = 1) +
  labs(
    y     = "# mRNAs",
    x     = "% of Single Cells",
    title = paste0("# Cells = ", total_cells)
  ) +
  theme_light() +
  theme(
    axis.title  = element_text(size = 22),
    axis.text   = element_text(size = 18),
    plot.title  = element_text(size = 22, hjust = 0.5)
  )

Data_Completeness_RNA_Dataset


# ── 11. ICA visualization (non-Dendritic cells only) ─────────────────────────────────

cells_no_dc    <- colnames(pbmc.so)[pbmc.so@active.ident != "Dendritic cells"]
scaled_mat_sub <- GetAssayData(pbmc.so, assay = "RNA", layer = "scale.data")[, cells_no_dc]

set.seed(123)
ica_res <- fastICA(t(scaled_mat_sub), n.comp = 2, alg.typ = "parallel")

ica_df <- data.frame(
  cell     = cells_no_dc,
  IC1      = ica_res$S[, 1],
  IC2      = ica_res$S[, 2],
  CellType = pbmc.so@active.ident[cells_no_dc]
)

celltype_colors_no_dc <- c(
  "CD4 T cells" = "#FFB347",
  "CD8 T cells" = "#5A8DEE",
  "B cells"     = "#9B3D6A",
  "NK cells"    = "#008B72",
  "Monocytes"   = "#D2BD92"
)


sct.ica <- ggplot(ica_df, aes(x = IC1, y = IC2, color = CellType)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = celltype_colors_no_dc) +
  facet_wrap(~ "mRNA") +   # creates a single facet strip labeled "mRNA"
  labs(x = "Independent Component 1", 
       y = "Independent Component 2") +
  theme_light(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.position = "none",
    strip.background = element_rect(fill = "grey60", color = NA),
    strip.text = element_text(color = "white", face = "bold", size = 8),
    plot.margin = margin(2, 2, 2, 2)
  )


sct.ica



# ── 12. Save outputs ───────────────────

saveRDS(pbmc.so, file = "/projects/slavov/LK/ImmPTR/data/final_data/UMI_EXON_ss3xpress_PBMC.rds")







############ SF1 Panel C — Canonical marker dotplot (mRNA side)
############ Uses pbmc.markers (already computed in Section 8) + custom ggplot

library(Seurat)
library(dplyr)
library(ggplot2)

### ========================================
### CONFIGURATION
### ========================================

# Cell type order with Dendritic cells at the bottom of the y-axis
celltype_order_rna_full <- c("Monocytes", "NK cells", "CD8 T cells",
                             "CD4 T cells", "B cells", "Dendritic cells")

# Same canonical marker vector used for the protein panel
canonical_markers <- c(
  "CD3D", "CD3E", "CD247", "ZAP70",   # T cell core
  "CD8A", "GZMA",                      # CD8 T
  "GNLY",                               # NK
  "MS4A1", "CD74", "CD37",             # B cell
  "CD14", "ITGAM", "FCER1G",           # Monocyte / myeloid
  "FCER1A", "CD1C"             # Dendritic cells (cDC2 + pDC)
)

### ========================================
### ANALYSIS BLOCK — reuse pbmc.markers
### ========================================

rna_marker_stats <- purrr::map_dfr(celltype_order_rna_full, function(ct) {
  message("mRNA marker stats for: ", ct)
  res <- FindMarkers(
    pbmc.so,
    ident.1         = ct,
    features        = canonical_markers,
    min.pct         = 0,
    logfc.threshold = 0,
    only.pos        = FALSE
  )
  res$gene    <- rownames(res)
  res$cluster <- ct
  res
})

cat("\nSF1 Panel C — mRNA marker stats computed:\n")
print(dplyr::count(rna_marker_stats, cluster))

### ========================================
### Order markers by cell-type specificity
### ========================================

# For each canonical marker, find the cell type where it has max avg_log2FC
marker_assignment_rna <- rna_marker_stats %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, assigned_cluster = cluster, top_log2FC = avg_log2FC)

marker_order_rna <- marker_assignment_rna %>%
  dplyr::mutate(assigned_cluster = factor(assigned_cluster,
                                          levels = celltype_order_rna_full)) %>%
  dplyr::arrange(assigned_cluster, dplyr::desc(top_log2FC)) %>%
  dplyr::pull(gene)

cat("\nSF1 Panel C — mRNA marker x-axis order:\n")
print(marker_order_rna)

dotplot_rna_df <- rna_marker_stats %>%
  dplyr::filter(gene %in% marker_order_rna) %>%
  dplyr::mutate(
    gene    = factor(gene,    levels = marker_order_rna),
    cluster = factor(cluster, levels = rev(celltype_order_rna_full)),
    avg_log2FC_clip = pmax(pmin(avg_log2FC, 2), -2)
  )

### ========================================
### PLOTTING BLOCK — custom ggplot dotplot
### ========================================

C_rna_dotplot <- ggplot(dotplot_rna_df,
                        aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct.1, fill = avg_log2FC_clip),
             shape = 21, color = "grey30", stroke = 0.2) +
  scale_size_continuous(
    name   = "% Detected",
    range  = c(0, 3),
    limits = c(0, 1),
    breaks = c(0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_fill_gradient2(
    name     = expression(log[2]~"FC"),
    low      = "darkblue",
    mid      = "white",
    high     = "#CE2029",
    midpoint = 0,
    limits   = c(-2, 2)
  ) +
  facet_wrap(~ "mRNA") +
  labs(x = NULL, y = NULL) +
  theme_light(base_size = 6.5) +
  theme(
    axis.text.x          = element_text(size = 5.5, angle = 45, hjust = 1),
    axis.text.y          = element_text(size = 5.5),
    strip.background     = element_rect(fill = "grey60", color = NA),
    strip.text           = element_text(color = "white", face = "bold", size = 6.5),
    legend.direction     = "horizontal",
    legend.text          = element_text(size = 5),
    legend.title         = element_text(size = 5.5),
    legend.key.size      = unit(2, "mm"),
    legend.margin        = margin(0, 0, 0, 0),
    legend.box.spacing   = unit(1, "mm"),
    panel.grid.minor     = element_blank(),
    plot.margin          = margin(2, 2, 2, 2)
  ) +
  guides(
    size = guide_legend(position = "top",                                    # CHANGED
                        override.aes = list(fill = "grey70")),
    fill = guide_colorbar(position = "bottom",                                # CHANGED
                          barwidth = unit(20, "mm"),
                          barheight = unit(2, "mm"))
  )

C_rna_dotplot




############ SF1 Panel I (mRNA) — Abundance rank plot with protein-detected genes highlighted
############ Uses pre-gene-filter mRNA data (cell QC only, no gene-level filter)
############ Mirrors the Panel I protein rank plot structure

library(Seurat)
library(tidyverse)
library(ggrepel)
library(ggplot2)

### ========================================
### CONFIGURATION
### ========================================

# Colors mirror the protein panel: bulk in grey, highlighted subset in protein modality color
modality_palette <- c("mRNA"          = "grey40",
                      "mRNA + protein" = "#D55E00")

### ========================================
### ANALYSIS BLOCK
### ========================================

## ---- 1. Regenerate pre-gene-filter Seurat object ----
##       Starts from rna.counts (built in Section 1 of the mRNA processing script).
##       Applies cell-level QC and log-normalization but NO gene-level filter.
cat("\n=== 1. Rebuilding pre-gene-filter Seurat object ===\n")

pbmc.so.prefilter <- CreateSeuratObject(
  counts       = rna.counts,
  min.cells    = 3,
  min.features = 200
)

# Key-based metadata attachment (mirrors Section 3)
prefilter_meta_idx <- match(colnames(pbmc.so.prefilter),
                            barcode_donor_meta_labels_full$XC)
stopifnot(!anyNA(prefilter_meta_idx))

pbmc.so.prefilter <- AddMetaData(pbmc.so.prefilter,
                                 metadata = barcode_donor_meta_labels_full$donor[prefilter_meta_idx],
                                 col.name = "donor")
pbmc.so.prefilter <- AddMetaData(pbmc.so.prefilter,
                                 metadata = barcode_donor_meta_labels_full$PlateName[prefilter_meta_idx],
                                 col.name = "PlateName")
pbmc.so.prefilter <- AddMetaData(pbmc.so.prefilter,
                                 metadata = barcode_donor_meta_labels_full$lonza_donorid[prefilter_meta_idx],
                                 col.name = "lonza_donorid")
pbmc.so.prefilter <- AddMetaData(pbmc.so.prefilter,
                                 metadata = barcode_donor_meta_labels_full$XC_NovaSeq[prefilter_meta_idx],
                                 col.name = "XC_NovaSeq")
pbmc.so.prefilter$PlateName <- gsub("Slavov_", "", pbmc.so.prefilter$PlateName)

pbmc.so.prefilter[["percent.mt"]] <- PercentageFeatureSet(
  pbmc.so.prefilter, pattern = "^MT-"
)

# Cell-level QC identical to Section 3 of the mRNA script
pbmc.so.prefilter <- subset(
  pbmc.so.prefilter,
  subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10
)

# Log-normalize exactly as in the main script
pbmc.so.prefilter <- NormalizeData(
  pbmc.so.prefilter,
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)

cat(sprintf("  Cells in pre-filter object: %d\n", ncol(pbmc.so.prefilter)))
cat(sprintf("  Genes in pre-filter object: %d\n", nrow(pbmc.so.prefilter)))

## ---- 2. Select non-DC cells (using the annotated pbmc.so's cell IDs) ----
cat("\n=== 2. Selecting non-DC cells ===\n")

non_dc_cells <- colnames(pbmc.so)[pbmc.so@active.ident != "Dendritic cells"]
cells_to_use <- intersect(non_dc_cells, colnames(pbmc.so.prefilter))

cat(sprintf("  Non-DC cells in annotated pbmc.so: %d\n", length(non_dc_cells)))
cat(sprintf("  Cells used for rank analysis: %d\n", length(cells_to_use)))

## ---- 3. Compute mean log2 abundance per gene ----
##       Seurat's "data" layer is ln(x+1) normalized; convert to log2 by dividing by ln(2).
##       Vectorized via rowMeans.
cat("\n=== 3. Computing per-gene mean log2 abundance ===\n")

rna_norm <- GetAssayData(pbmc.so.prefilter, layer = "data", assay = "RNA")
rna_mat_log2 <- as.matrix(rna_norm[, cells_to_use] / log(2))

gene_mean <- rowMeans(rna_mat_log2, na.rm = TRUE)
gene_n_detected <- rowSums(rna_mat_log2 > 0)

## ---- 4. Build rank table with protein-detection flag ----
protein_genes <- rownames(sc.batch_cor.ab) # sc.batch_cor.ab from protein data

rna_abundance <- data.frame(
  Gene           = rownames(rna_mat_log2),
  mean_abundance = gene_mean,
  n_detected     = gene_n_detected,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(n_detected > 0) %>%
  dplyr::arrange(dplyr::desc(mean_abundance)) %>%
  dplyr::mutate(
    rank                 = dplyr::row_number(),
    detected_in_protein  = Gene %in% protein_genes
  )

cat(sprintf("  Total genes ranked: %d\n", nrow(rna_abundance)))
cat(sprintf("  Genes also detected in protein dataset: %d\n",
            sum(rna_abundance$detected_in_protein)))

## ---- 5. Identify top 5 and bottom 5 genes WITHIN the protein-detected subset ----
protein_detected_ranked <- rna_abundance %>%
  dplyr::filter(detected_in_protein) %>%
  dplyr::arrange(rank)

top5_genes <- protein_detected_ranked %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::pull(Gene)

bottom5_genes <- protein_detected_ranked %>%
  dplyr::slice_tail(n = 5) %>%
  dplyr::pull(Gene)

genes_to_label <- c(top5_genes, bottom5_genes)

cat("\n=== Top 5 most abundant protein-detected mRNAs ===\n")
print(protein_detected_ranked %>% dplyr::slice_head(n = 5) %>%
        dplyr::select(rank, Gene, mean_abundance))
cat("\n=== Bottom 5 least abundant protein-detected mRNAs ===\n")
print(protein_detected_ranked %>% dplyr::slice_tail(n = 5) %>%
        dplyr::select(rank, Gene, mean_abundance))

rna_abundance <- rna_abundance %>%
  dplyr::mutate(label = ifelse(Gene %in% genes_to_label, Gene, ""))

## ---- 6. Summary statistics (for captions / methods) ----
cat("\n=== 6. Summary ===\n")
cat(sprintf("  Highest:  %s (log2 = %.2f)\n",
            rna_abundance$Gene[1], rna_abundance$mean_abundance[1]))
cat(sprintf("  Lowest:   %s (log2 = %.2f)\n",
            rna_abundance$Gene[nrow(rna_abundance)],
            rna_abundance$mean_abundance[nrow(rna_abundance)]))
cat(sprintf("  Dynamic range: %.2f log2 units\n",
            diff(range(rna_abundance$mean_abundance))))

prot_detect_stats <- rna_abundance %>%
  dplyr::filter(detected_in_protein) %>%
  dplyr::summarise(
    n          = dplyr::n(),
    med_abund  = median(mean_abundance),
    med_rank   = median(rank),
    min_rank   = min(rank),
    max_rank   = max(rank)
  )
cat(sprintf(
  "  Protein-detected genes: n = %d | median rank = %d (of %d) | range %d–%d\n",
  prot_detect_stats$n, round(prot_detect_stats$med_rank),
  nrow(rna_abundance),
  prot_detect_stats$min_rank, prot_detect_stats$max_rank
))

### ========================================
### PLOTTING BLOCK — SF1 Panel I (mRNA)
### ========================================

# Order so protein-detected genes plot on top of the bulk cloud
rna_abundance <- rna_abundance %>%
  dplyr::arrange(detected_in_protein)

I_rna_rank <- ggplot(rna_abundance,
                     aes(x = rank, y = mean_abundance)) +
  geom_point(
    data = dplyr::filter(rna_abundance, !detected_in_protein),
    aes(color = "mRNA"),
    size = 2.5, alpha = 0.5
  ) +
  geom_point(
    data = dplyr::filter(rna_abundance, detected_in_protein),
    aes(color = "mRNA + protein"),
    size = 1, alpha = 0.15
  ) +
  scale_color_manual(
    name   = "Detection",
    values = modality_palette,
    breaks = c("mRNA", "mRNA + protein")
  )+
  facet_wrap(~ "mRNA") +
  xlim(0, 15500) +
  labs(
    x = "Gene Rank",
    y = expression("Mean log"[2]~"(abundance)")
  ) +
  theme_light(base_size = 7) +
  theme(
    axis.title           = element_text(size = 7.5),
    axis.text            = element_text(size = 6.5),
    strip.background     = element_rect(fill = "grey60", color = NA),
    strip.text           = element_text(color = "white", face = "bold", size = 7),
    legend.position      = c(0.5, 0.97),
    legend.justification = c(0.5, 1),
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 6.5),
    legend.text          = element_text(size = 6),
    legend.key.size      = unit(2.5, "mm"),
    legend.background    = element_rect(fill = NA, color = NA),
    legend.margin        = margin(1, 2, 1, 2),
    panel.grid.minor     = element_blank(),
    plot.margin          = margin(2, 2, 2, 2)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

I_rna_rank

saveRDS(pbmc.so.prefilter, "/projects/slavov/LK/ImmPTR/data/final_data/ss3xpress_prefilter.rds")






############ SF1 Panel — FACS marker abundance (mRNA-defined cell types)


library(dplyr)
library(tidyr)
library(ggplot2)

### ========================================
### CONFIGURATION — FACS marker panel
### ========================================

# All 17 immune-relevant FACS marker channels in the panel (FSC.A and SSC.A
# excluded — scatter, not markers).
facs_markers <- c(
  "HLA.DR.APC", "CD16.A7", "CD3.A700", "CD8.BUV395", "CD4.BUV661",
  "CD19.BUV737", "CCR7.BV421", "CD11c.BV605", "FCER1a.BV650",
  "CD123.BV711", "CD45RA.BV786", "CD25.VioBright.FITC", "CX3CR1.PE",
  "CD14.PC5", "CD127.PC7", "CD56.PE.CF594", "CD94.BB700"
)

# Pretty marker labels (channel suffix stripped) for plotting
facs_marker_labels <- c(
  "HLA.DR.APC"          = "HLA-DR",
  "CD16.A7"             = "CD16",
  "CD3.A700"            = "CD3",
  "CD8.BUV395"          = "CD8",
  "CD4.BUV661"          = "CD4",
  "CD19.BUV737"         = "CD19",
  "CCR7.BV421"          = "CCR7",
  "CD11c.BV605"         = "CD11c",
  "FCER1a.BV650"        = "FceR1a",
  "CD123.BV711"         = "CD123",
  "CD45RA.BV786"        = "CD45RA",
  "CD25.VioBright.FITC" = "CD25",
  "CX3CR1.PE"           = "CX3CR1",
  "CD14.PC5"            = "CD14",
  "CD127.PC7"           = "CD127",
  "CD56.PE.CF594"       = "CD56",
  "CD94.BB700"          = "CD94"
)

# Cell type order — match the mRNA/protein panels.
# Keep Dendritic cells at the bottom for consistency with Panel C.
celltype_order_facs <- c("Monocytes", "NK cells", "CD8 T cells",
                         "CD4 T cells", "B cells", "Dendritic cells")

# arcsinh cofactor — 150 is a standard default for fluorescence flow data.
arcsinh_cofactor <- 150

### ========================================
### ANALYSIS BLOCK — Step 1: join FACS to mRNA cells
### ========================================

## ---- 1. Pull mRNA cell-level metadata: barcode + cell type ----
mrna_cell_meta <- data.frame(
  XC       = colnames(pbmc.so),
  CellType = as.character(pbmc.so@active.ident),
  stringsAsFactors = FALSE
)

cat(sprintf("\nmRNA-annotated cells (pbmc.so): %d\n", nrow(mrna_cell_meta)))

## ---- 2. Subset FACS-bearing rows from the full barcode metadata table ----
facs_raw <- barcode_donor_meta_labels_full %>%
  dplyr::select(XC, XC_NovaSeq, dplyr::all_of(facs_markers))

facs_complete_mask <- complete.cases(facs_raw[, facs_markers, drop = FALSE])
facs_raw <- facs_raw[facs_complete_mask, , drop = FALSE]

cat(sprintf("Barcodes with complete FACS data (in metadata file): %d\n",
            nrow(facs_raw)))

## ---- 3. Join FACS to mRNA cells on XC ----
mrna_facs   <- dplyr::inner_join(mrna_cell_meta, facs_raw, by = "XC")
n_facs_match <- nrow(mrna_facs)

cat(sprintf("mRNA cells with paired FACS: %d / %d (%.1f%%)\n",
            n_facs_match, nrow(mrna_cell_meta),
            100 * n_facs_match / nrow(mrna_cell_meta)))

# Per-cell-type breakdown of FACS coverage
coverage_by_celltype <- mrna_cell_meta %>%
  dplyr::mutate(has_facs = XC %in% mrna_facs$XC) %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(
    n_total     = dplyr::n(),
    n_with_facs = sum(has_facs),
    pct_facs    = 100 * n_with_facs / n_total,
    .groups = "drop"
  )

cat("\n=== FACS coverage per cell type ===\n")
print(coverage_by_celltype)

## ---- 4. Reshape to long format for downstream thresholding & plotting ----
mrna_facs_long <- mrna_facs %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(facs_markers),
    names_to  = "marker_channel",
    values_to = "intensity_raw"
  ) %>%
  dplyr::mutate(
    marker_label = facs_marker_labels[marker_channel]
  )

cat(sprintf("\nLong-format FACS table: %d rows (%d cells × %d markers)\n",
            nrow(mrna_facs_long), n_facs_match, length(facs_markers)))

# Sanity check: every cell × marker combination present exactly once
stopifnot(
  nrow(mrna_facs_long) == n_facs_match * length(facs_markers),
  !anyNA(mrna_facs_long$intensity_raw)
)

### ========================================
### ANALYSIS BLOCK — Step 2: transform, threshold, summarize
### ========================================

## ---- 1. arcsinh transform ----
mrna_facs_long <- mrna_facs_long %>%
  dplyr::mutate(
    intensity_asinh = asinh(intensity_raw / arcsinh_cofactor)
  )

## ---- 2. Per-marker z-score across cells ----
mrna_facs_long <- mrna_facs_long %>%
  dplyr::group_by(marker_channel) %>%
  dplyr::mutate(
    intensity_z = as.numeric(scale(intensity_asinh))
  ) %>%
  dplyr::ungroup()

## ---- 3. Define "positive" cells for the % detected dot size ----
marker_thresholds <- mrna_facs_long %>%
  dplyr::group_by(marker_channel) %>%
  dplyr::summarise(
    threshold = median(intensity_asinh) + sd(intensity_asinh),
    .groups   = "drop"
  )

mrna_facs_long <- mrna_facs_long %>%
  dplyr::left_join(marker_thresholds, by = "marker_channel") %>%
  dplyr::mutate(is_positive = intensity_asinh >= threshold)

## ---- 4. Per-cell-type × marker summary for dotplot ----
facs_dotplot_df <- mrna_facs_long %>%
  dplyr::group_by(CellType, marker_channel, marker_label) %>%
  dplyr::summarise(
    mean_z       = mean(intensity_z, na.rm = TRUE),
    pct_positive = mean(is_positive, na.rm = TRUE),
    n_cells      = dplyr::n(),
    .groups      = "drop"
  )

## ---- 5. Order markers by cell-type specificity ----
marker_assignment_facs <- facs_dotplot_df %>%
  dplyr::group_by(marker_label) %>%
  dplyr::slice_max(mean_z, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(marker_label,
                assigned_cluster = CellType,
                top_z = mean_z)

marker_order_facs <- marker_assignment_facs %>%
  dplyr::mutate(assigned_cluster = factor(assigned_cluster,
                                          levels = celltype_order_facs)) %>%
  dplyr::arrange(assigned_cluster, dplyr::desc(top_z)) %>%
  dplyr::pull(marker_label)

cat("\nSF1 Panel — FACS marker x-axis order:\n")
print(marker_order_facs)

facs_dotplot_df <- facs_dotplot_df %>%
  dplyr::mutate(
    marker_label = factor(marker_label, levels = marker_order_facs),
    CellType     = factor(CellType,     levels = rev(celltype_order_facs)),
    mean_z_clip  = pmax(pmin(mean_z, 2), -2)
  )

### ========================================
### PLOTTING BLOCK
### ========================================

C_facs_dotplot <- ggplot(facs_dotplot_df,
                         aes(x = marker_label, y = CellType)) +
  geom_point(aes(size = pct_positive, fill = mean_z_clip),
             shape = 21, color = "grey30", stroke = 0.2) +
  scale_size_continuous(
    name   = "% Positive",
    range  = c(0, 3),
    limits = c(0, 1),
    breaks = c(0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_fill_gradient2(
    name     = "Mean z-score",
    low      = "darkblue",
    mid      = "white",
    high     = "#CE2029",
    midpoint = 0,
    limits   = c(-2, 2)
  ) +
  facet_wrap(~ "FACS (paired RNA cells)") +
  labs(x = NULL, y = NULL) +
  theme_light(base_size = 6.5) +
  theme(
    axis.text.x          = element_text(size = 5.5, angle = 45, hjust = 1),
    axis.text.y          = element_blank(),
    axis.ticks.y         = element_blank(),
    strip.background     = element_rect(fill = "grey60", color = NA),
    strip.text           = element_text(color = "white", face = "bold", size = 6.5),
    legend.direction     = "horizontal",
    legend.text          = element_text(size = 5),
    legend.title         = element_text(size = 5.5),
    legend.key.size      = unit(2, "mm"),
    legend.margin        = margin(0, 0, 0, 0),
    legend.box.spacing   = unit(1, "mm"),
    panel.grid.minor     = element_blank(),
    plot.margin          = margin(2, 2, 2, 2)
  ) +
  guides(
    size = guide_legend(position = "top",                                    # CHANGED
                        override.aes = list(fill = "grey70")),
    fill = guide_colorbar(position = "bottom",                                # CHANGED
                          barwidth = unit(20, "mm"),
                          barheight = unit(2, "mm"))
  )

C_facs_dotplot



