########################################################################################
## Integration of protein data: single-cells
########################################################################################


# ── 1. Load per-batch identification tables ───────────────────────────────────────────

base_id <- "/projects/slavov/LK/ImmPTR/results/final_scripts_output/"

ID_df_n2p1 <- read.csv(paste0(base_id, "IdentificationsList_n2p1.csv"))
ID_df_n2p2 <- read.csv(paste0(base_id, "IdentificationsList_n2p2.csv"))
ID_df_n3p1 <- read.csv(paste0(base_id, "IdentificationsList_n3p1.csv"))
ID_df_n6p1 <- read.csv(paste0(base_id, "IdentificationsList_n6p1.csv"))
ID_df_n6p2 <- read.csv(paste0(base_id, "IdentificationsList_n6p2.csv"))

# Append batch suffix so File.Name is globally unique across batches
ID_df_n2p1$File.Name <- paste0(ID_df_n2p1$File.Name, "_n2p1")
ID_df_n2p2$File.Name <- paste0(ID_df_n2p2$File.Name, "_n2p2")
ID_df_n3p1$File.Name <- paste0(ID_df_n3p1$File.Name, "_n3p1")
ID_df_n6p1$File.Name <- paste0(ID_df_n6p1$File.Name, "_n6p1")
ID_df_n6p2$File.Name <- paste0(ID_df_n6p2$File.Name, "_n6p2")

ID_df <- rbind(ID_df_n2p1, ID_df_n2p2, ID_df_n3p1, ID_df_n6p1, ID_df_n6p2)


# ── 2. Load per-batch quantification tables and build protein matrix ──────────────────

base_prot <- "/projects/slavov/LK/ImmPTR/results/final_scripts_output/"

prot_df_n2p1 <- read.csv(paste0(base_prot, "Genes_n2p1_for_prot_integration.csv"))
prot_df_n2p2 <- read.csv(paste0(base_prot, "Genes_n2p2_for_prot_integration.csv"))
prot_df_n3p1 <- read.csv(paste0(base_prot, "Genes_n3p1_for_prot_integration.csv"))
prot_df_n6p1 <- read.csv(paste0(base_prot, "Genes_n6p1_for_prot_integration.csv"))
prot_df_n6p2 <- read.csv(paste0(base_prot, "Genes_n6p2_for_prot_integration.csv"))

prot_df <- rbind(prot_df_n2p1, prot_df_n2p2, prot_df_n3p1, prot_df_n6p1, prot_df_n6p2)

# Harmonise donor labels
prot_df$sample[prot_df$sample == "Donor_04"] <- "Donor04"
prot_df$sample[prot_df$sample == "Donor_05"] <- "Donor05"

# Flag batches with corrected run order
prot_df <- prot_df %>%
  mutate(runorder_correct = if_else(dataset_batch %in% c("n2p1", "n6p1"), "True", "False"))

# Wide protein × cell matrix
prot_mat <- prot_df %>%
  dplyr::select(Genes, ID, Prot_FC_log2) %>%
  pivot_wider(names_from = ID, values_from = Prot_FC_log2) %>%
  as.data.frame()
rownames(prot_mat) <- prot_mat$Genes
prot_mat$Genes <- NULL
prot_mat <- as.matrix(prot_mat)

# Initial normalization (Normalize_saad_log defined in functions.R)
prot_mat <- Normalize_saad_log(prot_mat)


base_prot_ab <- "/projects/slavov/LK/ImmPTR/results/final_scripts_output/"

prot_df_ab_n2p1 <- read.csv(paste0(base_prot_ab, "Genes_n2p1_for_abundance_analysis.csv"))
prot_df_ab_n2p2 <- read.csv(paste0(base_prot_ab, "Genes_n2p2_for_abundance_analysis.csv"))
prot_df_ab_n3p1 <- read.csv(paste0(base_prot_ab, "Genes_n3p1_for_abundance_analysis.csv"))
prot_df_ab_n6p1 <- read.csv(paste0(base_prot_ab, "Genes_n6p1_for_abundance_analysis.csv"))
prot_df_ab_n6p2 <- read.csv(paste0(base_prot_ab, "Genes_n6p2_for_abundance_analysis.csv"))

prot_df_ab <- rbind(prot_df_ab_n2p1, prot_df_ab_n2p2, prot_df_ab_n3p1,
                    prot_df_ab_n6p1, prot_df_ab_n6p2)

prot_df_ab$sample[prot_df_ab$sample == "Donor_04"] <- "Donor04"
prot_df_ab$sample[prot_df_ab$sample == "Donor_05"] <- "Donor05"

prot_df_ab <- prot_df_ab %>%
  mutate(runorder_correct = if_else(dataset_batch %in% c("n2p1", "n6p1"), "True", "False"))

prot_mat_ab <- prot_df_ab %>%
  dplyr::select(Genes, ID, Prot_FC_log2) %>%
  pivot_wider(names_from = ID, values_from = Prot_FC_log2) %>%
  as.data.frame()
rownames(prot_mat_ab) <- prot_mat_ab$Genes
prot_mat_ab$Genes <- NULL
prot_mat_ab <- as.matrix(prot_mat_ab)

# Column-only normalization — no row centering at any step
prot_mat_ab <- Normalize_saad_log_column_only(prot_mat_ab)


# Batch label lookup used by ComBat and Seurat metadata
prot_df_dataset_labels <- prot_df %>%
  dplyr::select(ID, dataset_batch, sample, runorder_correct) %>%
  distinct()

batch_label <- as.data.frame(colnames(prot_mat))
colnames(batch_label)[1] <- "ID"
batch_label <- batch_label %>%
  left_join(prot_df_dataset_labels, by = "ID") %>%
  dplyr::select(ID, sample, dataset_batch) %>%
  distinct()


# ── 3. kNN imputation → ComBat batch correction → re-normalization ────────────────────

# Absolute and relative abundances

sc.imp.ab <- hknn(prot_mat_ab, 3)
sc.imp.ab <- Normalize_saad_log_column_only(sc.imp.ab)

sc.batch_cor.ab <- ComBat(sc.imp.ab, batch = factor(batch_label$dataset_batch))
sc.batch_cor.ab <- Normalize_saad_log_column_only(sc.batch_cor.ab)

sc.batch_cor.ab.imp <- sc.batch_cor.ab

# Restore original missing values
sc.batch_cor.ab[is.na(prot_mat_ab)] <- NA

# hknn defined in functions.R
sc.imp <- hknn(prot_mat, 3)
sc.imp <- Normalize_saad_log(sc.imp)

sc.batch_cor <- ComBat(sc.imp, batch = factor(batch_label$dataset_batch))
sc.batch_cor <- Normalize_saad_log(sc.batch_cor)

# Fully imputed version retained for ICA input (used in section 9)
sc.batch_cor.imp <- sc.batch_cor

# Restore original missing values in the primary analysis matrices
sc.batch_cor[is.na(prot_mat)] <- NA


# Save sc.batch_cor.ab (column-normalized absolute abundance)
fwrite(data.table(Genes = rownames(sc.batch_cor.ab), sc.batch_cor.ab),
       file = "/projects/slavov/LK/ImmPTR/data/final_data/Protein_x_Cell_Matrix_Column_Normalized.csv")

# Save sc.batch_cor (row-centered for covariation)
fwrite(data.table(Genes = rownames(sc.batch_cor), sc.batch_cor),
       file = "/projects/slavov/LK/ImmPTR/data/final_data/Protein_x_Cell_Matrix_Normalized_Centered.csv")

# Save sc.batch_cor.ab.imp (imputed, column-normalized absolute abundance)
fwrite(data.table(Genes = rownames(sc.batch_cor.ab.imp), sc.batch_cor.ab.imp),
       file = "/projects/slavov/LK/ImmPTR/data/final_data/Protein_x_Cell_Matrix_Imputed_Column_Normalized.csv")

# Save sc.batch_cor.imp (imputed, row-centered for covariation)
fwrite(data.table(Genes = rownames(sc.batch_cor.imp), sc.batch_cor.imp),
       file = "/projects/slavov/LK/ImmPTR/data/final_data/Protein_x_Cell_Matrix_Imputed_Normalized_Centered.csv")


###### Alternatively, load protein matrices directly:

# Protein matrices (NA-restored, gene x cell)
sc.batch_cor <- read.csv(
  "/projects/slavov/LK/ImmPTR/data/final_data/Protein_x_Cell_Matrix_Normalized_Centered.csv"
)
rownames(sc.batch_cor) <- sc.batch_cor$Genes
sc.batch_cor$Genes     <- NULL
sc.batch_cor          <- as.matrix(sc.batch_cor)

sc.batch_cor.ab <- read.csv(
  "/projects/slavov/LK/ImmPTR/data/final_data/Protein_x_Cell_Matrix_Column_Normalized.csv"
)
rownames(sc.batch_cor.ab) <- sc.batch_cor.ab$Genes
sc.batch_cor.ab$Genes     <- NULL
sc.batch_cor.ab          <- as.matrix(sc.batch_cor.ab)



# ── 4. Seurat object construction and dimensionality reduction ────────────────────────

# Median-impute any remaining NAs row-wise for Seurat input
prot_mat_imputed <- apply(sc.batch_cor, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
})
prot_mat_imputed <- t(prot_mat_imputed)   # apply() transposes; restore orientation


pbmc.prot <- CreateSeuratObject(
  counts       = prot_mat_imputed,
  assay        = "Proteomics",
  min.cells    = 10,
  min.features = 10
)
LayerData(pbmc.prot[["Proteomics"]], "data") <- LayerData(pbmc.prot[["Proteomics"]], "counts")

# Attach batch / donor metadata
meta_to_add <- prot_df_dataset_labels %>%
  filter(ID %in% colnames(pbmc.prot)) %>%
  column_to_rownames("ID")
pbmc.prot <- AddMetaData(pbmc.prot, metadata = meta_to_add)

# Variable features, scale (regressing donor), PCA, clustering, UMAP
vars <- apply(prot_mat_imputed, 1, var, na.rm = TRUE)
top_features <- names(sort(vars, decreasing = TRUE))[1:2000]
VariableFeatures(pbmc.prot) <- top_features
pbmc.prot <- ScaleData(pbmc.prot, vars.to.regress = "sample")
pbmc.prot <- RunPCA(pbmc.prot, features = VariableFeatures(pbmc.prot))
ElbowPlot(pbmc.prot, ndims = 35)

pbmc.prot <- FindNeighbors(pbmc.prot, dims = 1:4)
pbmc.prot <- FindClusters(pbmc.prot, resolution = 0.3, random.seed = 0)
pbmc.prot <- RunUMAP(pbmc.prot, dims = 1:4, seed.use = 42)

# Diagnostic DimPlots — batch, donor, run-order
DimPlot(pbmc.prot, reduction = "umap", label = TRUE, pt.size = 1, alpha = 0.2)
DimPlot(pbmc.prot, reduction = "pca", label = TRUE, pt.size = 1, alpha = 0.4)
DimPlot(pbmc.prot, reduction = "pca",  group.by = "sample",           label = FALSE, pt.size = 1, alpha = 0.4)
DimPlot(pbmc.prot, reduction = "pca",  group.by = "dataset_batch",    label = FALSE, pt.size = 1, alpha = 0.4)
DimPlot(pbmc.prot, reduction = "pca",  group.by = "runorder_correct", label = FALSE, pt.size = 1, alpha = 0.4)
DimPlot(pbmc.prot, reduction = "umap", group.by = "sample",           label = FALSE, pt.size = 1, alpha = 0.2)
DimPlot(pbmc.prot, reduction = "umap", group.by = "dataset_batch",    label = FALSE, pt.size = 1, alpha = 0.2)
DimPlot(pbmc.prot, reduction = "umap", group.by = "runorder_correct", label = FALSE, pt.size = 1, alpha = 0.2)


# ── 5. NA-aware differential expression for cell-type marker discovery ────────────────
############ SF1 Panel C — Canonical marker dotplots (protein side)
############ NA-aware DE + custom ggplot for style consistency

### ========================================
### CONFIGURATION
### ========================================

# Canonical markers used for cell-type annotation
canonical_markers <- c(
  "CD3D", "CD3E", "CD247", "ZAP70",   # T cell 
  "CD8A", "GZMA",             # CD8 T
  "GNLY",                      # NK
  "MS4A1", "CD74", "CD37",    # B cell
  "CD14", "ITGAM", "FCER1G"   # Monocyte / myeloid
)

# Cell-type row ordering for the dotplot (y-axis)
celltype_order <- c("Monocytes", "NK cells", "CD8 T cells", "CD4 T cells", "B cells")

### ========================================
### ANALYSIS BLOCK — NA-aware DE (vectorized)
### ========================================

run_na_aware_de <- function(expr_mat, clusters) {
  cluster_levels <- levels(clusters)
  results <- vector("list", length(cluster_levels))
  names(results) <- cluster_levels
  
  for (cl in cluster_levels) {
    message("DE for cluster: ", cl)
    in_clust  <- which(clusters == cl)
    out_clust <- which(clusters != cl)
    
    v1_mat <- expr_mat[, in_clust,  drop = FALSE]
    v2_mat <- expr_mat[, out_clust, drop = FALSE]
    
    # Per-gene summary statistics, fully vectorized
    n1    <- rowSums(!is.na(v1_mat))
    n2    <- rowSums(!is.na(v2_mat))
    m1    <- rowMeans(v1_mat, na.rm = TRUE)
    m2    <- rowMeans(v2_mat, na.rm = TRUE)
    # rowVars via rowSums on squared residuals (avoids matrixStats dep)
    s1_sq <- rowSums((v1_mat - m1)^2, na.rm = TRUE) / pmax(n1 - 1, 1)
    s2_sq <- rowSums((v2_mat - m2)^2, na.rm = TRUE) / pmax(n2 - 1, 1)
    
    # Welch's t and df
    se     <- sqrt(s1_sq / n1 + s2_sq / n2)
    t_stat <- (m1 - m2) / se
    df     <- (s1_sq / n1 + s2_sq / n2)^2 /
      ((s1_sq / n1)^2 / pmax(n1 - 1, 1) +
         (s2_sq / n2)^2 / pmax(n2 - 1, 1))
    
    # Two-sided p-value; NA where either group has <3 valid
    valid <- n1 >= 3 & n2 >= 3 & is.finite(se) & se > 0
    p_val <- rep(NA_real_, length(t_stat))
    p_val[valid] <- 2 * pt(-abs(t_stat[valid]), df = df[valid])
    
    results[[cl]] <- data.frame(
      gene       = rownames(expr_mat),
      cluster    = cl,
      avg_log2FC = m1 - m2,
      p_val      = p_val,
      pct.1      = n1 / ncol(v1_mat),
      pct.2      = n2 / ncol(v2_mat),
      stringsAsFactors = FALSE
    )
  }
  
  markers <- dplyr::bind_rows(results)
  markers$p_val_adj <- p.adjust(markers$p_val, method = "BH")
  markers
}

## ---- Run DE and store ----
expr_mat <- sc.batch_cor   # NA-preserved matrix used for DE
clusters <- Idents(pbmc.prot)

prot.markers <- run_na_aware_de(expr_mat, clusters)
significant_prot_markers <- prot.markers %>%
  dplyr::filter(!is.na(p_val_adj) & p_val_adj < 0.05)

pbmc.prot@misc$prot_markers     <- prot.markers
pbmc.prot@misc$prot_markers_sig <- significant_prot_markers

cat(sprintf("SF1 Panel C (protein) — total markers: %d; significant (FDR<0.05): %d\n",
            nrow(prot.markers), nrow(significant_prot_markers)))


GeneMarkersUsed_numbered <- significant_prot_markers %>%
  dplyr::filter(grepl(
    "CD4|CD3E|CD3D|IL7R|CCR7|SELL|GATA3|CD8A|CD8B|GNLY|GZMB|S1PR5|NKG7|MS4A1|CD79A|FCGR3A|MS4A7|ITGAX|CD14|S100A8|S100A9|IL1B|CCR1|FCER1A|CD1C|SPIB|JCHAIN",
    gene))

# Inspect DotPlot 
DotPlot(
  pbmc.prot,
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



# ── 6. Cell-type annotation ───────────────────────────────────────────────────────────

# Build the mapping using the numeric levels
new.cluster.ids <- c("Monocytes", "CD4 T cells", "CD8 T cells", "B cells", "NK cells")
names(new.cluster.ids) <- levels(pbmc.prot)   # while still "0".."4"

# Relabel prot.markers using the numeric-named mapping
prot.markers$cluster <- new.cluster.ids[as.character(prot.markers$cluster)]

# Now rename the Seurat object's idents
pbmc.prot <- RenameIdents(pbmc.prot, new.cluster.ids)

# Refresh cached copies
significant_prot_markers <- prot.markers %>%
  dplyr::filter(!is.na(p_val_adj) & p_val_adj < 0.05)
pbmc.prot@misc$prot_markers     <- prot.markers
pbmc.prot@misc$prot_markers_sig <- significant_prot_markers

celltype_colors <- c(
  "CD4 T cells" = "#FFB347",
  "CD8 T cells" = "#5A8DEE",
  "B cells"     = "#9B3D6A",
  "NK cells"    = "#008B72",
  "Monocytes"   = "#D2BD92"
)

celltype_order <- c("Monocytes", "NK cells", "CD8 T cells", "CD4 T cells", "B cells")

# Shared cell-type lookup used by ID depth plots and ICA below
cell_types_df <- pbmc.prot@active.ident %>%
  as.data.frame() %>%
  rownames_to_column("Cell") %>%
  dplyr::rename(CellType = ".")


### ========================================
### Subset + order markers by cell-type specificity
### ========================================

# For each canonical marker, find the cell type where it has max avg_log2FC.
marker_assignment <- prot.markers %>%
  dplyr::filter(gene %in% canonical_markers,
                cluster %in% celltype_order) %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, assigned_cluster = cluster, top_log2FC = avg_log2FC)

# Order markers: by assigned cell type (in celltype_order), then by
# decreasing avg_log2FC within each assigned cell type block.
marker_order_prot <- marker_assignment %>%
  dplyr::mutate(assigned_cluster = factor(assigned_cluster, levels = celltype_order)) %>%
  dplyr::arrange(assigned_cluster, dplyr::desc(top_log2FC)) %>%
  dplyr::pull(gene)

cat("SF1 Panel C — protein marker x-axis order:\n")
print(marker_order_prot)

dotplot_prot_df <- prot.markers %>%
  dplyr::filter(gene %in% marker_order_prot,
                cluster %in% celltype_order) %>%
  dplyr::mutate(
    gene    = factor(gene,    levels = marker_order_prot),
    cluster = factor(cluster, levels = rev(celltype_order)),
    avg_log2FC_clip = pmax(pmin(avg_log2FC, 2), -2)
  )

### ========================================
### PLOTTING BLOCK — custom ggplot dotplot
### ========================================

C_prot_dotplot <- ggplot(dotplot_prot_df,
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
  facet_wrap(~ "Protein") +
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
    size = guide_legend(position = "top",                                    
                        override.aes = list(fill = "grey70")),
    fill = guide_colorbar(position = "bottom",                                
                          barwidth = unit(20, "mm"),
                          barheight = unit(2, "mm"))
  )

C_prot_dotplot


# ── 7. QC: Identification depth per cell ─────────────────────────────────────────────

# Attach annotated cell types to the combined identification table
ID_df <- ID_df %>% left_join(cell_types_df, by = c("File.Name" = "Cell"))

## 7a. Per-cell summary across all cells (cell-type agnostic) — 4-panel sensitivity plot
sum_IDs_all_cells <- ID_df %>%
  group_by(File.Name) %>%
  summarise(
    Num_Precursors    = n_distinct(seqcharge),
    Num_Peptides      = n_distinct(Stripped.Sequence),
    Num_ProteinGroups = n_distinct(Protein.Group),
    Num_Proteins      = n_distinct(Genes),
    .groups = "drop"
  )

## 7b. Per-cell summary stratified by cell type — cell-type violin
sum_IDs_by_celltype <- ID_df %>%
  group_by(CellType, File.Name) %>%
  summarise(
    Num_Precursors    = n_distinct(seqcharge),
    Num_Peptides      = n_distinct(Stripped.Sequence),
    Num_ProteinGroups = n_distinct(Protein.Group),
    Num_Proteins      = n_distinct(Genes),
    .groups = "drop"
  )
sum_IDs_by_celltype$CellType <- factor(sum_IDs_by_celltype$CellType, levels = celltype_order)

## 7c. Proteins per cell by cell type (violin + boxplot)
SCP_IDs_Across_CellTypes <- ggplot(
  sum_IDs_by_celltype,
  aes(x = CellType, y = Num_Proteins, fill = CellType)
) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", alpha = 0.8) +
  scale_fill_manual(values = celltype_colors) +
  ylim(0, 2000) +
  labs(x = "", y = "# Proteins") +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 16, angle = 45, hjust = 1),
    axis.title.y    = element_text(size = 22),
    axis.text.y     = element_text(size = 16)
  )

SCP_IDs_Across_CellTypes

## 7d. Detection sensitivity: Precursors / Peptides / Proteins across ALL cells (1 × 3)

violin_theme <- theme_light(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.title.x    = element_blank(),
    axis.title.y    = element_text(size = 8),
    axis.text.y     = element_text(size = 7),
    plot.margin     = margin(2, 2, 2, 2)
  )

p_precursors <- ggplot(sum_IDs_all_cells, aes(x = "", y = Num_Precursors)) +
  geom_violin(fill = "grey40", trim = FALSE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black",
               fill = "white", alpha = 0.8, linewidth = 0.3) +
  coord_cartesian(ylim = c(5000, 20000)) +    # explicit, keeps existing range
  labs(y = "# Precursor Ions / Cell") +
  violin_theme

p_peptides <- ggplot(sum_IDs_all_cells, aes(x = "", y = Num_Peptides)) +
  geom_violin(fill = "grey40", trim = FALSE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black",
               fill = "white", alpha = 0.8, linewidth = 0.3) +
  coord_cartesian(ylim = c(2000, 10000)) +    # CHANGED
  labs(y = "# Peptides / Cell") +
  violin_theme

p_proteins <- ggplot(sum_IDs_all_cells, aes(x = "", y = Num_Proteins)) +
  geom_violin(fill = "#D55E00", trim = FALSE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black",
               fill = "white", alpha = 0.8, linewidth = 0.3) +
  coord_cartesian(ylim = c(1000, 2000)) +     # CHANGED
  labs(y = "# Proteins / Cell") +
  violin_theme

SCP_ID_Depth_panel <- p_precursors | p_peptides | p_proteins

SCP_ID_Depth_panel

## 7e. Detection sensitivity per cell type: Precursors / Peptides / Proteins (3 × 1)

# Shared theme for the stacked per-cell-type plots
ct_violin_theme_top <- theme_light(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.title.x    = element_blank()
  )

ct_violin_theme_bottom <- theme_light(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 13, angle = 45, hjust = 1),
    axis.title.x    = element_blank()
  )

ct_precursors <- ggplot(
  sum_IDs_by_celltype,
  aes(x = CellType, y = Num_Precursors, fill = CellType)
) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", alpha = 0.8) +
  scale_fill_manual(values = celltype_colors) +
  labs(y = "# Precursor Ions / Cell") +
  ct_violin_theme_top

ct_peptides <- ggplot(
  sum_IDs_by_celltype,
  aes(x = CellType, y = Num_Peptides, fill = CellType)
) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", alpha = 0.8) +
  scale_fill_manual(values = celltype_colors) +
  labs(y = "# Peptides / Cell") +
  ct_violin_theme_top

ct_proteins <- ggplot(
  sum_IDs_by_celltype,
  aes(x = CellType, y = Num_Proteins, fill = CellType)
) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", alpha = 0.8) +
  scale_fill_manual(values = celltype_colors) +
  labs(y = "# Proteins / Cell") +
  ct_violin_theme_bottom

SCP_ID_Depth_CellType_panel <- ct_precursors / ct_peptides / ct_proteins

SCP_ID_Depth_CellType_panel


# ── 8. UMAP visualization (annotated cell types) ─────────────────────────────────────

umap_df           <- as.data.frame(Embeddings(pbmc.prot, "umap"))
umap_df$celltype  <- pbmc.prot@active.ident

centroids <- umap_df %>%
  group_by(celltype) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

scp.umap <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 2.5, alpha = 0.45) +
  scale_color_manual(values = celltype_colors) +
  geom_text_repel(
    data         = centroids,
    aes(label    = celltype),
    size         = 4,
    color        = "black",
    box.padding  = 0.1,
    max.overlaps = Inf
  ) +
  labs(title = "Protein", x = "UMAP 1", y = "UMAP 2") +
  theme_light(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "none"
  )

scp.umap


# ── 9. ICA visualization (fully imputed input) ───────────────────────────────────────

sc.batch_cor.imp <- read.csv(
  "/projects/slavov/LK/ImmPTR/data/final_data/Protein_x_Cell_Matrix_Imputed_Normalized_Centered.csv"
)
rownames(sc.batch_cor.imp) <- sc.batch_cor.imp$Genes
sc.batch_cor.imp$Genes     <- NULL
sc.batch_cor.imp          <- as.matrix(sc.batch_cor.imp)

set.seed(123)
ica_res <- fastICA(t(sc.batch_cor.imp), n.comp = 3, alg.typ = "deflation")

ica_df     <- as.data.frame(ica_res$S)   # cells × independent components
ica_df$ID  <- rownames(ica_df)

plot_df <- cell_types_df %>%
  inner_join(ica_df, by = c("Cell" = "ID"))

ICA_kNN_BatchCor_Norm <- ggplot(plot_df, aes(x = V1, y = V2, color = CellType)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = celltype_colors) +
  facet_wrap(~ "Protein") +
  labs(x = "Independent Component 1",
       y = "Independent Component 2") +
  theme_light(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 7),
    legend.position = "none",
    strip.background = element_rect(fill = "grey60", color = NA),
    strip.text = element_text(color = "white", face = "bold", size = 8),
    plot.margin = margin(2, 2, 2, 2)
  )

ICA_kNN_BatchCor_Norm


# ── 10. Save outputs ───────────────────

saveRDS(pbmc.prot,
        file = "/projects/slavov/LK/ImmPTR/data/final_data/SCP_PBMC.rds")


# ── 11. ICA visualization (complete proteins only, no imputation) ─────────────────────

complete_proteins <- rownames(sc.batch_cor)[
  rowSums(is.na(sc.batch_cor)) == 0
]
sc.complete <- sc.batch_cor[complete_proteins, ]

set.seed(123)
ica_res <- fastICA(t(sc.complete), n.comp = 3, alg.typ = "deflation")

ica_df     <- as.data.frame(ica_res$S)
ica_df$ID  <- rownames(ica_df)

plot_df <- cell_types_df %>%
  dplyr::inner_join(ica_df, by = c("Cell" = "ID"))

ICA_Complete_BatchCor <- ggplot(plot_df, aes(x = V1, y = V2, color = CellType)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = celltype_colors) +
  facet_wrap(~ "Protein") +
  annotate("text",
           x = -Inf, y = Inf,
           label = paste0("n = ", length(complete_proteins), " Proteins"),
           hjust = -0.1, vjust = 1.4,
           size = 1.8, family = "sans") +
  labs(x = "Independent Component 1",
       y = "Independent Component 2") +
  theme_light(base_size = 7) +
  theme(
    aspect.ratio         = 1,
    axis.title           = element_text(size = 6.5),
    axis.text            = element_text(size = 5.5),
    legend.position      = "bottom",                                    
    legend.justification = "center",                                    
    legend.direction     = "horizontal",                                
    legend.title         = element_blank(),                             
    legend.text          = element_text(size = 4.5),                     
    legend.key.size      = unit(1.8, "mm"),                             
    legend.background    = element_rect(fill = NA, color = NA),         
    legend.margin        = margin(0, 0, 0, 0),
    legend.box.spacing   = unit(1, "mm"),                               
    strip.background     = element_rect(fill = "grey60", color = NA),
    strip.text           = element_text(color = "white", face = "bold", size = 6.5),
    plot.margin          = margin(2, 2, 2, 2)
  ) +
  guides(color = guide_legend(nrow = 2,
                              override.aes = list(size = 2, alpha = 1)))

ICA_Complete_BatchCor
cat("SF1 Panel A — Proteins at 100% completeness:", length(complete_proteins), "\n")




########################################################################################
## Integration of peptide data: single-cells (parallel to protein integration)
########################################################################################

# ── P1. Load per-batch peptide quantification + mapping tables ────────────────────────

base_pep <- "/projects/slavov/LK/ImmPTR/results/final_scripts_output/"

pep_df_n2p1 <- read.csv(paste0(base_pep, "Peptides_n2p1_for_integration.csv"))
pep_df_n2p2 <- read.csv(paste0(base_pep, "Peptides_n2p2_for_integration.csv"))
pep_df_n3p1 <- read.csv(paste0(base_pep, "Peptides_n3p1_for_integration.csv"))
pep_df_n6p1 <- read.csv(paste0(base_pep, "Peptides_n6p1_for_integration.csv"))
pep_df_n6p2 <- read.csv(paste0(base_pep, "Peptides_n6p2_for_integration.csv"))

pep_df <- rbind(pep_df_n2p1, pep_df_n2p2, pep_df_n3p1, pep_df_n6p1, pep_df_n6p2)
pep_df$sample[pep_df$sample == "Donor_04"] <- "Donor04"
pep_df$sample[pep_df$sample == "Donor_05"] <- "Donor05"
pep_df <- pep_df %>%
  mutate(runorder_correct = if_else(dataset_batch %in% c("n2p1", "n6p1"), "True", "False"))

pep_df_ab_n2p1 <- read.csv(paste0(base_pep, "Peptides_n2p1_for_abundance_analysis.csv"))
pep_df_ab_n2p2 <- read.csv(paste0(base_pep, "Peptides_n2p2_for_abundance_analysis.csv"))
pep_df_ab_n3p1 <- read.csv(paste0(base_pep, "Peptides_n3p1_for_abundance_analysis.csv"))
pep_df_ab_n6p1 <- read.csv(paste0(base_pep, "Peptides_n6p1_for_abundance_analysis.csv"))
pep_df_ab_n6p2 <- read.csv(paste0(base_pep, "Peptides_n6p2_for_abundance_analysis.csv"))

pep_df_ab <- rbind(pep_df_ab_n2p1, pep_df_ab_n2p2, pep_df_ab_n3p1, pep_df_ab_n6p1, pep_df_ab_n6p2)
pep_df_ab$sample[pep_df_ab$sample == "Donor_04"] <- "Donor04"
pep_df_ab$sample[pep_df_ab$sample == "Donor_05"] <- "Donor05"
pep_df_ab <- pep_df_ab %>%
  mutate(runorder_correct = if_else(dataset_batch %in% c("n2p1", "n6p1"), "True", "False"))

ppm_n2p1 <- read.csv(paste0(base_pep, "PepProtMap_n2p1.csv"))
ppm_n2p2 <- read.csv(paste0(base_pep, "PepProtMap_n2p2.csv"))
ppm_n3p1 <- read.csv(paste0(base_pep, "PepProtMap_n3p1.csv"))
ppm_n6p1 <- read.csv(paste0(base_pep, "PepProtMap_n6p1.csv"))
ppm_n6p2 <- read.csv(paste0(base_pep, "PepProtMap_n6p2.csv"))

ppm_all <- rbind(ppm_n2p1, ppm_n2p2, ppm_n3p1, ppm_n6p1, ppm_n6p2)


# ── P2. Build integrated peptide-to-protein mapping with strict proteotypic filter ────

pep_prot_map_integrated <- ppm_all %>%
  group_by(Stripped.Sequence) %>%
  summarise(
    is_proteotypic_strict = all(is_proteotypic),
    n_datasets_detected   = n_distinct(dataset_batch),
    Genes                 = paste(unique(Genes), collapse = "|"),
    Protein.Group         = paste(unique(Protein.Group), collapse = "|"),
    .groups = "drop"
  ) %>%
  mutate(genes_consistent = !grepl("\\|", Genes))

cat(sprintf("Integrated peptides: %d total, %d proteotypic (strict), %d gene-consistent\n",
            nrow(pep_prot_map_integrated),
            sum(pep_prot_map_integrated$is_proteotypic_strict),
            sum(pep_prot_map_integrated$genes_consistent)))


# ── P3. Build wide peptide × cell matrices (relative and absolute) ─────────────────────

pep_mat <- pep_df %>%
  dplyr::select(Stripped.Sequence, ID, Peptide_FC_log2) %>%
  pivot_wider(names_from = ID, values_from = Peptide_FC_log2) %>%
  as.data.frame()
rownames(pep_mat) <- pep_mat$Stripped.Sequence
pep_mat$Stripped.Sequence <- NULL
pep_mat <- as.matrix(pep_mat)
pep_mat <- Normalize_saad_log(pep_mat)

pep_mat_ab <- pep_df_ab %>%
  dplyr::select(Stripped.Sequence, ID, Peptide_FC_log2) %>%
  pivot_wider(names_from = ID, values_from = Peptide_FC_log2) %>%
  as.data.frame()
rownames(pep_mat_ab) <- pep_mat_ab$Stripped.Sequence
pep_mat_ab$Stripped.Sequence <- NULL
pep_mat_ab <- as.matrix(pep_mat_ab)
pep_mat_ab <- Normalize_saad_log_column_only(pep_mat_ab)


# ── P4. Cross-dataset ComBat integration on peptide matrices ──────────────────────────

pep_df_dataset_labels <- pep_df %>%
  dplyr::select(ID, dataset_batch, sample, runorder_correct) %>%
  distinct()

batch_label_pep <- as.data.frame(colnames(pep_mat))
colnames(batch_label_pep)[1] <- "ID"
batch_label_pep <- batch_label_pep %>%
  left_join(pep_df_dataset_labels, by = "ID") %>%
  dplyr::select(ID, sample, dataset_batch) %>%
  distinct()

# Relative
sc.imp.pep <- hknn(pep_mat, 3)
sc.imp.pep <- Normalize_saad_log(sc.imp.pep)
sc.batch_cor.pep <- ComBat(sc.imp.pep, batch = factor(batch_label_pep$dataset_batch))
sc.batch_cor.pep <- Normalize_saad_log(sc.batch_cor.pep)
sc.batch_cor.pep[is.na(pep_mat)] <- NA

# Absolute
sc.imp.pep.ab <- hknn(pep_mat_ab, 3)
sc.imp.pep.ab <- Normalize_saad_log_column_only(sc.imp.pep.ab)
sc.batch_cor.pep.ab <- ComBat(sc.imp.pep.ab, batch = factor(batch_label_pep$dataset_batch))
sc.batch_cor.pep.ab <- Normalize_saad_log_column_only(sc.batch_cor.pep.ab)
sc.batch_cor.pep.ab[is.na(pep_mat_ab)] <- NA

fwrite(data.table(Stripped.Sequence = rownames(sc.batch_cor.pep), sc.batch_cor.pep),
       file = "/projects/slavov/LK/ImmPTR/data/final_data/Peptide_x_Cell_Matrix_Normalized_Centered.csv")
fwrite(data.table(Stripped.Sequence = rownames(sc.batch_cor.pep.ab), sc.batch_cor.pep.ab),
       file = "/projects/slavov/LK/ImmPTR/data/final_data/Peptide_x_Cell_Matrix_Column_Normalized.csv")
fwrite(pep_prot_map_integrated,
       file = "/projects/slavov/LK/ImmPTR/data/final_data/PepProtMap_integrated.csv")



