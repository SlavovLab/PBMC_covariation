########################################################################################
## Processing raw data and quality control for proteomic dataset: n3p1
########################################################################################

stats <- diann_load("/projects/slavov/LK/ImmPTR/data/raw_data/n3p1/diann_b17/report.stats.tsv")
linker <- stats %>% dplyr::select(File.Name)
linker <- linker %>% filter(!grepl("d8only", File.Name))
linker <- linker %>%
  mutate(
    Run = gsub(".*\\\\|\\.d$", "", File.Name),
    Well = gsub(".*NewSource_R([A-Z0-9]+)_.*", "\\1", Run)
  )
linker$plate <- "1"
linker$plate <- as.numeric(linker$plate)
linker <- linker %>% dplyr::select(Run, Well, plate)
linker <- linker %>%
  mutate(Order = as.integer(str_sub(Run, -4))) %>%
  arrange(Order) %>%
  mutate(Order = rank(Order, ties.method = "first"))


### Read in the path to DIA-NN file
data_path <- '/projects/slavov/LK/ImmPTR/data/raw_data/n3p1/diann_b17/report.tsv'
link_path <- '/projects/slavov/LK/PBMC_PTR/datasets/n3p1/meta_data/meta_n3p1.csv'
## Read in cell isolation files from CellenONE 
two <- '/projects/slavov/LK/ImmPTR/data/meta_data/n3p1/don4sorted2_use_isolated.xls'
one <- '/projects/slavov/LK/ImmPTR/data/meta_data/n3p1/don5sorted1_use_isolated.xls'
all_cells <- list(Donor05 = one,
                  Donor04 = two)

# Use QQC to load raw data, generate cell isolation data, and map labels and wells 
PBMC <- DIANN_to_QQC2(data_path,link_path, plex = 2, carrier = T) ### create QQC object
PBMC <- cellXpeptide2(PBMC, TQVal = 1, chQVal = 1) ### don't use the pep/prot matrices from this, just for mapping sorting info in next line
PBMC <- link_cellenONE_Raw(PBMC,all_cells)
PlotSlideLayout_celltype(PBMC)
PlotSlideLayout_label(PBMC)
rep.full <- PBMC@raw_data
meta.data <- PBMC@meta.data


rep <- rep.full
rep <- rep %>% filter(Lib.Q.Value <= 0.01, Lib.PG.Q.Value <= 0.05, Ms1.Area > 0, !grepl("MBR", Run), grepl("HUMAN", Protein.Names)) 
sel.ID.diameters <- meta.data %>% dplyr::select(diameter, ID)
rep <- rep %>% left_join(sel.ID.diameters, by = c("ID"="ID"))
rep$SampleType <- ifelse(rep$plex == 8 & is.na(rep$diameter), "carrier",
                         ifelse(is.na(rep$diameter) & rep$plex != 8, "neg",
                                ifelse(rep$plex %in% c(0, 4) & !is.na(rep$diameter), "sc", NA)))
meta.data.sample.ID <- meta.data %>% dplyr::select(ID, sample, diameter)
sample_map <- rep %>% dplyr::select(File.Name, SampleType, plex, Order, ID, Well) %>% distinct()
sample_map <- meta.data.sample.ID %>% full_join(sample_map, by=c("ID"="ID"))
sample_map <- sample_map %>%
  mutate(Cell_Volume = (4/3) * pi * (diameter / 2)^3)
seqcharge_PG_map <- rep %>% dplyr::select(seqcharge, Protein.Group) %>% distinct()
pep_PG_map <- rep %>% dplyr::select(Stripped.Sequence, seqcharge, Protein.Group) %>% distinct()


### Data filtering
CTQ_filtered_rep <- rep %>% filter(Channel.Q.Value <= 0.1, Translated.Q.Value <= 0.1)
CTQ_filtered_raw_data_stats <- CTQ_filtered_rep
CTQ_filtered_raw_data_stats <- collect_raw_statistics(CTQ_filtered_raw_data_stats, sample_map)
CTQ_filtered_raw_data_stats_NOCARRIER <- CTQ_filtered_raw_data_stats %>% filter(!grepl("carrier", SampleType))

### Determine Median Absolute Deviation thresholds for IDs and signal metrics
MAD_threshold_IDs <- calculate_MAD_threshold(CTQ_filtered_raw_data_stats_NOCARRIER, "precursor_counts", thresh = 1)
MAD_threshold_Signal <- calculate_MAD_threshold(CTQ_filtered_raw_data_stats_NOCARRIER, "total_MS1_signal", thresh = 1)
PrecursorIDThreshold <- MAD_threshold_IDs$MAD_LowerThreshold
PrecursorMS1Threshold <- MAD_threshold_Signal$MAD_LowerThreshold

### Plot cell selection criteria
filtered_IDsSignal_CellSelection <- ggplot(CTQ_filtered_raw_data_stats_NOCARRIER, aes(y = log10(total_MS1_signal), x = precursor_counts)) +
  geom_rect(aes(xmin = PrecursorIDThreshold, xmax = Inf, ymin = log10(PrecursorMS1Threshold), ymax = Inf), 
            fill = "#ffe4ed", alpha = 0.3) +
  geom_point(aes(color = ifelse(SampleType == "sc" & (precursor_counts < PrecursorIDThreshold | log10(total_MS1_signal) < log10(PrecursorMS1Threshold)), 
                                "sc removed", SampleType)), size = 4, alpha = 0.5) +
  scale_color_manual(values = c("sc" = "#009E73", "neg" = "#D55E00", "other_sample_type" = "#009E73", "sc removed" = "purple")) +
  labs(title = "Single Cell Selection by IDs and Signal (n3p1)",
       x="# Precursors", y=expression("Total MS1 Signal, Log"[10])) +
  geom_vline(xintercept = PrecursorIDThreshold) +
  geom_hline(yintercept = log10(PrecursorMS1Threshold)) +
  theme_light() +
  guides(color = guide_legend(title = "SampleType")) +  # Adjust the legend title
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20, hjust = 0.5))


### Select cells for downstream analysis by IDs and signal criteria
### Remove low quality samples, check negative controls are removed
CellsToKeep_IDsSignal <- CTQ_filtered_raw_data_stats %>%
  filter(SampleType=="sc", precursor_counts > PrecursorIDThreshold, total_MS1_signal > PrecursorMS1Threshold) %>%
  dplyr::select(ID)

### Second filter by median protein CV per single cell
### Column and row normalize precursor MS1 Areas before calculating median protein CV per cell
rep.colrownorm <- rep %>% inner_join(CellsToKeep_IDsSignal, by=c("ID"="ID"))
rep.colrownorm <- rep.colrownorm %>% dplyr::select(seqcharge, ID, Ms1.Area)
data_d04 <- rep.colrownorm %>% filter(grepl("0$|4$", ID)) %>% distinct()
data_d04 <- data_d04 %>% pivot_wider(names_from = ID, values_from = Ms1.Area) %>% as.data.frame()
rownames(data_d04) <- data_d04$seqcharge
data_d04$seqcharge <- NULL
colrownorm_data_d04 <- Normalize_saad(data_d04) %>% as.data.frame()
colrownorm_data_d04$seqcharge <- rownames(colrownorm_data_d04)
rownames(colrownorm_data_d04) <- NULL
normalized_data <- colrownorm_data_d04
normalized_data <- reshape2::melt(normalized_data)
colnames(normalized_data) <- c("seqcharge", "ID", "Ms1.Area")
normalized_data <- normalized_data %>% left_join(seqcharge_PG_map, by=c("seqcharge"="seqcharge"))

### Calculate median protein CV per cell
median_CV_df <- calculate_median_CV(normalized_data)
median_CV_df <- median_CV_df %>% left_join(sample_map, by=c("ID"="ID"))

### Determine Median Absolute Deviation thresholds for CV metric
MAD_CV_df <- median_CV_df %>% filter(SampleType == "sc")
MAD_thresholds <- calculate_MAD_threshold(MAD_CV_df, "Median_CV", thresh = 1)
CV_thresh <- MAD_thresholds$MAD_UpperThreshold

### Plot cell selection criteria
filtered_plot_CV_CellSelection <- ggplot(median_CV_df, aes(x = Median_CV, fill = SampleType)) +
  geom_rect(aes(xmin = -Inf, xmax = CV_thresh, ymin = -Inf, ymax = Inf),
            fill = "#ffe4ed", alpha = 0.3) +
  geom_histogram(position = "identity", bins = 50, alpha = 0.5) +
  geom_vline(aes(xintercept = CV_thresh), color = "grey40", size = 1) + 
  scale_fill_manual(values = c("#009E73")) +
  labs(title = "Single Cell Selection by CV (n3p1)",
       x = "Median Protein CV per Cell",
       y = "Count") +
  theme_light() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20, hjust = 0.5))

CellsToKeep_CV <- median_CV_df %>%
  filter(SampleType=="sc", Median_CV < CV_thresh) %>%
  dplyr::select(ID)

##### Combine filtering steps and apply to original report data frame:
CellsToKeep <- CellsToKeep_CV
rep.selected <- CellsToKeep %>% left_join(rep, by=c("ID"="ID"))

############ Collapse to gene level abundance using maxLFQ algorithm in DIA-NN R package for selected single cells

# ### Select for only single-cells passing filtering criteria
rep.selected.d04 <- rep.selected %>% filter(grepl("0|4", plex))
rep.selected.d04 <- rep.selected.d04 %>% dplyr::select(seqcharge, File.Name, Ms1.Area, Stripped.Sequence, Genes, Protein.Group)
rep.selected.d04 <- na.omit(rep.selected.d04)

### Use MS1 Area for gene level quantification
prot.sc <- diann_maxlfq(
  rep.selected.d04,
  sample.header = "File.Name",
  group.header = "Genes",
  id.header = "seqcharge",
  quantity.header = "Ms1.Area")

# Normalize matrix and log2 transform
prot.sc.norm.log <- prot.sc %>% Normalize_saad(log = "yes")

prot.sc.norm.log.colonly <- prot.sc %>% Normalize_saad_column_only(log = "yes")

############ Evaluate protein group data-completeness across single cells

### Tabulate data
prot.sc.norm.log.df <- reshape2::melt(prot.sc.norm.log)
colnames(prot.sc.norm.log.df) <- c("Genes", "File.Name", "PG_FC_log2")
prot.sc.norm.log.df <- prot.sc.norm.log.df %>% left_join(sample_map, by=c("File.Name"="File.Name"))

### Calculate total signal and number of protein groups per cell
summary_df_prot.sc <- prot.sc.norm.log.df %>%
  group_by(File.Name) %>%
  summarise(Total_Protein_Signal = sum(2^(PG_FC_log2), na.rm = TRUE),
            Proteins_Quantified = n_distinct(Genes[!is.na(PG_FC_log2)]))
summary_df_prot.sc <- summary_df_prot.sc %>% left_join(sample_map, by=c("File.Name"="File.Name"))

### completeness calculation
summary_df_sampcomp.sc <- prot.sc.norm.log.df %>%
  group_by(Genes) %>%
  summarise(Samples_Detected = n_distinct(File.Name[!is.na(PG_FC_log2)]),
            Total_Samples = n_distinct(File.Name),
            Dataset_Percentage = (Samples_Detected / Total_Samples)*100)

plot_scProtCompleteness <- ggplot(summary_df_sampcomp.sc, aes(x=Dataset_Percentage))+
  geom_histogram(bins = 50, position = "identity", alpha = 0.8)+
  labs(y = "# Proteins",
       x = "% of Single Cells",
       title = paste0('# Cells = ', length(CellsToKeep$ID))) +
  theme_light() +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 22),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20))

############ Impute single cell x protein group matrix with hknn
############ Then batch correct for mTRAQ label with ComBat
############ Then remove imputed values from matrix and process matrix for PCA

### To protein matrix that is normalized and log2 transformed, map the mTRAQ labels used for batch correction
batch_label  <- as.data.frame(colnames(prot.sc.norm.log))
colnames(batch_label)[1] <- 'ID'
batch_label <- batch_label %>% left_join(sample_map,by = c('ID'))
batch_label <- batch_label %>% dplyr::select(ID, plex)
batch_label <- unique(batch_label)

### Impute kNN, use hknn function
sc.imp <- hknn(prot.sc.norm.log,3)
sc.imp <- Normalize_saad_log(sc.imp)

sc.imp.colonly <- hknn(prot.sc.norm.log.colonly, 3)
sc.imp.colonly <- Normalize_saad_log_column_only(sc.imp.colonly)

### Batch correct for mTRAQ label bias using ComBat
sc.batch_cor <- ComBat(sc.imp, batch=factor(batch_label$plex))
sc.batch_cor <- Normalize_saad_log(sc.batch_cor)

sc.batch_cor.colonly <- ComBat(sc.imp.colonly, batch = factor(batch_label$plex))
sc.batch_cor.colonly <- Normalize_saad_log_column_only(sc.batch_cor.colonly)

### NAing imputed values post batch correction
sc.batch_cor[is.na(prot.sc.norm.log)==T] <- NA

sc.batch_cor.colonly[is.na(prot.sc.norm.log.colonly)] <- NA



### Steps for PCA
cor_mat <- cor(sc.batch_cor,use = 'pairwise.complete.obs')
sc.pca <- eigen(cor_mat)
scx<-as.data.frame(sc.pca$vectors)
percent_var <- sc.pca$values/sum(sc.pca$values)*100
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$ID <-colnames(sc.batch_cor)

### Map the data to calculated quantities and cell meta data
scx <- scx %>% left_join(batch_label,by = c('ID'))
scx <- scx %>% dplyr::select(!plex)
scx <- scx %>% left_join(summary_df_prot.sc, by=c("ID"="ID"))
sc.batch_cor.df <- as.data.frame(sc.batch_cor)
sc.batch_cor.df$Genes <- rownames(sc.batch_cor.df)
sc.batch_cor.df$Genes <- gsub(";", "", sc.batch_cor.df$Genes)
sc.batch_cor.df$Genes <- gsub("-", "", sc.batch_cor.df$Genes)
rownames(sc.batch_cor.df) <- sc.batch_cor.df$Genes
sc.batch_cor.df <- sc.batch_cor.df %>% dplyr::select(!Genes)
sc.batch_cor.t <- t(sc.batch_cor.df)
sc.batch_cor.t <- as.data.frame(sc.batch_cor.t)
sc.batch_cor.t$ID <- rownames(sc.batch_cor.t)
scx.pg <- scx %>% inner_join(sc.batch_cor.t, by=c("ID"="ID"))




### Visualize batch stats/protein abundances mapped onto PCA
# Mass Tag
ggplot(scx,aes(x = PC1,y = PC2, color = plex))+
  geom_point(size = 5, alpha = 0.5) +
  theme_classic() +
  ggtitle('Label Bias')
# Total protein abundance
ggplot(scx,aes(x = PC1,y = PC2, color = Total_Protein_Signal)) +
  geom_point(size = 5, alpha = 1) +
  theme_classic()+
  ggtitle('Cell Intensity')+
  scale_color_gradient2(midpoint = median(scx$Total_Protein_Signal,na.rm = T), low = "blue", mid = "white",high = "red",name = '')
# Sample type
ggplot(scx,aes(x = PC1,y = PC2, color = sample)) +
  geom_point(size = 5, alpha = 0.5) +
  theme_classic()+
  ggtitle('Cell Type') +
  xlab(paste0('PC1 ',round(percent_var[1]),'%'))+ylab(paste0('PC2 ',round(percent_var[2]),'%'))
# Run order
ggplot(scx,aes(x = PC1,y = PC2, color = Order)) +
  geom_point(size = 5, alpha = 0.5) +
  theme_classic()+
  ggtitle('Run Order') +
  xlab(paste0('PC1 ',round(percent_var[1]),'%'))+ylab(paste0('PC2 ',round(percent_var[2]),'%'))
# All cells
ggplot(scx.pg, aes(x=PC1, y=PC2)) +
  geom_point(size = 5,  alpha = 0.9) +
  ggtitle('') +
  xlab(paste0('PC1 ',round(percent_var[1]),'%'))+
  ylab(paste0('PC2 ',round(percent_var[2]),'%')) +
  theme_light()+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        plot.title = element_text(size = 22, hjust = 0.5),
        legend.text = element_text(size = 22),
        legend.title = element_blank())
# Marker protein - CD14
ggplot(scx.pg, aes(x=PC1, y=PC2, color = CD14)) +
  geom_point(size = 5,  alpha = 0.9, aes(colour = CD14)) +
  ggtitle('CD14; Monocyte Differentiation Antigen') +
  xlab(paste0('PC1 ',round(percent_var[1]),'%')) +
  ylab(paste0('PC2 ',round(percent_var[2]),'%')) +
  scale_colour_gradient2(low = "navy", mid = "white", high = "firebrick3")+
  labs(color = "Protein FC") +
  theme_light()+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        plot.title = element_text(size = 22, hjust = 0.5),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18))


### Want the protein matrix (mTRAQ tag batch corrected)
### And save it with corresponding Genes for downstream integration
prot_matrix_n3p1 <- sc.batch_cor
prot_df_n3p1 <- reshape2::melt(prot_matrix_n3p1)
colnames(prot_df_n3p1) <- c("Genes", "ID", "Prot_FC_log2")
cellID_donorID <- meta.data %>% dplyr::select(sample, ID) %>% distinct()
prot_df_n3p1 <- prot_df_n3p1 %>% left_join(cellID_donorID, by=c("ID"="ID"))
prot_df_n3p1$dataset_batch <- "n3p1"
prot_df_n3p1$ID <- paste0(prot_df_n3p1$ID, "_", prot_df_n3p1$dataset_batch)


### Save a list of selected cell IDs that passed quality control

fwrite(prot_df_n3p1,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/Genes_n3p1_for_prot_integration.csv",
       row.names = FALSE)

# Save parallel output CSV
prot_matrix_n3p1.colonly <- sc.batch_cor.colonly
prot_df_n3p1.colonly <- reshape2::melt(prot_matrix_n3p1.colonly)
colnames(prot_df_n3p1.colonly) <- c("Genes", "ID", "Prot_FC_log2")
prot_df_n3p1.colonly <- prot_df_n3p1.colonly %>% left_join(cellID_donorID, by = c("ID" = "ID"))
prot_df_n3p1.colonly$dataset_batch <- "n3p1"
prot_df_n3p1.colonly$ID <- paste0(prot_df_n3p1.colonly$ID, "_", prot_df_n3p1.colonly$dataset_batch)
fwrite(prot_df_n3p1.colonly,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/Genes_n3p1_for_abundance_analysis.csv",
       row.names = FALSE)


### Peptide level

### Use MS1 Area for peptide level quantification (same input as protein-level maxLFQ)
peptide.sc <- diann_maxlfq(
  rep.selected.d04,
  sample.header  = "File.Name",
  group.header   = "Stripped.Sequence",
  id.header      = "seqcharge",
  quantity.header= "Ms1.Area")

### Normalize: relative (row+column) and absolute-friendly (column-only), both log2
peptide.sc.norm.log         <- peptide.sc %>% Normalize_saad(log = "yes")
peptide.sc.norm.log.colonly <- peptide.sc %>% Normalize_saad_column_only(log = "yes")

### Batch label (plex), derived from the peptide matrix columns
batch_label_pep <- as.data.frame(colnames(peptide.sc.norm.log))
colnames(batch_label_pep)[1] <- "ID"
batch_label_pep <- batch_label_pep %>% left_join(sample_map, by = "ID")
batch_label_pep <- batch_label_pep %>% dplyr::select(ID, plex) %>% unique()

### kNN imputation
sc.imp.pep         <- hknn(peptide.sc.norm.log, 3)
sc.imp.pep         <- Normalize_saad_log(sc.imp.pep)

sc.imp.pep.colonly <- hknn(peptide.sc.norm.log.colonly, 3)
sc.imp.pep.colonly <- Normalize_saad_log_column_only(sc.imp.pep.colonly)

### ComBat batch correction on mTRAQ plex
sc.batch_cor.pep         <- ComBat(sc.imp.pep,         batch = factor(batch_label_pep$plex))
sc.batch_cor.pep         <- Normalize_saad_log(sc.batch_cor.pep)

sc.batch_cor.pep.colonly <- ComBat(sc.imp.pep.colonly, batch = factor(batch_label_pep$plex))
sc.batch_cor.pep.colonly <- Normalize_saad_log_column_only(sc.batch_cor.pep.colonly)

### NA imputed values post batch correction (mirrors protein pipeline)
sc.batch_cor.pep[is.na(peptide.sc.norm.log) == TRUE]                 <- NA
sc.batch_cor.pep.colonly[is.na(peptide.sc.norm.log.colonly) == TRUE] <- NA

### Peptide -> protein mapping table (filtered to peptides in quantified matrix)
pep_prot_map_n3p1 <- rep.selected.d04 %>%
  dplyr::select(Stripped.Sequence, Genes, Protein.Group) %>%
  distinct()

pep_prot_map_n3p1 <- pep_prot_map_n3p1 %>%
  filter(Stripped.Sequence %in% rownames(sc.batch_cor.pep))

### Identify proteotypic peptides: unambiguously map to one gene and one canonical
### protein group, with no isoform suffix and no semicolon-separated multi-assignment
proteotypic_flag <- pep_prot_map_n3p1 %>%
  group_by(Stripped.Sequence) %>%
  summarise(
    n_genes            = n_distinct(Genes),
    n_protein_groups   = n_distinct(Protein.Group),
    multi_gene_string  = any(grepl(";", Genes)),
    multi_pg_string    = any(grepl(";", Protein.Group)),
    isoform_assignment = any(grepl("-[0-9]+$", Protein.Group)),
    .groups = "drop"
  ) %>%
  mutate(is_proteotypic = (n_genes == 1) & (n_protein_groups == 1) &
           !multi_gene_string & !multi_pg_string & !isoform_assignment)

pep_prot_map_n3p1 <- pep_prot_map_n3p1 %>%
  left_join(proteotypic_flag, by = "Stripped.Sequence")
pep_prot_map_n3p1$dataset_batch <- "n3p1"

### Long-format peptide outputs, parallel to prot_df_n3p1 / prot_df_n3p1.colonly
cellID_donorID_pep <- sample_map %>% dplyr::select(sample, ID) %>% distinct()

# Relative (row+column normalized + ComBat plex)
peptide_df_n3p1 <- reshape2::melt(sc.batch_cor.pep)
colnames(peptide_df_n3p1) <- c("Stripped.Sequence", "ID", "Peptide_FC_log2")
peptide_df_n3p1 <- peptide_df_n3p1 %>% left_join(cellID_donorID_pep, by = "ID")
peptide_df_n3p1$dataset_batch <- "n3p1"
peptide_df_n3p1$ID <- paste0(peptide_df_n3p1$ID, "_", peptide_df_n3p1$dataset_batch)

# Absolute-friendly (column-only normalized + ComBat plex)
peptide_df_n3p1.colonly <- reshape2::melt(sc.batch_cor.pep.colonly)
colnames(peptide_df_n3p1.colonly) <- c("Stripped.Sequence", "ID", "Peptide_FC_log2")
peptide_df_n3p1.colonly <- peptide_df_n3p1.colonly %>% left_join(cellID_donorID_pep, by = "ID")
peptide_df_n3p1.colonly$dataset_batch <- "n3p1"
peptide_df_n3p1.colonly$ID <- paste0(peptide_df_n3p1.colonly$ID, "_", peptide_df_n3p1.colonly$dataset_batch)

### Save outputs
fwrite(peptide_df_n3p1,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/Peptides_n3p1_for_integration.csv",
       row.names = FALSE)
fwrite(peptide_df_n3p1.colonly,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/Peptides_n3p1_for_abundance_analysis.csv",
       row.names = FALSE)
fwrite(pep_prot_map_n3p1,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/PepProtMap_n3p1.csv",
       row.names = FALSE)



### Save a list of selected cell IDs that passed quality control


IdentificationsList_n3p1 <- rep.selected.d04 %>% distinct()

fwrite(IdentificationsList_n3p1,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/IdentificationsList_n3p1.csv",
       row.names = FALSE)





############ Measure accuracy and precision of synthetic yeast spike-ins

### Load spike in map
spike_in_map <- read.csv("/projects/slavov/LK/ImmPTR/data/meta_data/SyntheticYeast_SpikeIn_Map.csv")
spike_in_map$spike_ID <- paste0(spike_in_map$Well, "_", spike_in_map$Channel)

### Gather raw data, filter
### Initial light filter
### Remove extra bulk samples included for improving spectral library second-pass DIA-NN search
Raw_data_with_spike_in <- rep.full %>% filter(Lib.Q.Value <= 0.01, Lib.PG.Q.Value <= 0.05, Ms1.Area > 0, !grepl("MBR", Run)) 
sel.ID.diameters <- meta.data %>% dplyr::select(diameter, ID)
Raw_data_with_spike_in <- Raw_data_with_spike_in %>% left_join(sel.ID.diameters, by = c("ID"="ID"))
Raw_data_with_spike_in$SampleType <- ifelse(Raw_data_with_spike_in$plex == 8 & is.na(Raw_data_with_spike_in$diameter), "carrier",
                                            ifelse(is.na(Raw_data_with_spike_in$diameter) & Raw_data_with_spike_in$plex != 8, "neg",
                                                   ifelse(Raw_data_with_spike_in$plex %in% c(0, 4) & !is.na(Raw_data_with_spike_in$diameter), "sc", NA)))

############ Begin processing fully mapped raw data with initial FDR filters

### Filter using Q-value metrics from DIA-NN output
### Second-pass search: Run-specific Precursor FDR 1%
### Second-pass search: Run-specific Protein Group FDR 5%
### Keep precursors with non-zero MS1 signal
### Assign carrier vs single-cell vs negative control sample type information
### Load the saved list of selected cells passing QC, analyze those cells only
rep_with_spike_in <- Raw_data_with_spike_in
rep_analyze_spike_ins <- rep_with_spike_in
rep_analyze_spike_ins$spike_ID <- paste0(rep_analyze_spike_ins$Well, "_", rep_analyze_spike_ins$plex)
rep_analyze_spike_ins <- rep_analyze_spike_ins %>% filter(grepl("sc|neg", SampleType))
rep_analyze_spike_ins <- rep_analyze_spike_ins %>% dplyr::select(Ms1.Area, spike_ID, seqcharge, Protein.Names, File.Name)
rep_analyze_spike_ins <- rep_analyze_spike_ins %>% left_join(spike_in_map, by=c("spike_ID"="spike_ID"))

### Add Species column based on conditions in Protein.Names
rep_analyze_spike_ins <- rep_analyze_spike_ins %>%
  mutate(Analyte_Class = case_when(
    grepl("HUMAN", Protein.Names) ~ "Single-Cell",
    grepl("synYeast", Protein.Names) ~ "Spike-In",
    TRUE ~ "Other"  # Optional: Handle unmatched cases
  ))

### Perform regression analysis of spike-ins to evaluate precision and accuracy
### Start by filtering for data-completeness

spike_ins_seqcharges <- rep_analyze_spike_ins %>% filter(Analyte_Class == "Spike-In")
spike_ins_seqcharges <- spike_ins_seqcharges %>% dplyr::select(spike_ID, seqcharge, Ms1.Area)
spike_ins_seqcharges <- spike_ins_seqcharges %>% pivot_wider(names_from = spike_ID, values_from = Ms1.Area) %>% as.data.frame()
rownames(spike_ins_seqcharges) <- spike_ins_seqcharges$seqcharge
spike_ins_seqcharges$seqcharge <- NULL
spike_ins_seqcharges$completeness <- rowMeans(!is.na(spike_ins_seqcharges)) * 100
spike_ins_seqcharges$seqcharge <- rownames(spike_ins_seqcharges)
spike_ins_seqcharges <- spike_ins_seqcharges %>% filter(completeness > 70)

seqcharges_quantified <- as.data.frame(spike_ins_seqcharges$seqcharge)
colnames(seqcharges_quantified) <- c("seqcharge")

### Then, map the spike in levels to each sample 
spike_ins_seqcharges <- spike_ins_seqcharges %>% dplyr::select(-completeness)
rownames(spike_ins_seqcharges) <- NULL
spike_ins_seqcharges <- reshape2::melt(spike_ins_seqcharges) %>% na.omit()
colnames(spike_ins_seqcharges) <- c("seqcharge","spike_ID","Ms1.Area")
spike_ins <- rep_analyze_spike_ins %>% filter(Analyte_Class == "Spike-In") %>% dplyr::select(spike_ID, SpikeIn_Level)
spike_ins_seqcharges <- spike_ins_seqcharges %>% left_join(spike_ins, by=c("spike_ID"="spike_ID")) %>% distinct() 

### Per precursor and per level, remove observations that are below the median value across samples - 1 * median absolute deviation 
### Remove noisy, low quality signal, accounting for issues in accurate dispensing of spike-ins during sample prep
spike_ins <- spike_ins_seqcharges %>%
  group_by(seqcharge, SpikeIn_Level) %>%
  mutate(mad_value = mad(Ms1.Area, constant = 1.4826, na.rm = TRUE),
         median_value = median(Ms1.Area, na.rm = TRUE)) %>%
  filter(Ms1.Area >= median_value - 1 * mad_value) %>%
  dplyr::select(-mad_value, -median_value) %>%  # Remove temporary MAD/median columns
  ungroup()

### Calculate CVs per precursor per level after mean normalization per sample
### Then take the mean precursor CV per level
sum_CV_spike_ins <- spike_ins %>% dplyr::select(-SpikeIn_Level) %>% 
  pivot_wider(names_from = spike_ID, values_from = Ms1.Area) %>% as.data.frame()
rownames(sum_CV_spike_ins) <- sum_CV_spike_ins$seqcharge
sum_CV_spike_ins <- sweep(sum_CV_spike_ins %>% dplyr::select(-seqcharge), 
                          2, 
                          colMeans(sum_CV_spike_ins %>% dplyr::select(-seqcharge), na.rm = TRUE), 
                          FUN = "/") %>% as.data.frame()
sum_CV_spike_ins$seqcharge <- rownames(sum_CV_spike_ins)
rownames(sum_CV_spike_ins) <- NULL
sum_CV_spike_ins <- reshape2::melt(sum_CV_spike_ins)
colnames(sum_CV_spike_ins) <- c("seqcharge","spike_ID", "Ms1.Area")
sum_CV_spike_ins <- sum_CV_spike_ins %>% left_join(spike_in_map, by=c("spike_ID"="spike_ID")) %>% distinct()
sum_CV_spike_ins <- sum_CV_spike_ins %>% group_by(seqcharge, SpikeIn_Level) %>%
  summarise(CV_seqcharge = sd(Ms1.Area, na.rm = TRUE) / 
              mean(Ms1.Area, na.rm = TRUE)) %>% ungroup()
sum_CV_spike_ins <- sum_CV_spike_ins %>% group_by(SpikeIn_Level) %>%
  summarise(CV_level = mean(CV_seqcharge))

### Calculate a reference vector for the precursors by taking their median values across all samples with level 1 spike in
ref_median_spike_in_level_1 <- spike_ins %>% filter(SpikeIn_Level == 1) %>%
  group_by(seqcharge) %>%
  summarise(ref_median = median(Ms1.Area, na.rm=TRUE))

### Normalize each sample to this reference vector to get measurements relative to the lowest level
spike_ins <- spike_ins %>% dplyr::select(spike_ID, Ms1.Area, seqcharge) %>%
  pivot_wider(names_from = spike_ID, values_from = Ms1.Area)
spike_ins_normalized <- spike_ins %>%
  left_join(ref_median_spike_in_level_1, by = "seqcharge") %>%
  mutate(across(
    where(is.numeric) & !all_of("ref_median"),
    ~ ifelse(is.na(.), NA, . / ref_median)
  )) %>%
  dplyr::select(-ref_median)
spike_ins_normalized <- reshape2::melt(spike_ins_normalized)
colnames(spike_ins_normalized) <- c("seqcharge","spike_ID","Ms1.Area.Ref.Med.Norm")
spike_ins_normalized$sequence <- gsub('[[:digit:]]+', '', spike_ins_normalized$seqcharge)

spike_ins_normalized_ListSeqCharge <- spike_ins_normalized %>% dplyr::select(seqcharge) %>% distinct()

spike_ins_normalized <- spike_ins_normalized %>% group_by(spike_ID, sequence) %>%
  summarise(Ms1.Area.Ref.Med.Norm.Per.Seq = median(Ms1.Area.Ref.Med.Norm, na.rm = TRUE)) %>%
  ungroup()
spike_ins_normalized <- spike_ins_normalized %>% left_join(spike_in_map, by=c("spike_ID"="spike_ID")) %>%
  dplyr::select(Ms1.Area.Ref.Med.Norm.Per.Seq, sequence, spike_ID, SpikeIn_Level) %>% distinct()

### Collect results, count number of datapoints per level
sum_stats_spike_ins <- spike_ins_normalized %>% 
  group_by(SpikeIn_Level) %>%
  summarise(
    Median_spike_in_measured_per_level = median(Ms1.Area.Ref.Med.Norm.Per.Seq, na.rm = TRUE),
    Num_data_points = n_distinct(spike_ID))

### Join with raw precursor median CV per level
sum_stats_spike_ins <- sum_stats_spike_ins %>% left_join(sum_CV_spike_ins, by=c("SpikeIn_Level"="SpikeIn_Level"))

### Begin the robust regression analysis and plot the data
# Add log2-transformed columns
sum_stats_spike_ins <- sum_stats_spike_ins %>% 
  mutate(log2_SpikeIn_Level = log2(SpikeIn_Level),
         log2_Median_spike_in_measured_per_level = log2(Median_spike_in_measured_per_level))
# Perform robust regression using rlm()
robust_model <- rlm(log2_Median_spike_in_measured_per_level ~ log2_SpikeIn_Level,
                    data = sum_stats_spike_ins)
# Extract slope (coefficient)
slope <- coef(robust_model)[2]  # Second coefficient is the slope
predicted <- predict(robust_model, newdata = sum_stats_spike_ins)
r2 <- cor(sum_stats_spike_ins$log2_Median_spike_in_measured_per_level, predicted)^2
# Add predicted values for plotting
sum_stats_spike_ins$predicted <- predicted
# Plot with robust regression line and annotations

### Collect the set of spike in precursors that made it passed filtering steps and were analyzed by robust regression
### Plot their log2 transformed MS1 Areas compared to the rest of the precursors measured per single-cell
spike_ins_normalized_ListIDs <- spike_ins_normalized %>% dplyr::select(spike_ID) %>% distinct()
rep_analyze_spike_ins <- rep_analyze_spike_ins[rep_analyze_spike_ins$spike_ID %in% spike_ins_normalized_ListIDs$spike_ID, ]
rep_analyze_spike_ins_HUMAN <- rep_analyze_spike_ins %>% filter(grepl("HUMAN", Protein.Names))
rep_analyze_spike_ins_SPIKE <- rep_analyze_spike_ins %>% filter(!grepl("HUMAN", Protein.Names))

spike_ins_normalized_ListSeq <- spike_ins_normalized_ListSeqCharge %>% dplyr::select(seqcharge) %>% distinct()

rep_analyze_spike_ins_SPIKE <- rep_analyze_spike_ins_SPIKE %>%
  filter(seqcharge %in% spike_ins_normalized_ListSeq$seqcharge)

rep_analyze_spike_ins_SPIKE <- rep_analyze_spike_ins_SPIKE %>% dplyr::select(Ms1.Area, Analyte_Class)
rep_analyze_spike_ins_HUMAN <- rep_analyze_spike_ins_HUMAN %>% dplyr::select(Ms1.Area, Analyte_Class)
rep_analyze_spike_ins <- rbind(rep_analyze_spike_ins_SPIKE, rep_analyze_spike_ins_HUMAN)


RobustRegression_SpikeIns_n3p1 <- ggplot(sum_stats_spike_ins, 
                                         aes(x = log2_SpikeIn_Level, 
                                             y = log2_Median_spike_in_measured_per_level)) +
  geom_point(aes(color = CV_level), size = 2.2) +
  geom_line(aes(y = predicted), color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_text(aes(x = log2_SpikeIn_Level, y = -0.15,
                label = Num_data_points),
            family = "sans", vjust = 1, size = 1.8) +
  scale_color_gradient(low = "lightblue", high = "blue4", name = "Mean CV") +
  xlim(-0.1, 3.2) +
  ylim(-0.3, 3.2) +
  labs(x = expression(log[2]~"(Spike-In Amount)"),
       y = expression(log[2]~"(Measured Amount)")) +
  annotate("text", x = 0.1, y = 3.1,
           label = paste0("Slope: ", round(slope, 2), "\nR²: ", round(r2, 3)),
           hjust = 0, vjust = 1, size = 2.2, color = "black") +
  theme_light(base_size = 8) +
  theme(
    aspect.ratio = 1,
    plot.title = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  )



SingleCell_Versus_SpikeIn_Abundances_n3p1 <- ggplot(rep_analyze_spike_ins, 
       aes(y = log2(Ms1.Area), x = Analyte_Class)) +
  rasterise(geom_violin(width = 1, linewidth = 0.3, fill = "grey85"), dpi = 600) +
  rasterise(geom_boxplot(width = 0.1, color = "black", alpha = 0.6, 
                         linewidth = 0.3, outlier.size = 0.4), dpi = 600) +
  labs(y = expression(log[2]~"(Precursor Abundance)")) +
  theme_light(base_size = 8) +
  theme(
    aspect.ratio = 1,
    plot.title = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    plot.margin = margin(2, 2, 2, 2)
  )