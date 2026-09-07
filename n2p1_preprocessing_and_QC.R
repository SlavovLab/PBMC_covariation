########################################################################################
## Processing raw data and quality control for proteomic dataset: n2p1
########################################################################################





############ Below are dataset-specific functions required for mapping cell meta data from an earlier version of nPOP via the CellenONE platform
############ Load these functions for analysis of this dataset only, otherwise they are not required
############ Additionally, load all functions first from: library>functions>Functions_for_mapping_and_analysis.R 
############ As they are also required for this analysis, then run this script





### Load functions specific to this dataset that were generated using earlier version of nPOP/QuantQC mapping
### This function performs mapping of cell sorting meta data onto the raw data
analyzeCellenONE_mTRAQ <- function(allDays,pickupPath2,pickupPath1,labelPath){
  allDays[grepl("Transmission",allDays$X),]$X <- NA
  allDays <- allDays %>% fill(2:7, .direction = "up") %>% drop_na(XPos)
  #### Labelling and Pickup Field Files
  ## labelling file
  con_lab <-file(labelPath)
  lines_lab<- readLines(con_lab)
  close(con_lab)
  slines_lab <- strsplit(lines_lab,"\t")
  colCount_lab <- max(unlist(lapply(slines_lab, length)))
  label <- read.table(labelPath, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)
  fieldRaw_label <- label[grepl("\\[\\d",label$position),]$position
  fieldBrack_label <- gsub("\\[|\\]", "", fieldRaw_label)
  fieldComma_label <- strsplit(fieldBrack_label, "," )
  fieldNum_label <- unlist(lapply(fieldComma_label, `[[`, 1))
  label$field[grepl("\\[\\d",label$position)] <- fieldNum_label
  label <- label %>% fill(field, .direction = "down")
  label <-  label[!label$well=="",]
  label$field <- as.numeric(label$field) + 1
  labelxyPos <- strsplit(label$position, "\\/")
  label$yPos <- unlist(lapply(labelxyPos, '[[', 1))
  label$xPos <- unlist(lapply(labelxyPos, '[[', 2))
  ## sample pickup file
  con_pickup <-file(pickupPath1)
  lines_pickup<- readLines(con_pickup)
  close(con_pickup)
  slines_pickup <- strsplit(lines_pickup,"\t")
  colCount_pickup <- max(unlist(lapply(slines_pickup, length)))
  pickup <- read.table(pickupPath1, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)
  con_pickup <-file(pickupPath2)
  lines_pickup<- readLines(con_pickup)
  close(con_pickup)
  slines_pickup <- strsplit(lines_pickup,"\t")
  colCount_pickup <- max(unlist(lapply(slines_pickup, length)))
  pickup2 <- read.table(pickupPath2, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 27)
  pickup <- rbind(pickup,pickup2)
  fieldRaw_pickup <- pickup[grepl("\\[\\d",pickup$position),]$position
  fieldBrack_pickup <- gsub("\\[|\\]", "",fieldRaw_pickup)
  fieldComma_pickup <- strsplit(fieldBrack_pickup, "," )
  fieldNum_pickup <- unlist(lapply(fieldComma_pickup, `[[`, 1))
  pickup$field[grepl("\\[\\d",pickup$position)] <- fieldNum_pickup
  pickup <- pickup %>% fill(field, .direction = "down")
  pickup <-  pickup[!pickup$well=="",]
  pickup$field <- as.numeric(pickup$field) + 1
  pickupxyPos <- strsplit(pickup$position, "\\/")
  pickup$yPos <- unlist(lapply(pickupxyPos, '[[', 1))
  pickup$xPos <- unlist(lapply(pickupxyPos, '[[', 2))
  ### making sure that label "slide" refers to fields or actual slides.
  ## There is a sequence of 108 x,y and there are 108 spots per slide
  # order is y/x positions
  LabelxyPos <- strsplit(label$position, "\\/")
  label$yPos <- unlist(lapply(LabelxyPos, '[[', 1))
  label$xPos <- unlist(lapply(LabelxyPos, '[[', 2))
  label$well <- substring(label$well, 3)
  label$well <- as.numeric(gsub(",","",label$well))
  matchTMTSCP <- paste0("TMT", 1:length(unique(label$well)))
  for (i in 1:length(matchTMTSCP)) {
    label[which(label$well == i),]$well <- matchTMTSCP[i]
  }
  label$well <- substring(label$well, 4)
  #### Trying to map label to cell
  ###  lets try and keep all three important in one
  allDays <- transform(allDays, xyf = paste0(XPos, YPos, Field))
  label <- transform(label, xyf = paste0(xPos, yPos, field))
  ### Isolation and Label merged
  isoLab <- allDays %>% group_by(Target) %>% left_join(label, by = 'xyf')
  labelCount <- isoLab %>% group_by(well) %>% dplyr::summarize(count=n())
  ## trying to get pickup working in the same way
  # fixing pickup file
  pickup$target <- 1
  pickup$well <- substring(pickup$well, 2)
  pickup$well <- gsub(",","",pickup$well)
  slidesUsedForPrep <- unique(isoLab$Target)
  
  ### clustering together pickup points with sample points using ANN
  ## super convolutedly/idiotically
  isoLab$ann <- NA
  isoLab$pickupX <- NA
  isoLab$pickupY <- NA
  ann_123 <- ann(ref = as.matrix(unique(pickup[(-which(pickup$field == 4)) , c("xPos","yPos")])),  target = as.matrix(isoLab[(-which(isoLab$Field == 4)) , c("xPos","yPos")]), k=1)
  
  ann_4 <- ann(ref = as.matrix(unique(pickup[(which(pickup$field == 4)) , c("xPos","yPos")])),  target = as.matrix(isoLab[(which(isoLab$Field == 4)) , c("xPos","yPos")]), k=1)
  isoLab[(-which(isoLab$Field == 4)),]$ann <-  ann_123$knnIndexDist[,1]
  isoLab[(which(isoLab$Field == 4)),]$ann <-  ann_4$knnIndexDist[,1]
  ## split - combine
  isoLab_123 <- isoLab[-which(isoLab$Field == 4),]
  isoLab_4 <- isoLab[which(isoLab$Field == 4),]
  notFieldFourPickUnique <- unique(pickup[-(which(pickup$field == 4)) , c("xPos","yPos")])
  fieldFourPickUnique <- unique(pickup[(which(pickup$field == 4)) , c("xPos","yPos")])
  isoLab_123$pickupX <- notFieldFourPickUnique[isoLab_123$ann,]$xPos
  isoLab_123$pickupY <- notFieldFourPickUnique[isoLab_123$ann,]$yPos
  isoLab_4$pickupX <- fieldFourPickUnique[isoLab_4$ann,]$xPos
  isoLab_4$pickupY <- fieldFourPickUnique[isoLab_4$ann,]$yPos
  isoLab_bound <- rbind(isoLab_123, isoLab_4)
  ### Merge pickup and isoLab
  isoLab_bound <- transform(isoLab_bound, xyft = paste0(pickupX, pickupY, Field, Target))
  pickup <-  transform(pickup, xyft = paste0(xPos, yPos, field, target))
  
  isoLab_final <- isoLab_bound %>% left_join(pickup, by = 'xyft')
  wellCount <- isoLab_final %>% group_by(well.y) %>% dplyr::summarize(count=n())
  ### Clean up to yield final dataframe
  cellenOne_data <- data.frame (sample = isoLab_final$condition, isoTime = isoLab_final$Time, diameter = isoLab_final$Diameter, elongation = isoLab_final$Elongation, slide = isoLab_final$Target, field = isoLab_final$Field, dropXPos = isoLab_final$XPos, dropYPos = isoLab_final$YPos, label = isoLab_final$well.x, pickupXPos = isoLab_final$pickupX, pickupYPos = isoLab_final$pickupY, injectWell = isoLab_final$well.y)
  
  cellenOne_data$wellAlph <- substring(cellenOne_data$injectWell, 1, 1)
  cellenOne_data$wellNum <- substring(cellenOne_data$injectWell, 2)
  cellenOne_data <- cellenOne_data %>% arrange(wellAlph, as.numeric(wellNum), as.numeric(label))
  
  cellenOne_data$pickupXPos_numb <- as.numeric(cellenOne_data$pickupXPos)
  cellenOne_data$pickupYPos_numb  <- as.numeric(cellenOne_data$pickupYPos)
  return(cellenOne_data)
}


### This function reads in the DIANN output containing raw data and extracts mTRAW mass tag information
Read.DIANN <- function(path){
  columns_to_read <-c("Genes", "Run", "Q.Value", "PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "RT", 
                      "Precursor.Id", "Stripped.Sequence", "Precursor.Mz", 
                      "Precursor.Charge", "Precursor.Quantity", "Ms1.Area", 
                      "Protein.Group", "Translated.Q.Value", "Channel.Q.Value")
  Raw_data <- data.table::fread(path,select = columns_to_read)
  Raw_data <- as.data.frame(Raw_data)
  Raw_data <- Raw_data %>% filter(Lib.PG.Q.Value < .05)
  Raw_data <- Raw_data %>% filter(Run %in% linker$Run)
  Raw_data <- Raw_data %>% left_join(linker, by = c('Run'))
  #Unique precursor ID
  Raw_data$seqcharge <- paste0(Raw_data$Stripped.Sequence,Raw_data$Precursor.Charge)
  Raw_data <- Raw_data %>% filter(Protein.Group != '')
  # this grabs the mTRAQ tag used (may need to be adjusted if using different multiplexing tag)
  Raw_data$plex <- substr(Raw_data$Precursor.Id[1:nrow(Raw_data)], 10, 10)
  # Unique cell ID
  Raw_data$ID <- paste0(Raw_data$Well,Raw_data$Plate,'.',Raw_data$plex)
  Raw_data$File.Name <- Raw_data$ID
  #Remove redundant data points
  Raw_data$uq <- paste0(Raw_data$File.Name,Raw_data$Protein.Group,Raw_data$seqcharge)
  Raw_data <- Raw_data %>% distinct(uq,.keep_all = T)
  Raw_data$uq <- NULL
  return(Raw_data)
}





############ Begin loading, mapping, and processing the raw data

### Read in the path to DIA-NN file
data_path <- '/projects/slavov/LK/ImmPTR/data/raw_data/n2p1/diann_b17/report.tsv'

### List cell types
cell_types <- c('Donor_04','Donor_05')

### Read in cell isolation files from CellenONE and assign cell type
c1 <- read.table(file = "/projects/slavov/LK/ImmPTR/data/meta_data/n2p1/Donor04_2plex_isolated_t2.csv", sep = ",", header = TRUE)
c1$condition <- "Donor_04"
c2 <- read.table(file = "/projects/slavov/LK/ImmPTR/data/meta_data/n2p1/Donor05_2plex_isolated.csv", sep = ",", header = TRUE)
c2$condition <- "Donor_05"
all_cells <- rbind(c1,c2)
pickupPath2 <- '/projects/slavov/LK/ImmPTR/data/meta_data/n2p1/Pickup_2_mock_jason.fld'
pickupPath1 <- '/projects/slavov/LK/ImmPTR/data/meta_data/n2p1/Pickup_1_mock_jason.fld'
labelPath <- '/projects/slavov/LK/ImmPTR/data/meta_data/n2p1/Labels_fromJDused.fld'





############ Load the DIA-NN stats output file to generate a linker data frame

### Extract sample and batch information to generate a linker used to integrate with above sample prep data
stats <- diann_load("/projects/slavov/LK/ImmPTR/data/raw_data/n2p1/diann_b17/report.stats.tsv")

linker <- stats %>%
  dplyr::select(File.Name) %>%
  filter(grepl("n2p1", File.Name), !grepl("MBR", File.Name)) %>%
  mutate(
    Well = str_extract(File.Name, "(?<=_B-)[A-Z]\\d{1,2}(?=_)"),
    Run = gsub(".*\\\\", "", File.Name),
    Run = gsub("\\.d$", "", Run),
    Plate = "1"
  ) %>%
  dplyr::select(Run, Well, Plate) %>%
  mutate(
    Order = as.integer(str_sub(Run, -3)),
    Order = rank(Order, ties.method = "first")
  )





############ Load the full raw data from the DIA-NN output

### Load raw data
Raw_data_full <- data.table::fread(data_path,select = c("Genes", "Run", "Q.Value", "PG.Q.Value", "Lib.Q.Value", 
                                                        "Lib.PG.Q.Value", "RT", "Precursor.Id", "Stripped.Sequence", 
                                                        "Precursor.Mz", "Precursor.Charge", "Precursor.Quantity", 
                                                        "Ms1.Area", "Protein.Group", "Translated.Q.Value", 
                                                        "Channel.Q.Value", "Protein.Names"))

### Initial light filter
### Remove extra bulk samples included for improving spectral library second-pass DIA-NN search
Raw_data <- Raw_data_full %>% filter(Lib.PG.Q.Value < .05)
Raw_data <- Raw_data %>% left_join(linker, by = c('Run'))
Raw_data <- Raw_data %>% filter(!grepl("MBR", Run), grepl("HUMAN", Protein.Names))

ChQval <- 1
CV_thresh <- 1

carrier <- T
if(carrier == T){
  carrier_CH <- c('8')
}
plex_sc <- c('0','4') # here put the mTRAQ channels that were used for single cells
plex_used <- c(plex_sc,carrier_CH) # here put the mTRAQ channels that were used 





############ Mapping cell data to raw data

### Call mapping function
cellenOne_data <- analyzeCellenONE_mTRAQ(all_cells,pickupPath2,pickupPath1,labelPath)

### If using plex/scopeDIA samples are injected into 2 different plates,
### There are 16 fields and samples from 1-8 go into plate 1, 9-16 go into plate 2
cellenOne_data$plate <- NA
cellenOne_data$plate[cellenOne_data$field > 8] <- 2
cellenOne_data$plate[cellenOne_data$field <= 8] <- 1

### Assigning mTRAQ tags to wells where label was picked up out of during prep
### Each tag was dispensed multiple times so maps to multiple wells of plate
cellenOne_data$tag <- NA

if(length(plex_sc == 2)){
  cellenOne_data$tag[cellenOne_data$label %in% c('1','3','5','7')] <- '0'
  cellenOne_data$tag[cellenOne_data$label %in% c('2','4','6','8')] <- '4'
}

### Create cell ID to match DIA-NN report
cellenOne_data$ID <- paste0(cellenOne_data$injectWell,cellenOne_data$plate,cellenOne_data$tag)
cellenOne_data_small <- cellenOne_data %>% dplyr::select(ID,diameter,sample,tag,plate)

### Unique precursor ID
Raw_data <- as.data.frame(Raw_data)
Raw_data$seqcharge <- paste0(Raw_data$Stripped.Sequence,Raw_data$Precursor.Charge)
Raw_data <- Raw_data %>% filter(Protein.Group != '')

### This grabs the mTRAQ tag used (may need to be adjusted if using different multiplexing tag)
Raw_data$plex <- substr(Raw_data$Precursor.Id[1:nrow(Raw_data)], 10, 10)

### Filter report for relevant channels
### Create data frame that maps all cell types from cellenONE data to IDs from DIA-NN and 
### Also specifies which channels are negative controls (i.e. no match to cell sorting data)

if(carrier == F){
  Raw_data <- Raw_data %>% filter(plex %in% plex_sc)
  # Analogous cell ID made from DIA-NN report and linker file
  Raw_data$ID <- paste0(Raw_data$Well,Raw_data$Plate,Raw_data$plex)
  Raw_data$File.Name <- Raw_data$ID
  cellID <- unique(Raw_data$ID)
  cellID <- as.data.frame(cellID)
  colnames(cellID) <- 'ID'
  cellID <- cellID %>% left_join(cellenOne_data_small,by = c('ID'))
  cellID$sample[is.na(cellID$sample)==T] <- 'neg'
  
}else{
  # Analogous cell ID made from DIA-NN report and linker file
  Raw_data$ID <- paste0(Raw_data$Well,Raw_data$Plate,Raw_data$plex)
  Raw_data$Run <- paste0(Raw_data$Well,Raw_data$Plate)
  Raw_data$File.Name <- Raw_data$ID
  Raw_data_cell <- Raw_data %>% filter(plex %in% plex_sc)
  cellID <-as.data.frame(unique(Raw_data_cell$ID))
  colnames(cellID) <- 'ID'
  cellID <- cellID %>% left_join(cellenOne_data_small,by = c('ID'))
  cellID$sample[is.na(cellID$sample)==T ] <- 'neg'
}

### Match IDs from raw data to meta data
Raw_data <- Raw_data %>% left_join(cellID, by=c("ID"="ID"))





############ Begin processing fully mapped raw data with initial FDR filters

### Filter using Q-value metrics from DIA-NN output
### Second-pass search: Run-specific Precursor FDR 1%
### Second-pass search: Run-specific Protein Group FDR 5%
### Keep precursors with non-zero MS1 signal
### Assign carrier vs single-cell vs negative control sample type information
rep <- Raw_data %>% filter(Lib.Q.Value <= 0.01, Lib.PG.Q.Value <= 0.05, Ms1.Area > 0)
rep$SampleType <- ifelse(rep$plex == 8 & is.na(rep$diameter), "carrier",
                         ifelse(is.na(rep$diameter) & rep$plex != 8, "neg",
                                ifelse(rep$plex %in% c(0, 4) & !is.na(rep$diameter), "sc", NA)))


### Generate a sample map and a protein/gene/unique precursor sequence map
sample_map <- rep %>% dplyr::select(File.Name, SampleType, plex, Order, ID, Well, diameter, sample) %>% distinct()
sample_map <- sample_map %>%
  mutate(Cell_Volume = (4/3) * pi * (diameter / 2)^3)
seqcharge_PG_map <- rep %>% dplyr::select(seqcharge, Protein.Group) %>% distinct()
pep_PG_map <- rep %>% dplyr::select(Stripped.Sequence, seqcharge, Protein.Group) %>% distinct()





############ Next, collect batch statistics to filter on cell quality
############ Filtering data by selection criteria: IDs, signal, CV

### First filter by precursor IDs and total signal with Channel.Q.Value and Translated.Q.Value 10% FDR
### Apply filter, collect batch statistics
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
  labs(title = "Single Cell Selection by IDs and Signal (n2p1)",
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
  labs(title = "Single Cell Selection by CV (n2p1)",
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
batch_label <- batch_label %>% dplyr::select(ID, plex, Order)
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

# --- Batch correction for continuous Order using limma
batch_label <- batch_label[match(colnames(sc.batch_cor), batch_label$ID), ]
design <- model.matrix(~ Order, data=batch_label)
sc.batch_cor.limma <- removeBatchEffect(sc.batch_cor, covariates=design[, "Order", drop=FALSE])

# --- Normalize + restore NAs
sc.batch_cor.limma <- Normalize_saad_log(sc.batch_cor.limma)
sc.batch_cor.limma[is.na(prot.sc.norm.log)] <- NA

sc.batch_cor <- sc.batch_cor.limma


sc.batch_cor.limma.colonly <- removeBatchEffect(sc.batch_cor.colonly, covariates = design[, "Order", drop = FALSE])
sc.batch_cor.limma.colonly <- Normalize_saad_log_column_only(sc.batch_cor.limma.colonly)
sc.batch_cor.limma.colonly[is.na(prot.sc.norm.log.colonly)] <- NA
sc.batch_cor.colonly <- sc.batch_cor.limma.colonly

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
ggplot(scx.pg,aes(x = PC1,y = PC2, color = sample)) +
  geom_point(size = 5, alpha = 0.5) +
  theme_classic()+
  ggtitle('Cell Type') +
  xlab(paste0('PC1 ',round(percent_var[1]),'%'))+ylab(paste0('PC2 ',round(percent_var[2]),'%'))
# Run order
ggplot(scx,aes(x = PC1,y = PC2, color = Order.x)) +
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
prot_matrix_n2p1 <- sc.batch_cor
prot_df_n2p1 <- reshape2::melt(prot_matrix_n2p1)
colnames(prot_df_n2p1) <- c("Genes", "ID", "Prot_FC_log2")
cellID_donorID <- Raw_data %>% dplyr::select(sample, ID) %>% distinct()
prot_df_n2p1 <- prot_df_n2p1 %>% left_join(cellID_donorID, by=c("ID"="ID"))
prot_df_n2p1$dataset_batch <- "n2p1"
prot_df_n2p1$ID <- paste0(prot_df_n2p1$ID, "_", prot_df_n2p1$dataset_batch)


### Save a list of selected cell IDs that passed quality control

fwrite(prot_df_n2p1,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/Genes_n2p1_for_prot_integration.csv",
       row.names = FALSE)

prot_matrix_n2p1.colonly <- sc.batch_cor.colonly
prot_df_n2p1.colonly <- reshape2::melt(prot_matrix_n2p1.colonly)
colnames(prot_df_n2p1.colonly) <- c("Genes", "ID", "Prot_FC_log2")
prot_df_n2p1.colonly <- prot_df_n2p1.colonly %>% left_join(cellID_donorID, by = c("ID" = "ID"))
prot_df_n2p1.colonly$dataset_batch <- "n2p1"
prot_df_n2p1.colonly$ID <- paste0(prot_df_n2p1.colonly$ID, "_", prot_df_n2p1.colonly$dataset_batch)
fwrite(prot_df_n2p1.colonly,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/Genes_n2p1_for_abundance_analysis.csv",
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

### Batch label (plex + Order), derived from the peptide matrix columns
batch_label_pep <- as.data.frame(colnames(peptide.sc.norm.log))
colnames(batch_label_pep)[1] <- "ID"
batch_label_pep <- batch_label_pep %>% left_join(sample_map, by = "ID")
batch_label_pep <- batch_label_pep %>% dplyr::select(ID, plex, Order) %>% unique()

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

### --- limma correction for continuous Order (mirrors protein pipeline)
batch_label_pep <- batch_label_pep[match(colnames(sc.batch_cor.pep), batch_label_pep$ID), ]
design_pep <- model.matrix(~ Order, data = batch_label_pep)

sc.batch_cor.pep.limma         <- removeBatchEffect(sc.batch_cor.pep,
                                                    covariates = design_pep[, "Order", drop = FALSE])
sc.batch_cor.pep.limma         <- Normalize_saad_log(sc.batch_cor.pep.limma)
sc.batch_cor.pep.limma[is.na(peptide.sc.norm.log)] <- NA
sc.batch_cor.pep <- sc.batch_cor.pep.limma

sc.batch_cor.pep.limma.colonly <- removeBatchEffect(sc.batch_cor.pep.colonly,
                                                    covariates = design_pep[, "Order", drop = FALSE])
sc.batch_cor.pep.limma.colonly <- Normalize_saad_log_column_only(sc.batch_cor.pep.limma.colonly)
sc.batch_cor.pep.limma.colonly[is.na(peptide.sc.norm.log.colonly)] <- NA
sc.batch_cor.pep.colonly <- sc.batch_cor.pep.limma.colonly

### Peptide -> protein mapping table (filtered to peptides in quantified matrix)
pep_prot_map_n2p1 <- rep.selected.d04 %>%
  dplyr::select(Stripped.Sequence, Genes, Protein.Group) %>%
  distinct()

pep_prot_map_n2p1 <- pep_prot_map_n2p1 %>%
  filter(Stripped.Sequence %in% rownames(sc.batch_cor.pep))

### Identify proteotypic peptides: unambiguously map to one gene and one canonical
### protein group, with no isoform suffix and no semicolon-separated multi-assignment
proteotypic_flag <- pep_prot_map_n2p1 %>%
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

pep_prot_map_n2p1 <- pep_prot_map_n2p1 %>%
  left_join(proteotypic_flag, by = "Stripped.Sequence")
pep_prot_map_n2p1$dataset_batch <- "n2p1"

### Long-format peptide outputs, parallel to prot_df_n2p1 / prot_df_n2p1.colonly
cellID_donorID_pep <- sample_map %>% dplyr::select(sample, ID) %>% distinct()

# Relative (row+column normalized, ComBat plex + limma Order corrected)
peptide_df_n2p1 <- reshape2::melt(sc.batch_cor.pep)
colnames(peptide_df_n2p1) <- c("Stripped.Sequence", "ID", "Peptide_FC_log2")
peptide_df_n2p1 <- peptide_df_n2p1 %>% left_join(cellID_donorID_pep, by = "ID")
peptide_df_n2p1$dataset_batch <- "n2p1"
peptide_df_n2p1$ID <- paste0(peptide_df_n2p1$ID, "_", peptide_df_n2p1$dataset_batch)

# Absolute-friendly (column-only normalized, ComBat plex + limma Order corrected)
peptide_df_n2p1.colonly <- reshape2::melt(sc.batch_cor.pep.colonly)
colnames(peptide_df_n2p1.colonly) <- c("Stripped.Sequence", "ID", "Peptide_FC_log2")
peptide_df_n2p1.colonly <- peptide_df_n2p1.colonly %>% left_join(cellID_donorID_pep, by = "ID")
peptide_df_n2p1.colonly$dataset_batch <- "n2p1"
peptide_df_n2p1.colonly$ID <- paste0(peptide_df_n2p1.colonly$ID, "_", peptide_df_n2p1.colonly$dataset_batch)

### Save outputs
fwrite(peptide_df_n2p1,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/Peptides_n2p1_for_integration.csv",
       row.names = FALSE)
fwrite(peptide_df_n2p1.colonly,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/Peptides_n2p1_for_abundance_analysis.csv",
       row.names = FALSE)
fwrite(pep_prot_map_n2p1,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/PepProtMap_n2p1.csv",
       row.names = FALSE)



### Save a list of selected cell IDs that passed quality control


IdentificationsList_n2p1 <- rep.selected.d04 %>% distinct()

fwrite(IdentificationsList_n2p1,
       file = "/projects/slavov/LK/ImmPTR/results/final_scripts_output/IdentificationsList_n2p1.csv",
       row.names = FALSE)







