########################################################################################
## Load these functions prior to processing raw data or analyzing processed data
########################################################################################





############ Load functions:

## Note:
# Mapping single-cell meta data acquired from the CellenONE during sample prep requires the QuantQC package
# Certain functions from the QuantQC package have been edited to accomodate differences in the data structure for analysis here versus data structures created using newer versions of QuantQC/nPOP
# Mapping functions from QuantQC have evolved over time, and previous versions of nPOP require distinct mapping functions as the sample prep has been improved
# Certain datasets analyzed here are thus mapped using slightly different code
# For the "n2p1 and n2p2" datasets, specific mapping functions will be included with those scripts alone, which should be run for processing
# I suggest the user loads all of these functions if the intent is to reproduce part of all of the analysis
# Load these functions prior to processing each of the raw datasets
# Comments will describe the purpose of each function

## Important: certain mapping functions are exclusive to certain single-cell proteomic datasets that were prepared using earlier versions of nPOP
## They are instead included with the preprocessing R scripts specific to those datasets (n2p1, n2p2)





### This function generates a QuantQC S4 object but does not apply any filtering of the data yet
### A dataframe containing all the required raw data will be extracted from generated object in a following step of analysis
DIANN_to_QQC2 <- function (data_path, linker_path, plex, carrier = F) 
{
  linker <- read.csv(linker_path)
  # linker$Order <- 1:nrow(linker)
  columns_to_read <- c("Genes", "Run", "Q.Value", "PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "RT", 
                       "Precursor.Id", "Stripped.Sequence", "Precursor.Mz", 
                       "Precursor.Charge", "Precursor.Quantity", "Ms1.Area", 
                       "Protein.Group", "Protein.Names", "Translated.Q.Value", "Channel.Q.Value")
  if (dir.exists(data_path) == F) {
    Raw_data <- data.table::fread(data_path, select = columns_to_read)
  }
  if (dir.exists(data_path)) {
    if (str_sub(data_path, -1) != "/") {
      data_path <- paste0(data_path, "/")
    }
    raw_files <- list.files(data_path)
    for (i in 1:length(raw_files)) {
      if (i == 1) {
        Raw_data <- data.table::fread(paste0(data_path, 
                                             raw_files[i]), select = columns_to_read)
      }
      if (i > 1) {
        data_next <- data.table::fread(paste0(data_path, 
                                              raw_files[i]), select = columns_to_read)
        data_list <- list(Raw_data, data_next)
        Raw_data <- data.table::rbindlist(data_list)
      }
    }
  }
  Raw_data <- as.data.frame(Raw_data)
  # Raw_data <- Raw_data %>% filter(Lib.PG.Q.Value < 0.01)
  Raw_data <- Raw_data %>% filter(Run %in% linker$Run)
  Raw_data <- Raw_data %>% left_join(linker, by = c("Run"))
  Raw_data$seqcharge <- paste0(Raw_data$Stripped.Sequence, 
                               Raw_data$Precursor.Charge)
  Raw_data <- Raw_data %>% filter(Protein.Group != "")
  Raw_data$plex <- substr(Raw_data$Precursor.Id[1:nrow(Raw_data)], 
                          10, 10)
  Raw_data$ID <- paste0(Raw_data$Well, Raw_data$plate, ".", 
                        Raw_data$plex)
  Raw_data$File.Name <- Raw_data$ID
  Raw_data$uq <- paste0(Raw_data$File.Name, Raw_data$Protein.Group, 
                        Raw_data$seqcharge)
  Raw_data <- Raw_data %>% distinct(uq, .keep_all = T)
  Raw_data$uq <- NULL
  if (carrier == F) {
    QQC <- new("QQC", raw_data = Raw_data, ms_type = "DIA", 
               meta.data = linker, misc = list(plex = plex))
  }
  if (carrier == T) {
    QQC <- new("QQC", raw_data = Raw_data, meta.data = linker, 
               ms_type = "DIA_C", misc = list(plex = plex))
  }
  return(QQC)
}





### This function maps mTRAQ mass tag chanels onto the raw single-cell data
cellXpeptide_DIA2 <- function (QQC, TQVal, chQVal) 
{
  Raw_data <- QQC@raw_data
  plex <- QQC@misc[["plex"]]
  type <- QQC@ms_type
  if (plex == 2 & type == "DIA_C") {
    plex_used <- c(0, 4, 8)
  }
  else if (plex == 2 & type == "DIA") {
    plex_used <- c(0, 4)
  }
  else if (plex == 3) {
    plex_used <- c(0, 4, 8)
  }
  else {
    return("plex not valid")
  }
  Raw_data <- Raw_data %>% filter(plex %in% plex_used)
  Raw_data_filt <- Raw_data %>% filter(Channel.Q.Value <= chQVal)
  Raw_data_filt <- Raw_data_filt %>% filter(Translated.Q.Value <= 
                                              TQVal)
  Raw_data_lim_filt <- Raw_data_filt %>% dplyr::select(Protein.Group, 
                                                       seqcharge, Ms1.Area, File.Name)
  Raw_data_lim.d_filt <- reshape2::dcast(Raw_data_lim_filt, 
                                         Protein.Group + seqcharge ~ File.Name, value.var = "Ms1.Area")
  Raw_data_lim.d_filt[Raw_data_lim.d_filt == 0] <- NA
  Raw_data_lim_NF <- Raw_data %>% dplyr::select(Protein.Group, 
                                                seqcharge, Ms1.Area, File.Name)
  Raw_data_lim.d_NF <- reshape2::dcast(Raw_data_lim_NF, Protein.Group + 
                                         seqcharge ~ File.Name, value.var = "Ms1.Area")
  Raw_data_lim.d_NF <- Raw_data_lim.d_NF %>% filter(seqcharge %in% 
                                                      Raw_data_lim.d_filt$seqcharge)
  if (type == "DIA_C") {
    Raw_data_lim.d_filt <- DIA_carrier_norm(Raw_data_lim.d_filt, 
                                            8, plex_used)
    Raw_data_lim.d_NF <- Raw_data_lim.d_NF %>% dplyr::select(colnames(Raw_data_lim.d_filt))
  }
  peptide_protein_map <- Raw_data_lim.d_filt %>% dplyr::select(seqcharge, 
                                                               Protein.Group)
  colnames(peptide_protein_map) <- c("seqcharge", "Protein")
  Raw_data_lim.d_filt <- as.matrix(Raw_data_lim.d_filt[, 3:ncol(Raw_data_lim.d_filt)])
  Raw_data_lim.d_NF <- as.matrix(Raw_data_lim.d_NF[, 3:ncol(Raw_data_lim.d_NF)])
  Raw_data_lim.d_NF[Raw_data_lim.d_NF == 0] <- NA
  pep_mask <- is.na(Raw_data_lim.d_NF) == F
  data_matricies <- new("matricies_DIA", peptide = Raw_data_lim.d_filt, 
                        peptide_mask = pep_mask, peptide_protein_map = peptide_protein_map)
  QQC@matricies <- data_matricies
  QQC@misc[["ChQ"]] <- chQVal
  return(QQC)
}





### This function uses the above function "cellXpeptide_DIA2" to complete raw data mapping, but does not apply any filters yet
cellXpeptide2 <- function (QQC, TQVal = 1, chQVal = 1) 
{
  if (QQC@ms_type == "DIA" | QQC@ms_type == "DIA_C") {
    QQC <- cellXpeptide_DIA2(QQC, TQVal, chQVal)
  }
  if (QQC@ms_type == "DDA") {
    QQC <- TMT_Reference_channel_norm(QQC)
  }
  return(QQC)
}





### This function maps the raw data to the single-cell meta data from the CellenONE aquired during sample prep
link_cellenONE_Raw <- function(QQC,cells_file){
  linker <- QQC@meta.data
  if(QQC@ms_type == 'DDA'){
    cellenOne_data <- analyzeCellenONE_TMT(cells_file,QQC@misc[['plex']])
  }
  if(QQC@ms_type == 'DIA' | QQC@ms_type =='DIA_C'){
    cellenOne_data <- analyzeCellenONE_mTRAQ(cells_file,QQC@misc[['plex']])
  }
  peptide_data <- QQC@matricies@peptide
  # Get list of unique cell IDs
  cellID <- colnames(peptide_data)
  cellID <- as.data.frame(cellID)
  colnames(cellID) <- 'ID'
  if(sum(colnames(cellenOne_data)=='Intensity') == 0){
    cellenOne_data_small <- cellenOne_data %>% dplyr::select(any_of(c('ID','diameter','sample','label','injectWell','plate')))
    cellenOne_data_small <- as.data.frame(cellenOne_data_small)
  }
  if(sum(colnames(cellenOne_data)=='Intensity') == 1){
    cellenOne_data_small <- cellenOne_data %>% dplyr::select(any_of(c('ID','diameter','sample','label','injectWell','plate','Intensity','Stain_Diameter')))
    cellenOne_data_small <- as.data.frame(cellenOne_data_small)
  }
  #cellenOne_data_small <- cellenOne_data_small %>% filter(injectWell %in%  QQC@raw_data$Well)
  #cellenOne_data_small <- cellenOne_data_small %>% filter(plate %in%  QQC@raw_data$plate)
  cellID <- cellID %>% dplyr::left_join(cellenOne_data_small,by = c('ID'))
  cellID$sample[is.na(cellID$sample)==T] <- 'neg'
  # cellID$prot_total <- log2(colSums(peptide_data[,1:ncol(peptide_data)],na.rm = T))
  QQC@cellenONE.meta <- cellenOne_data
  cellID$WP <- paste0(cellID$plate,cellID$injectWell)
  cellID$plate <- NULL
  linker$WP <- paste0(linker$plate,linker$Well)
  cellID <- cellID %>% left_join(linker, by = c('WP'))
  cellID$WP <- NULL
  QQC@meta.data <- cellID
  return(QQC)
}





### This function collects relevant batch statistics to evaluate general LC-MS performance and raw data quality
collect_raw_statistics <- function(data, map) {
  raw_statistics <- data %>%
  dplyr::group_by(File.Name) %>%
  summarise(total_MS1_signal = sum(Ms1.Area, na.rm = TRUE),
            median_MS1_signal = median(Ms1.Area, na.rm = TRUE),
            total_MS2_signal = sum(Precursor.Quantity, na.rm = TRUE),
            total_MS1MS2_signal_ratio = (sum(Ms1.Area, na.rm = TRUE))/(sum(Precursor.Quantity, na.rm = TRUE)),
            precursor_counts = n_distinct(seqcharge),
            proteingroup_counts = n_distinct(Protein.Group)) %>%
  dplyr::left_join(map, by=c("File.Name"="File.Name"))
return(raw_statistics)
}





### This function performs column and then row normalization on a data matrix
### It scales columns to the row medians and rows to have a mean of 1 for better comparability across samples
### It optionally performs log2 transformation
Normalize_saad <- function(dat,log = 'no'){
  dat <- as.matrix(dat)
  dat[dat==0] <- NA
  refVec <- matrixStats::rowMedians(x = dat, cols = 1:ncol(dat), na.rm = T)
  for(k in 1:ncol(dat)){
    dat[,k]<-dat[,k] * median(refVec/dat[,k], na.rm = T)
  }
  for(k in 1:nrow(dat)){
    dat[k,]<-dat[k,]/mean(dat[k,], na.rm = T)
  }
  if(log == 'yes'){
    dat <- log2(dat)
  }
  return(dat)
}





### This function performs column-only normalization on a data matrix
### It scales columns to the row medians for robustness
Normalize_saad_column_only <- function(dat, log = 'no') {
  dat <- as.matrix(dat)
  dat[dat == 0] <- NA  # Replace 0s with NA to avoid skewing medians
  refVec <- matrixStats::rowMedians(dat, na.rm = TRUE)
  for (k in 1:ncol(dat)) {
    dat[, k] <- dat[, k] * median(refVec / dat[, k], na.rm = TRUE)
  }
  if (log == 'yes') {
    dat <- log2(dat)
  }
  return(dat)
}





### This function performs column and then row normalization on a data matrix
### It scales columns to the row medians and rows to have a mean of 1 for better comparability across samples
### This is used when the data is already log2 transformed
Normalize_saad_log <- function(dat){
  refVec <- matrixStats::rowMedians(x = dat, cols = 1:ncol(dat), na.rm = T)
  for(k in 1:ncol(dat)){
    dat[,k]<-dat[,k] + median(refVec - dat[,k], na.rm = T)
  }
  for(k in 1:nrow(dat)){
    dat[k,]<-dat[k,]-mean(dat[k,], na.rm = T)
  }
  return(dat)
}



### This function performs column-only normalization on a data matrix
### It scales columns to the row medians for robustness
### This is used when the data is already log2 transformed
Normalize_saad_log_column_only <- function(dat){
  refVec <- matrixStats::rowMedians(dat, na.rm = TRUE)
  for(k in 1:ncol(dat)){
    dat[,k] <- dat[,k] + median(refVec - dat[,k], na.rm = TRUE)
  }
  return(dat)
}





### This function is used for calculating the median protein CV per single cell
### For each protein in each cell, it calculates the CV of normalized precursor intensities and takes the median across all proteins
### This is used for filtering out low quality cells
calculate_median_CV <- function(normalized_data) {
  # Create an empty data frame to store the results
  results <- data.frame(ID = character(), Median_CV = numeric(), stringsAsFactors = FALSE)
  # Loop through each unique entry in normalized_data$ID
  for (id in unique(normalized_data$ID)) {
    # Subset the data for the current ID
    id_data <- normalized_data[normalized_data$ID == id, ]
    # Initialize a vector to store the CVs for each Protein.Group
    cv_list <- numeric()
    # Loop through each unique Protein.Group within the current ID
    for (protein in unique(id_data$Protein.Group)) {
      # Subset the data for the current Protein.Group
      protein_data <- id_data[id_data$Protein.Group == protein, ]
      # Filter out rows with NA in Ms1.Area
      valid_data <- protein_data$Ms1.Area[!is.na(protein_data$Ms1.Area)]
      # Check if there are at least 3 observations
      if (length(valid_data) >= 3) {
        # Calculate the standard deviation and mean
        sd_value <- sd(valid_data, na.rm = TRUE)
        mean_value <- mean(valid_data, na.rm = TRUE)
        # Calculate the CV (coefficient of variation)
        cv <- sd_value / mean_value
        # Add the CV to the list
        cv_list <- c(cv_list, cv)
      }
    }
    # Calculate the median CV for the current ID
    if (length(cv_list) > 0) {
      median_cv <- median(cv_list, na.rm = TRUE)
    } else {
      median_cv <- NA  # If no CVs were calculated, return NA
    }
    # Append the result to the results data frame
    results <- rbind(results, data.frame(ID = id, Median_CV = median_cv, stringsAsFactors = FALSE))
  }
  return(results)
}





### This function calculates the median absolute deviation using a threshold of 1 x the MAD as an upper or lower boundary
### This is used to systematically filtering out low quality cells 
calculate_MAD_threshold <- function(df, col, thresh = 1) {
  # Extract the column of interest
  values <- df[[col]]
  # Calculate the median and MAD
  median_val <- median(values, na.rm = TRUE)
  mad_val <- mad(values, constant = 1.4826, na.rm = TRUE)
  # Define the MAD-based thresholds
  MAD_UpperThreshold <- median_val + (thresh * mad_val)
  MAD_LowerThreshold <- median_val - (thresh * mad_val)
  # Return the thresholds
  return(list(MAD_UpperThreshold = MAD_UpperThreshold, MAD_LowerThreshold = MAD_LowerThreshold))
}





#### This function performs imputation by kNN
hknn<-function(dat, k){
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<- 1-as.matrix(cor((dat), use="pairwise.complete.obs"))
  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  # For each column in the data... 
  for(X in cnames){
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    # Take the values in the column of interest
    vec<-dat[, X]
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      #print(length(closest.columns))
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
      }
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ] )
      }
    }
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
  }
  return(dat.imp)
}




