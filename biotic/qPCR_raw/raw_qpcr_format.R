# Raw qPCR data output to formatted file for filtering
# MGE 8/05/2025, last updated 5/11/26

plate_positions <- read.csv("/Users/madelineeppley/Desktop/cviqpcr/qPCR_runs - qpcr_plate_postitions_all.csv")
individuals <- read.csv("/Users/madelineeppley/Desktop/cviqpcr/seascape_samp_full_subset.csv")
extractions <- read.csv("/Users/madelineeppley/Desktop/cviqpcr/SeascapeSamples Extractions - plate_extraction.csv")

extractions$Vial_Label <- extractions$Vial_Label_Check
dim(individuals)
dim(extractions)
str(individuals)
str(extractions)

extraction_counts <- extractions %>%
  count(Vial_Label) %>%
  mutate(is_duplicate = n > 1)

# join this information with extractions
extractions_processed <- extractions %>%
  left_join(extraction_counts, by = "Vial_Label") %>%
  # for dups, set mg_tissue to 20
  mutate(mg_tissue = ifelse(is_duplicate, 20, mg_tissue)) %>%
  # keep 1 entry per Vial_Label
  distinct(Vial_Label, .keep_all = TRUE) %>%
  select(-n, -is_duplicate)

# now join with individuals
inds_weights <- individuals %>%
  left_join(extractions_processed, by = "Vial_Label") %>%
  # these are the only columns i want for final dataset
  select(Vial_Label, ID_SiteDate, clean_ID, mg_tissue)

inds_weights <- inds_weights %>%
  mutate(mg_tissue = ifelse(is.na(mg_tissue), 20, mg_tissue))

str(plate_positions)

inds_weights_plates <- inds_weights %>%
  left_join(plate_positions, by = "clean_ID")

id_counts <- inds_weights_plates %>%
  count(clean_ID, sort = TRUE)

# most frequent samples by clean_ID
head(id_counts, 20) # we want all individuals to appear twice because there are two wells assigned for each qPCR run (duplicates)

na_count <- sum(is.na(inds_weights_plates$qPCR_plate))
na_count #13 inds - these are individuals that were extracted but we don't have qPCR data for yet 
# not all of these individuals will get qPCR data since there was a subset that failed sequencing

# inds that appear only once (n=1) in id_counts
singles_count <- sum(id_counts$n == 1)
singles_count #13 inds, matches the data above 

write.csv(inds_weights_plates, "/Users/madelineeppley/Desktop/cviqpcr/all_sequenced_forqpcr.csv")

# qpcr files added here in order 
qpcr_files <- c(
  "/Users/madelineeppley/Desktop/cviqpcr/20250417_091935_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 1
  "/Users/madelineeppley/Desktop/cviqpcr/20250422_091919_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 2
  "/Users/madelineeppley/Desktop/cviqpcr/20250424_093750_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 3
  "/Users/madelineeppley/Desktop/cviqpcr/20250507_084825_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 4
  "/Users/madelineeppley/Desktop/cviqpcr/20250508_083618_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 5
  "/Users/madelineeppley/Desktop/cviqpcr/20250702_124802_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 6
  "/Users/madelineeppley/Desktop/cviqpcr/20250707_084915_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 7
  "/Users/madelineeppley/Desktop/cviqpcr/20250707_130017_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 8
  "/Users/madelineeppley/Desktop/cviqpcr/20250709_090150_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 9
  "/Users/madelineeppley/Desktop/cviqpcr/20250709_122810_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 10
  "/Users/madelineeppley/Desktop/cviqpcr/20250710_091457_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 11
  "/Users/madelineeppley/Desktop/cviqpcr/20250710_140319_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 12
  "/Users/madelineeppley/Desktop/cviqpcr/20250714_095328_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 13
  "/Users/madelineeppley/Desktop/cviqpcr/20250715_114917_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 14
  "/Users/madelineeppley/Desktop/cviqpcr/20250716_085638_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 15
  "/Users/madelineeppley/Desktop/cviqpcr/20250716_124553_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 16
  "/Users/madelineeppley/Desktop/cviqpcr/20250717_121536_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 17
  "/Users/madelineeppley/Desktop/cviqpcr/20250721_091720_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 18
  "/Users/madelineeppley/Desktop/cviqpcr/20250721_133749_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 19
  "/Users/madelineeppley/Desktop/cviqpcr/20250722_085441_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv", # plate 20
  "/Users/madelineeppley/Desktop/cviqpcr/20250722_115147_CT023639_OYSTER -  Quantification Cq Results.xlsx - 0.csv"  # plate 21
  )

plate_numbers <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

# file paths
#plate_positions_file <- "/Users/madelineeppley/Desktop/cviqpcr/qPCR_runs - qpcr_plate_postitions_all.csv"
#individual_file <- "/Users/madelineeppley/Desktop/cviqpcr/samp_full_subset.csv"
#extractions_file <- "/Users/madelineeppley/Desktop/cviqpcr/SeascapeSamples Extractions - plate_extraction.csv"
# now we will replace all three of the above files with the single file that contains all of the sample IDs, tissue weights, and populations
inds <- "/Users/madelineeppley/Desktop/cviqpcr/all_sequenced_forqpcr.csv"
output_directory <- "/Users/madelineeppley/Desktop/cviqpcr/results"

# constants
cq_thresh <- 2  # threshold for Cq difference
ext_weight <- 20  # 20mg default extraction weight if we can't find the measurement from extraction

# get qpcr data from csv
process_qpcr_data <- function(file_path, plate_number) {
  qpcr_data <- read.csv(file_path)
  # get rid of standards, positive controls, and non-template controls
  filter_pattern <- "Std|NTC|Pos"
  control_rows <- grep(filter_pattern, qpcr_data$Content)
  if (length(control_rows) > 0) {
    qpcr_data <- qpcr_data[-control_rows, ]}
  
  results <- data.frame(
    Plate = plate_number,
    Well = qpcr_data$Well,
    Fluor = qpcr_data$Fluor,
    Content = qpcr_data$Content,
    Cq = as.numeric(qpcr_data$Cq),
    Cq_Mean = as.numeric(qpcr_data$`Cq.Mean`),
    SQ_Mean = as.numeric(qpcr_data$`SQ.Mean`),
    stringsAsFactors = FALSE)
  results$Pathogen <- ifelse(
    grepl("FAM", results$Fluor),
    "Dermo",
    ifelse(grepl("VIC|HEX", results$Fluor), "MSX", "Unknown"))
  return(results)}

map_wells_to_cleanid <- function(qpcr_data, inds_dataframe, plate_number) {
  plate_name <- paste0("qpcr_pl", plate_number)
  plate_positions <- inds_dataframe[inds_dataframe$qPCR_plate == plate_name, ]
  qpcr_data$clean_ID <- NA
  qpcr_data$Vial_Label <- NA
  qpcr_data$ID_SiteDate <- NA
  qpcr_data$mg_tissue <- NA
  mapped_count <- 0
  for (i in 1:nrow(qpcr_data)) {
    well <- qpcr_data$Well[i]
    match_id <- which(plate_positions$cq_well_position == well)
    if (length(match_id) > 0) {
      qpcr_data$clean_ID[i] <- plate_positions$clean_ID[match_id[1]]
      qpcr_data$Vial_Label[i] <- plate_positions$Vial_Label[match_id[1]]
      qpcr_data$ID_SiteDate[i] <- plate_positions$ID_SiteDate[match_id[1]]
      qpcr_data$mg_tissue[i] <- plate_positions$mg_tissue[match_id[1]]
      mapped_count <- mapped_count + 1}}
  
  no_weight <- is.na(qpcr_data$mg_tissue)
    qpcr_data$mg_tissue[no_weight] <- ext_weight}
  return(qpcr_data)}

# look at the replicates - each sample was run 2x
analyze_replicates <- function(qpcr_data, cq_thresh) {
  qpcr_data$sample_id_for_grouping <- ifelse(
    is.na(qpcr_data$clean_ID), 
    paste0("row_", 1:nrow(qpcr_data)),
    qpcr_data$clean_ID)
  # include plate number to make sure replicates come from the same plate
  qpcr_data$sample_pathogen <- paste(qpcr_data$Plate, qpcr_data$sample_id_for_grouping, qpcr_data$Pathogen, sep = "_")
  qpcr_data$replicate <- NA
  qpcr_data$cq_difference <- NA
  qpcr_data$status <- NA
  qpcr_data$need_rerun <- FALSE
  unique_samples <- unique(qpcr_data$sample_pathogen)
  positive_count <- 0
  negative_count <- 0
  inconclusive_count <- 0
  
  for (sample in unique_samples) {
    rows <- which(qpcr_data$sample_pathogen == sample)
    
    if (length(rows) == 0) next
    
    qpcr_data$replicate[rows] <- 1:length(rows)
    
    # calculate the Cq difference and use it to determine inf status
    valid_cq <- !is.na(qpcr_data$Cq[rows]) & qpcr_data$Cq[rows] > 0
    
    if (sum(valid_cq) > 1) {
      # calculate the difference between Cq values
      cq_values <- qpcr_data$Cq[rows][valid_cq]
      cq_diff <- max(cq_values) - min(cq_values)
      qpcr_data$cq_difference[rows] <- cq_diff
      
      # mark as needing a rerun if difference ≥ cq threshold
      qpcr_data$need_rerun[rows] <- cq_diff >= cq_thresh
      
      if (cq_diff < cq_thresh) {
        qpcr_data$status[rows] <- "Positive"
        positive_count <- positive_count + 1
      } else {
        qpcr_data$status[rows] <- "Inconclusive"
        inconclusive_count <- inconclusive_count + 1
      }
    } else if (sum(valid_cq) == 1) {
      qpcr_data$need_rerun[rows] <- TRUE
      qpcr_data$status[rows] <- "Inconclusive"
      inconclusive_count <- inconclusive_count + 1
    } else {
      qpcr_data$need_rerun[rows] <- FALSE
      qpcr_data$status[rows] <- "Negative"
      negative_count <- negative_count + 1
    }
  }
  
  print(positive_count)
  print(negative_count)
  print(inconclusive_count)
  print(sum(qpcr_data$need_rerun))
  return(qpcr_data)}

# function to calculate intensity
calculate_intensity <- function(qpcr_data) {
  qpcr_data$weighted_sq <- NA
  dermo_rows <- qpcr_data$Pathogen == "Dermo"
  msx_rows <- qpcr_data$Pathogen == "MSX"
  if (sum(dermo_rows) > 0) {
    qpcr_data$weighted_sq[dermo_rows] <- qpcr_data$SQ_Mean[dermo_rows] * 6973476000 / qpcr_data$mg_tissue[dermo_rows]
  }
  if (sum(msx_rows) > 0) {
    qpcr_data$weighted_sq[msx_rows] <- qpcr_data$SQ_Mean[msx_rows] * 7810534000 / qpcr_data$mg_tissue[msx_rows]
  }
  qpcr_data$intensity <- 0
  for (i in 1:nrow(qpcr_data)) {
    # apply default 0 or  intensity only to positive samples
    if (is.na(qpcr_data$status[i]) || qpcr_data$status[i] != "Positive") {
      qpcr_data$intensity[i] <- 0
      next
    }
    wsq <- qpcr_data$weighted_sq[i]
    if (is.na(wsq) || wsq == 0) {
      qpcr_data$intensity[i] <- 0
    } else if (wsq < 100) {
      qpcr_data$intensity[i] <- 0.5
    } else if (wsq < 1000) {
      qpcr_data$intensity[i] <- 1
    } else if (wsq < 10000) {
      qpcr_data$intensity[i] <- 2
    } else {
      qpcr_data$intensity[i] <- 3
    }
  }
  
  positive_samples <- sum(qpcr_data$status == "Positive", na.rm = TRUE)
  if (positive_samples > 0) {
    intensity_counts <- table(qpcr_data$intensity[qpcr_data$status == "Positive"])
    print(intensity_counts)
  } else {
    print("no pos. samples found for intensity\n")
  }
  
  return(qpcr_data)}

generate_rerun_list <- function(qpcr_data) {
  need_rerun <- which(qpcr_data$need_rerun)
  if (length(need_rerun) == 0) {
    print("no reruns\n")
    return(data.frame())}
  
  rerun_ids <- unique(qpcr_data$clean_ID[need_rerun])
  rerun_ids <- rerun_ids[!is.na(rerun_ids)]
  
  if (length(rerun_ids) == 0) {
    print("no reruns\n")
    return(data.frame())}
  # make rerun list
  rerun_list <- data.frame(
    clean_ID = character(),
    Vial_Label = character(),
    Plate = integer(),
    Wells = character(),
    Dermo_Status = character(),
    MSX_Status = character(),
    Dermo_Cq = character(),
    MSX_Cq = character(),
    Dermo_Diff = numeric(),
    MSX_Diff = numeric(),
    Rerun_Reason = character(),
    stringsAsFactors = FALSE)
  
  for (id in rerun_ids) {
    id_rows <- which(qpcr_data$clean_ID == id & qpcr_data$need_rerun)
    if (length(id_rows) == 0) next
    unique_plates <- unique(qpcr_data$Plate[id_rows])
    
    for (plate in unique_plates) {
      plate_id_rows <- id_rows[qpcr_data$Plate[id_rows] == plate]
      all_plate_id_rows <- which(qpcr_data$clean_ID == id & qpcr_data$Plate == plate)
      vial_label <- qpcr_data$Vial_Label[all_plate_id_rows[1]]
      wells <- paste(unique(qpcr_data$Well[plate_id_rows]), collapse = ", ")
      
      dermo_rows <- all_plate_id_rows[qpcr_data$Pathogen[all_plate_id_rows] == "Dermo"]
      dermo_status <- if (length(dermo_rows) > 0) qpcr_data$status[dermo_rows[1]] else "Unknown"
      dermo_cq <- if (length(dermo_rows) > 0) 
        paste(round(qpcr_data$Cq[dermo_rows], 2), collapse = ", ") else "NA"
      dermo_diff <- if (length(dermo_rows) > 0 && !is.na(qpcr_data$cq_difference[dermo_rows[1]]))
        qpcr_data$cq_difference[dermo_rows[1]] else NA
      
      msx_rows <- all_plate_id_rows[qpcr_data$Pathogen[all_plate_id_rows] == "MSX"]
      msx_status <- if (length(msx_rows) > 0) qpcr_data$status[msx_rows[1]] else "Unknown"
      msx_cq <- if (length(msx_rows) > 0)
        paste(round(qpcr_data$Cq[msx_rows], 2), collapse = ", ") else "NA"
      msx_diff <- if (length(msx_rows) > 0 && !is.na(qpcr_data$cq_difference[msx_rows[1]]))
        qpcr_data$cq_difference[msx_rows[1]] else NA
      
      rerun_reason <- NA
      if (dermo_status == "Inconclusive" && msx_status == "Inconclusive") {
        rerun_reason <- "Both pathogens inconclusive"
      } else if (dermo_status == "Inconclusive") {
        rerun_reason <- paste("Dermo inconclusive, MSX", tolower(msx_status))
      } else if (msx_status == "Inconclusive") {
        rerun_reason <- paste("MSX inconclusive, Dermo", tolower(dermo_status))
      }
      
      rerun_list <- rbind(rerun_list,
                          data.frame(
                            clean_ID = id,
                            Vial_Label = vial_label,
                            Plate = plate,
                            Wells = wells,
                            Dermo_Status = dermo_status,
                            MSX_Status = msx_status,
                            Dermo_Cq = dermo_cq,
                            MSX_Cq = msx_cq,
                            Dermo_Diff = dermo_diff,
                            MSX_Diff = msx_diff,
                            Rerun_Reason = rerun_reason,
                            stringsAsFactors = FALSE))}}
  return(rerun_list)}

summarize_by_sample <- function(qpcr_data) {
  unique_ids <- unique(qpcr_data$clean_ID)
  unique_ids <- unique_ids[!is.na(unique_ids)]
  
  # results data frame
  results <- data.frame(
    clean_ID = character(),
    Vial_Label = character(),
    ID_SiteDate = character(),
    Plate = integer(),
    Dermo_Status = character(),
    MSX_Status = character(),
    Dermo_Intensity = numeric(),
    MSX_Intensity = numeric(),
    Dermo_Weighted_SQ = numeric(),
    MSX_Weighted_SQ = numeric(),
    mg_tissue = numeric(),
    Result_Combination = character(),
    stringsAsFactors = FALSE)
  
  for (id in unique_ids) {
    id_rows <- which(qpcr_data$clean_ID == id)
    unique_plates <- unique(qpcr_data$Plate[id_rows])
    
    for (plate in unique_plates) {
      plate_id_rows <- id_rows[qpcr_data$Plate[id_rows] == plate]
      
      vial_label <- qpcr_data$Vial_Label[plate_id_rows[1]]
      id_sitedate <- qpcr_data$ID_SiteDate[plate_id_rows[1]]
      mg_tissue_val <- qpcr_data$mg_tissue[plate_id_rows[1]]
      
      dermo_rows <- plate_id_rows[qpcr_data$Pathogen[plate_id_rows] == "Dermo"]
      dermo_status <- if (length(dermo_rows) > 0) qpcr_data$status[dermo_rows[1]] else "Unknown"
      dermo_intensity <- if (length(dermo_rows) > 0) qpcr_data$intensity[dermo_rows[1]] else 0
      dermo_wsq <- if (length(dermo_rows) > 0) qpcr_data$weighted_sq[dermo_rows[1]] else NA
      
      msx_rows <- plate_id_rows[qpcr_data$Pathogen[plate_id_rows] == "MSX"]
      msx_status <- if (length(msx_rows) > 0) qpcr_data$status[msx_rows[1]] else "Unknown"
      msx_intensity <- if (length(msx_rows) > 0) qpcr_data$intensity[msx_rows[1]] else 0
      msx_wsq <- if (length(msx_rows) > 0) qpcr_data$weighted_sq[msx_rows[1]] else NA
      combination <- paste0(dermo_status, "/", msx_status)
      
      results <- rbind(results,
                       data.frame(
                         clean_ID = id,
                         Vial_Label = vial_label,
                         ID_SiteDate = id_sitedate,
                         Plate = plate,
                         Dermo_Status = dermo_status,
                         MSX_Status = msx_status,
                         Dermo_Intensity = dermo_intensity,
                         MSX_Intensity = msx_intensity,
                         Dermo_Weighted_SQ = dermo_wsq,
                         MSX_Weighted_SQ = msx_wsq,
                         mg_tissue = mg_tissue_val,
                         Result_Combination = combination,
                         stringsAsFactors = FALSE))}}
  
  combo_counts <- table(results$Result_Combination)
  print(combo_counts)
  
  return(results)}


calculate_population_stats <- function(sample_summary) {
  populations <- unique(sample_summary$ID_SiteDate)
  populations <- populations[!is.na(populations)]
  
  population_stats <- data.frame(
    Population = character(),
    Pathogen = character(),
    Total_Samples = integer(),
    Positive_Samples = integer(),
    Negative_Samples = integer(),
    Inconclusive_Samples = integer(),
    Prevalence = numeric(),
    Total_Intensity = numeric(),
    Average_Intensity = numeric(),
    Weighted_Prevalence = numeric(),
    stringsAsFactors = FALSE)
  
  for (pop in populations) {
    pop_samples <- sample_summary[sample_summary$ID_SiteDate == pop, ]
    
    for (pathogen in c("Dermo", "MSX")) {
      status_col <- paste0(pathogen, "_Status")
      intensity_col <- paste0(pathogen, "_Intensity")
      
      total <- nrow(pop_samples)
      positives <- sum(pop_samples[[status_col]] == "Positive", na.rm = TRUE)
      negatives <- sum(pop_samples[[status_col]] == "Negative", na.rm = TRUE)
      inconclusives <- sum(pop_samples[[status_col]] == "Inconclusive", na.rm = TRUE)
      
      # calculate prevalence (positives / [positives + negatives])
      conclusives <- positives + negatives
      prevalence <- if (conclusives > 0) positives / conclusives else 0
      
      # calculate intensity metrics
      total_intensity <- sum(pop_samples[[intensity_col]], na.rm = TRUE)
      avg_intensity <- if (positives > 0) total_intensity / positives else 0
      weighted_prev <- if (conclusives > 0) total_intensity / conclusives else 0
      
      population_stats <- rbind(population_stats,
                                data.frame(
                                  Population = pop,
                                  Pathogen = pathogen,
                                  Total_Samples = total,
                                  Positive_Samples = positives,
                                  Negative_Samples = negatives,
                                  Inconclusive_Samples = inconclusives,
                                  Prevalence = prevalence,
                                  Total_Intensity = total_intensity,
                                  Average_Intensity = avg_intensity,
                                  Weighted_Prevalence = weighted_prev,
                                  stringsAsFactors = FALSE))}}
  return(population_stats)}

calculate_overall_stats <- function(sample_summary) {
  overall_stats <- data.frame(
    Pathogen = c("Dermo", "MSX"),
    Total_Samples = 0,
    Positive_Samples = 0,
    Negative_Samples = 0,
    Inconclusive_Samples = 0,
    Prevalence = 0,
    Total_Intensity = 0,
    Average_Intensity = 0,
    Weighted_Prevalence = 0,
    stringsAsFactors = FALSE)
  
  for (i in 1:nrow(overall_stats)) {
    pathogen <- overall_stats$Pathogen[i]
    
    status_col <- paste0(pathogen, "_Status")
    intensity_col <- paste0(pathogen, "_Intensity")
    total <- nrow(sample_summary)
    positives <- sum(sample_summary[[status_col]] == "Positive", na.rm = TRUE)
    negatives <- sum(sample_summary[[status_col]] == "Negative", na.rm = TRUE)
    inconclusives <- sum(sample_summary[[status_col]] == "Inconclusive", na.rm = TRUE)
    conclusives <- positives + negatives
    prevalence <- if (conclusives > 0) positives / conclusives else 0
    total_intensity <- sum(sample_summary[[intensity_col]], na.rm = TRUE)
    avg_intensity <- if (positives > 0) total_intensity / positives else 0
    weighted_prev <- if (conclusives > 0) total_intensity / conclusives else 0
    overall_stats$Total_Samples[i] <- total
    overall_stats$Positive_Samples[i] <- positives
    overall_stats$Negative_Samples[i] <- negatives
    overall_stats$Inconclusive_Samples[i] <- inconclusives
    overall_stats$Prevalence[i] <- prevalence
    overall_stats$Total_Intensity[i] <- total_intensity
    overall_stats$Average_Intensity[i] <- avg_intensity
    overall_stats$Weighted_Prevalence[i] <- weighted_prev}
  return(overall_stats)}

process_plate <- function(file_path, plate_number, individuals_df) {
  cat("\n qpcr plate", plate_number, ":", basename(file_path), "\n")
  
  # load qPCR data
  qpcr_data <- process_qpcr_data(file_path, plate_number)
  # get wells
  qpcr_data <- map_wells_to_cleanid(qpcr_data, individuals_df, plate_number)
  # replicates
  qpcr_data <- analyze_replicates(qpcr_data, cq_thresh)
  # intensity
  qpcr_data <- calculate_intensity(qpcr_data)
  # reruns
  rerun_list <- generate_rerun_list(qpcr_data)
  # summary
  sample_summary <- summarize_by_sample(qpcr_data)
  # results
  results <- list(
    qpcr_data = qpcr_data,
    rerun_list = rerun_list,
    sample_summary = sample_summary)
  return(results)}


main_multi_plates <- function(qpcr_files, plate_numbers, individuals_file, output_directory) {
  individuals_df <- read.csv(individuals_file)
  print(nrow(individuals_df))
  
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)}
  
  all_qpcr_data <- list()
  all_rerun_lists <- list()
  all_sample_summaries <- list()
  
  # run through for each plate
  for (i in 1:length(qpcr_files)) {
    plate_results <- process_plate(qpcr_files[i], plate_numbers[i], individuals_df)
    all_qpcr_data[[i]] <- plate_results$qpcr_data
    all_rerun_lists[[i]] <- plate_results$rerun_list
    all_sample_summaries[[i]] <- plate_results$sample_summary
    write.csv(plate_results$qpcr_data,
              file.path(output_directory, paste0("plate", plate_numbers[i], "_results.csv")),
              row.names = FALSE)
    if (nrow(plate_results$rerun_list) > 0) {
      write.csv(plate_results$rerun_list,
                file.path(output_directory, paste0("plate", plate_numbers[i], "_rerun_list.csv")),
                row.names = FALSE)}
    write.csv(plate_results$sample_summary,
              file.path(output_directory, paste0("plate", plate_numbers[i], "_sample_summary.csv")),
              row.names = FALSE)}
  
  # combine all data
  combined_qpcr_data <- do.call(rbind, all_qpcr_data)
  write.csv(combined_qpcr_data,
            file.path(output_directory, "all_plates_qpcr_data.csv"),
            row.names = FALSE)
  
  # combine rerun lists
  all_reruns <- all_rerun_lists[sapply(all_rerun_lists, nrow) > 0]
  if (length(all_reruns) > 0) {
    combined_reruns <- do.call(rbind, all_reruns)
    write.csv(combined_reruns,
              file.path(output_directory, "all_plates_rerun_list.csv"),
              row.names = FALSE)
  } else {
    combined_reruns <- data.frame()}
  
  # combine sample summaries
  combined_sample_summary <- do.call(rbind, all_sample_summaries)
  write.csv(combined_sample_summary,
            file.path(output_directory, "all_plates_sample_summary.csv"),
            row.names = FALSE)
  
  # calculate population-level results
  population_stats <- calculate_population_stats(combined_sample_summary)
  write.csv(population_stats,
            file.path(output_directory, "all_plates_population_statistics.csv"),
            row.names = FALSE)
  
  # calculate overall statistics on combined data
  overall_stats <- calculate_overall_stats(combined_sample_summary)
  write.csv(overall_stats,
            file.path(output_directory, "all_plates_overall_statistics.csv"),
            row.names = FALSE)
  
  # return all results
  results <- list(
    all_qpcr_data = all_qpcr_data,
    combined_qpcr_data = combined_qpcr_data,
    all_rerun_lists = all_rerun_lists,
    combined_reruns = combined_reruns,
    all_sample_summaries = all_sample_summaries,
    combined_sample_summary = combined_sample_summary,
    population_stats = population_stats,
    overall_stats = overall_stats)
  return(results)}

# all sequenced data
individuals_file <- "/Users/madelineeppley/Desktop/cviqpcr/all_sequenced_forqpcr.csv"

# output directory - can update each time with the new date to avoid overwriting
output_directory <- "/Users/madelineeppley/Desktop/cviqpcr/20250806_results"

# run
results <- main_multi_plates(
  qpcr_files = qpcr_files,
  plate_numbers = plate_numbers,
  individuals_file = individuals_file,
  output_directory = output_directory)

# summary
# overall stats
print(results$overall_stats)

# pop level stats
if (nrow(results$population_stats) > 0) {
  pop_summary <- aggregate(
    cbind(Prevalence, Average_Intensity, Weighted_Prevalence) ~ Pathogen, 
    data = results$population_stats, 
    FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))
  print(pop_summary)
  
  # some stats look at dermo pops
  dermo_pops <- results$population_stats[results$population_stats$Pathogen == "Dermo", ]
  dermo_pops <- dermo_pops[order(dermo_pops$Prevalence, decreasing = TRUE), ]
  print(head(dermo_pops[, c("Population", "Prevalence", "Average_Intensity", "Weighted_Prevalence")], 5))
  msx_pops <- results$population_stats[results$population_stats$Pathogen == "MSX", ]
  msx_pops <- msx_pops[order(msx_pops$Prevalence, decreasing = TRUE), ]
  print(head(msx_pops[, c("Population", "Prevalence", "Average_Intensity", "Weighted_Prevalence")], 5))}

# sample summary
print(nrow(results$combined_sample_summary)

# samples with disease combinations
status_counts <- table(results$combined_sample_summary$Result_Combination)
print(status_counts)

# intensity distributions
dermo_intensity <- table(results$combined_sample_summary$Dermo_Intensity)
msx_intensity <- table(results$combined_sample_summary$MSX_Intensity)
print(dermo_intensity)
print(msx_intensity)

# rerun sample information
print(nrow(results$combined_reruns))
