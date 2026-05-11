### MGE original version 1/26/26
### this script will now output data with a hard cutoff of 37 for Cq values with adjusted intensity scores. 
## update 4/27/26 will be using the cutoff of 37 for Dermo and 27 for MSX.

library(dplyr)
setwd("/Users/madelineeppley/Desktop/seadispub")

cq_cutoff <- 37  # samples with Cq > this value will be considered NEG
# qpcr data subset from raw outputs
qpcr <- read.csv("/Users/madelineeppley/Desktop/cviqpcr/20250806_results/all_plates_qpcr_data.csv")

# subset into msx and dermo
msx_sub <- qpcr[qpcr$Fluor == "HEX",]
dermo_sub <- qpcr[qpcr$Fluor == "FAM",]

# filter out inconclusives
msx_conc <- msx_sub[msx_sub$status != "Inconclusive",]
dermo_conc <- dermo_sub[dermo_sub$status != "Inconclusive",]

# keep only rep 1
msx_rep1 <- msx_conc[msx_conc$replicate == 1,]
dermo_rep1 <- dermo_conc[dermo_conc$replicate == 1,]

# filter out rows without clean_ID
msx <- msx_rep1[!is.na(msx_rep1$clean_ID),]
dermo <- dermo_rep1[!is.na(dermo_rep1$clean_ID),]

# for dermo apply cq cutoff for inf status
dermo$status_Cq37 <- ifelse(
  dermo$Cq_Mean > 0 & dermo$Cq_Mean <= cq_cutoff, 
  "Positive", 
  "Negative")

# intensity scale
dermo$intensity_Cq37 <- 0  # default to neg
dermo$intensity_Cq37[dermo$Cq_Mean >= 15 & dermo$Cq_Mean <= 23] <- 3 # high
dermo$intensity_Cq37[dermo$Cq_Mean > 23 & dermo$Cq_Mean <= 25] <- 2 # moderate
dermo$intensity_Cq37[dermo$Cq_Mean > 25 & dermo$Cq_Mean <= 37] <- 1 # low
# everything else (Cq > 37 or < 1) stays at 0 (negative)

# msx cq cutoff
msx$status_Cq37 <- ifelse(
  msx$Cq_Mean > 0 & msx$Cq_Mean <= cq_cutoff, 
  "Positive", 
  "Negative")

# intensity scale
msx$intensity_Cq37 <- 0  # Default to negative
msx$intensity_Cq37[msx$Cq_Mean >= 15 & msx$Cq_Mean <= 23] <- 3 # high
msx$intensity_Cq37[msx$Cq_Mean > 23 & msx$Cq_Mean <= 25] <- 2 # moderate
msx$intensity_Cq37[msx$Cq_Mean > 25 & msx$Cq_Mean <= 37] <- 1 # low

# dermo dups
dermo_dup_ids <- unique(dermo$clean_ID[duplicated(dermo$clean_ID)])
dermo_filtered <- dermo

for (id in dermo_dup_ids) {
  id_rows <- dermo[dermo$clean_ID == id, ]
  if ("Positive" %in% id_rows$status_Cq37) {
    keep_row <- id_rows[id_rows$status_Cq37 == "Positive", ][1, ]
  } else {
    keep_row <- id_rows[which.min(id_rows$cq_difference), ]}
  dermo_filtered <- dermo_filtered[dermo_filtered$clean_ID != id, ]
  dermo_filtered <- rbind(dermo_filtered, keep_row)}

# msx dups
msx_dup_ids <- unique(msx$clean_ID[duplicated(msx$clean_ID)])
msx_filtered <- msx

for (id in msx_dup_ids) {
  id_rows <- msx[msx$clean_ID == id, ]
  if ("Positive" %in% id_rows$status_Cq37) {
    keep_row <- id_rows[id_rows$status_Cq37 == "Positive", ][1, ]
  } else {
    keep_row <- id_rows[which.min(id_rows$cq_difference), ]}
  msx_filtered <- msx_filtered[msx_filtered$clean_ID != id, ]
  msx_filtered <- rbind(msx_filtered, keep_row)}

nrow(dermo_filtered) #646 samples
nrow(msx_filtered) #640 samples

# pop stats
calculate_pop_stats_Cq37 <- function(data, pathogen_name, status_col, intensity_col) {
  populations <- unique(data$ID_SiteDate)
  populations <- populations[!is.na(populations)]
  
  pop_stats <- data.frame(
    Population = character(),
    Pathogen = character(),
    Total_Samples = integer(),
    Positive_Samples = integer(),
    Negative_Samples = integer(),
    Prevalence = numeric(),
    Total_Intensity = numeric(),
    Average_Intensity = numeric(),
    Weighted_Prevalence = numeric(),
    stringsAsFactors = FALSE)
  
  for (pop in populations) {
    pop_data <- data[data$ID_SiteDate == pop, ]
    
    total <- nrow(pop_data)
    positives <- sum(pop_data[[status_col]] == "Positive", na.rm = TRUE)
    negatives <- sum(pop_data[[status_col]] == "Negative", na.rm = TRUE)
    
    conclusives <- positives + negatives
    prevalence <- if (conclusives > 0) positives / conclusives else 0
    
    total_intensity <- sum(pop_data[[intensity_col]], na.rm = TRUE)
    avg_intensity <- if (positives > 0) total_intensity / positives else 0
    weighted_prev <- if (conclusives > 0) total_intensity / conclusives else 0
    
    pop_stats <- rbind(pop_stats, data.frame(
      Population = pop,
      Pathogen = pathogen_name,
      Total_Samples = total,
      Positive_Samples = positives,
      Negative_Samples = negatives,
      Prevalence = round(prevalence, 4),
      Total_Intensity = total_intensity,
      Average_Intensity = round(avg_intensity, 4),
      Weighted_Prevalence = round(weighted_prev, 4),
      stringsAsFactors = FALSE
    ))}
  return(pop_stats)}

# now calculate for both pathogens
dermo_pop_stats <- calculate_pop_stats_Cq37(dermo_filtered, "Dermo", "status_Cq37", "intensity_Cq37")
msx_pop_stats <- calculate_pop_stats_Cq37(msx_filtered, "MSX", "status_Cq37", "intensity_Cq37")

# rbind
population_stats_Cq37 <- rbind(dermo_pop_stats, msx_pop_stats)

# now overall stats
overall_stats_Cq37 <- data.frame(
  Pathogen = c("Dermo", "MSX"),
  Total_Samples = c(nrow(dermo_filtered), nrow(msx_filtered)),
  Positive_Samples = c(
    sum(dermo_filtered$status_Cq37 == "Positive"),
    sum(msx_filtered$status_Cq37 == "Positive")),
  Negative_Samples = c(
    sum(dermo_filtered$status_Cq37 == "Negative"),
    sum(msx_filtered$status_Cq37 == "Negative")),
  stringsAsFactors = FALSE)

overall_stats_Cq37$Prevalence <- round(
  overall_stats_Cq37$Positive_Samples / 
    (overall_stats_Cq37$Positive_Samples + overall_stats_Cq37$Negative_Samples), 4)

# check
print(overall_stats_Cq37)

# save output files
write.csv(population_stats_Cq37, "/Users/madelineeppley/Desktop/cviqpcr/final/population_statistics_Cq37cutoff.csv", row.names = FALSE)
write.csv(overall_stats_Cq37, "/Users/madelineeppley/Desktop/cviqpcr/final/overall_statistics_Cq37cutoff.csv", row.names = FALSE)
write.csv(dermo_filtered, "/Users/madelineeppley/Desktop/cviqpcr/final/dermo_final_Cq37cutoff.csv", row.names = FALSE)

# msx - all southeast
msx_sorted <- msx_pop_stats[order(-msx_pop_stats$Prevalence), ]
print(head(msx_sorted[, c("Population", "Prevalence", "Positive_Samples", "Total_Samples")], 10))

# dermo - mostly atlantic
dermo_sorted <- dermo_pop_stats[order(-dermo_pop_stats$Prevalence), ]
print(head(dermo_sorted[, c("Population", "Prevalence", "Positive_Samples", "Total_Samples")], 10))
