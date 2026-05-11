## MGE filter qpcr data to produce population-level stats
## last updated 5/11/26

# load in all qpcr data from processing script
qpcr <- read.csv("/Users/madelineeppley/Desktop/cviqpcr/20250806_results/all_plates_qpcr_data.csv")

# subset into msx and dermo
msx_sub <- qpcr[qpcr$Fluor == "HEX",]
dermo_sub <- qpcr[qpcr$Fluor == "FAM",]

# filter out any inconclusive results
inconc <- c("Inconclusive")
msx_conc <- msx_sub[!(msx_sub$status %in% inconc),]
dermo_conc <- dermo_sub[!(dermo_sub$status %in% inconc),]

# now i have 2 replicates for each sample, need to remove one replicate first
rep_2 <- c("2")
msx_rep1 <- msx_conc[!(msx_conc$replicate %in% rep_2),]
dermo_rep1 <- dermo_conc[!(dermo_conc$replicate %in% rep_2),]

msx <- msx_rep1 %>% filter(!is.na(clean_ID))
dermo <- dermo_rep1 %>% filter(!is.na(clean_ID))

dim(msx) #682
dim(dermo) #656

# deal with dermo first, make Cq intensity scale & filter out any replicates
hist(dermo$Cq_Mean)

# make an intensity scale based off of the Cq values - we can connect this back to copy # of parasite DNA
dermo$cq_intensity[dermo$Cq_Mean >= 35] <- 0.5
dermo$cq_intensity[dermo$Cq_Mean >= 29 & dermo$Cq_Mean <= 35] <- 1
dermo$cq_intensity[dermo$Cq_Mean >= 23 & dermo$Cq_Mean <= 29] <- 2
dermo$cq_intensity[dermo$Cq_Mean <= 23 & dermo$Cq_Mean >= 15] <- 3
dermo$cq_intensity[dermo$Cq_Mean < 1] <- 0

hist(dermo$cq_intensity)

# get our repeated IDs - these are samples that had >1 definitive result bc of re-runs/tech reps 
dermo_dups <- dermo[duplicated(dermo$clean_ID), ] # 10 samples re-run

# let's drop the repeat that is negative, if one is negative and one is positive
# if both positive, let's drop the repeat that has the higher difference in Cq between reps
dim(dermo) #656

dermo_filtered <- dermo
dermo_dup_ids <- unique(dermo_dups$clean_ID)

for (id in dermo_dup_ids) {
  id_rows <- dermo[dermo$clean_ID == id, ]
  if ("Positive" %in% id_rows$status) {
    keep_row <- id_rows[id_rows$status == "Positive", ][1, ]
  } else {
    keep_row <- id_rows[which.min(id_rows$cq_difference), ]}
  dermo_filtered <- dermo_filtered[dermo_filtered$clean_ID != id, ]
  dermo_filtered <- rbind(dermo_filtered, keep_row)}

dim(dermo_filtered) #646, successfully removed the dups


# now MSX make Cq intensity scale & filter out any replicates
hist(msx$Cq_Mean)

# next make an intensity scale based off of the Cq values to connect this back to copy # of parasite DNA
msx$cq_intensity[msx$Cq_Mean >= 35] <- 0.5
msx$cq_intensity[msx$Cq_Mean >= 29 & msx$Cq_Mean <= 35] <- 1
msx$cq_intensity[msx$Cq_Mean >= 23 & msx$Cq_Mean <= 29] <- 2
msx$cq_intensity[msx$Cq_Mean <= 23 & msx$Cq_Mean >= 15] <- 3
msx$cq_intensity[msx$Cq_Mean < 1] <- 0

hist(msx$cq_intensity)

# get repeated IDs, these are samples that had >1 definitive result bc of re-runs/tech reps 
msx_dups <- msx[duplicated(msx$clean_ID), ] # 28, but some of these were run multiple times (look at rep # - we actually have 3 extra samples)

msx_filtered <- msx
msx_dup_ids <- unique(msx_dups$clean_ID) # actually 25 samples unique

dim(msx_filtered) # 682

for (id in msx_dup_ids) {
  id_rows <- msx[msx$clean_ID == id, ]
  if ("Positive" %in% id_rows$status) {
    keep_row <- id_rows[id_rows$status == "Positive", ][1, ]
  } else {
    keep_row <- id_rows[which.min(id_rows$cq_difference), ]}
  msx_filtered <- msx_filtered[msx_filtered$clean_ID != id, ]
  msx_filtered <- rbind(msx_filtered, keep_row)}
dim(msx_filtered) #639

###################################################
# now export data and join to see which samples we have in common
# i will just keep the data that we want, and re-name everything to keep the disease info

# the data that we want is: qPCR plate, Pathogen, Status, clean_ID, Cq_Mean, ID_SiteDate, Cq_Intensity
dermo_final <- dermo_filtered[, c("clean_ID", "ID_SiteDate", "Pathogen", "status", "cq_intensity", "Cq_Mean", "SQ_Mean")]
msx_final <- msx_filtered[, c("clean_ID", "ID_SiteDate", "Pathogen", "status", "cq_intensity", "Cq_Mean", "SQ_Mean")]

names(dermo_final) <- c("clean_ID", "ID_SiteDate", "Pathogen_Dermo", "status_dermo", "intensity_dermo", "Cq_Mean_Dermo", "SQ_Mean_Dermo")
names(msx_final) <- c("clean_ID", "ID_SiteDate", "Pathogen_MSX", "status_msx", "intensity_msx", "Cq_Mean_MSX", "SQ_Mean_MSX")

write.csv(dermo_final, "/Users/madelineeppley/Desktop/cviqpcr/final/dermo_final.csv")
write.csv(msx_final, "/Users/madelineeppley/Desktop/cviqpcr/final/msx_final.csv")

# merge together to see which samples we have complete data for
qpcr_merge <- merge(dermo_final, msx_final)
write.csv(qpcr_merge, "/Users/madelineeppley/Desktop/cviqpcr/final/complete_samples.csv")

#### read back in for some further processing
dermo_final <- read.csv("/Users/madelineeppley/Desktop/cviqpcr/final/dermo_final.csv")
msx_final <- read.csv("/Users/madelineeppley/Desktop/cviqpcr/final/msx_final.csv")

# merge together to see which samples we have complete data for
qpcr_merge <- read.csv ("/Users/madelineeppley/Desktop/cviqpcr/final/complete_samples.csv")


# remove laguna madre for pop counts
clean_sites <- qpcr_merge[!grepl("LM_TX_2022-7-10", qpcr_merge$ID_SiteDate), ]
site_counts <- table(clean_sites$ID_SiteDate)
range(site_counts) # minimum 10, max 19
quantile(site_counts) # 10, 12, 15, 17, 19
mean(site_counts) # 14.53
hist(site_counts)

clean_sites_msx <- msx_final[!grepl("LM_TX_2022-7-10", msx_final$ID_SiteDate), ]
site_counts_msx <- table(clean_sites_msx$ID_SiteDate)
range(site_counts_msx) # minimum 12, max 20
quantile(site_counts_msx) # 12, 15, 16, 18, 20
mean(site_counts_msx) # 16.35
hist(site_counts_msx)

clean_sites_dermo <- dermo_final[!grepl("LM_TX_2022-7-10", dermo_final$ID_SiteDate), ]
site_counts_dermo <- table(clean_sites_dermo$ID_SiteDate)
range(site_counts_dermo) # minimum 12, max 20
quantile(site_counts_dermo) # 12, 14.5, 17, 18.5, 20
mean(site_counts_dermo) # 16.53
hist(site_counts_dermo)


dermo_pop1 <- aggregate(intensity_dermo ~ ID_SiteDate, FUN = mean, data = clean_sites_dermo)
dermo_pop2 <- aggregate(Cq_Mean_Dermo ~ ID_SiteDate, FUN = mean, data = clean_sites_dermo)
dermo_pop3 <- aggregate(SQ_Mean_Dermo ~ ID_SiteDate, FUN = mean, data = clean_sites_dermo)

dermo_pop_final <- cbind(dermo_pop1, dermo_pop2$Cq_Mean_Dermo, dermo_pop3$SQ_Mean_Dermo)
colnames(dermo_pop_final) <- c("ID_SiteDate", "intensity_dermo", "Cq_Mean_Dermo", "SQ_Mean_Dermo")
write.csv(dermo_pop_final, "/Users/madelineeppley/Desktop/cviqpcr/final/dermo_pop_final.csv")


msx_pop1 <- aggregate(intensity_msx ~ ID_SiteDate, FUN = mean, data = clean_sites_msx)
msx_pop2 <- aggregate(Cq_Mean_MSX ~ ID_SiteDate, FUN = mean, data = clean_sites_msx)
msx_pop3 <- aggregate(SQ_Mean_MSX ~ ID_SiteDate, FUN = mean, data = clean_sites_msx)

msx_pop_final <- cbind(msx_pop1, msx_pop2$Cq_Mean_MSX, msx_pop3$SQ_Mean_MSX)
colnames(msx_pop_final) <- c("ID_SiteDate", "intensity_msx", "Cq_Mean_MSX", "SQ_Mean_MSX")
write.csv(msx_pop_final, "/Users/madelineeppley/Desktop/cviqpcr/final/msx_pop_final.csv")
