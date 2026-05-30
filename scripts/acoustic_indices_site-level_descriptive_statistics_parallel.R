# Purpose: Calculates descriptive statistics on acoustic indices
# Matthew Clark, March 27, 2026

library(dplyr)
library(doParallel)
library(foreach)
library(tidyr)
library(e1071)

# Input directory
inDir <- "G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/github"

# Output file
outFile <- file.path(inDir, "data/acoustic_indices_minutes_dawn-dusk_AGI90_descriptive_statistics.csv")

# Desired time period
# NA = use all; "PredawnQuiet", "DawnChorus","Midday", "DuskChorus", "Nighttime"

#period <- NA # none
periods <- c("DawnChorus","DuskChorus") # dawn & dusk chorus only

# Load acoustic indices by campaign
df1 <- read.csv(file.path(inDir, "data/bioscape_acoustic_indices_campaign1_231022.csv"))
df2 <- read.csv(file.path(inDir, "data/bioscape_acoustic_indices_campaign2_240203.csv"))
df <- rbind(df1,df2) %>% select(-c(path,YYYY,MM,DD,hh,mm,DDhh))
remove(df1,df2)

# Filter out minutes with high level of geophony, interference or anthropophony?
agi_filter <- 1 # 1 = Yes, 0 = No
agi_thresh <- 90 # minutes with over this amount of AGI combined will be filtered
agbi_minutes <- file.path(inDir, "data/wildmon_minute-level_species_taxon_abgi.csv") # file with minutes summary

# Descriptive statistics function
getStats <- function(dd) {
  stats <- c(
    mean = mean(dd, na.rm = TRUE),
    median = median(dd, na.rm = TRUE),
    sd = sd(dd, na.rm = TRUE),
    iqr = IQR(dd, na.rm = TRUE),
    q25 = as.numeric(quantile(dd, probs = 0.25, na.rm = TRUE)),
    q75 = as.numeric(quantile(dd, probs = 0.75, na.rm = TRUE)),
    skewness = skewness(dd, na.rm = TRUE, type = 2),
    kurtosis = kurtosis(dd, na.rm = TRUE, type = 2)
  )
  return(stats)
}

# Custom time binning
binTimeOfDay <- function(dd){
  hrb <- ifelse(dd >= 000 & dd < 500, "PredawnQuiet",  # Pre-dawn/Quiet
                ifelse(dd >= 500 & dd < 900, "DawnChorus",  # Dawn chorus
                       ifelse(dd >= 900 & dd < 1600, "Midday",  # Midday
                              ifelse(dd >= 1600 & dd < 2100, "DuskChorus", # Dusk chorus
                                     ifelse(dd >= 2100 & dd < 2400, "Nighttime", NA)))))  # Night
  return(hrb)
}

# Extract site ID
extract_site <- function(filename) substr(basename(filename), 1, 15)

# Get files list
files <- df$file

# List of acoustic indices to summarize
cn <- c("ACI", "ADI", "AEI", "NDSI", "NDSI_A", "NDSI_B",
        "BI", "H", "Ht", "Hs", "M", "R", "sfm", "rugo", "zcr_mean")

# Sites
sites <- unique(sapply(files, extract_site))

# Read in AGI minutes data (1 = Yes)
if (agi_filter == 1){
  df_agi <- read.csv(agbi_minutes) %>% select(file, Geophony, Anthropophony, Interference)
  df_agi <- df_agi %>%
    mutate(across(c(Geophony, Anthropophony, Interference), ~replace_na(as.integer(.), 0)),
           AGI_percent = pmin(59, Geophony + Anthropophony + Interference) / 59 * 100)
  df_agi <- df_agi %>% select(file, AGI_percent)
} 

# Parallel setup
#numCores <- 30
numCores <- parallel::detectCores() - 2
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Main loop per site

stats_vectors <- list()

stats_vectors <- foreach(s = sites, .packages = c("dplyr", "tidyr", "e1071")) %dopar% {

  tryCatch({
    siteFiles <- files[grep(s, files)]
    if (length(siteFiles) == 0) return(NULL)
    df_site <- df %>% filter(file %in% siteFiles)
    
    # Filter high AGI minutes (1 = Yes)
    if (agi_filter == 1){
      
      # Get file name
      df_site$file <- paste0(df_site$file,".WAV")
      
      # Combine with AGI percent data
      df_site <- left_join(df_site,df_agi, by="file")
      
      # Filter those minutes with AGI percent less than the desired threshold
      df_site <- df_site %>% filter(AGI_percent <= agi_thresh)   
      
    }
    
    # Filter to desired time period
    df_site$time_period <- binTimeOfDay(df_site$hhmm)
    if (sum(is.na(periods)) == 0) {
      df_site <- df_site %>% filter(time_period %in% periods)
    }

    # Combine into matrix
    df_site <- df_site[,cn]
    if (nrow(df_site) < 10) return(NULL)
    
    # Calculate statistics
    stats <- apply(df_site, 2, getStats)
    
    # Format into a data frame
    stats_df <- as.data.frame(t(stats))
    stats_df$feature <- rownames(stats_df)
    stats_df$feature <- cn
    rownames(stats_df) <- NULL
    
    # Reshape: make feature_stat columns
    stats_wide <- stats_df %>%
      pivot_longer(cols = -feature, names_to = "stat", values_to = "value") %>%
      mutate(feature_stat = paste0(feature, "_", stat)) %>%
      select(feature_stat, value) %>%
      pivot_wider(names_from = feature_stat, values_from = value)
    
    # Add site name column
    stats_wide$SiteID <- s
    
    rm(df_site, stats, stats_df)
    stats_wide
  }, error = function(e) {
    message(paste0("Error at site: ", s, " --> ", e$message))
  })
}

# Stop parallel cluster
stopCluster(cl)

# Write out combined statistics file
stats_vectors <- stats_vectors[!sapply(stats_vectors, is.null)]
df_final <- do.call(rbind, stats_vectors)
write.csv(df_final, outFile, row.names = FALSE)
