# Purpose: Calculates harmonic statistics on acoustic indices
# Matthew Clark, March 27, 2026

# ---- Harmonic Feature Summary  ----
# For each acoustic feature, we fit a two-term harmonic model:
#   y(t) = a0 + a1*cos(2 pi t/24) + b1*sin(2 pi t/24) + a2*cos(4 pi t/24) + b2*sin(4 pi t/24)
#
# From this model, we extract the following metrics per feature:
#
# 1. a0  - Baseline activity level (mean across the 24-hour cycle)
#    . Ecological interpretation: Higher a0 may indicate consistently high acoustic presence (e.g., overall vocal activity).
#
# 2. amp1 = sqrt(a1^2 + b1^2) - Amplitude of the primary (24h) cycle
#    . Reflects the strength of daily (diel) rhythmicity.
#    . Higher amp1 may indicate stronger diurnal structure (e.g., consistent dawn chorus).
#
# 3. amp2 = sqrt(a2^2 + b2^2) - Amplitude of the secondary (12h) cycle
#    . Captures bimodal patterns (e.g., dawn and dusk).
#    . Useful for identifying species with multiple activity peaks.
#
# 4. amp_ratio = amp2 / amp1
#    . Relative importance of the secondary cycle.
#    . High values may indicate crepuscular (twilight-active) or complex temporal structure.
#
# 5. sin_phase1 = sin(phase1), cos_phase1 = cos(phase1)
#    . phase1 = atan2(-b1, a1) - timing of the primary (24h) peak, in radians
#    . Converted to sine/cosine to handle circularity for machine learning models.
#    . Indicates whether peak activity is centered in the morning, afternoon, evening, etc.
#
# 6. sin_phase2 = sin(phase2), cos_phase2 = cos(phase2)
#    . phase2 = atan2(-b2, a2) - timing of the secondary (12h) peak, in radians
#    . Also transformed to sine/cosine for compatibility with ML models.
#    . Useful for detecting timing of additional peaks (e.g., dusk activity following a dawn peak).
#
# These features serve as interpretable, compact summaries of daily activity structure
# that can be used in Random Forest models to predict ecological response variables like species richness.

library(dplyr)
library(doParallel)
library(foreach)
library(tidyr)

# Input directory
inDir <- "G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/github"

# Output file
outFile <- file.path(inDir, "data/acoustic_indices_minutes_AGI90_2phase-harmonic_statistics.csv")

# Load acoustic indices by campaign
df1 <- read.csv(file.path(inDir, "data/bioscape_acoustic_indices_campaign1_231022.csv"))
df2 <- read.csv(file.path(inDir, "data/bioscape_acoustic_indices_campaign2_240203.csv"))
df <- rbind(df1,df2) %>% select(-c(path,YYYY,MM,DD,hh,mm,DDhh))
remove(df1,df2)

# Filter out minutes with high level of geophony, interference or anthropophony?
agi_filter <- 1 # 1 = Yes, 0 = No
agi_thresh <- 90 # minutes with over this amount of AGI combined will be filtered
agbi_minutes <- file.path(inDir, "data/wildmon_minute-level_species_taxon_abgi.csv") # file with minutes summary

# Two-term harmonic fit function
fit_harmonic_two_term <- function(t, y) {
  t_rad <- 2 * pi * t / 24
  model <- tryCatch({
    lm(y ~ cos(t_rad) + sin(t_rad) + cos(2 * t_rad) + sin(2 * t_rad))
  }, error = function(e) return(rep(NA, 10)))  # updated length to 10
  
  if (is.null(model)) return(rep(NA, 10))
  
  coefs <- coef(model)
  a0 <- coefs[1]
  a1 <- coefs[2]
  b1 <- coefs[3]
  a2 <- coefs[4]
  b2 <- coefs[5]
  
  amp1 <- sqrt(a1^2 + b1^2)
  amp2 <- sqrt(a2^2 + b2^2)
  amp_ratio <- ifelse(amp1 > 0, amp2 / amp1, NA)
  
  # Phase 1
  phase1_rad <- atan2(-b1, a1)
  phase1_hour <- (phase1_rad * 24 / (2 * pi)) %% 24
  sin_phase1 <- sin(2 * pi * phase1_hour / 24)
  cos_phase1 <- cos(2 * pi * phase1_hour / 24)
  
  # Phase 2
  phase2_rad <- atan2(-b2, a2)
  phase2_hour <- (phase2_rad * 24 / (2 * pi)) %% 24
  sin_phase2 <- sin(2 * pi * phase2_hour / 24)
  cos_phase2 <- cos(2 * pi * phase2_hour / 24)
  
  return(c(a0 = a0, amp1 = amp1, amp2 = amp2, amp_ratio = amp_ratio,
           sin_phase1 = sin_phase1, cos_phase1 = cos_phase1,
           sin_phase2 = sin_phase2, cos_phase2 = cos_phase2))
}

# Extract site ID
extract_site <- function(filename) substr(basename(filename), 1, 15)

# Extract time-of-day as fractional hour
extract_time <- function(t) {
  hours <- t %/% 100
  minutes <- t %% 100
  hours + minutes / 60
}

# Get files list
files <- df$file
sites <- unique(sapply(files, extract_site))

# List of acoustic indices to summarize
cn <- c("ACI", "ADI", "AEI", "NDSI", "NDSI_A", "NDSI_B",
        "BI", "H", "Ht", "Hs", "M", "R", "sfm", "rugo", "zcr_mean")

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

stats_vectors <- foreach(s = sites, .packages = c("dplyr", "tidyr")) %dopar% {

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
    
    # Combine into matrix
    X <- df_site[,cn]
    if (nrow(X) < 10) return(NULL)
    times <- extract_time(df_site$hhmm)
    
    # Fit harmonic model to each feature
    harmonic_results <- apply(X, 2, function(y) fit_harmonic_two_term(times, y))
    harmonic_df <- as.data.frame(t(harmonic_results))
    harmonic_df$feature <- cn
    rownames(harmonic_df) <- NULL
    harmonic_df$SiteID <- s
    
    colnames(harmonic_df) <- c("a0", "amp1", "amp2", "amp_ratio", "sin_phase1", "cos_phase1", "sin_phase2", "cos_phase2", "feature","SiteID")
    
    harmonic_wide <- harmonic_df %>%
      pivot_longer(cols = -c(SiteID, feature), names_to = "stat", values_to = "value") %>%
      mutate(feature_stat = paste0(feature, "_", stat)) %>%
      select(SiteID, feature_stat, value) %>%
      pivot_wider(names_from = feature_stat, values_from = value)
    
    rm(df_site, stats, stats_df)
    harmonic_wide
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
