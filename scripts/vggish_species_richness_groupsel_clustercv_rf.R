# Purpose: Random forest modeling of species richness from VGGish embedding statistics
# Matthew Clark, March 27, 2026

library(dplyr)
library(tidyr)
library(tibble)
library(ranger)
library(doParallel)
library(foreach)
library(caret)

inDir <- "G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/github"

source(file.path(inDir,"scripts/fit_grouped_rfe_cluster_rf.R"))


predictor_files <- c(
  
  # vggish embedding descriptive statistics
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_descriptive_statistics.csv",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_descriptive_statistics.csv",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_descriptive_statistics.csv",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_descriptive_statistics.csv",
  
  # vggish embedding harmonic statistics
  "vggish_embeddings_minutes_nonzero_mean_AGI90_2phase-harmonic_statistics.csv",
  "vggish_embeddings_minutes_nonzero_mean_AGI75_2phase-harmonic_statistics.csv",
  "vggish_embeddings_minutes_nonzero_mean_AGI50_2phase-harmonic_statistics.csv",
  "vggish_embeddings_minutes_nonzero_mean_2phase-harmonic_statistics.csv"
  
)

# Top number of variables to select per group prior to RFE
top_n_per_group <- 50

# Number of threads for parallel processing
num.threads = 2

# Perform RFE? (1=Yes)
rfe_on <- 1

# Load cluster information
siteInfo <- read.csv(file.path(inDir,"data/biosoundscape_sites_daac_250507_dbscan_cluster_2000m.csv"))
siteInfo$LogRecNum <- log1p(siteInfo$RecordingN)
siteInfo <- siteInfo %>% select(SiteID,CLUSTER_ID,LogRecNum)

# Load richness data
species_pc <- read.csv(file.path(inDir,"data/bioscape_point_count_richness_v20250602.csv"))
species_pc <- species_pc %>% rename(PCrichness = PC_Richness) %>%
  select(SiteID,PCrichness, Campaign)

species_ai <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) 
species_ai <- species_ai %>% rename(AIrichness = wm_richness_ge3, SiteID = siteid) %>%
  select(SiteID,AIrichness)

species_df <- left_join(species_ai,species_pc,by="SiteID")

diversity_vars <- c("PCrichness", "AIrichness")

# Parallel setup
cores <- parallel::detectCores()
cl <- makeCluster(cores - 2)
registerDoParallel(cl)
message("Parallel RFE using ", cores - 2, " workers")

# Loop through predictor files
for (predictor_file in predictor_files){
  
  df <- read.csv(file.path(inDir,"data",predictor_file))
  
  species_filtered <- left_join(species_df, df, by = "SiteID")
  species_filtered <- left_join(species_filtered, siteInfo, by = "SiteID")
  
  # Prepare data
  dry_data <- species_filtered %>%
    filter(Campaign == "Dry season") %>%
    select(PCrichness, AIrichness, CLUSTER_ID, LogRecNum, matches("^fea"))
  wet_data <- species_filtered %>%
    filter(Campaign == "Wet season") %>%
    select(PCrichness, AIrichness, CLUSTER_ID, LogRecNum, matches("^fea"))
  drywet_data <- species_filtered %>%
    select(PCrichness, AIrichness, CLUSTER_ID, LogRecNum, matches("^fea"))
  
  predictor_vars <- dry_data %>% select(matches("fea")) %>% colnames()
  
  # Run models
  dry_models <- list()
  wet_models <- list()
  drywet_models <- list()
  
  for (metric in diversity_vars) {
    cat("\nRunning dry season model for:", metric, "\n")
    result <- fit_grouped_rfe_cluster_rf(dry_data, predictor_vars, metric, num.threads = num.threads, 
                                         top_n_per_group = top_n_per_group, rfe_on = rfe_on, importance = "permutation")
    dry_models[[metric]] <- result
    
    cat("\nRunning wet season model for:", metric, "\n")
    result <- fit_grouped_rfe_cluster_rf(wet_data, predictor_vars, metric, num.threads = num.threads, 
                                         top_n_per_group = top_n_per_group, rfe_on = rfe_on, importance = "permutation")
    wet_models[[metric]] <- result
    
    cat("\nRunning dry-wet season model for:", metric, "\n")
    result <- fit_grouped_rfe_cluster_rf(drywet_data, predictor_vars, metric, num.threads = num.threads, 
                                         top_n_per_group = top_n_per_group, rfe_on = rfe_on, importance = "permutation")
    drywet_models[[metric]] <- result
  }
  
  # Save results
  outfile <- sub("\\.csv$", "_dry_season.RData", predictor_file)
  save(dry_models, file = file.path(inDir,"models/rf_vggish",outfile))
  
  outfile <- sub("\\.csv$", "_wet_season.RData", predictor_file)
  save(wet_models, file = file.path(inDir,"models/rf_vggish",outfile))
  
  outfile <- sub("\\.csv$", "_dry-wet_season.RData", predictor_file)
  save(drywet_models, file = file.path(inDir,"models/rf_vggish",outfile))
  
}

# Cleanup parallel
stopCluster(cl)
registerDoSEQ()