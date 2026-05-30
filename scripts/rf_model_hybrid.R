library(dplyr)
library(tidyr)
library(tibble)
library(ranger)
library(doParallel)
library(foreach)
library(caret)

inDir <- "G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/github"

source(file.path(inDir,"scripts/fit_rfe_cluster_rf.R"))

# Number of threads for parallel processing
num.threads = 4

# Perform RFE? (1=Yes)
rfe_on <- 1

model_files <- data.frame(
  descriptive_models = c(
    
    "rf_birdnet/birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_descriptive_statistics",
    "rf_birdnet/birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_descriptive_statistics",
    "rf_birdnet/birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_descriptive_statistics",
    "rf_birdnet/birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_descriptive_statistics",

    "rf_vggish/vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_descriptive_statistics",
    "rf_vggish/vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_descriptive_statistics",
    "rf_vggish/vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_descriptive_statistics",
    "rf_vggish/vggish_embeddings_minutes_nonzero_mean_dawn-dusk_descriptive_statistics",

    "rf_acoustic_indices/acoustic_indices_minutes_dawn-dusk_AGI90_descriptive_statistics",
    "rf_acoustic_indices/acoustic_indices_minutes_dawn-dusk_AGI75_descriptive_statistics",
    "rf_acoustic_indices/acoustic_indices_minutes_dawn-dusk_AGI50_descriptive_statistics",
    "rf_acoustic_indices/acoustic_indices_minutes_dawn-dusk_descriptive_statistics",

    "rf_aves/aves_embeddings_minutes_nonzero_mean_dawn-dusk_descriptive_statistics",
    "rf_aves/aves_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_descriptive_statistics",
    "rf_aves/aves_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_descriptive_statistics",
    "rf_aves/aves_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_descriptive_statistics"
  ),
  harmonic_models = c(

    "rf_birdnet/birdnet_embeddings_minutes_nonzero_mean_AGI90_2phase-harmonic_statistics",
    "rf_birdnet/birdnet_embeddings_minutes_nonzero_mean_AGI75_2phase-harmonic_statistics",
    "rf_birdnet/birdnet_embeddings_minutes_nonzero_mean_AGI50_2phase-harmonic_statistics",
    "rf_birdnet/birdnet_embeddings_minutes_nonzero_mean_2phase-harmonic_statistics",

    "rf_vggish/vggish_embeddings_minutes_nonzero_mean_AGI90_2phase-harmonic_statistics",
    "rf_vggish/vggish_embeddings_minutes_nonzero_mean_AGI75_2phase-harmonic_statistics",
    "rf_vggish/vggish_embeddings_minutes_nonzero_mean_AGI50_2phase-harmonic_statistics",
    "rf_vggish/vggish_embeddings_minutes_nonzero_mean_2phase-harmonic_statistics",

    "rf_acoustic_indices/acoustic_indices_minutes_AGI90_2phase-harmonic_statistics",
    "rf_acoustic_indices/acoustic_indices_minutes_AGI75_2phase-harmonic_statistics",
    "rf_acoustic_indices/acoustic_indices_minutes_AGI50_2phase-harmonic_statistics",
    "rf_acoustic_indices/acoustic_indices_minutes_2phase-harmonic_statistics",

    "rf_aves/aves_embeddings_minutes_nonzero_mean_2phase-harmonic_statistics",
    "rf_aves/aves_embeddings_minutes_nonzero_mean_AGI50_2phase-harmonic_statistics",
    "rf_aves/aves_embeddings_minutes_nonzero_mean_AGI75_2phase-harmonic_statistics",
    "rf_aves/aves_embeddings_minutes_nonzero_mean_AGI90_2phase-harmonic_statistics"

  )
)

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

species_filtered <- left_join(species_df, siteInfo, by = "SiteID")

# Parallel setup
cores <- parallel::detectCores()
cl <- makeCluster(cores - 2)
registerDoParallel(cl)
message("Parallel RFE using ", cores - 2, " workers")

# Loop through predictor files
for (i in 1:dim(model_files)[1]){
  
  models <- model_files[i,]
  
  load(file.path(inDir, paste0("models/",models$descriptive_models, "_dry_season.RData")))
  df_dry_pc_desc <- dry_models$PCrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  df_dry_ai_desc <- dry_models$AIrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  load(file.path(inDir, paste0("models/",models$harmonic_models, "_dry_season.RData")))
  df_dry_pc_harm <- dry_models$PCrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  df_dry_ai_harm <- dry_models$AIrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  
  load(file.path(inDir, paste0("models/",models$descriptive_models, "_wet_season.RData")))
  df_wet_pc_desc <- wet_models$PCrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  df_wet_ai_desc <- wet_models$AIrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  load(file.path(inDir, paste0("models/",models$harmonic_models, "_wet_season.RData")))
  df_wet_pc_harm <- wet_models$PCrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  df_wet_ai_harm <- wet_models$AIrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  
  load(file.path(inDir, paste0("models/",models$descriptive_models, "_dry-wet_season.RData")))
  df_drywet_pc_desc <- drywet_models$PCrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  df_drywet_ai_desc <- drywet_models$AIrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  load(file.path(inDir, paste0("models/",models$harmonic_models, "_dry-wet_season.RData")))
  df_drywet_pc_harm <- drywet_models$PCrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  df_drywet_ai_harm <- drywet_models$AIrichness$final_data %>% select(-response,-any_of("LogRecNum"))
  

  df_dry_pc <- cbind(df_dry_pc_desc,df_dry_pc_harm)
  remove(df_dry_pc_desc,df_dry_pc_harm)
  df_wet_pc <- cbind(df_wet_pc_desc,df_wet_pc_harm)
  remove(df_wet_pc_desc,df_wet_pc_harm)
  df_drywet_pc <- cbind(df_drywet_pc_desc,df_drywet_pc_harm)
  remove(df_drywet_pc_desc,df_drywet_pc_harm)
  
  df_dry_ai <- cbind(df_dry_ai_desc,df_dry_ai_harm) %>% select(-any_of("LogRecNum"))
  remove(df_dry_ai_desc,df_dry_ai_harm)
  df_wet_ai <- cbind(df_wet_ai_desc,df_wet_ai_harm) %>% select(-any_of("LogRecNum"))
  remove(df_wet_ai_desc,df_wet_ai_harm)
  df_drywet_ai <- cbind(df_drywet_ai_desc,df_drywet_ai_harm) %>% select(-any_of("LogRecNum"))
  remove(df_drywet_ai_desc,df_drywet_ai_harm)
  
  # Prepare data
  dry_data_pc <- species_filtered %>%
    filter(Campaign == "Dry season" & !is.na(PCrichness)) %>%
    select(CLUSTER_ID, LogRecNum, PCrichness) %>% cbind(df_dry_pc)
  wet_data_pc <- species_filtered %>%
    filter(Campaign == "Wet season" & !is.na(PCrichness)) %>%
    select(CLUSTER_ID, LogRecNum, PCrichness) %>% cbind(df_wet_pc)
  drywet_data_pc <- species_filtered %>%
    filter(!is.na(PCrichness)) %>%
    select(CLUSTER_ID, LogRecNum, PCrichness) %>% cbind(df_drywet_pc)
   
  dry_data_ai <- species_filtered %>%
    filter(Campaign == "Dry season" & !is.na(AIrichness)) %>%
    select(CLUSTER_ID, LogRecNum, AIrichness) %>% cbind(df_dry_ai)
  wet_data_ai <- species_filtered %>%
    filter(Campaign == "Wet season" & !is.na(AIrichness)) %>%
    select(CLUSTER_ID, LogRecNum, AIrichness) %>% cbind(df_wet_ai)
  drywet_data_ai <- species_filtered %>%
    filter(!is.na(AIrichness)) %>%
    select(CLUSTER_ID, LogRecNum, AIrichness) %>% cbind(df_drywet_ai)

  # Run models
  dry_models <- list()
  wet_models <- list()
  drywet_models <- list()
  
  # RF models for PC richness data
  metric <- "PCrichness"
  
  cat("\nRunning dry season model for:", metric, "\n")
  predictor_vars <- dry_data_pc %>% select(-PCrichness,-CLUSTER_ID,-LogRecNum) %>% colnames()
  result <- fit_rfe_cluster_rf(dry_data_pc, predictor_vars, metric, num.threads = num.threads, 
                               rfe_on = rfe_on, importance = "permutation")
  dry_models[[metric]] <- result
  
  cat("\nRunning wet season model for:", metric, "\n")
  predictor_vars <- wet_data_pc %>% select(-PCrichness,-CLUSTER_ID,-LogRecNum) %>% colnames()
  result <- fit_rfe_cluster_rf(wet_data_pc, predictor_vars, metric, num.threads = num.threads, 
                               rfe_on = rfe_on, importance = "permutation")
  wet_models[[metric]] <- result
  
  cat("\nRunning dry-wet season model for:", metric, "\n")
  predictor_vars <- drywet_data_pc %>% select(-PCrichness,-CLUSTER_ID,-LogRecNum) %>% colnames()
  result <- fit_rfe_cluster_rf(drywet_data_pc, predictor_vars, metric, num.threads = num.threads, 
                               rfe_on = rfe_on, importance = "permutation")
  drywet_models[[metric]] <- result
  
  # RF models for AI richness data
  metric <- "AIrichness"

  cat("\nRunning dry season model for:", metric, "\n")
  predictor_vars <- dry_data_ai %>% select(-AIrichness,-CLUSTER_ID,-LogRecNum) %>% colnames()
  result <- fit_rfe_cluster_rf(dry_data_ai, predictor_vars, metric, num.threads = num.threads, 
                               rfe_on = rfe_on, importance = "permutation")
  dry_models[[metric]] <- result
  
  cat("\nRunning wet season model for:", metric, "\n")
  predictor_vars <- wet_data_ai %>% select(-AIrichness,-CLUSTER_ID,-LogRecNum) %>% colnames()
  result <- fit_rfe_cluster_rf(wet_data_ai, predictor_vars, metric, num.threads = num.threads, 
                               rfe_on = rfe_on, importance = "permutation")
  wet_models[[metric]] <- result
  
  cat("\nRunning dry-wet season model for:", metric, "\n")
  predictor_vars <- drywet_data_ai %>% select(-AIrichness,-CLUSTER_ID,-LogRecNum) %>% colnames()
  result <- fit_rfe_cluster_rf(drywet_data_ai, predictor_vars, metric, num.threads = num.threads, 
                               rfe_on = rfe_on, importance = "permutation")
  drywet_models[[metric]] <- result
  
  # Save results
  predictor_file <- sub("_descriptive_statistics$", "",basename(models$descriptive_models)) 
  outfile <- paste0(predictor_file, "_hybrid_dry_season.RData")
  save(dry_models, file = file.path(inDir,"models/rf_hybrid",outfile))
  
  outfile <- paste0(predictor_file, "_hybrid_wet_season.RData")
  save(wet_models, file = file.path(inDir,"models/rf_hybrid",outfile))
  
  outfile <- paste0(predictor_file, "_hybrid_dry-wet_season.RData")
  save(drywet_models, file = file.path(inDir,"models/rf_hybrid",outfile))
  
}

# Cleanup parallel
stopCluster(cl)
registerDoSEQ()