library(dplyr)
library(tidyr)
library(tibble)
library(ranger)
library(doParallel)
library(foreach)
library(caret)

inDir <- "G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/github"
outDir <- file.path(inDir, "models", "rf_hybrid")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

source(file.path(inDir, "scripts/fit_rfe_cluster_rf_v2.R"))

# -----------------------------
# User settings
# -----------------------------
num.threads <- 1
rfe_on <- 1
n_workers <- max(1, parallel::detectCores(logical = FALSE) - 1)

# -----------------------------
# Model file table
# -----------------------------
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
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# Shared input data
# -----------------------------
siteInfo <- read.csv(file.path(inDir, "data/biosoundscape_sites_daac_250507_dbscan_cluster_2000m.csv"))
siteInfo$LogRecNum <- log1p(siteInfo$RecordingN)
siteInfo <- siteInfo %>% select(SiteID, CLUSTER_ID, LogRecNum)

species_pc <- read.csv(file.path(inDir, "data/bioscape_point_count_richness_v20250602.csv"))
species_pc <- species_pc %>%
  rename(PCrichness = PC_Richness) %>%
  select(SiteID, PCrichness, Campaign)

species_ai <- read.csv(file.path(inDir, "data/wildmon_site-level_species_250907.csv"))
species_ai <- species_ai %>%
  rename(AIrichness = wm_richness_ge3, SiteID = siteid) %>%
  select(SiteID, AIrichness)

species_df <- left_join(species_ai, species_pc, by = "SiteID")
species_filtered <- left_join(species_df, siteInfo, by = "SiteID")

# -----------------------------
# Helpers
# -----------------------------
safe_load_object <- function(filepath, object_name) {
  e <- new.env(parent = emptyenv())
  load(filepath, envir = e)
  if (!exists(object_name, envir = e, inherits = FALSE)) {
    stop("Object '", object_name, "' not found in ", filepath)
  }
  get(object_name, envir = e, inherits = FALSE)
}

extract_final_predictors <- function(model_list, metric_name) {
  model_list[[metric_name]]$final_data %>%
    select(-response, -any_of("LogRecNum"))
}

build_combined_predictors <- function(desc_models, harm_models, metric_name, drop_logrecnum_again = FALSE) {
  out <- cbind(
    extract_final_predictors(desc_models, metric_name),
    extract_final_predictors(harm_models, metric_name)
  )
  out <- out[, !duplicated(colnames(out)), drop = FALSE]
  if (drop_logrecnum_again) {
    out <- out %>% select(-any_of("LogRecNum"))
  }
  out
}

get_predictor_file_base <- function(descriptive_stub) {
  sub("_descriptive_statistics$", "", basename(descriptive_stub))
}

get_output_file <- function(predictor_file, season) {
  file.path(outDir, paste0(predictor_file, "_hybrid_", season, ".RData"))
}

prepare_model_group_data <- function(i) {
  models <- model_files[i, ]
  predictor_file <- get_predictor_file_base(models$descriptive_models)
  
  dry_desc_file <- file.path(inDir, paste0("models/", models$descriptive_models, "_dry_season.RData"))
  dry_harm_file <- file.path(inDir, paste0("models/", models$harmonic_models, "_dry_season.RData"))
  wet_desc_file <- file.path(inDir, paste0("models/", models$descriptive_models, "_wet_season.RData"))
  wet_harm_file <- file.path(inDir, paste0("models/", models$harmonic_models, "_wet_season.RData"))
  drywet_desc_file <- file.path(inDir, paste0("models/", models$descriptive_models, "_dry-wet_season.RData"))
  drywet_harm_file <- file.path(inDir, paste0("models/", models$harmonic_models, "_dry-wet_season.RData"))
  
  dry_desc <- safe_load_object(dry_desc_file, "dry_models")
  dry_harm <- safe_load_object(dry_harm_file, "dry_models")
  wet_desc <- safe_load_object(wet_desc_file, "wet_models")
  wet_harm <- safe_load_object(wet_harm_file, "wet_models")
  drywet_desc <- safe_load_object(drywet_desc_file, "drywet_models")
  drywet_harm <- safe_load_object(drywet_harm_file, "drywet_models")
  
  df_dry_pc <- build_combined_predictors(dry_desc, dry_harm, "PCrichness", FALSE)
  df_wet_pc <- build_combined_predictors(wet_desc, wet_harm, "PCrichness", FALSE)
  df_drywet_pc <- build_combined_predictors(drywet_desc, drywet_harm, "PCrichness", FALSE)
  
  df_dry_ai <- build_combined_predictors(dry_desc, dry_harm, "AIrichness", TRUE)
  df_wet_ai <- build_combined_predictors(wet_desc, wet_harm, "AIrichness", TRUE)
  df_drywet_ai <- build_combined_predictors(drywet_desc, drywet_harm, "AIrichness", TRUE)
  
  dry_data_pc <- species_filtered %>%
    filter(Campaign == "Dry season", !is.na(PCrichness)) %>%
    select(CLUSTER_ID, LogRecNum, PCrichness) %>%
    cbind(df_dry_pc)
  
  wet_data_pc <- species_filtered %>%
    filter(Campaign == "Wet season", !is.na(PCrichness)) %>%
    select(CLUSTER_ID, LogRecNum, PCrichness) %>%
    cbind(df_wet_pc)
  
  drywet_data_pc <- species_filtered %>%
    filter(!is.na(PCrichness)) %>%
    select(CLUSTER_ID, LogRecNum, PCrichness) %>%
    cbind(df_drywet_pc)
  
  dry_data_ai <- species_filtered %>%
    filter(Campaign == "Dry season", !is.na(AIrichness)) %>%
    select(CLUSTER_ID, LogRecNum, AIrichness) %>%
    cbind(df_dry_ai)
  
  wet_data_ai <- species_filtered %>%
    filter(Campaign == "Wet season", !is.na(AIrichness)) %>%
    select(CLUSTER_ID, LogRecNum, AIrichness) %>%
    cbind(df_wet_ai)
  
  drywet_data_ai <- species_filtered %>%
    filter(!is.na(AIrichness)) %>%
    select(CLUSTER_ID, LogRecNum, AIrichness) %>%
    cbind(df_drywet_ai)
  
  list(
    predictor_file = predictor_file,
    dry_data_pc = dry_data_pc,
    wet_data_pc = wet_data_pc,
    drywet_data_pc = drywet_data_pc,
    dry_data_ai = dry_data_ai,
    wet_data_ai = wet_data_ai,
    drywet_data_ai = drywet_data_ai
  )
}

build_task_table <- function() {
  tasks <- vector("list", length = nrow(model_files) * 6)
  k <- 1
  
  for (i in seq_len(nrow(model_files))) {
    predictor_file <- get_predictor_file_base(model_files$descriptive_models[i])
    
    tasks[[k]] <- data.frame(model_index = i, predictor_file = predictor_file, season = "dry_season", metric = "PCrichness", stringsAsFactors = FALSE); k <- k + 1
    tasks[[k]] <- data.frame(model_index = i, predictor_file = predictor_file, season = "wet_season", metric = "PCrichness", stringsAsFactors = FALSE); k <- k + 1
    tasks[[k]] <- data.frame(model_index = i, predictor_file = predictor_file, season = "dry-wet_season", metric = "PCrichness", stringsAsFactors = FALSE); k <- k + 1
    tasks[[k]] <- data.frame(model_index = i, predictor_file = predictor_file, season = "dry_season", metric = "AIrichness", stringsAsFactors = FALSE); k <- k + 1
    tasks[[k]] <- data.frame(model_index = i, predictor_file = predictor_file, season = "wet_season", metric = "AIrichness", stringsAsFactors = FALSE); k <- k + 1
    tasks[[k]] <- data.frame(model_index = i, predictor_file = predictor_file, season = "dry-wet_season", metric = "AIrichness", stringsAsFactors = FALSE); k <- k + 1
  }
  
  bind_rows(tasks)
}

fit_task <- function(task_row) {
  i <- task_row$model_index
  predictor_file <- task_row$predictor_file
  season <- task_row$season
  metric <- task_row$metric
  
  message("[START] ", predictor_file, " | ", season, " | ", metric)
  
  dat_list <- prepare_model_group_data(i)
  
  dat <- switch(
    paste(season, metric, sep = "|"),
    "dry_season|PCrichness" = dat_list$dry_data_pc,
    "wet_season|PCrichness" = dat_list$wet_data_pc,
    "dry-wet_season|PCrichness" = dat_list$drywet_data_pc,
    "dry_season|AIrichness" = dat_list$dry_data_ai,
    "wet_season|AIrichness" = dat_list$wet_data_ai,
    "dry-wet_season|AIrichness" = dat_list$drywet_data_ai,
    stop("Unknown season/metric combination")
  )
  
  predictor_vars <- dat %>%
    select(-all_of(metric), -CLUSTER_ID, -LogRecNum) %>%
    colnames()
  
  result <- fit_rfe_cluster_rf(
    cor_data = dat,
    predictor_vars = predictor_vars,
    metric = metric,
    num.threads = num.threads,
    rfe_on = rfe_on,
    importance = "permutation"
  )
  
  rm(dat_list, dat)
  gc()
  
  message("[DONE] ", predictor_file, " | ", season, " | ", metric)
  
  list(
    model_index = i,
    predictor_file = predictor_file,
    season = season,
    metric = metric,
    result = result,
    status = "completed"
  )
}

save_group_results <- function(group_results, predictor_file) {
  dry_models <- list()
  wet_models <- list()
  drywet_models <- list()
  
  for (x in group_results) {
    if (x$season == "dry_season") {
      dry_models[[x$metric]] <- x$result
    } else if (x$season == "wet_season") {
      wet_models[[x$metric]] <- x$result
    } else if (x$season == "dry-wet_season") {
      drywet_models[[x$metric]] <- x$result
    }
  }
  
  dry_file <- get_output_file(predictor_file, "dry_season")
  wet_file <- get_output_file(predictor_file, "wet_season")
  drywet_file <- get_output_file(predictor_file, "dry-wet_season")
  
  if (length(dry_models) > 0) save(dry_models, file = dry_file)
  if (length(wet_models) > 0) save(wet_models, file = wet_file)
  if (length(drywet_models) > 0) save(drywet_models, file = drywet_file)
}

# -----------------------------
# Build task table
# -----------------------------
tasks <- build_task_table()

# Skip tasks where final season file already exists and likely contains both metrics
tasks <- tasks %>%
  rowwise() %>%
  mutate(
    outfile = get_output_file(predictor_file, season),
    skip_existing = file.exists(outfile)
  ) %>%
  ungroup()

message("Total tasks: ", nrow(tasks))
message("Tasks skipped due to existing output file: ", sum(tasks$skip_existing))
message("Tasks to run: ", sum(!tasks$skip_existing))

tasks_to_run <- tasks %>% filter(!skip_existing)

# -----------------------------
# Parallel execution
# -----------------------------
cl <- makeCluster(n_workers)
registerDoParallel(cl)

message("Running ", nrow(tasks_to_run), " model fits with ", n_workers, " workers")
message("Each task uses num.threads = ", num.threads)

results_list <- foreach(
  j = seq_len(nrow(tasks_to_run)),
  .packages = c("dplyr", "tidyr", "tibble", "ranger", "caret"),
  .export = c(
    "tasks_to_run", "model_files", "species_filtered", "inDir", "outDir",
    "num.threads", "rfe_on",
    "safe_load_object", "extract_final_predictors", "build_combined_predictors",
    "get_predictor_file_base", "get_output_file", "prepare_model_group_data",
    "build_task_table", "fit_task", "save_group_results", "fit_rfe_cluster_rf"
  ),
  .errorhandling = "pass",
  .options.snow = list(preschedule = FALSE)
) %dopar% {
  task_row <- tasks_to_run[j, ]
  tryCatch(
    fit_task(task_row),
    error = function(e) {
      list(
        model_index = task_row$model_index,
        predictor_file = task_row$predictor_file,
        season = task_row$season,
        metric = task_row$metric,
        result = NULL,
        status = paste("error:", conditionMessage(e))
      )
    }
  )
}

stopCluster(cl)
registerDoSEQ()

# -----------------------------
# Save completed results grouped back into season files
# -----------------------------
successful_results <- Filter(function(x) is.list(x) && identical(x$status, "completed"), results_list)
failed_results <- Filter(function(x) is.list(x) && !identical(x$status, "completed"), results_list)

if (length(successful_results) > 0) {
  split_keys <- vapply(successful_results, function(x) x$predictor_file, character(1))
  grouped_results <- split(successful_results, split_keys)
  
  for (nm in names(grouped_results)) {
    save_group_results(grouped_results[[nm]], nm)
  }
}

# -----------------------------
# Write run log
# -----------------------------
result_log <- bind_rows(lapply(results_list, function(x) {
  data.frame(
    model_index = x$model_index,
    predictor_file = x$predictor_file,
    season = x$season,
    metric = x$metric,
    status = x$status,
    stringsAsFactors = FALSE
  )
}))

if (nrow(tasks %>% filter(skip_existing)) > 0) {
  skipped_log <- tasks %>%
    filter(skip_existing) %>%
    transmute(
      model_index,
      predictor_file,
      season,
      metric,
      status = "skipped_existing"
    )
  result_log <- bind_rows(result_log, skipped_log)
}

result_log <- result_log %>%
  arrange(model_index, predictor_file, season, metric)

print(result_log)
write.csv(result_log, file.path(outDir, "rf_hybrid_parallel_fit_task_status.csv"), row.names = FALSE)