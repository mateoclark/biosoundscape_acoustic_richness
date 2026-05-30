# Purpose: Summarize statistics from random forest modeling of species richness 
# Matthew Clark, April 25, 2026

# Load libraries
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(tibble)

baseDir <- "G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/github"

# Set up parallel backend
future::plan(multisession, workers = parallel::detectCores() - 2)

# ---- Summary helper function ----
summarize_model <- function(df) {
  summarise(df,
            RMSE_mean = mean(RMSE),
            RMSE_se   = sd(RMSE) / sqrt(n()),
            Rsq_mean  = mean(Rsquared),
            Rsq_se    = sd(Rsquared) / sqrt(n())
  )
}

# ---- One-file processor ----
process_model_file <- function(m, inDir, richness_types = c("AIrichness", "PCrichness")) {
  load(file.path(inDir, paste0(m, "_dry_season.RData")))
  load(file.path(inDir, paste0(m, "_wet_season.RData")))
  load(file.path(inDir, paste0(m, "_dry-wet_season.RData")))
  
  model_data <- setNames(
    lapply(richness_types, function(rt) {
      list(
        dry    = dry_models[[rt]]$fold_results,
        wet    = wet_models[[rt]]$fold_results,
        drywet = drywet_models[[rt]]$fold_results
      )
    }),
    richness_types
  )
  
  summary_stats <- imap_dfr(
    model_data,
    ~ map_dfr(.x, summarize_model, .id = "season") %>%
      mutate(richness_type = .y),
    .id = NULL
  ) %>%
    relocate(richness_type, season) %>%
    mutate(
      file_name = m,
      predictors = case_when(
        grepl("2phase-harmonic", m) ~ "2phase-harmonic",
        grepl("descriptive", m)     ~ "descriptive",
        grepl("hybrid", m)          ~ "hybrid",
        TRUE                        ~ NA_character_
      ),
      periods = case_when(
        grepl("dawn-dusk", m) ~ "Dawn & Dusk",
        TRUE                  ~ "None"
      ),
      AGIfilter = case_when(
        grepl("AGI90", m) ~ "AGI90",
        grepl("AGI75", m) ~ "AGI75",
        grepl("AGI50", m) ~ "AGI50",
        TRUE              ~ "None"
      )
    )
  
  fold_df_all <- imap_dfr(
    model_data,
    ~ bind_rows(
      mutate(.x$dry, season = "dry"),
      mutate(.x$wet, season = "wet"),
      mutate(.x$drywet, season = "drywet")
    ) %>%
      mutate(
        richness_type = .y,
        file_name = m,
        predictors = case_when(
          grepl("2phase-harmonic", m) ~ "2phase-harmonic",
          grepl("descriptive", m)     ~ "descriptive",
          grepl("hybrid", m)          ~ "hybrid",
          TRUE                        ~ NA_character_
        ),
        periods = case_when(
          grepl("dawn-dusk", m) ~ "Dawn & Dusk",
          TRUE                  ~ "None"
        ),
        AGIfilter = case_when(
          grepl("AGI90", m) ~ "AGI90",
          grepl("AGI75", m) ~ "AGI75",
          grepl("AGI50", m) ~ "AGI50",
          TRUE              ~ "None"
        )
      ) %>%
      relocate(richness_type, season, file_name, predictors, periods, AGIfilter, Resample),
    .id = NULL
  )
  
  list(summary = summary_stats, folds = fold_df_all)
}

# ---- Main function ----
getStats <- function(model_files, inDir, richness_types = c("AIrichness", "PCrichness")) {
  result_list <- future_map(
    model_files,
    ~ process_model_file(.x, inDir,richness_types),
    .options = furrr_options(seed = TRUE)
  )
  
  summary_df <- bind_rows(map(result_list, "summary"))
  fold_df    <- bind_rows(map(result_list, "folds"))
  
  comparisons <- c("AGI50", "AGI75", "AGI90")
  
  stats_table <- fold_df %>%
    filter((predictors == "descriptive" & periods == "Dawn & Dusk") |
             (predictors %in% c("2phase-harmonic","hybrid"))) %>%
    group_by(richness_type, season, predictors) %>%
    group_split() %>%
    map_dfr(function(df_group) {
      map_dfr(comparisons, function(comp) {
        if (all(c("None", comp) %in% unique(df_group$AGIfilter))) {
          df_comp <- df_group %>% filter(AGIfilter %in% c("None", comp))
          
          df_rs <- df_comp %>%
            group_by(AGIfilter, Resample) %>%
            summarise(Rsquared = mean(Rsquared), RMSE = mean(RMSE), .groups = "drop")
          
          rsq_wide <- df_rs %>%
            select(AGIfilter, Resample, Rsquared) %>%
            pivot_wider(names_from = AGIfilter, values_from = Rsquared) %>%
            drop_na()
          
          rmse_wide <- df_rs %>%
            select(AGIfilter, Resample, RMSE) %>%
            pivot_wider(names_from = AGIfilter, values_from = RMSE) %>%
            drop_na()
          
          tibble(
            richness_type = unique(df_group$richness_type),
            season = unique(df_group$season),
            predictors = unique(df_group$predictors),
            comparison = comp,
            Rsq_p = if (nrow(rsq_wide) > 1) wilcox.test(rsq_wide[[comp]], rsq_wide$None, paired = TRUE, exact = FALSE)$p.value else NA_real_,
            RMSE_p = if (nrow(rmse_wide) > 1) wilcox.test(rmse_wide[[comp]], rmse_wide$None, paired = TRUE, exact = FALSE)$p.value else NA_real_,
            Rsq_diff = if ("None" %in% names(rsq_wide) && comp %in% names(rsq_wide)) mean(rsq_wide[[comp]] - rsq_wide$None, na.rm = TRUE) else NA_real_,
            RMSE_diff = if ("None" %in% names(rmse_wide) && comp %in% names(rmse_wide)) mean(rmse_wide[[comp]] - rmse_wide$None, na.rm = TRUE) else NA_real_,
            n_folds = nrow(rsq_wide)
          )
  
        } else {
          tibble(
            richness_type = unique(df_group$richness_type),
            season = unique(df_group$season),
            predictors = unique(df_group$predictors),
            comparison = comp,
            Rsq_p = NA_real_,
            RMSE_p = NA_real_,
            Rsq_diff = NA_real_,
            RMSE_diff = NA_real_,
            n_folds = 0
          )
        }
      })
    })
  
  list(summary = summary_df, folds = fold_df, wilcox = stats_table)
}

##################################
# BirdNET statistics
##################################
inDir <- file.path(baseDir,"models/rf_birdnet")

model_files <- c(
  "birdnet_embeddings_minutes_nonzero_mean_AGI90_2phase-harmonic_statistics",
  "birdnet_embeddings_minutes_nonzero_mean_AGI75_2phase-harmonic_statistics",
  "birdnet_embeddings_minutes_nonzero_mean_AGI50_2phase-harmonic_statistics",
  "birdnet_embeddings_minutes_nonzero_mean_2phase-harmonic_statistics",
  "birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_descriptive_statistics",
  "birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_descriptive_statistics",
  "birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_descriptive_statistics",
  "birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_descriptive_statistics",
  "birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_hybrid",
  "birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_hybrid",
  "birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_hybrid",
  "birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_hybrid"
)

# Execute and store result
birdnet_results <- getStats(model_files, inDir)
birdnet_results$folds$features <- "BirdNET embeddings"
birdnet_stats <- birdnet_results$summary
birdnet_stats$features <- "BirdNET embeddings"

# Access Wilcoxon tests
birdnet_wilcox <- birdnet_results$wilcox
birdnet_wilcox$features <- "BirdNET embeddings"

##################################
# VGGish statistics
##################################
inDir <- file.path(baseDir,"models/rf_vggish")

model_files <- c(
  "vggish_embeddings_minutes_nonzero_mean_AGI90_2phase-harmonic_statistics",
  "vggish_embeddings_minutes_nonzero_mean_AGI75_2phase-harmonic_statistics",
  "vggish_embeddings_minutes_nonzero_mean_AGI50_2phase-harmonic_statistics",
  "vggish_embeddings_minutes_nonzero_mean_2phase-harmonic_statistics",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_descriptive_statistics",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_descriptive_statistics",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_descriptive_statistics",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_descriptive_statistics",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_hybrid",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_hybrid",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_hybrid",
  "vggish_embeddings_minutes_nonzero_mean_dawn-dusk_hybrid"
)

# Execute and store result
vggish_results <- getStats(model_files, inDir)
vggish_results$folds$features <- "VGGish embeddings"
vggish_stats <- vggish_results$summary
vggish_stats$features <- "VGGish embeddings"

# Access Wilcoxon tests
vggish_wilcox <- vggish_results$wilcox
vggish_wilcox$features <- "VGGish embeddings"

##################################
# AVES statistics
##################################
inDir <- file.path(baseDir,"models/rf_aves")

model_files <- c(
  "aves_embeddings_minutes_nonzero_mean_AGI90_2phase-harmonic_statistics",
  "aves_embeddings_minutes_nonzero_mean_AGI75_2phase-harmonic_statistics",
  "aves_embeddings_minutes_nonzero_mean_AGI50_2phase-harmonic_statistics",
  "aves_embeddings_minutes_nonzero_mean_2phase-harmonic_statistics",
  "aves_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_descriptive_statistics",
  "aves_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_descriptive_statistics",
  "aves_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_descriptive_statistics",
  "aves_embeddings_minutes_nonzero_mean_dawn-dusk_descriptive_statistics",
  "aves_embeddings_minutes_nonzero_mean_dawn-dusk_AGI90_hybrid",
  "aves_embeddings_minutes_nonzero_mean_dawn-dusk_AGI75_hybrid",
  "aves_embeddings_minutes_nonzero_mean_dawn-dusk_AGI50_hybrid",
  "aves_embeddings_minutes_nonzero_mean_dawn-dusk_hybrid"
)

# Execute and store result
aves_results <- getStats(model_files, inDir)
aves_results$folds$features <- "AVES embeddings"
aves_stats <- aves_results$summary
aves_stats$features <- "AVES embeddings"

# Access Wilcoxon tests
aves_wilcox <- aves_results$wilcox
aves_wilcox$features <- "AVES embeddings"

##################################
# Acoustic indices statistics
##################################
inDir <- file.path(baseDir,"models/rf_acoustic_indices")

model_files <- c(
  "acoustic_indices_minutes_AGI90_2phase-harmonic_statistics",
  "acoustic_indices_minutes_AGI75_2phase-harmonic_statistics",
  "acoustic_indices_minutes_AGI50_2phase-harmonic_statistics",
  "acoustic_indices_minutes_2phase-harmonic_statistics",
  "acoustic_indices_minutes_dawn-dusk_AGI90_descriptive_statistics",
  "acoustic_indices_minutes_dawn-dusk_AGI75_descriptive_statistics",
  "acoustic_indices_minutes_dawn-dusk_AGI50_descriptive_statistics",
  "acoustic_indices_minutes_dawn-dusk_descriptive_statistics",
  "acoustic_indices_minutes_dawn-dusk_AGI90_hybrid",
  "acoustic_indices_minutes_dawn-dusk_AGI75_hybrid",
  "acoustic_indices_minutes_dawn-dusk_AGI50_hybrid",
  "acoustic_indices_minutes_dawn-dusk_hybrid"
)

# Execute and store result
acoustic_indices_results <- getStats(model_files, inDir)
acoustic_indices_results$folds$features <- "Acoustic Indices"
acoustic_indices_stats <- acoustic_indices_results$summary
acoustic_indices_stats$features <- "Acoustic Indices"

# Access Wilcoxon tests
acoustic_indices_wilcox <- acoustic_indices_results$wilcox
acoustic_indices_wilcox$features <- "Acoustic Indices"

##################################
# Combine and export
##################################

final_stats <- rbind(birdnet_stats,
                     vggish_stats,
                     aves_stats,
                     acoustic_indices_stats)
outFile <-  file.path(baseDir,"data/rf_model_statistics_260410.csv")
write.csv(final_stats, outFile, row.names = FALSE)

final_wilcox <- rbind(birdnet_wilcox,
                      vggish_wilcox,
                      aves_wilcox,
                      acoustic_indices_wilcox)
outFile <- file.path(baseDir,"data/rf_model_wilcox_260410.csv")
write.csv(final_wilcox, outFile, row.names = FALSE)

##################################
# BirdNET PCA, UMAP, TPD statistics
##################################

inDir <- file.path(baseDir,"models/rf_pca")

# PCA results - descriptive statistics
model_files <- "birdnet_pca_minutes_nonzero_mean_dawn-dusk_descriptive_statistics" # PCA models
pca_results_desc <- getStats(model_files, inDir)
pca_results_desc$folds$features <- "PCA"
pca_stats_desc <- pca_results_desc$summary
pca_stats_desc$features <- "PCA"

# PCA results - harmonic statistics
model_files <- "birdnet_pca_minutes_nonzero_mean_2phase-harmonic_statistics" # PCA models
pca_results_harm <- getStats(model_files, inDir)
pca_results_harm$folds$features <- "PCA"
pca_stats_harm <- pca_results_harm$summary
pca_stats_harm$features <- "PCA"

param_cols <- data.frame(n_neighbors = rep("NA",12), min_dist = rep("NA",12))

# UMAP input directory
inDir <- file.path(baseDir,"models/rf_umap")

# UMAP parameter grid
param_grid <- expand.grid(
  n_neighbors = c(5, 20, 50, 100, 200, 500, 1000),
  min_dist = c(0.001, 0.1, 0.5, 0.8, 0.9, 0.99),
  stringsAsFactors = FALSE
)

# UMAP results - descriptive statistics
model_files <- NULL
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  nn <- params$n_neighbors
  md <- params$min_dist
  infile <- paste0("birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_nn",nn,"_mdist",md,"_descriptive_statistics")
  model_files <- c(model_files,infile)
  df_param <- data.frame(n_neighbors = rep(nn,3), min_dist = rep(md,3))
  param_cols <- rbind(param_cols,df_param)
}

# UMAP results
umap_results_desc <- getStats(model_files, inDir, richness_types = "AIrichness")
umap_results_desc$folds$features <- "UMAP"
umap_stats_desc <- umap_results_desc$summary
umap_stats_desc$features <- "UMAP"

# UMAP results - harmonic statistics
model_files <- NULL
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  nn <- params$n_neighbors
  md <- params$min_dist
  infile <- paste0("birdnet_embeddings_minute_nonzero_mean_dry-wet_pca_umap_nn",nn,"_mdist",md,"_2phase-harmonic_statistics")
  model_files <- c(model_files,infile)
  df_param <- data.frame(n_neighbors = rep(nn,3), min_dist = rep(md,3))
  param_cols <- rbind(param_cols,df_param)
}

# UMAP results
umap_results_harm <- getStats(model_files, inDir, richness_types = "AIrichness")
umap_results_harm$folds$features <- "UMAP"
umap_stats_harm <- umap_results_harm$summary
umap_stats_harm$features <- "UMAP"

##################################
# Combine and export
##################################

final_stats <- rbind(pca_stats_desc,pca_stats_harm,
                     umap_stats_desc, umap_stats_harm)
final_stats <- cbind(final_stats,param_cols)
  
outFile <- file.path(baseDir,"/data/rf_model_pca_umap_tpd_statistics_260410.csv")
write.csv(final_stats, outFile, row.names = FALSE)

# Write out folds detailed results
final_folds <- rbind(birdnet_results$folds,
                     vggish_results$folds,
                     aves_results$folds,
                     acoustic_indices_results$folds,
                     pca_results_desc$folds,
                     pca_results_harm$folds,
                     umap_results_desc$folds,
                     umap_results_harm$folds)
outFile <-  file.path(baseDir,"data/rf_model_statistics_folds_260410.csv")
write.csv(final_folds, outFile, row.names = FALSE)

