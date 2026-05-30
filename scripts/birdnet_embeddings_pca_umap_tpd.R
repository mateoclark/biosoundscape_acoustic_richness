# Purpose: TPD traits calculated on UMAP embeddings
# Matthew Clark, April 11, 2026

library(dplyr)
library(tidyr)
library(e1071)
library(readr)
library(TPD)

inDir <- "G:/Shared drives/BioSoundScape/Paper development/RQ1.1 Animal-acoustic diversity relationships"

# TPD alpha; A number between 0 and 1, indicating the proportion of the probability density function of each population to include
alpha_level <- 0.95

# Time periods to keep
periods <- c("DawnChorus","DuskChorus")

# UMAP parameter grid
param_grid <- expand.grid(
  #n_neighbors = c(5, 20, 50, 100, 200, 500, 1000),
  #min_dist     = c(0.001, 0.1, 0.5, 0.8, 0.9, 0.99),
  n_neighbors = c(20), # From best RF dry-wet model
  min_dist     = c(0.9), # From best RF dry-wet model
  stringsAsFactors = FALSE
)

# Custom time binning
binTimeOfDay <- function(dd){
  ifelse(dd >=   0 & dd <  500, "PredawnQuiet",
         ifelse(dd >= 500 & dd <  900, "DawnChorus",
                ifelse(dd >= 900 & dd < 1600, "Midday",
                       ifelse(dd >=1600 & dd < 2100, "DuskChorus",
                              ifelse(dd >=2100 & dd < 2400, "Nighttime", NA)))))
}

# List of UMAP variables for TPD
umap_cols <- c("UMAP5", "UMAP1", "UMAP3", "UMAP2") # From RF order of importance from RF dry-wet model

outDir <- paste0(inDir,"/umap_tpd_modeling/BirdNET_pca_umap_tpd_stats")

for (i in seq_len(nrow(param_grid))) {
  
  nn <- param_grid$n_neighbors[i]
  md <- param_grid$min_dist[i]
  message("Processing: n_neighbors=", nn, "  min_dist=", md)
  
  # Build input path
  infile <- sprintf(
    "birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_nn%d_mdist%g.csv",
    nn, md
  )
  inpath <- file.path(
    inDir,
    "umap_tpd_modeling",
    "BirdNET_pca_umap_embeddings",
    infile
  )
  
  # Load embeddings
  df <- read_csv(inpath, show_col_types = FALSE)
  
  # Filter to dawn/dusk chorus
  df <- df %>%
    mutate(time_period = binTimeOfDay(time_code)) %>%
    filter(time_period %in% periods)
  
  # Compute TPDs and REND metrics
  tpd_list     <- TPDs(species = df$SiteID, traits = df[, umap_cols], alpha = alpha_level)
  rend_metrics <- REND(TPDs = tpd_list)
  
  # Tidy output
  td_df <- tibble(
    SiteID      = names(rend_metrics$species$FRichness),
    FRichness   = rend_metrics$species$FRichness,
    FEvenness   = rend_metrics$species$FEvenness,
    FDivergence = rend_metrics$species$FDivergence
  )
  
  # Write CSV
  umap <- paste(umap_cols, collapse = "-")
  outfile <- sprintf(
    "birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_tpd_nn%d_mdist%g_%s.csv",
    nn, md,umap
  )
  write_csv(td_df, file.path(outDir, outfile))
}

message("All parameter combinations processed.")
