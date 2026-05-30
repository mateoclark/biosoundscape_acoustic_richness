# Purpose: Load non-zero minute-level BirdNET embeddings, reduce via PCA, then UMAP across all sites.
# Usage: adjust file paths as needed.
# Matthew Clark, April 11, 2026

# 1. Load libraries
library(data.table)
library(irlba)
library(uwot)
library(dplyr)
library(future.apply)

# UMAP parameter grid
param_grid <- expand.grid(
  n_neighbors = c(5, 20, 50, 100, 200, 500, 1000),
  min_dist = c(0.001, 0.1, 0.5, 0.8, 0.9, 0.99),
  stringsAsFactors = FALSE
)


# Custom time binning
binTimeOfDay <- function(dd){
  hrb <- ifelse(dd >= 000 & dd < 500, "PredawnQuiet",  # Pre-dawn/Quiet
                ifelse(dd >= 500 & dd < 900, "DawnChorus",  # Dawn chorus
                       ifelse(dd >= 900 & dd < 1600, "Midday",  # Midday
                              ifelse(dd >= 1600 & dd < 2100, "DuskChorus", # Dusk chorus
                                     ifelse(dd >= 2100 & dd < 2400, "Nighttime", NA)))))  # Night
  return(hrb)
}


# Set up parallel plan (adjust workers to match your system)
plan(multisession, workers = parallel::detectCores() - 1)

# Get list of CSV file paths
csv_files <- list.files("E:/active/project/birdnet/BirdNET_embeddings_minutes_nonzero_mean", pattern = "\\.csv$", full.names = TRUE)

# Read in parallel with fread
csv_list <- future_lapply(csv_files, function(file) {
  site_id <- substr(basename(file), 1, 15)
  tryCatch({
    dt <- fread(file)
    dt[, SiteID := site_id]
    dt
  }, error = function(e) {
    message("Failed: ", file)
    NULL
  })
})

# Desired time period
# NA = use all; "PredawnQuiet", "DawnChorus","Midday", "DuskChorus", "Nighttime"

#period <- NA # none
periods <- c("DawnChorus","DuskChorus") # dawn & dusk chorus only

# Combine all into a single data.table
dt_min <- rbindlist(csv_list, use.names = TRUE, fill = TRUE)

# Load site info
daac <- fread("G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/data/biosoundscape_sites_daac_250507.csv")
daac <- daac %>% select(SiteID,Campaign)

# Filter to sites in DAAC file
dt_min <- dt_min %>% filter(dt_min$SiteID %in% daac$SiteID)
dt_min <- left_join(dt_min, daac, by = "SiteID")
dt_min$period <- binTimeOfDay(dt_min$time_code)
dt_min <- dt_min %>% filter(period %in% periods)

# 3. Identify embedding columns
embed_cols <- dt_min %>% select(matches("fea")) %>% colnames()

# 7. PCA: reduce 1024 embeddings to 20 principal components
pc_n <- 20
mat <- as.matrix(dt_min[, ..embed_cols])
mat_scaled <- scale(mat)
pca_res <- prcomp_irlba(mat_scaled, n = pc_n, center = FALSE, scale. = FALSE)
pcs <- as.data.frame(pca_res$x)
setnames(pcs, paste0("PC", 1:pc_n))
pcs$SiteID <- dt_min$SiteID
pcs$date <- dt_min$date
pcs$time_code <- dt_min$time_code

# Save PCs
fwrite(pcs, "E:/active/project/birdnet/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca.csv")
#fwrite(pcs, "E:/active/project/birdnet/birdnet_embeddings_minute_nonzero_mean_dry-wet_pca.csv")

# Save PCA to RData
save(pca_res,file = "E:/active/project/birdnet/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_20pcs.RData")
#save(pca_res,file = "E:/active/project/birdnet/birdnet_embeddings_minute_nonzero_mean_dry-wet_pca_20pcs.RData")

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]

  nn <- params$n_neighbors
  md <- params$min_dist
  
  print(paste0("NN:",nn," Min Dist:",md))
  
  # 8. UMAP on PCA scores
  set.seed(42)
  umap_n <- 5 # number of UMAP dimensions
  umap_out <- uwot::umap(
    pcs[, paste0("PC", 1:pc_n)],
    n_neighbors   = nn,
    min_dist      = md,
    metric        = "euclidean",
    n_components  = umap_n,
    n_threads = 30
  )
  
  # Build column names dynamically for UMAP1 through UMAPn
  umap_dt <- data.table(SiteID = pcs$SiteID, date = pcs$date, time_code = pcs$time_code)
  umap_cols <- paste0("UMAP", 1:umap_n)
  umap_dt[, (umap_cols) := as.data.table(umap_out)]
  
  # Write out UMAP
  #outfile <- paste0("E:/active/project/birdnet/BirdNET_pca_umap_embeddings/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_nn",nn,"_mdist",md,".csv")
  outfile <- paste0("E:/active/project/birdnet/BirdNET_pca_umap_embeddings/birdnet_embeddings_minute_nonzero_mean_dry-wet_pca_umap_nn",nn,"_mdist",md,".csv")
  fwrite(umap_dt, outfile)
  
}


