library(dplyr)
library(parallel)
library(tidyr)
library(lubridate)
library(stringr)
library(reticulate)

# Path setup
inDir <- "Y:/aves_campaign"
files <- list.files(path = inDir, pattern = "\\.npy$", full.names = TRUE)

# Extract site ID
extract_site <- function(filename) substr(basename(filename), 1, 15)

# Extract time information from file
extract_time <- function(file) {
  # Extract datetime from filename (e.g., s2lam165_230725_2023-08-01_11-10)
  f <- basename(file)
  datetime_str <- stringr::str_extract(f, "\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}")
  datetime <- lubridate::ymd_hm(datetime_str)
  date <- as.Date(datetime)
  hh <- hour(datetime)
  mm <- minute(datetime)
  time_code <- hh * 100 + mm
  d <- data.frame(time_code=time_code, date=date)
  return(d)
}

# Feature names
cn <- paste0("fea", 1:1024)
sites <- unique(sapply(files, extract_site))

# Output directory
outDir <- paste0(inDir, "_minutes_mean/")
if (!dir.exists(outDir)) dir.create(outDir)

# List of site IDs or file groups
#n_cores <- detectCores() - 1
n_cores <- 30
cl <- makeCluster(n_cores)

np <- import("numpy")  # Load in main session, but re-import inside workers too

# Loop through sites
for (s in sites) {
  
  outFile <- file.path(outDir, paste0(s, "_aves_embeddings_minutes_mean.csv"))
  
  if (!file.exists(outFile)){
    cat("Processing site:", s, "\n")
    tryCatch({
      siteFiles <- files[grep(s, files)]
      if (length(siteFiles) == 0) next
      
      # Setup cluster
      cl <- makeCluster(n_cores)
      clusterExport(cl, varlist = c("siteFiles", "extract_time"), envir = environment())
      
      # Load required packages in each worker
      clusterEvalQ(cl, {
        library(reticulate)
        library(lubridate)   # <- This line fixes the hour() error
        np <- import("numpy")
        return(TRUE)
      })
      
      # Process each file in parallel
      results <- parLapply(cl, siteFiles, function(f) {
        np <- import("numpy")
        extract_time <- get("extract_time", envir = .GlobalEnv)
        
        t <- extract_time(f)
        
        array <- np$load(f)
        
        r_array <- py_to_r(array)
        if (length(r_array) == 0 || any(dim(r_array) == 0)) return(NULL)
        df <- data.frame(r_array)
        
        # Mean of values
        minute_mean <- colMeans(df, na.rm = TRUE)
        
        list(mean = minute_mean, time = t)
      })
      
      
      stopCluster(cl)
      
      # Combine results
      minute_vectors <- lapply(results, `[[`, "mean")
      times <- do.call(rbind, lapply(results, `[[`, "time"))
      
      X <- do.call(rbind, minute_vectors)
      colnames(X) <- cn
      final_df <- cbind(times, X)
      
      # Save CSV
      write.csv(final_df, outFile, row.names = FALSE)
      
      rm(minute_vectors, final_df, times, results)
      gc()
      
    }, error = function(e) {
      message(sprintf("Error at site %s: %s", s, e$message))
    })
  } else {
    cat("Skipping site:", s, "\n")
  }
}

