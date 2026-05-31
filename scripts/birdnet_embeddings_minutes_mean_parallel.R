library(dplyr)
library(readr)
library(doParallel)
library(foreach)
library(tidyr)
library(lubridate)
library(stringr)

# Path setup
inDir <- "Z:/BirdNET_embeddings"
files <- list.files(path = inDir, pattern = "\\.txt$", full.names = TRUE)

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

# Parallel setup
numCores <- 20
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Output directory
outDir <- paste0(inDir, "_minutes_mean/")
if (!dir.exists(outDir)) dir.create(outDir)

# Main parallel loop
foreach(s = sites, .packages = c("readr", "dplyr", "tidyr","lubridate","stringr")) %dopar% {
#for (s in sites){
  print(s)
  tryCatch({
    siteFiles <- files[grep(s, files)]
    if (length(siteFiles) == 0) return(NULL)
    
    minute_vectors <- list()
    times <- numeric()
    
    for (f in siteFiles) {
      
      # Get time and date information from file name
      t <- extract_time(f)
      
      # Read embeddings
      fread <- gsub("\t", ",", readLines(f))
      conn <- textConnection(fread)
      df <- read.csv(conn, header = FALSE, stringsAsFactors = FALSE)
      close(conn)
      colnames(df) <- c("start", "stop", cn)
      
      # Compute mean of values for each embedding feature (including zeros)
      minute_mean <- colMeans(df[, 3:1026], na.rm = TRUE)
      
      minute_vectors[[length(minute_vectors) + 1]] <- minute_mean
      times <- rbind(times,t)
    }
    
    X <- do.call(rbind, minute_vectors)
    colnames(X) <- cn
    final_df <- cbind(times,X)
    
    # Write output CSV
    outFile <- paste0(outDir, s, "_birdnet_embeddings_minutes_mean.csv")
    write.csv(final_df, outFile, row.names = FALSE)
    
    rm(minute_vectors, final_df, times)
    gc()
  }, error = function(e) {
    message(paste0("Error at site: ", s, " --> ", e$message))
  })
}

# Stop parallel cluster
stopCluster(cl)
