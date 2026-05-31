library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(stringr)
library(reticulate)

# Path setup
inDir <- "X:"

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

# read pickle function
readPickle <- function(file){
  pickleData <- data.frame(read_pickle_file(file)$raw_audioset_feats_960ms)
  return(pickleData)
}

# Load pickle reader
use_python("C:/ProgramData/anaconda3/python.exe", required = TRUE)
source_python(paste0(inDir,"/pickle_reader.py"))

# Input directories
inDir1 <- paste0(inDir,"/VGGishFeaturesCampaign1")
inDir2 <- paste0(inDir,"/VGGishFeaturesCampaign2")

# Get pickle files
files1 <- list.files(path = inDir1, pattern = "\\.pickle$", full.names = TRUE)
files2 <- list.files(path = inDir2, pattern = "\\.pickle$", full.names = TRUE)
files <- c(files1, files2)

# Feature names
cn <- paste0("fea", 1:128)
sites <- unique(sapply(files, extract_site))

# Output directory
outDir <- paste0(inDir, "/vggish_minutes_mean/")
if (!dir.exists(outDir)) dir.create(outDir)

# Main loop
for (s in sites){
  print(s)
  tryCatch({
    siteFiles <- files[grep(s, files)]
    if (length(siteFiles) == 0) return(NULL)
    
    minute_vectors <- list()
    times <- numeric()
    df <- numeric()
    
    outFile <- paste0(outDir, s, "_vggish_embeddings_minutes_mean.csv")
    
    if (!file.exists(outFile)){
      for (f in siteFiles) {
        
        # Get time and date information from file name
        t <- extract_time(f)
        
        # Read embeddings
        data <- readPickle(f)
        df <- rbind(df, data)
        
        # Compute mean of values for each embedding feature (including zeros)
        minute_mean <- colMeans(df, na.rm = TRUE)
        
        minute_vectors[[length(minute_vectors) + 1]] <- minute_mean
        times <- rbind(times,t)
      }
      
      X <- do.call(rbind, minute_vectors)
      colnames(X) <- cn
      final_df <- cbind(times,X)
      
      # Write output CSV
      write.csv(final_df, outFile, row.names = FALSE)
      
      rm(minute_vectors, final_df, times)
      gc()
    }
  }, error = function(e) {
    message(paste0("Error at site: ", s, " --> ", e$message))
  })
}

