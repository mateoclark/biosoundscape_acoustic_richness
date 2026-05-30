# Acoustic Index calculations
# Created: 29/Jan/19
# Updated: 27/Oct/2022
# Author: Colin Quinn
# Project: Soundscapes 2 Landscapes

# Modified from original SLURM version for single parallel processing on a single computer
# Matthew Clark, 10/27/2022

##########################################
# Load libraries
library(seewave)
library(soundecology)
library(tuneR)
library(crayon)
library(dplyr)
library(tidyr)
library(reshape2)
library(doParallel)
library(parallelly)

# Set parameters

#ncores <- detectCores()-2 # number of parallel processing cores to use
ncores <- 15

wd_out <- "E:/archive/project/biosoundscape/AcousticIndicesCampaign2/"    # results location
filenames_csv <- "E:/archive/project/biosoundscape/biosoundscape_campaign2_site_recordings_fulllist_240102.csv"  # csv listing toDo wavs
file_dir = "E:/archive/project/biosoundscape/to_process"  # directory with .wav files

print(paste("Results dir:", wd_out))
print(paste("CSV with ToDo files = ",filenames_csv))
print(paste("Directory with .wav files = ",file_dir))

# use csv to get list of filenames to do
file_names <- read.csv(filenames_csv)$files
file_list <- as.vector(paste0(file_dir,"/",file_names))
n <- length(file_list)

print("\nLoaded audio files list...")
print(paste0("Length of audio file list:", n))
print(paste0("First file = ",file_list[1]))
print(paste0("Last file = ",file_list[n]))


##########################################
# ACOUSTIC INDEX CALCULATIONS #
# Read in wav files and string soundecology and seewave metrics with parameters influenced by Eldridge et al. 2018
  # outputs a csv with acoustic indices for each wav file
  # this combination of AI was chosen based on Eldridge et al. 2018, Mammides et al. 2017, Moreno-Gomez et al. 2019, Sueur 2018
  
index_list <- list()
error_log <- paste0("E:/archive/project/biosoundscape/error_log/error_log-acoustic_indices.txt")


workers = rep("localhost",ncores)
port = seq(from=9000, to = 9000+ncores)
cl1 <- parallelly::makeClusterPSOCK(workers, port=port, autoStop = TRUE)
registerDoParallel(cl1)
foreach(i = 1:n,.packages=c("seewave","soundecology","tuneR","crayon",'dplyr','tidyr','reshape2')) %dopar% {
  temp_wav <- file_list[i]
  temp_wav_sans_ext <- tools::file_path_sans_ext(temp_wav) # remove extension	
  temp_out_wav <- strsplit(temp_wav_sans_ext, '/')[[1]] # split file path
  temp_out_wav <- tail(temp_out_wav, n = 1) # file name only
  out_file <- paste0(wd_out, temp_out_wav,'.csv') # csv output
  
  # check if wav indices have been computed already
  if (file.exists(out_file)){  
    cat('File exists: skipping calculations')
    
  }else{  	  
    tryCatch({ 
      # current file
      cat(paste0(i,"- Working on:", temp_out_wav))
      
      # read in wav files from wd
      soundfile <- readWave(temp_wav)
      samp_freq <- soundfile@samp.rate # 44100 for LG; 48000 for AM
      
      # calculate each acoustic index 
      # Soundecology derived indices
      ndsi_indices <- soundecology::ndsi(soundfile, fft_w = 1024, anthro_min = 1000, anthro_max = 2000, bio_min = 2000, bio_max = 10000)
      
      ndsi <- ndsi_indices$ndsi_left
      
      anthro <- ndsi_indices$anthrophony_left
      
      bio <- ndsi_indices$biophony_left
      
      bi <- soundecology::bioacoustic_index(soundfile, min_freq = 2000, max_freq = 10000)$left_area 
      
      aei <- soundecology::acoustic_evenness(soundfile, max_freq = 10000, db_threshold = -50, freq_step = 1000)$aei_left
      
      adi <- soundecology::acoustic_diversity(soundfile, max_freq = 10000, db_threshold = -50, freq_step = 1000, shannon = TRUE)$adi_left
      
      
      # SEEWAVE derived indices
      spec = spec(soundfile, f = samp_freq, plot = FALSE)
      
      env <- seewave::env(soundfile, plot = FALSE)
      
      aci <- seewave::ACI(soundfile, flim = c(1,10))
      
      h <- seewave::H(soundfile, wl = 512)
      
      ht <- seewave::th(env)
      
      sh <- seewave::sh(spec)
      
      rms <- seewave::rms(env)
      
      zcr <- seewave::zcr(soundfile, plot = FALSE)
      
      M <- seewave::M(soundfile)
      
      rough <- seewave::roughness(meanspec(soundfile, plot = FALSE)[,2])
      
      sfm <- seewave::sfm(spec)
      
      rugo <- seewave::rugo(soundfile@left/max(soundfile@left))
      
      ##########################################
      # Extract datestamp inputs for dataframe
      sitename <- strsplit(temp_out_wav,'_')[[1]][1]
      YYMMDD <- strsplit(temp_out_wav,'_')[[1]][2]
      YY <- substr(YYMMDD,1,2)
      YYYY <- paste0('20',YY)		  
      MM <- substr(YYMMDD,3,4)
      DD <- substr(YYMMDD,5,6)
      hh_mm <- strsplit(temp_out_wav,'_')[[1]][4]
      hh <- substr(hh_mm,1,2)
      mm <- substr(hh_mm,4,5)
      DDhh <- paste0(DD,hh)
      hhmm <- paste0(hh,mm)
      
      ##########################################
      # organize derived products for each file in loop
      temp_df <- data.frame('path' = file_list[i], 
                            'file' = temp_out_wav,
                            
                            'ACI' = aci,
                            'ADI' = adi,
                            'AEI' = aei,
                            'NDSI' = ndsi,
                            'NDSI_A' = anthro,
                            'NDSI_B' = bio,
                            'BI' = bi,
                            'H' = h,
                            'Ht' = ht,
                            'Hs' = sh, # signal noisiness/pureness (similar to sfm)
                            'M' = M,
                            'R' = rough,
                            'sfm' = sfm, # spectral flatness measure (noisy -> 1; pure tone -> 0)
                            'rugo' = rugo, # noisy is higher
                            'zcr_mean' = mean(zcr[,2]),
                            
                            'YYYY' = YYYY, 
                            'MM' = MM, 
                            'DD' = DD, 
                            'hh' = hh, 
                            'mm' = mm, 
                            'DDhh' = DDhh, 
                            'hhmm' = hhmm)
      
      # Data Output #########################################
      # save the csv
      print(paste("SAVING OUTPUT FILE:", out_file))
      write.csv(temp_df, file = out_file, row.names = FALSE)  
      
    }, error = function(e){
      cat("ERROR :",conditionMessage(e), "\n")
      write(toString(file_list[i]), error_log, append=TRUE)
    })  
  }
  } # End foreach loop
parallelly::autoStopCluster(cl1)
gc()
