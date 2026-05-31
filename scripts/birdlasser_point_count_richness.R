library(dplyr)

df <- read.csv("G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/2024-03-28 All BioSCape Point Counts Spatial Attr.csv")

richness <- df %>%
  group_by(Location_ID) %>%
  filter(Species != "Unidentified") %>%
  summarise(PCrichness = n_distinct(Species)) %>%
  rename(SiteID = Location_ID)

write.csv(richness,"G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/birdlasser_point_count_richness_240328.csv", row.names = F)