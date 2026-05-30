# Purpose: Makes graphs and some light statistics found in the paper, including additional exploratory graphs.
# Matthew Clark, May 29, 2026

library(ggplot2)
library(ggpattern)
library(ggpubr)
library(dplyr)
library(readr)
library(forcats)
library(tidyr)
library(hexbin)
library(viridis)
library(RColorBrewer)
library(lubridate)
library(tidyverse)
library(mgcv)
library(nlme)
library(rcompanion)
library(Ternary)
library(viridisLite)
library(scales)
library(grid)
library(cowplot)
library(patchwork)
library(emmeans)

# input directory (change to local version)
inDir <- "G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships/github"

# ========================================================================
# Point Count richness R2 & RMSE - Features grouped - predictors & seasons
# ========================================================================

# Input file
inFile <- file.path(inDir,"data/rf_model_statistics_260410.csv")

# Read CSV with statistics
df <- read.csv(inFile)

# Prepare data
df_plot <- df %>%
  filter(
    richness_type == "PCrichness",
    AGIfilter == "None",
    !(predictors == "descriptive" & periods == "None")
  ) %>%
  mutate(
    season_predictors = paste(season, predictors, sep = " | "),
    season_predictors = recode(season_predictors,
                               "wet | descriptive" = "Wet | Descriptive",
                               "wet | 2phase-harmonic" = "Wet | Harmonics",
                               "wet | hybrid" = "Wet | Hybrid",
                               "dry | descriptive" = "Dry | Descriptive",
                               "dry | 2phase-harmonic" = "Dry | Harmonics",
                               "dry | hybrid" = "Dry | Hybrid",
                               "drywet | descriptive" = "Wet-Dry | Descriptive",
                               "drywet | 2phase-harmonic" = "Wet-Dry | Harmonics",
                               "drywet | hybrid" = "Wet-Dry | Hybrid"
    ),
    features = factor(features, 
                      levels = c("Acoustic Indices", "VGGish embeddings", "BirdNET embeddings", "AVES embeddings"),
                      labels = c("Indices", "VGGish", "BirdNET", "BirdAVES")),
    season_predictors = factor(season_predictors, levels = c(
      "Wet | Descriptive", "Wet | Harmonics", "Wet | Hybrid",
      "Dry | Descriptive", "Dry | Harmonics", "Dry | Hybrid",
      "Wet-Dry | Descriptive", "Wet-Dry | Harmonics", "Wet-Dry | Hybrid"
    ))
  )

# Color palette with final labels
pal_wet <- RColorBrewer::brewer.pal(8, "Blues")
pal_dry <- RColorBrewer::brewer.pal(8, "Oranges")
pal_wetdry <- RColorBrewer::brewer.pal(8, "Greens")
palette <- c(
  "Wet | Descriptive" = pal_wet[4],
  "Wet | Harmonics" = pal_wet[5],
  "Wet | Hybrid" = pal_wet[6],
  "Dry | Descriptive" = pal_dry[4],
  "Dry | Harmonics" = pal_dry[5],
  "Dry | Hybrid" = pal_dry[6],
  "Wet-Dry | Descriptive" = pal_wetdry[4],
  "Wet-Dry | Harmonics" = pal_wetdry[5],
  "Wet-Dry | Hybrid" = pal_wetdry[6]
)

# Plot A: R2
plot_r2 <- ggplot(df_plot, aes(x = features, y = Rsq_mean, fill = season_predictors)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Rsq_mean - Rsq_se, ymax = Rsq_mean + Rsq_se),
                position = position_dodge(0.8), width = 0.2) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  scale_fill_manual(values = palette, guide = guide_legend(title = "Predictors")) +
  labs(
    x = NULL,
    y = expression("CV R"^2)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_blank()
  )

# Plot B: RMSE
plot_rmse <- ggplot(df_plot, aes(x = features, y = RMSE_mean, fill = season_predictors)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = RMSE_mean - RMSE_se, ymax = RMSE_mean + RMSE_se),
                position = position_dodge(0.8), width = 0.2) +
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) +
  scale_fill_manual(values = palette, guide = "none") +
  labs(
    x = NULL,
    y = "CV RMSE"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_blank()
  )

# Combine plots with equal height
combined_plot <- ggarrange(
  plot_r2, plot_rmse,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "right",
  heights = c(1, 1)  # equal vertical space
)

combined_plot <- annotate_figure(
  combined_plot,
  top = text_grob("Point Count Richness", face = "bold", size = 14)
)

# Show and save plot
dev.new(); print(combined_plot)

outFile <- file.path(inDir, "figures/PCrichness_features_r2_rmse.png")
ggsave(outFile, combined_plot, width = 6, height = 4, units = "in", dpi = 600, bg = "white")

# highest R2
df_plot %>%
  slice_max(Rsq_mean, n = 1)

# lowest R2
df_plot %>%
  slice_min(Rsq_mean, n = 1)

df_plot %>% filter(features == "Indices")
df_plot %>% filter(features == "Indices" & season == "dry")

# BirdNET summary
birdnet_summary <- df_plot %>%
  filter(features == "BirdNET") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
birdnet_summary

diffs_season <- df_plot %>%
  #filter(features == "BirdNET") %>%
  select(season, Rsq_mean, predictors, periods, file_name) %>%
  pivot_wider(
    names_from = season,
    values_from = Rsq_mean
  ) %>%
  mutate(
    diff_dry_wet    = dry - wet,
    diff_dry_drywet = dry - drywet
  )
diffs_season

diffs_predictors <- df_plot %>%
  select(richness_type, features, season, predictors, Rsq_mean) %>%
  pivot_wider(
    names_from = predictors,
    values_from = Rsq_mean
  ) %>%
  mutate(
    diff_desc_harm = descriptive - `2phase-harmonic`,
    diff_desc_hybrid = descriptive - hybrid
  )
diffs_predictors

# VGGish summary
vggish_summary <- df_plot %>%
  filter(features == "VGGish") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
vggish_summary

# BirdAVES summary
aves_summary <- df_plot %>%
  filter(features == "BirdAVES") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
aves_summary

# Acoustic Indices summary
acoustic_indices_summary <- df_plot %>%
  filter(features == "Indices") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
acoustic_indices_summary

# Summary of differences in R2 from acoustic indices vs CNN models
df_comp <- df_plot %>%
  select(richness_type, season, predictors, periods, AGIfilter,
         features, Rsq_mean)

df_wide <- df_comp %>%
  pivot_wider(
    names_from = features,
    values_from = Rsq_mean
  )

df_diff <- df_wide %>%
  mutate(
    diff_Indices_BirdNET = Indices - BirdNET,
    diff_Indices_BirdAVES = Indices - BirdAVES,
    diff_Indices_VGGish = Indices - VGGish
  )

range_summary <- df_diff %>%
  summarise(
    min_Indices_BirdNET = min(diff_Indices_BirdNET, na.rm = TRUE),
    max_Indices_BirdNET = max(diff_Indices_BirdNET, na.rm = TRUE),
    
    min_Indices_BirdAVES = min(diff_Indices_BirdAVES, na.rm = TRUE),
    max_Indices_BirdAVES = max(diff_Indices_BirdAVES, na.rm = TRUE),
    
    min_Indices_VGGish = min(diff_Indices_VGGish, na.rm = TRUE),
    max_Indices_VGGish = max(diff_Indices_VGGish, na.rm = TRUE)
  )

range_summary

# Test statistical differences
inFile <- file.path(inDir,"data/rf_model_statistics_folds_260410.csv")
df <- read.csv(inFile) %>% filter(richness_type == "PCrichness" & 
                                    AGIfilter == "None" & 
                                    !features %in% c("PCA","UMAP"))

ggplot(df, aes(x = features, y = Rsquared, fill = features)) +
  geom_boxplot() +
  facet_wrap(~ season) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# ANOVA of effect of features and season on R2
model <- lm(Rsquared ~ features * season, data = df)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ features)
emmeans(model, pairwise ~ season)

# ANOVA of effect of predictor on R2
model <- lm(Rsquared ~ predictors, data = df)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ predictors)

# ANOVA of effect of predictor on R2, just embeddings
df_embed <- df %>% filter(features != "Acoustic Indices")
model <- lm(Rsquared ~ predictors, data = df_embed)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ predictors)


# ========================================================================
# AI richness R2 & RMSE - Features grouped - predictors & seasons
# ========================================================================

# Input file
inFile <- file.path(inDir,"data/rf_model_statistics_260410.csv")

# Read CSV with statistics
df <- read.csv(inFile)

# Prepare data
df_plot <- df %>%
  filter(
    richness_type == "AIrichness",
    AGIfilter == "None",
    !(predictors == "descriptive" & periods == "None")
  ) %>%
  mutate(
    season_predictors = paste(season, predictors, sep = " | "),
    season_predictors = recode(season_predictors,
                               "wet | descriptive" = "Wet | Descriptive",
                               "wet | 2phase-harmonic" = "Wet | Harmonics",
                               "wet | hybrid" = "Wet | Hybrid",
                               "dry | descriptive" = "Dry | Descriptive",
                               "dry | 2phase-harmonic" = "Dry | Harmonics",
                               "dry | hybrid" = "Dry | Hybrid",
                               "drywet | descriptive" = "Wet-Dry | Descriptive",
                               "drywet | 2phase-harmonic" = "Wet-Dry | Harmonics",
                               "drywet | hybrid" = "Wet-Dry | Hybrid"
    ),
    features = factor(features, 
                      levels = c("Acoustic Indices", "VGGish embeddings", "BirdNET embeddings", "AVES embeddings"),
                      labels = c("Indices", "VGGish", "BirdNET", "BirdAVES")),
    season_predictors = factor(season_predictors, levels = c(
      "Wet | Descriptive", "Wet | Harmonics", "Wet | Hybrid",
      "Dry | Descriptive", "Dry | Harmonics", "Dry | Hybrid",
      "Wet-Dry | Descriptive", "Wet-Dry | Harmonics", "Wet-Dry | Hybrid"
    ))
  )


# Color palette with final labels
pal_wet <- RColorBrewer::brewer.pal(8, "Blues")
pal_dry <- RColorBrewer::brewer.pal(8, "Oranges")
pal_wetdry <- RColorBrewer::brewer.pal(8, "Greens")
palette <- c(
  "Wet | Descriptive" = pal_wet[4],
  "Wet | Harmonics" = pal_wet[5],
  "Wet | Hybrid" = pal_wet[6],
  "Dry | Descriptive" = pal_dry[4],
  "Dry | Harmonics" = pal_dry[5],
  "Dry | Hybrid" = pal_dry[6],
  "Wet-Dry | Descriptive" = pal_wetdry[4],
  "Wet-Dry | Harmonics" = pal_wetdry[5],
  "Wet-Dry | Hybrid" = pal_wetdry[6]
)

# Plot A: R2
plot_r2 <- ggplot(df_plot, aes(x = features, y = Rsq_mean, fill = season_predictors)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Rsq_mean - Rsq_se, ymax = Rsq_mean + Rsq_se),
                position = position_dodge(0.8), width = 0.2) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.1)) +
  scale_fill_manual(values = palette, guide = guide_legend(title = "Predictors")) +
  labs(
    x = NULL,
    y = expression("CV R"^2)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_blank()
  )

# Plot B: RMSE
plot_rmse <- ggplot(df_plot, aes(x = features, y = RMSE_mean, fill = season_predictors)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = RMSE_mean - RMSE_se, ymax = RMSE_mean + RMSE_se),
                position = position_dodge(0.8), width = 0.2) +
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) +
  scale_fill_manual(values = palette, guide = "none") +
  labs(
    x = NULL,
    y = "CV RMSE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_blank()
  )

# Combine plots with equal height
combined_plot <- ggarrange(
  plot_r2, plot_rmse,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "right",
  heights = c(1, 1)  # equal vertical space
)

combined_plot <- annotate_figure(
  combined_plot,
  top = text_grob("AI-based Richness", face = "bold", size = 14)
)


# Show and save plot
dev.new(); print(combined_plot)

outFile <- file.path(inDir, "figures/AIrichness_features_r2_rmse.png")
ggsave(outFile, combined_plot, width = 6, height = 4, units = "in", dpi = 600, bg = "white")

# highest R2
df_plot %>%
  slice_max(Rsq_mean, n = 1)

# lowest R2
df_plot %>%
  slice_min(Rsq_mean, n = 1)

# BirdNET summary
birdnet_summary <- df_plot %>%
  filter(features == "BirdNET") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
birdnet_summary

diffs_season <- df_plot %>%
  #filter(features == "BirdNET") %>%
  select(season, Rsq_mean, predictors, periods, file_name) %>%
  pivot_wider(
    names_from = season,
    values_from = Rsq_mean
  ) %>%
  mutate(
    diff_dry_wet    = dry - wet,
    diff_dry_drywet = dry - drywet
  )
diffs_season

diffs_predictors <- df_plot %>%
  select(richness_type, features, season, predictors, Rsq_mean) %>%
  pivot_wider(
    names_from = predictors,
    values_from = Rsq_mean
  ) %>%
  mutate(
    diff_desc_harm = descriptive - `2phase-harmonic`,
    diff_desc_hybrid = descriptive - hybrid
  )
diffs_predictors

# VGGish summary
vggish_summary <- df_plot %>%
  filter(features == "VGGish") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
vggish_summary

# BirdAVES summary
aves_summary <- df_plot %>%
  filter(features == "BirdAVES") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
aves_summary

# Acoustic Indices summary
acoustic_indices_summary <- df_plot %>%
  filter(features == "Indices") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
acoustic_indices_summary

# Summary of differences in R2 from acoustic indices vs CNN models
df_comp <- df_plot %>%
  select(richness_type, season, predictors, periods, AGIfilter,
         features, Rsq_mean)

df_wide <- df_comp %>%
  pivot_wider(
    names_from = features,
    values_from = Rsq_mean
  )

df_diff <- df_wide %>%
  mutate(
    diff_Indices_BirdNET = Indices - BirdNET,
    diff_Indices_BirdAVES = Indices - BirdAVES,
    diff_Indices_VGGish = Indices - VGGish
  )

range_summary <- df_diff %>%
  summarise(
    min_Indices_BirdNET = min(diff_Indices_BirdNET, na.rm = TRUE),
    max_Indices_BirdNET = max(diff_Indices_BirdNET, na.rm = TRUE),
    
    min_Indices_BirdAVES = min(diff_Indices_BirdAVES, na.rm = TRUE),
    max_Indices_BirdAVES = max(diff_Indices_BirdAVES, na.rm = TRUE),
    
    min_Indices_VGGish = min(diff_Indices_VGGish, na.rm = TRUE),
    max_Indices_VGGish = max(diff_Indices_VGGish, na.rm = TRUE)
  )

range_summary

# Test statistical differences
inFile <- file.path(inDir,"data/rf_model_statistics_folds_260410.csv")
df <- read.csv(inFile) %>% filter(richness_type == "AIrichness" & 
                                    AGIfilter == "None" & 
                                    !features %in% c("PCA","UMAP"))

ggplot(df, aes(x = features, y = Rsquared, fill = features)) +
  geom_boxplot() +
  facet_wrap(~ season) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# ANOVA of effect of features and season on R2
model <- lm(Rsquared ~ features * season, data = df)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ features)
emmeans(model, pairwise ~ season)

# ANOVA of effect of predictor on R2
model <- lm(Rsquared ~ predictors, data = df)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ predictors)

# ANOVA of effect of predictor on R2, just embeddings
df_embed <- df %>% filter(features != "Acoustic Indices")
model <- lm(Rsquared ~ predictors, data = df_embed)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ predictors)

# ========================================================================
# AI richness Diff AGI R2 & RMSE - Features grouped - predictors & seasons
# ========================================================================

## Load t-test output
df <- read_csv(paste0(inDir,"/data/rf_model_wilcox_260410.csv"))

# Prepare data
df_plot <- df %>%
  filter(richness_type == "AIrichness") %>%
  mutate(
    predictors = factor(predictors, levels = c("descriptive", "2phase-harmonic", "hybrid"),
                        labels = c("Descriptive", "Harmonics", "Hybrid")),
    season = factor(season, levels = c("wet", "dry", "drywet"),
                    labels = c("Wet", "Dry", "Wet-Dry")),
    features = factor(features, levels = c("Acoustic Indices", "VGGish embeddings", "BirdNET embeddings", "AVES embeddings"),
                                labels = c("Indices", "VGGish", "BirdNET", "BirdAVES")),
    AGIfilter = factor(comparison, levels = c("AGI50", "AGI75", "AGI90"),
                       labels = c("≥50%", "≥75%", "≥90%")),
    pattern = ifelse(Rsq_p < 0.05, "stripe", "none"),
    pattern = factor(pattern, levels = c("none", "stripe"))
  )

# Define fill colors
pal <- RColorBrewer::brewer.pal(8, "Set2")
fill_colors <- c("Descriptive" = pal[1], "Harmonics" = pal[2], "Hybrid" = pal[3])

# Create plot with ggpattern
p <- ggplot(df_plot, aes(x = AGIfilter, y = Rsq_diff,
                         fill = predictors, pattern = pattern)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = NA,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.5
  ) +
  scale_y_continuous(
    limits = c(-0.11, 0.11),
    breaks = seq(-0.1, 0.1, 0.05),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  facet_grid(rows = vars(features), cols = vars(season)) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"), guide = "none") +
  labs(
    x = "AGI Filter",
    y = expression(Delta~CV~R^2)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "gray90", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.8)
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/AIrichness_r2_agi_filter_difference.png")
ggsave(outFile, p, width = 6, height = 6, units = "in", dpi = 600, bg = "white")

# Statistics
summary(df_plot$Rsq_diff)
summary(df_plot$RMSE_diff)


# Test statistical differences
inFile <- file.path(inDir,"data/rf_model_statistics_folds_260410.csv")
df <- read.csv(inFile) %>% filter(richness_type == "AIrichness" &
                                  !features %in% c("PCA","UMAP"))

# ANOVA of effect of features and AGI filter on R2
model <- lm(Rsquared ~ features * AGIfilter, data = df)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ AGIfilter)
emmeans(model, pairwise ~ AGIfilter | features)

# ========================================================================
# PC richness Diff AGI R2 & RMSE - Features grouped - predictors & seasons
# ========================================================================

## Load t-test output
df <- read_csv(file.path(inDir,"data/rf_model_wilcox_260410.csv"))

# Prepare data
df_plot <- df %>%
  filter(richness_type == "PCrichness") %>%
  mutate(
    predictors = factor(predictors, levels = c("descriptive", "2phase-harmonic", "hybrid"),
                        labels = c("Descriptive", "Harmonics", "Hybrid")),
    season = factor(season, levels = c("wet", "dry", "drywet"),
                    labels = c("Wet", "Dry", "Wet-Dry")),
    features = factor(features, levels = c("Acoustic Indices", "VGGish embeddings", "BirdNET embeddings", "AVES embeddings"),
                      labels = c("Indices", "VGGish", "BirdNET", "BirdAVES")),
    AGIfilter = factor(comparison, levels = c("AGI50", "AGI75", "AGI90"),
                       labels = c("≥50%", "≥75%", "≥90%")),
    pattern = ifelse(Rsq_p < 0.05, "stripe", "none"),
    pattern = factor(pattern, levels = c("none", "stripe"))
  )

# Define fill colors
pal <- RColorBrewer::brewer.pal(8, "Set2")
fill_colors <- c("Descriptive" = pal[1], "Harmonics" = pal[2], "Hybrid" = pal[3])

# Create plot with ggpattern
p <- ggplot(df_plot, aes(x = AGIfilter, y = Rsq_diff,
                         fill = predictors, pattern = pattern)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = NA,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.5
  ) +
  scale_y_continuous(
    limits = c(-0.20, 0.1),
    breaks = seq(-0.20, 0.1, 0.05),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  facet_grid(rows = vars(features), cols = vars(season)) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"), guide = "none") +
  labs(
    x = "AGI Filter",
    y = expression(Delta~CV~R^2)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "gray90", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.8)
  )
# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir,"figures/PCrichness_r2_agi_filter_difference.png")
ggsave(outFile, p, width = 6, height = 6, units = "in", dpi = 600, bg = "white")

# ========================================================================
# AI richness PCA - descriptive + harmonic statistics
# ========================================================================

# Load and filter
df_pca <- read_csv(file.path(inDir,"data/rf_model_pca_umap_tpd_statistics_260410.csv")) %>%
  filter(
    richness_type == "AIrichness",
    features != "UMAP",
  ) %>%
  mutate(
    season = factor(season, levels = c("wet", "dry", "drywet")),
    predictors = factor(predictors, levels = c("descriptive", "2phase-harmonic"))
  )

# Set up season colors
pal <- brewer.pal(3, "Set2")
season_colors <- c(
  "wet" = pal[3],
  "dry" = pal[2],
  "drywet" = pal[1]
)

# Plot
plot_pca_r2 <- ggplot(df_pca, aes(
  x = predictors,
  y = Rsq_mean,
  fill = season
)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = Rsq_mean - Rsq_se,
      ymax = Rsq_mean + Rsq_se
    ),
    position = position_dodge(0.8),
    width = 0.2
  ) +
  scale_fill_manual(
    values = season_colors,
    labels = c("Wet", "Dry", "Wet-Dry")
  ) +
  scale_x_discrete(
    labels = c(
      "descriptive" = "Descriptive",
      "2phase-harmonic" = "Harmonics"
    )
  ) +
  labs(
    y = expression("Mean CV " * R^2),
    fill = "Season",
    x = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom"
  )
# Show plot
#dev.new(); print(plot_pca_r2)

plot_pca_rmse <- ggplot(df_pca, aes(
  x = predictors,
  y = RMSE_mean,
  fill = season
)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = RMSE_mean - RMSE_se,
      ymax = RMSE_mean + RMSE_se
    ),
    position = position_dodge(0.8),
    width = 0.2
  ) +
  scale_fill_manual(
    values = season_colors,
    labels = c("Wet", "Dry", "Wet-Dry")
  ) +
  scale_x_discrete(
    labels = c(
      "descriptive" = "Descriptive",
      "2phase-harmonic" = "Harmonics"
    )
  ) +
  labs(
    y = "Mean CV RMSE",
    fill = "Season",
    x = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom"
  )
# Show plot
#dev.new(); print(plot_pca_rmse)

# Combine plots with equal height
combined_pca_plot <- ggarrange(
  plot_pca_r2, plot_pca_rmse,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "right",
  heights = c(1, 1)  # equal vertical space
)

# Show and save plot
dev.new(); print(combined_pca_plot)

outFile <- file.path(inDir,"figures/AIrichness_pca_descriptive+harmonic_statistics.png")
ggsave(outFile, combined_pca_plot, width = 5, height = 5, units = "in", dpi = 600, bg = "white")

# Statistics
harm <- df_pca %>% filter(predictors == "2phase-harmonic")
harm <- harm %>% select(season,Rsq_mean)
summary(harm$Rsq_mean)

desc <- df_pca %>% filter(predictors == "descriptive")
desc <- desc %>% select(season,Rsq_mean)
summary(desc$Rsq_mean)

# Get raw embeddings model statistics for comparison
df_raw <- read_csv(file.path(inDir,"data/rf_model_statistics_260410.csv"))

df_raw <- df_raw %>%
  filter(
    richness_type == "AIrichness",
    AGIfilter == "None",
    !(predictors == "descriptive" & periods == "None"),
    !predictors == "hybrid",
    features == "BirdNET embeddings",
  )
desc_raw <- df_raw %>% filter(predictors == "descriptive")
desc <- desc %>% mutate(diff = Rsq_mean - desc_raw$Rsq_mean)
harm_raw <- df_raw %>% filter(predictors == "2phase-harmonic")
harm <- harm %>% mutate(diff = Rsq_mean - harm_raw$Rsq_mean)

# Test statistical differences
inFile <- file.path(inDir,"data/rf_model_statistics_folds_260410.csv")
df <- read.csv(inFile) %>% filter(richness_type == "AIrichness" &
                                    predictors != "hybrid" &
                                    features %in% c("PCA","BirdNET embeddings"))

# ANOVA of effect of features and season on R2
model <- lm(Rsquared ~ features * season, data = df)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ features)
emmeans(model, pairwise ~ season | features)

# ANOVA of effect of features and predictors on R2
model <- lm(Rsquared ~ features * predictors, data = df)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ features)
emmeans(model, pairwise ~ predictors | features)

# ========================================================================
# AI richness UMAP - descriptive + harmonic statistics
# ========================================================================

# Load and filter
df_umap <- read_csv(file.path(inDir,"data/rf_model_pca_umap_tpd_statistics_260410.csv")) %>%
  filter(
    richness_type == "AIrichness",
    features != "PCA",
    !is.na(n_neighbors),
    !is.na(min_dist)
  ) %>%
  mutate(
    # combine and order the param factor
    param = paste0(n_neighbors, "_", min_dist),
    param = fct_reorder(
      param,
      # reorder by numeric neighbors then numeric min_dist
      as.numeric(n_neighbors) * 1e3 + as.numeric(min_dist)
    ),
    season = factor(season, levels = c("wet", "dry", "drywet")),
    predictors = factor(predictors, levels = c("descriptive", "2phase-harmonic"))
  )

# Set up season colors
pal <- brewer.pal(3, "Set2")
season_colors <- c(
  "wet" = pal[3],
  "dry" = pal[2],
  "drywet" = pal[1]
)
param_levels <- levels(factor(df_umap$param))
x_pos <- which(param_levels == "20_0.9")

# Plot: consistent color by season, dashed lines for descriptive
plot_umap_r2 <- ggplot(df_umap, aes(
  x = param,
  y = Rsq_mean,
  color = season,
  linetype = predictors,
  shape = predictors, 
  group = interaction(season, predictors)
)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  scale_y_continuous(
    limits = c(0.46, 0.62),
    breaks = seq(0.46, 0.62, 0.02),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_color_manual(
    values = season_colors,
    labels = c("Wet", "Dry", "Wet-Dry")
  ) +
  scale_shape_manual(
    values = c("descriptive" = 16, "2phase-harmonic" = 1)  # solid circle, open circle
  ) +
  scale_linetype_manual(
    values = c("descriptive" = "solid", "2phase-harmonic" = "twodash"),
    labels = c("Descriptive", "Harmonics")
  ) +
  labs(
    x = "UMAP parameters (nearest neighbors_minimum distance)",
    y = expression("Mean CV " * R^2),
    color = "Season",
    linetype = "Statistics"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 6),
    legend.position = "top"
  ) +
  guides(shape = "none") +
  geom_vline(xintercept = x_pos, color = "black", linetype = "solid", linewidth = 0.5)


# Show and save plot
#dev.new(); print(p_umap_r2)

# Plot: consistent color by season, dashed lines for descriptive
plot_umap_rmse <- ggplot(df_umap, aes(
  x = param,
  y = RMSE_mean,
  color = season,
  linetype = predictors,
  shape = predictors, 
  group = interaction(season, predictors)
)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  scale_y_continuous(
    limits = c(3.99, 5.2),
    breaks = seq(4, 5.2, 0.2),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_color_manual(
    values = season_colors,
    labels = c("Wet", "Dry", "Wet-Dry")
  ) +
  scale_shape_manual(
    values = c("descriptive" = 16, "2phase-harmonic" = 1)  # solid circle, open circle
  ) +
  scale_linetype_manual(
    values = c("descriptive" = "solid", "2phase-harmonic" = "twodash"),
    labels = c("Descriptive", "Harmonics")
  ) +
  labs(
    x = "UMAP parameters (nearest neighbors_minimum distance)",
    y = expression("Mean RMSE"),
    color = "Season",
    linetype = "Statistics"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 6),
    legend.position = "top"
  ) +
  guides(shape = "none") +
  geom_vline(xintercept = x_pos, color = "black", linetype = "solid", linewidth = 0.5)


# Show and save plot
#dev.new(); print(p_umap_rmse)

# Combine plots with equal height
combined_umap_plot <- ggarrange(
  plot_umap_r2, plot_umap_rmse,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "right",
  heights = c(1, 1)  # equal vertical space
)

# Show and save plot
dev.new(); print(combined_umap_plot)

outFile <- file.path(inDir,"figures/AIrichness_umap_descriptive+harmonic_statistics.png")
ggsave(outFile, combined_umap_plot, width = 7, height = 7, units = "in", dpi = 600, bg = "white")

# statistics of best model
best_model_stats <- df_umap %>% filter(n_neighbors == 20 & min_dist == 0.9) %>% select(season, predictors, Rsq_mean, RMSE_mean) 

best <- best_model_stats %>% filter(predictors == "descriptive")
best <- best %>% mutate(diff = Rsq_mean - desc$Rsq_mean)

# Statistics
umap_harm <- df_umap %>%
  filter(predictors == "2phase-harmonic") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
umap_harm

umap_desc <- df_umap %>%
  filter(predictors == "descriptive") %>%
  summarise(
    Rsq_mean_min = min(Rsq_mean, na.rm = TRUE),
    Rsq_mean_max = max(Rsq_mean, na.rm = TRUE),
    Rsq_mean_avg = mean(Rsq_mean, na.rm = TRUE),
    
    RMSE_mean_min   = min(RMSE_mean, na.rm = TRUE),
    RMSE_mean_max   = max(RMSE_mean, na.rm = TRUE),
    RMSE_mean_avg   = mean(RMSE_mean, na.rm = TRUE)
  )
umap_desc

# Test statistical differences
inFile <- file.path(inDir,"data/rf_model_statistics_folds_260410.csv")
df <- read.csv(inFile) %>%
  filter(
    richness_type == "AIrichness",
    predictors != "hybrid",
    features %in% c("PCA", "UMAP"),
    file_name %in% c("birdnet_pca_minutes_nonzero_mean_dawn-dusk_descriptive_statistics",
                     "birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_nn20_mdist0.9_descriptive_statistics")
  )

# ANOVA of effect of features and season on R2
model <- lm(Rsquared ~ features * season, data = df)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ features)
emmeans(model, pairwise ~ season | features)

# ANOVA of effect of features and predictors on R2
model <- lm(Rsquared ~ features * predictors, data = df)
anova(model)

# paired comparisons
emmeans(model, pairwise ~ features)
emmeans(model, pairwise ~ predictors | features)

# ========================================================================
# AI richness UMAP scatterplot
# ========================================================================

# Read embeddings file
df <- read_csv(file.path(inDir,
  "data/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_nn20_mdist0.9.csv"),
  col_types = cols_only(
    SiteID    = col_character(),
    date      = col_character(),
    time_code = col_integer(),
    UMAP5     = col_double(),
    UMAP1     = col_double()
  )
)

# Read site attributes 
sites <- read_csv(paste0(inDir,"/data/biosoundscape_sites_daac_250507.csv")) %>%
  select(SiteID, Campaign, LandCoverClass, ElevationClass)

# Bin time of day
binTimeOfDay <- function(dd){
  hrb <- ifelse(dd >= 000 & dd < 500, "PredawnQuiet",  # Pre-dawn/Quiet
                ifelse(dd >= 500 & dd < 900, "DawnChorus",  # Dawn chorus
                       ifelse(dd >= 900 & dd < 1600, "Midday",  # Midday
                              ifelse(dd >= 1600 & dd < 2100, "DuskChorus", # Dusk chorus
                                     ifelse(dd >= 2100 & dd < 2400, "Nighttime", NA)))))  # Night
  return(hrb)
}

# Merge and filter
plot_df <- df %>%
  mutate(time_period = binTimeOfDay(time_code)) %>%
  filter(time_period %in% c("DawnChorus", "DuskChorus")) %>%
  left_join(sites, by = "SiteID")

# Load species richness 
richness <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) %>%
  rename(SiteID = siteid)

# Join species to plot data
umap_df <- plot_df %>%
  inner_join(richness, by = "SiteID")

# base aes
base_aes <- aes(x = UMAP1, y = UMAP5, z = wm_richness_ge3)

# Classify sites into low, med and high richness
site_classes <- umap_df %>% 
  distinct(SiteID, wm_richness_ge3) %>%
  mutate(
    RichnessClass = cut(
      wm_richness_ge3,
      breaks = quantile(wm_richness_ge3, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      labels = c("Low Richness", "Moderate Richness", "High Richness"),
      include.lowest = TRUE
    )
  ) %>% select(-wm_richness_ge3)

# Join back to every minute record
plot_df <- umap_df %>%
  left_join(site_classes, by = "SiteID") %>%
  filter(!is.na(RichnessClass), wm_richness_ge3 >= 0)

# global x and y limits based on the union of both datasets
xlims <- c(-10,10)
ylims <- c(-10,10)

# Plot: minute counts per hex, faceted by richness class
p <- ggplot(plot_df, base_aes) +
  stat_binhex(bins = 50, aes(fill = ..count..), color = NA) +
  stat_density_2d(
    aes(x = UMAP1, y = UMAP5, group = 1),  # or group = RichnessClass
    geom   = "contour",
    bins   = 5,
    colour = "red",
    linewidth   = 0.4
  ) +
  scale_fill_viridis_c(name = "Temporal\ndensity\n(# of min)", option = "D") +
  coord_cartesian(xlim = xlims, ylim = ylims) + 
  facet_wrap(~ RichnessClass, ncol = 3) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5, linetype = "22") +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5, linetype = "22") +
  labs(
    x = "UMAP1",
    y = "UMAP5"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text      = element_text(face = "bold"),
    legend.position = "right"
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/AIrichness_umap_scatterplot.png")
ggsave(outFile, p, width = 7, height = 3, units = "in", dpi = 600, bg = "white")

# Taxa UMAP
taxa <- read.csv(file.path(inDir,"data/wildmon_minute-level_species_taxon_abgi.csv"))
taxa_select <- taxa %>% filter(TimePeriod %in% c("DawnChorus", "DuskChorus")) %>% 
  select(siteid,Bird_percent, Frog_percent, Insect_percent,Biophony_percent, Geophony_percent,Anthropophony_percent, Interference_percent, year,month,day,time_code) %>% 
  rename(SiteID = siteid)
threshold <- 10
taxa_select <- taxa_select  %>%
  mutate(
    TaxaClass = case_when(
      Bird_percent   >= threshold & Frog_percent   < threshold & Insect_percent < threshold ~ "Birds",
      Frog_percent   >= threshold & Bird_percent   < threshold & Insect_percent < threshold ~ "Amphibians",
      Insect_percent >= threshold & Bird_percent   < threshold & Frog_percent   < threshold ~ "Insects",
      Bird_percent   <= 1 & Frog_percent  <= 1 & Insect_percent <= 1 ~ "Quiet",
      TRUE ~ "Mixed Taxa"
    )
  )
taxa_select <- taxa_select  %>%
  mutate(
    ABGIClass = case_when(
      Geophony_percent >= threshold & Anthropophony_percent < threshold & Interference_percent < threshold ~ "Geophony",
      Anthropophony_percent   >= threshold & Geophony_percent   < threshold & Interference_percent < threshold ~ "Anthropophony",
      Interference_percent >= threshold & Geophony_percent   < threshold & Anthropophony_percent   < threshold ~ "Interference",
      Geophony_percent   <= 1 & Anthropophony_percent  <= 1 & Interference_percent <= 1 & Biophony_percent <= 1 ~ "Quiet",
      TRUE ~ "Mixed BAGI"
    )
  )
taxa_select$date <- with(taxa_select, 
                         paste(
                           year,
                           sprintf("%02d", month),
                           sprintf("%02d", day),
                           sep = "-"
                         )
)
taxa_select$key <- paste0(taxa_select$SiteID,"_",taxa_select$date,"_",taxa_select$time_code)
taxa_select <- taxa_select %>% select(-SiteID,date,time_code)
plot_df$key <- paste0(plot_df$SiteID,"_",plot_df$date,"_",plot_df$time_code)
plot_taxa_agi <- left_join(plot_df,taxa_select,by="key") 

# Filter out unwanted classes
plot_df7 <- plot_taxa_agi %>%
  filter(!TaxaClass  %in% c("Quiet")     & !is.na(TaxaClass))

plot_df8 <- plot_taxa_agi %>%
  filter(!ABGIClass %in% c("Mixed BAGI") & !is.na(ABGIClass))

# Base plot A (by TaxaClass)
pA <- ggplot(plot_df7, aes(x = UMAP1, y = UMAP5)) +
  stat_binhex(bins = 50, aes(fill = ..count..), color = NA) +
  stat_density_2d(aes(group = 1),
                  geom   = "contour",
                  bins   = 5,
                  colour = "red",
                  size   = 0.4) +
  facet_wrap(~ TaxaClass, ncol = 2) +
  labs(x = "UMAP1", y = "UMAP5") +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5, linetype = "22") +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5, linetype = "22") +
  theme_minimal(base_size = 12) +
  theme(strip.text      = element_text(face = "bold"),
        legend.position = "bottom") 

# Base plot B (by ABGIClass)
pB <- ggplot(plot_df8, aes(x = UMAP1, y = UMAP5)) +
  stat_binhex(bins = 50, aes(fill = ..count..), color = NA) +
  stat_density_2d(aes(group = 1),
                  geom   = "contour",
                  bins   = 5,
                  colour = "red",
                  size   = 0.4) +
  facet_wrap(~ ABGIClass, ncol = 2) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5, linetype = "22") +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5, linetype = "22") +
  labs(x = "UMAP1", y = "UMAP5") +
  theme_minimal(base_size = 12) +
  theme(strip.text      = element_text(face = "bold"),
        legend.position = "bottom")

# Compute global scale limits
# 1) global maximum count for color scale
bldA      <- ggplot_build(pA)
bldB      <- ggplot_build(pB)
maxA      <- max(bldA$data[[1]]$count, na.rm = TRUE)
maxB      <- max(bldB$data[[1]]$count, na.rm = TRUE)
global_max <- max(maxA, maxB)

# 2) global x and y limits based on the union of both datasets
#xlims <- range(c(plot_df7$UMAP1, plot_df8$UMAP1), na.rm = TRUE)
#ylims <- range(c(plot_df7$UMAP5, plot_df8$UMAP5), na.rm = TRUE)
xlims <- c(-10,10)
ylims <- c(-10,10)
# xlims <- c(-5,5)
# ylims <- c(-5,5)

# Shared fill scale
shared_fill <- scale_fill_viridis_c(
  name    = "Temporal\ndensity\n(# of min)",
  option  = "D",
  limits  = c(0, global_max),
  oob     = scales::squish
)

# Apply shared scales, manual limits, and legend settings
pA2 <- pA +
  shared_fill +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  guides(fill = guide_colorbar(title.position = "top")) +
  theme(legend.position = "right")

pB2 <- pB +
  shared_fill +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  guides(fill = guide_colorbar(title.position = "top")) +
  theme(legend.position = "right")

# Combine plots with a single legend, side by side
combined <- ggarrange(
  pA2, pB2,
  labels        = c("A", "B"),
  ncol          = 2,
  common.legend = TRUE,
  legend        = "right",
  align         = "hv"
)

# Show and save plot
dev.new(); print(combined)

outFile <- file.path(inDir, "figures/taxa_agi_umap_scatterplot.png")
ggsave(outFile, combined, width = 7, height = 4, units = "in", dpi = 600, bg = "white")

# ========================================================================
# AI richness log(FRichness) vs wm_richness_ge3 scatterplot
# ========================================================================

# Read in the data
richness <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) %>% 
  rename(SiteID = siteid) %>% select(SiteID, wm_richness_ge3)
embeddings <- read.csv(file.path(inDir,"data/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_tpd_nn20_mdist0.9_UMAP5-UMAP1-UMAP3-UMAP2.csv"))
sites <- read_csv(file.path(inDir,"data/biosoundscape_sites_daac_250507.csv"))

# Merge 
merged <- embeddings %>%
  inner_join(richness, by = "SiteID") 
merged <- left_join(merged, sites, by = "SiteID")

# Drop any NA FRichness and add elevation class
plot_data <- merged %>%
  filter(!is.na(FRichness)) %>% 
  drop_na(ElevationClass) %>%
  mutate(
    ElevationClass = fct_relevel(ElevationClass,
                                 "Low: 0-500 m",
                                 "Medium: 500-1000 m",
                                 "High: >1000 m"
    ))

# Log scale model

# Fit the model on log10(FRichness)
model_log <- lm(wm_richness_ge3 ~ log10(FRichness), data = plot_data)
r2_log   <- summary(model_log)$r.squared

p <- ggplot(plot_data, aes(x = log10(FRichness), y = wm_richness_ge3)) +
  geom_point(alpha = 0.6, color = "black") +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("R² = ", round(r2_log, 2)),
           hjust = 1.1, vjust = 1.5,
           size = 5) +
  scale_x_continuous(limits = c(-0.5, 3)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(
    x = "Acoustic Trait-based Richness (Log scale)",
    y = "AI-based Richness"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 12)
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/AIrichness_tpd_log_richness_scatterplot.png")
ggsave(outFile, p, width = 4, height = 4, units = "in", dpi = 600, bg = "white")

# Log model RGB version

# Read the input CSVs
abgi <- read_csv(file.path(inDir, "data/wildmon_hour-level_species_taxon_abgi.csv"))

# Filter minutes to dawn and dusk choruses
abgi <- abgi %>% filter((hour >= 5 & hour < 9) | (hour >= 16 & hour < 21))

# Total ABGI
abgi_site <- abgi %>%
  group_by(siteid) %>%
  summarise(
    Anthropophony_DawnDusk = sum(Anthropophony, na.rm = TRUE),
    Biophony_DawnDusk      = sum(Biophony, na.rm = TRUE),
    Geophony_DawnDusk      = sum(Geophony, na.rm = TRUE),
    Interference_DawnDusk  = sum(Interference, na.rm = TRUE),
    Bird_DawnDusk          = sum(Bird, na.rm = TRUE),
    Frog_DawnDusk          = sum(Frog, na.rm = TRUE),
    Insect_DawnDusk        = sum(Insect, na.rm = TRUE),
    .groups = "drop"
  ) %>% rename(SiteID = siteid)

# Merge datasets
plot_data <- plot_data %>%
  inner_join(abgi_site, by = "SiteID")

# Normalize RGB channels so color reflects relative composition of BGI
total_sound <- with(
  plot_data,
  Interference_DawnDusk + Biophony_DawnDusk + Geophony_DawnDusk
)

plot_data$rgb_col <- rgb(
  red   = ifelse(total_sound > 0, plot_data$Interference_DawnDusk / total_sound, 0),
  green = ifelse(total_sound > 0, plot_data$Biophony_DawnDusk / total_sound, 0),
  blue  = ifelse(total_sound > 0, plot_data$Geophony_DawnDusk / total_sound, 0)
)

# Manual RGB legend
legend_df <- data.frame(
  x = c(-0.45, -0.45, -0.45),
  y = c(38, 36, 34),
  label = c("Interference", "Biophony", "Geophony"),
  col = c(
    rgb(1, 0, 0),
    rgb(0, 1, 0),
    rgb(0, 0, 1)
  )
)


# Plot with ABGI as RBG
plot_rgb_abgi <- ggplot(plot_data, aes(x = log10(FRichness), y = wm_richness_ge3)) +
  geom_point(aes(color = rgb_col), alpha = 0.6, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  annotate(
    "text",
    x = 2.75, y = 40,
    label = paste0("R² = ", round(r2_log, 2)),
    hjust = 1, vjust = 1,
    size = 5
  ) +
  geom_point(
    data = legend_df,
    aes(x = x, y = y),
    color = legend_df$col,
    size = 1,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = legend_df,
    aes(x = x + 0.08, y = y, label = label),
    hjust = 0,
    size = 4,
    inherit.aes = FALSE
  ) +
  scale_color_identity() +
  scale_x_continuous(limits = c(-0.5, 3)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(
    x = "Acoustic Trait-based Richness (Log scale)",
    y = "AI-based Richness"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12)
  )

# Calculate percent Biophony
plot_data <- plot_data %>%
  mutate(Biophony_DawnDusk_Percent = Biophony_DawnDusk/(Interference_DawnDusk + Biophony_DawnDusk + Geophony_DawnDusk + Anthropophony_DawnDusk) 
)

# Plot with Biophony percent
plot_biophony_percent <- ggplot(plot_data, aes(
  x = log10(FRichness),
  y = wm_richness_ge3,
  color = Biophony_DawnDusk_Percent
)) +
  geom_point(alpha = 0.8, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  annotate(
    "text",
    x = 2.75, y = 40,
    label = paste0("R² = ", round(r2_log, 2)),
    hjust = 1, vjust = 1,
    size = 5
  ) +
  scale_color_viridis_c(option = "magma", name = "Biophony (%)") +
  scale_x_continuous(limits = c(-0.5, 3)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(
    x = "Acoustic Trait-based Richness (Log scale)",
    y = "AI-based Richness"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank()
  )

#dev.new(); print(plot_biophony_percent)

# Plot with taxa as RBG
total_sound <- with(
  plot_data,
  Bird_DawnDusk + Frog_DawnDusk + Insect_DawnDusk
)

plot_data$rgb_col <- rgb(
  red   = ifelse(total_sound > 0, plot_data$Insect_DawnDusk / total_sound, 0),
  green = ifelse(total_sound > 0, plot_data$Bird_DawnDusk / total_sound, 0),
  blue  = ifelse(total_sound > 0, plot_data$Frog_DawnDusk / total_sound, 0)
)

n <- 200
legend_rgb <- expand.grid(
  x = seq(0, 1, length.out = n),
  y = seq(0, 1, length.out = n)
) %>%
  mutate(z = 1 - x - y) %>%
  filter(z >= 0) %>%
  mutate(
    R = x,
    G = y,
    B = z,
    col = rgb(R, G, B)
  )

rgb_triangle <- ggplot(legend_rgb, aes(x, y)) +
  geom_raster(aes(fill = col)) +
  scale_fill_identity() +
  coord_equal(
    xlim = c(-0.08, 1.12),
    ylim = c(-0.08, 1.08),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(2, 2, 2, 2),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text",
    x = 1.02, y = 0.00,
    label = "Insects",
    color = "red",
    hjust = 0, vjust = 0,
    size = 4
  ) +
  annotate(
    "text",
    x = 0.00, y = 1.02,
    label = "Birds",
    color = "green4",
    hjust = 0, vjust = 0,
    size = 4
  ) +
  annotate(
    "text",
    x = 0.00, y = -0.04,
    label = "Amphibians",
    color = "blue",
    hjust = 0, vjust = 1,
    size = 4
  )

plot_rgb_taxa <- ggplot(plot_data, aes(x = log10(FRichness), y = wm_richness_ge3)) +
  geom_point(aes(color = rgb_col), alpha = 0.6, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  annotate(
    "text",
    x = 2.75, y = 40,
    label = paste0("R² = ", round(r2_log, 2)),
    hjust = 1, vjust = 1,
    size = 5
  ) +
  scale_color_identity() +
  scale_x_continuous(limits = c(-0.5, 3)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(
    x = "Acoustic Trait-based Richness (Log scale)",
    y = "AI-based Richness"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12)
  )

plot_rgb_taxa_inset <- plot_rgb_taxa +
  inset_element(
    rgb_triangle,
    left   = 0.02,
    bottom = 0.64,
    right  = 0.28,
    top    = 0.98,
    align_to = "panel",
    ignore_tag = TRUE
  )

#dev.new(); print(plot_rgb_taxa_inset)

# Plot with elevation as colors
plot_elev <- ggplot(plot_data, aes(x = log10(FRichness), y = wm_richness_ge3)) +
  geom_point(aes(color = ElevationClass), alpha = 0.6, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  annotate(
    "text",
    x = 2.75, y = 40,
    label = paste0("R² = ", round(r2_log, 2)),
    hjust = 1, vjust = 1,
    size = 5
  ) +
  scale_color_manual(
    values = c(
      "Low: 0-500 m" = "#4DAF4A",
      "Medium: 500-1000 m" = "#E41A1C",
      "High: >1000 m" = "#377EB8"
    ),
    name = "Elevation"
  ) +
  scale_x_continuous(limits = c(-0.5, 3)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(
    x = "Acoustic Trait-based Richness (Log scale)",
    y = "AI-based Richness"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    strip.text = element_text(size = 12)
  )

#dev.new(); print(plot_elev)

# Plot with UMAP5 statistics as color ramp
umap <- read_csv(file.path(inDir,"data/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_nn20_mdist0.9_descriptive_statistics.csv"))
plot_umap_data <- left_join(plot_data,umap,by="SiteID")
plot_umap <- ggplot(plot_umap_data, aes(
  x = log10(FRichness),
  y = wm_richness_ge3,
  color = UMAP5_q25
)) +
  geom_point(alpha = 0.6, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  annotate(
    "text",
    x = 2.75, y = 40,
    label = paste0("R² = ", round(r2_log, 2)),
    hjust = 1, vjust = 1,
    size = 5
  ) +
  scale_color_viridis_c(option = "magma", name = "UMAP5 Std Dev") +
  scale_x_continuous(limits = c(-0.5, 3)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(
    x = "Acoustic Trait-based Richness (Log scale)",
    y = "AI-based Richness"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank()
  )

#dev.new(); print(plot_umap)

# UMAP5 RGB version
plot_umap_data_rgb <- plot_umap_data %>%
  mutate(
    x_log = log10(FRichness),
    R = percent_rank(UMAP5_median)^0.8,
    G = percent_rank(UMAP5_iqr)^0.8,
    B = percent_rank(UMAP5_skewness)^0.8,
    rgb_col = rgb(R, G, B)
  )

n <- 200
legend_rgb <- expand.grid(
  x = seq(0, 1, length.out = n),
  y = seq(0, 1, length.out = n)
) %>%
  mutate(z = 1 - x - y) %>%
  filter(z >= 0) %>%
  mutate(
    R = x,
    G = y,
    B = z,
    col = rgb(R, G, B)
  )

rgb_triangle <- ggplot(legend_rgb, aes(x, y)) +
  geom_raster(aes(fill = col)) +
  scale_fill_identity() +
  coord_equal(
    xlim = c(-0.08, 1.12),
    ylim = c(-0.08, 1.08),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(2, 2, 2, 2),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text",
    x = 1.02, y = 0.00,
    label = "Median",
    color = "red",
    hjust = 0, vjust = 0,
    size = 4
  ) +
  annotate(
    "text",
    x = 0.00, y = 1.02,
    label = "IQR",
    color = "green4",
    hjust = 0, vjust = 0,
    size = 4
  ) +
  annotate(
    "text",
    x = 0.00, y = -0.04,
    label = "Skewness",
    color = "blue",
    hjust = 0, vjust = 1,
    size = 4
  )

plot_umap_rgb <- ggplot(plot_umap_data_rgb, aes(x = x_log, y = wm_richness_ge3)) +
  geom_point(aes(color = rgb_col), alpha = 0.8, size = 0.5) +
  scale_color_identity() +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  annotate(
    "text",
    x = 2.75, y = 40,
    label = paste0("R² = ", round(r2_log, 2)),
    hjust = 1, vjust = 1,
    size = 5
  ) +
  scale_x_continuous(limits = c(-0.5, 3)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(
    x = "Acoustic Trait-based Richness (Log scale)",
    y = "AI-based Richness"
  ) +
  theme_minimal(base_size = 12)

plot_umap_rgb_inset <- plot_umap_rgb +
  inset_element(
    rgb_triangle,
    left   = 0.02,
    bottom = 0.64,
    right  = 0.28,
    top    = 0.98,
    align_to = "panel",
    ignore_tag = TRUE
  )

combined <- plot_umap_rgb_inset + plot_rgb_taxa_inset +#plot_biophony_percent +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")

# Show and save plot
dev.new(); print(combined)

outFile <- file.path(inDir, "figures/umap_scatterplot.png")
ggsave(outFile, combined, width = 7.5, height = 4, units = "in", dpi = 600, bg = "white")

# ========================================================================
# AI richness vs TPD faceted scatterplot -- OLS, not scaled
# ========================================================================

# Read in the data
richness <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) %>%
  rename(SiteID = siteid)

embeddings <- read.csv(file.path(inDir,"data/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_tpd_nn20_mdist0.9_UMAP5-UMAP1-UMAP3-UMAP2.csv"))
sites <- read_csv(file.path(inDir,"data/biosoundscape_sites_daac_250507.csv"))

# Merge 
merged <- embeddings %>%
  inner_join(richness, by = "SiteID") 

merged <- left_join(merged, sites, by = "SiteID")

# Pivot the three functional metrics into long form
plot_data <- merged %>%
  dplyr::select(SiteID, wm_richness_ge3, FRichness, FEvenness, FDivergence, ElevationClass) %>%
  pivot_longer(
    cols = c(FRichness, FEvenness, FDivergence),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  mutate(Metric = factor(Metric,
                       levels = c("FRichness", "FEvenness", "FDivergence"),
                       labels = c("FRic",
                                  "FEve",
                                  "FDiv")
  ))

# Fit one linear model per facet (for annotation)
models <- plot_data %>%
  group_by(Metric) %>%
  do(model = lm(wm_richness_ge3 ~ Value, data = .)) %>%
  summarise(
    Metric,
    intercept = coef(model)[1],
    slope     = coef(model)[2],
    r2        = summary(model)$r.squared,
    pval      = summary(model)$coefficients[2, "Pr(>|t|)"],
  )
models <- models %>%
  mutate(
    stars = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ ""
    ),
    label = paste0("R^2 == ", round(r2, 2), " * '", stars, "'")
  )

# Join back the model statitics for annotation
plot_data <- plot_data %>%
  left_join(models, by = "Metric") %>%
  drop_na(ElevationClass) %>%
  mutate(
    ElevationClass = fct_relevel(ElevationClass,
                                 "Low: 0-500 m",
                                 "Medium: 500-1000 m",
                                 "High: >1000 m"
    )
  )

# Make the faceted plot
cols <- brewer.pal(8, "Set1")[c(5, 2, 1)]
p <- ggplot(plot_data, aes(x = Value, y = wm_richness_ge3)) +
  geom_point(aes(color = ElevationClass), alpha = 0.6, size = 0.5) +
  geom_abline(
    data = models,
    aes(intercept = intercept, slope = slope),
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Metric, scales = "free_x", ncol = 3) +
  geom_text(
    data = models,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1,
    vjust = 1.5,
    size = 4,
    parse = TRUE,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(expand = expansion(add = 0.1)) +
  scale_color_manual(
    values = c(
      "Low: 0-500 m" = cols[1],
      "Medium: 500-1000 m" = cols[2],
      "High: >1000 m" = cols[3]
    ),
    name = "Elevation"
  ) +
  labs(
    x = NULL,
    y = "AI-based Richness",
    color = "Elevation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 13),
    panel.spacing = unit(1, "lines")
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/AIrichness_tpd_scatterplot.png")
ggsave(outFile, p, width = 7, height = 4, units = "in", dpi = 600, bg = "white")

# ========================================================================
# AI richness vs TPD faceted scatterplot - GAM regression
# ========================================================================

# Read in the data
richness <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) %>%
  rename(SiteID = siteid)
embeddings <- read.csv(file.path(inDir,"data/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_tpd_nn20_mdist0.9_UMAP5-UMAP1-UMAP3-UMAP2.csv"))
sites <- read_csv(file.path(inDir,"data/biosoundscape_sites_daac_250507.csv"))
veg <- read_csv(file.path(inDir,"data/landcover_vegclass3_summary.csv"))

# Merge 
merged <- embeddings %>%
  inner_join(richness, by = "SiteID")

df <- left_join(merged, sites, by = "SiteID")
df <- df %>%
  mutate(LogRecordingNum = log(RecordingNum)) %>%
  drop_na(wm_richness_ge3, FRichness, FEvenness, FDivergence)

# GAM modeling with smooth spatial structure
gam_model_full <- gam(
  wm_richness_ge3 ~ log(FRichness) + FEvenness + FDivergence + 
    s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)
summary(gam_model_full)

gam_model_fric <- gam(
  wm_richness_ge3 ~ log(FRichness) + 
    s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)
summary(gam_model_fric)

gam_model_feve <- gam(
  wm_richness_ge3 ~ FEvenness + 
  s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)
summary(gam_model_feve)

gam_model_fdiv <- gam(
  wm_richness_ge3 ~ FDivergence + 
    s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)
summary(gam_model_fdiv)

gam_model_null <- gam(
  wm_richness_ge3 ~ 1 + s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)
summary(gam_model_null)


anova(gam_model_null, gam_model_full, test = "Chisq")
anova(gam_model_null, gam_model_fric, test = "Chisq")
anova(gam_model_null, gam_model_feve, test = "Chisq")
anova(gam_model_null, gam_model_fdiv, test = "Chisq")

# Predict from GAM models
df <- df %>%
  mutate(
    GAM_Full = predict(gam_model_full, type = "response"),
    GAM_FRic = predict(gam_model_fric, type = "response"),
    GAM_FEve = predict(gam_model_feve, type = "response"),
    GAM_FDiv = predict(gam_model_fdiv, type = "response")
  ) %>% left_join(veg, by = "SiteID")

# Convert to long format
df_long <- df %>%
  dplyr::select(wm_richness_ge3, GAM_Full, GAM_FRic, GAM_FEve, GAM_FDiv, MajorityClass3, HighVeg, LowVeg, Other, ElevationClass) %>%
  pivot_longer(cols = starts_with("GAM_"), names_to = "Model", values_to = "Predicted") %>%
  mutate(
    Model = recode(Model,
                   GAM_Full = "Full",
                   GAM_FRic = "Fric",
                   GAM_FEve = "Feve",
                   GAM_FDiv = "Fdiv"),
    Model = fct_drop(as.factor(Model))
  )

# Compute Nagelkerke R²
R2_nagelkerke <- function(model, null_model) {
  loglik_model <- logLik(model)
  loglik_null <- logLik(null_model)
  n <- nobs(model)
  
  R2 <- (1 - exp((2 / n) * (loglik_null - loglik_model))) /
    (1 - exp((2 / n) * loglik_null))
  return(as.numeric(R2))
}

# Calculate R² values per model
r2_values <- tibble(
  Model = c("Full", "Fric", "Feve", "Fdiv"),
  R2 = c(
    R2_nagelkerke(gam_model_full,gam_model_null),
    R2_nagelkerke(gam_model_fric,gam_model_null),
    R2_nagelkerke(gam_model_feve,gam_model_null),
    R2_nagelkerke(gam_model_fdiv,gam_model_null)
  )
) %>%
  mutate(label = paste0("Pseudo R² = ", round(R2, 2)))

df_long <- df_long %>%
  left_join(r2_values, by = "Model") %>%
  drop_na(HighVeg, LowVeg, Other) %>%
  mutate(
    r = LowVeg/100,
    g = HighVeg/100,
    b = Other/100,
    rgb_color = rgb(r, g, b)
  )

df_long <- df_long %>%
  drop_na(ElevationClass) %>%
  mutate(
    ElevationClass = fct_relevel(ElevationClass,
                                 "Low: 0-500 m",
                                 "Medium: 500-1000 m",
                                 "High: >1000 m"
    )
  )

# Plot with facets
cols <- brewer.pal(8, "Set1")[c(5, 2, 1)]
p <- ggplot(df_long, aes(x = wm_richness_ge3, y = Predicted)) +
  geom_point(aes(color = ElevationClass), alpha = 0.6, size = 0.5) +
  #geom_point(aes(color = MajorityClass3), alpha = 0.5, size = 0.5) +
  #geom_point(aes(color = rgb_color), alpha = 0.8, size = 0.7, show.legend = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "red", linewidth = 0.5) +
  facet_wrap(~Model) +
  geom_text(
    data = r2_values,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1,
    vjust = 2.0,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Low: 0-500 m" = cols[1],
      "Medium: 500-1000 m" = cols[2],
      "High: >1000 m" = cols[3]
    ),
    name = "Elevation"
  ) +
  labs(
    x = "Observed Richness",
    y = "Predicted Richness",
    color = "Elevation"
  ) +
  theme_minimal(base_size = 12) + 
  theme(
    legend.position = "bottom",
    strip.text       = element_text(size = 13),
    panel.spacing    = unit(1, "lines")
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/AIrichness_tpd_gam_regression_scatterplot.png")
ggsave(outFile, p, width = 6, height = 4, units = "in", dpi = 600, bg = "white")

# ========================================================================
# AI richness wm_richness_ge3 vs RF UMAP important variables scatterplot
# ========================================================================

# Read in the data
richness <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) %>%
  rename(SiteID = siteid)
embeddings <- read.csv(file.path(umapDir,"BirdNET_pca_umap_tpd_stats/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_tpd_nn20_mdist0.9_UMAP5-UMAP1-UMAP3-UMAP2.csv"))

# Merge 
merged <- embeddings %>%
  inner_join(richness, by = "SiteID")

load(file.path(inDir,"models/rf_birdnet/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_nn20_mdist0.9_descriptive_statistics_dry-wet_season.RData"))
imp <- drywet_models$AIrichness$final_model$variable.importance
s <- order(imp, decreasing = T)
imp <- imp[s]
vars <- names(imp[1:6])
vars <- vars[vars != "LogRecNum"]

# Pivot the three functional metrics into long form
plot_data <- merged %>%
  select(SiteID, wm_richness_ge3, vars) %>%
  pivot_longer(
    cols = vars,
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# Fit one linear model per facet (for annotation)
models <- plot_data %>%
  group_by(Metric) %>%
  do(
    model = lm(wm_richness_ge3 ~ Value, data = .)
  ) %>%
  summarise(
    Metric,
    intercept = coef(model)[1],
    slope     = coef(model)[2],
    r2        = summary(model)$r.squared
  )

# Join back the R2 values for annotation
plot_data <- plot_data %>%
  left_join(models, by = "Metric")

# Make the faceted plot
p <- ggplot(plot_data, aes(x = Value, y = wm_richness_ge3)) +
  geom_point(alpha = 0.5) +
  geom_abline(aes(intercept = intercept, slope = slope), data = models, color = "red", size = 0.8) +
  facet_wrap(~ Metric, scales = "free_x", ncol = 1) +
  geom_text(
    data = models,
    aes(
      x = Inf, y = Inf,
      label = paste0("R² = ", round(r2, 2))
    ),
    hjust = 1.1, vjust = 1.5, size = 4
  ) +
  scale_x_continuous(expand = expansion(add = 0.1)) +
  labs(
    x = NULL,
    y = "AI-based Richness",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text       = element_text(size = 13),
    panel.spacing    = unit(1, "lines")
  )

# Show plot
dev.new(); print(p)

# ========================================================================
# Histogram of point count sampling time
# ========================================================================

# Load data
data <- read_csv(file.path(inDir,"data/2024-03-28 All BioSCape Point Counts Spatial Attr.csv"))

# Convert Time to POSIXct and extract hour
data <- data %>%
  mutate(Time = as.character(Time),
         Time = parse_time(Time, format = "%H:%M:%S"),
         Hour = hour(Time)) %>%
  filter(!is.na(Hour))

# Get minimum hour per site
site_min_hour <- data %>%
  group_by(Location_ID) %>%
  summarize(MinHour = min(Hour, na.rm = TRUE)) %>%
  ungroup()

# Compute summary statistics on MinHour
mean_min_hour <- mean(site_min_hour$MinHour)
sd_min_hour <- sd(site_min_hour$MinHour)

# Define chorus hours
dawn_hours <- 5:9
dusk_hours <- 16:20
chorus_hours <- c(dawn_hours, dusk_hours)

# Percent of sites outside chorus hours
non_chorus_pct <- site_min_hour %>%
  filter(!MinHour %in% chorus_hours) %>%
  nrow() / nrow(site_min_hour) * 100

# Print summary
cat(sprintf("Mean of Minimum Hour per Site: %.2f\nSD: %.2f\nPercent Outside Chorus: %.1f%%\n",
            mean_min_hour, sd_min_hour, non_chorus_pct))

# Create dummy data for shading
shading_df <- tibble(
  xmin = c(5, 16),
  xmax = c(9, 21),
  ymin = 0,
  ymax = Inf,
  Chorus = c("Dawn Chorus (5–9 AM)", "Dusk Chorus (4–9 PM)")
)

# Plot
p <- ggplot(site_min_hour, aes(x = MinHour)) +
  # Shaded chorus periods
  geom_rect(data = shading_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Chorus),
            inherit.aes = FALSE, alpha = 0.3) +
  # Histogram
  geom_histogram(binwidth = 1, boundary = 0, fill = "gray70", color = "black") +
  # Mean of minimum hour per site
  geom_vline(aes(xintercept = mean_min_hour, color = "Mean Start Time"), linetype = "dashed", linewidth = 1) +
  scale_fill_manual(name = "Chorus Period", values = c("Dawn Chorus (5–9 AM)" = "orange", "Dusk Chorus (4–9 PM)" = "skyblue")) +
  scale_color_manual(name = "", values = c("Mean Start Time" = "red")) +
  labs(
       x = "Start hour per Site",
       y = "Number of Sites") +
  theme_minimal() +
  theme(legend.position = "right")

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/point_count_sampling_time_histogram.png")
ggsave(outFile, p, width = 6, height = 3, units = "in", dpi = 600, bg = "white")

# ========================================================================
# AI richness BirdNET importance of hybrid model
# ========================================================================
 
# Load model
load(file.path(inDir,"models/rf_birdnet/birdnet_embeddings_minutes_nonzero_mean_dawn-dusk_hybrid_dry-wet_season.RData"))

imp <- drywet_models$AIrichness$final_model$variable.importance
importance <- data.frame(variable = names(imp), importance = imp)

#write.csv(importance,"c:/temp/imp.csv", row.names = F)

# Sort by importance descending
df <- importance[order(-importance$importance),]

# Select top 20
df20 <- df[1:20,]

# Create a bar plot with gray bars, high to low (top to bottom)
p <- ggplot(df20, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "gray") +
  coord_flip() +
  labs(x = "Variable", y = "Importance") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10)
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/rf_variable_importance_birdnet_wet-dry_hybrid_statistics.png")
ggsave(outFile, p, width = 4, height = 4, units = "in", dpi = 600, bg = "white")

# Separate the 'variable' column into 'feature' and 'statistic'
df_sep <- df %>%
  separate(variable, into = c("feature", "statistic"), sep = "_", remove = FALSE)

# Count number of each feature type
feature_counts <- df_sep %>%
  count(feature, sort = TRUE)

# Count number of each statistic type
statistic_counts <- df_sep %>%
  count(statistic, sort = TRUE)

# View the results
print(feature_counts)
print(statistic_counts)

# ========================================================================
# AI richness VGGish importance of hybrid model
# ========================================================================

# Load model
load(file.path(inDir,"models/rf_vggish/vggish_embeddings_minutes_nonzero_mean_dawn-dusk_hybrid_dry-wet_season.RData"))

imp <- drywet_models$AIrichness$final_model$variable.importance
importance <- data.frame(variable = names(imp), importance = imp)

# Sort by importance descending
df <- importance[order(-importance$importance),]

# Select top 20
df20 <- df[1:20,]

# Create a bar plot with gray bars, high to low (top to bottom)
p <- ggplot(df20, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "gray") +
  coord_flip() +
  labs(x = "Variable", y = "Importance") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10)
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/rf_variable_importance_vggish_wet-dry_hybrid_statistics.png")
ggsave(outFile, p, width = 4, height = 4, units = "in", dpi = 600, bg = "white")

# Separate the 'variable' column into 'feature' and 'statistic'
df_sep <- df %>%
  separate(variable, into = c("feature", "statistic"), sep = "_", remove = FALSE)

# Count number of each feature type
feature_counts <- df_sep %>%
  count(feature, sort = TRUE)

# Count number of each statistic type
statistic_counts <- df_sep %>%
  count(statistic, sort = TRUE)

# View the results
print(feature_counts)
print(statistic_counts)

# ========================================================================
# AI richness AVES importance of hybrid model
# ========================================================================

# Load model
load(file.path(inDir,"models/rf_aves/aves_embeddings_minutes_nonzero_mean_dawn-dusk_hybrid_dry-wet_season.RData"))

imp <- drywet_models$AIrichness$final_model$variable.importance
importance <- data.frame(variable = names(imp), importance = imp)

# Sort by importance descending
df <- importance[order(-importance$importance),]

# Select top 20
df20 <- df[1:20,]

# Create a bar plot with gray bars, high to low (top to bottom)
p <- ggplot(df20, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "gray") +
  coord_flip() +
  labs(x = "Variable", y = "Importance") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10)
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/rf_variable_importance_aves_wet-dry_hybrid_statistics.png")
ggsave(outFile, p, width = 4, height = 4, units = "in", dpi = 600, bg = "white")

# Separate the 'variable' column into 'feature' and 'statistic'
df_sep <- df %>%
  separate(variable, into = c("feature", "statistic"), sep = "_", remove = FALSE)

# Count number of each feature type
feature_counts <- df_sep %>%
  count(feature, sort = TRUE)

# Count number of each statistic type
statistic_counts <- df_sep %>%
  count(statistic, sort = TRUE)

# View the results
print(feature_counts)
print(statistic_counts)

# ========================================================================
# AI richness Acoustic Indices importance of hybrid model
# ========================================================================

# Load model
load(file.path(inDir,"models/rf_acoustic_indices/acoustic_indices_minutes_dawn-dusk_hybrid_dry-wet_season.RData"))

imp <- drywet_models$AIrichness$final_model$variable.importance
importance <- data.frame(variable = names(imp), importance = imp)
importance$variable <- gsub("NDSI_A", "NDSI-A", importance$variable)
importance$variable <- gsub("NDSI_B", "NDSI-B", importance$variable)

#write.csv(importance,"c:/temp/imp.csv", row.names = F)

# Sort by importance descending
df <- importance[order(-importance$importance),]

# Select top 20
df20 <- df[1:20,]

# Create a bar plot with gray bars, high to low (top to bottom)
p <- ggplot(df20, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "gray") +
  coord_flip() +
  labs(x = "Variable", y = "Importance") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10)
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/rf_variable_importance_acoustic_indices_wet-dry_descriptive_statistics.png")
ggsave(outFile, p, width = 4, height = 4, units = "in", dpi = 600, bg = "white")

# Separate the 'variable' column into 'feature' and 'statistic'
df_sep <- df20 %>%
  separate(variable, into = c("feature", "statistic"), sep = "_", remove = FALSE)

# Count number of each feature type
feature_counts <- df_sep %>%
  count(feature, sort = TRUE)

# Count number of each statistic type
statistic_counts <- df_sep %>%
  count(statistic, sort = TRUE)

# View the results
print(feature_counts)
print(statistic_counts)

# ========================================================================
# UMAP Random Forest model importance
# ========================================================================

load(file.path(inDir,"models/rf_umap/birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_nn20_mdist0.9_descriptive_statistics_dry-wet_season.RData"))

imp <- drywet_models$AIrichness$final_model$variable.importance
importance <- data.frame(variable = names(imp), importance = imp)

#write.csv(importance,"c:/temp/imp.csv", row.names = F)

# Sort by importance descending
df <- importance[order(-importance$importance),]

# Create a bar plot with gray bars, high to low (top to bottom)
p <- ggplot(df, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "gray") +
  coord_flip() +
  labs(x = "Variable", y = "Importance") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10)
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/rf_variable_importance_umap_wet-dry_descriptive_statistics.png")
ggsave(outFile, p, width = 4, height = 4, units = "in", dpi = 600, bg = "white")

# ========================================================================
# Scree plot of PCA on embeddings
# ========================================================================

# Elbow detection (max distance to line method)
find_elbow <- function(x, y) {
  line_start <- c(x[1], y[1])
  line_end <- c(x[length(x)], y[length(y)])
  line_vec <- line_end - line_start
  line_vec_norm <- line_vec / sqrt(sum(line_vec^2))
  
  # Project each point onto the line and compute distance
  distances <- sapply(1:length(x), function(i) {
    point <- c(x[i], y[i])
    vec <- point - line_start
    proj_len <- sum(vec * line_vec_norm)
    proj_point <- line_start + proj_len * line_vec_norm
    dist <- sqrt(sum((point - proj_point)^2))
    return(dist)
  })
  
  which.max(distances)
}

# Input the variance explained by PCs
load(file.path(umapDir,"birdnet_embeddings_minute_nonzero_mean_dry-wet_pca_20pcs.RData"))
pc_sum <- t(summary(pca_res)$importance)

# Create data frame
df1 <- tibble(PC = seq(1:dim(pc_sum)[1]), Proportion = pc_sum[,2], Cumulative = pc_sum[,3])

elbow_index1 <- find_elbow(df1$PC, df1$Proportion)
elbow_pc1 <- df1$PC[elbow_index1]
remove(pc_sum)

# Add red point using a separate data frame with a single point
elbow_point1 <- data.frame(PC = elbow_pc1, Proportion = df1$Proportion[elbow_index1])

# Plot
pA3 <- ggplot(df1, aes(x = PC)) +
  geom_line(aes(y = Proportion), color = "black") +
  geom_point(aes(y = Proportion), color = "black") +
  geom_point(aes(y = Cumulative), color = "blue") +
  geom_line(aes(y = Cumulative), color = "blue") +
  geom_point(data = elbow_point1, aes(x = PC, y = Proportion), color = "red", size = 3) +
  # annotate("text", x = elbow_pc, y = df$Proportion[elbow_index] + 0.005,
  #          label = paste("Elbow: PC", elbow_pc), color = "red", hjust = 0) +
  labs(
    x = "Principal Component", y = "Variance Explained") +
  theme_minimal()

load(file.path(umapDir,"birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_20pcs.RData"))

pc_sum <- t(summary(pca_res)$importance)

# Create data frame
df2 <- tibble(PC = seq(1:dim(pc_sum)[1]), Proportion = pc_sum[,2], Cumulative = pc_sum[,3])

elbow_index2 <- find_elbow(df2$PC, df2$Proportion)
elbow_pc2 <- df2$PC[elbow_index2]

# Add red point using a separate data frame with a single point
elbow_point2 <- data.frame(PC = elbow_pc2, Proportion = df2$Proportion[elbow_index2])

# Plot
pB3 <- ggplot(df2, aes(x = PC)) +
  geom_line(aes(y = Proportion), color = "black") +
  geom_point(aes(y = Proportion), color = "black") +
  geom_point(aes(y = Cumulative), color = "blue") +
  geom_line(aes(y = Cumulative), color = "blue") +
  geom_point(data = elbow_point2, aes(x = PC, y = Proportion), color = "red", size = 3) +
  # annotate("text", x = elbow_pc, y = df$Proportion[elbow_index] + 0.005,
  #          label = paste("Elbow: PC", elbow_pc), color = "red", hjust = 0) +
  labs(
    x = "Principal Component", y = "Variance Explained") +
  theme_minimal()

# Combine plots with a single legend, side by side
combined <- ggarrange(
  pA3, pB3,
  labels        = c("A", "B"),
  ncol          = 2,
  common.legend = TRUE,
  legend        = FALSE,
  align         = "hv"
)

# Show and save plot
dev.new(); print(combined)

outFile <- file.path(inDir, "figures/birdnet_embeddings_pca_20pcs.png")

ggsave(outFile, combined, width = 7, height = 4, units = "in", dpi = 600, bg = "white")

# ========================================================================
# Box plot of richness by elevation
# ========================================================================

# Read in the data
species_pc <- read.csv(file.path(inDir,"data/bioscape_point_count_richness_v20250602.csv"))
species_pc <- species_pc %>% rename(PCrichness = PC_Richness) %>%
  select(SiteID,PCrichness, Campaign)

species_ai <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) 
species_ai <- species_ai %>% rename(AIrichness = wm_richness_ge3, SiteID = siteid) %>%
  select(SiteID,AIrichness)

richness <- left_join(species_ai,species_pc,by="SiteID")

sites <- read_csv(file.path(inDir,"data/biosoundscape_sites_daac_250507.csv"))

# Merge on SiteID and clean
merged <- richness %>%
  left_join(sites %>% select(SiteID, ElevationClass), by = "SiteID") %>%
  filter(!is.na(ElevationClass)) %>%
  mutate(Season = case_when(
    grepl("wet", Campaign, ignore.case = TRUE) ~ "Wet",
    grepl("dry", Campaign, ignore.case = TRUE) ~ "Dry",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Season))

# Reshape to long format for ggplot
long_data <- merged %>%
  select(SiteID, ElevationClass, Season, PCrichness, AIrichness) %>%
  pivot_longer(cols = c(PCrichness, AIrichness),
               names_to = "RichnessType",
               values_to = "Richness") %>%
  filter(!is.na(Richness)) %>%
  mutate(
    RichnessType = recode(RichnessType,
                          "PCrichness" = "Point Count",
                          "AIrichness" = "AI-based"),
    ElevationClass = factor(ElevationClass, 
                            levels = c("Low: 0-500 m", "Medium: 500-1000 m", "High: >1000 m"),
                            labels = c("Low", "Medium", "High")),
    Season = factor(Season, 
                    levels = c("Wet","Dry"),
                    labels = c("Wet Season", "Dry Season"))
  )

pal <- brewer.pal(6, "Set2")

# Updated plot
p <- ggplot(long_data, aes(x = ElevationClass, y = Richness, fill = RichnessType)) +
  geom_boxplot(outlier.color = "black", position = position_dodge(width = 0.75)) +
  facet_wrap(~Season) +
  scale_fill_manual(values = c("Point Count" = pal[1], "AI-based" = pal[2])) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Elevation Class",
    y = "Richness",
    fill = "Richness Type"
  )

# Show and save plot
dev.new(); print(p)

outFile <- file.path(inDir, "figures/species_richness_elevation.png")
ggsave(outFile, p, width = 6, height = 3, units = "in", dpi = 600, bg = "white")

# ========================================================================
# Ternary plot of bird, frog and insect percent by site richness 
# ========================================================================

# Read the input CSVs
abgi <- read_csv(file.path(inDir, "data/wildmon_hour-level_species_taxon_abgi.csv")) 
richness <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) %>% select(siteid, wm_richness_ge3)

# Merge and create richness class
merged <- abgi %>%
  inner_join(richness, by = "siteid") %>%
  mutate(RichnessClass = cut(
    wm_richness_ge3,
    breaks = quantile(wm_richness_ge3, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("Low Richness", "Moderate Richness", "High Richness"),
    include.lowest = TRUE
  ))

# Average percentages and normalize
site_percent <- merged %>%
  group_by(siteid, RichnessClass) %>%
  summarise(
    Bird = sum(Bird, na.rm = TRUE),
    Frog = sum(Frog, na.rm = TRUE),
    Insect = sum(Insect, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(total = Bird + Frog + Insect) %>%
  filter(!is.na(total) & total > 0) %>%
  mutate(
    Bird = Bird / total,
    Frog = Frog / total,
    Insect = Insect / total
  )

# Define color palette
class_colors <- setNames(brewer.pal(3, "Set2"),
                         c("Low Richness", "Moderate Richness", "High Richness"))

# Prepare matrix for ternary plot
coords <- as.matrix(site_percent[, c("Bird", "Frog", "Insect")])

# Begin ternary plot
TernaryPlot(alab = "Bird", blab = "Frog", clab = "Insect", grid.lines = 5)

# Plot points by RichnessClass
for (class in levels(site_percent$RichnessClass)) {
  points <- coords[site_percent$RichnessClass == class, ]
  alpha_color <- adjustcolor(class_colors[class], alpha.f = 0.5)
  TernaryPoints(points, pch = 19, col = alpha_color, cex = 0.7)
}

legend(
  "topright",                           # Position
  legend = levels(site_percent$RichnessClass),  # Labels
  pt.bg = adjustcolor(class_colors, alpha.f = 0.5),  # Fill color with transparency
  pch = 21,                             # Match point shape
  pt.cex = 1.2,                         # Point size in legend
  bty = "n",                            # No box around legend
  title = "Richness Class"
)


# Merge and summarize
merged <- abgi %>%
  inner_join(richness, by = "siteid") %>%
  group_by(siteid, wm_richness_ge3) %>%
  summarise(
    Bird = mean(Bird_percent, na.rm = TRUE),
    Frog = mean(Frog_percent, na.rm = TRUE),
    Insect = mean(Insect_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(total = Bird + Frog + Insect) %>%
  filter(!is.na(total) & total > 0) %>%
  mutate(
    Bird = Bird / total,
    Frog = Frog / total,
    Insect = Insect / total
  )

# Generate color ramp based on wm_richness_ge3
colors <- viridis(length(merged$wm_richness_ge3))
ranked <- rank(merged$wm_richness_ge3, ties.method = "first")
point_colors <- colors[ranked]

# Prepare ternary points
points_matrix <- as.matrix(merged[, c("Bird", "Frog", "Insect")])

# Plot
TernaryPlot(alab = "Bird", blab = "Frog", clab = "Insect", grid.lines = 5)
TernaryPoints(points_matrix, col = point_colors, pch = 19, cex = 0.8)

# Optional: Add color legend
legend_gradient <- colorRampPalette(viridis(100))
image(
  z = matrix(seq(min(merged$wm_richness_ge3), max(merged$wm_richness_ge3), length.out = 100), ncol = 1),
  col = legend_gradient(100),
  xaxt = "n", yaxt = "n",
  xlab = "", ylab = "",
  add = TRUE,
  zlim = range(merged$wm_richness_ge3),
  useRaster = TRUE
)

# ========================================================================
# Box plot of birds, frogs and insects by site richness class
# ========================================================================

# Read the input CSVs
abgi <- read_csv(file.path(inDir, "data/wildmon_hour-level_species_taxon_abgi.csv"))
richness <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) %>% select(siteid, wm_richness_ge3)
sites <- read_csv(file.path(inDir,"data/biosoundscape_sites_daac_250507.csv")) %>% rename(siteid = SiteID)

# Filter minutes to dawn and dusk choruses
abgi <- abgi %>% filter((hour >= 5 & hour < 9) | (hour >= 16 & hour < 21))

# Merge datasets
merged <- abgi %>%
  inner_join(richness, by = "siteid")
merged <- merged %>%
  inner_join(sites, by = "siteid")

# Create RichnessClass factor using tertiles of wm_richness_ge3
merged <- merged %>%
  mutate(RichnessClass = cut(
    wm_richness_ge3,
    breaks = quantile(wm_richness_ge3, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("Low", "Moderate", "High"),
    include.lowest = TRUE
  ))

# Summarize by RichnessClass and Campaign
summary_df <- merged %>%
  group_by(RichnessClass, Campaign) %>%
  summarise(
    Birds = sum(Bird, na.rm = TRUE),
    Amphibians = sum(Frog, na.rm = TRUE),
    Insects = sum(Insect, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Birds, Amphibians, Insects), names_to = "Taxon", values_to = "Count") %>%
  mutate(
    pattern = ifelse(Campaign == "Dry season", "stripe", "none")
  )

# Plot A: Detection counts with hatching for dry season
pA <- ggplot(summary_df, aes(x = RichnessClass, y = Count, fill = Taxon, pattern = pattern)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_stack(),
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.3,
    pattern_spacing = 0.05
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::label_comma()) + 
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe")) +
  labs(y = "Detection Count", fill = "Taxon") +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank()) + 
  guides(pattern = "none") 

# Plot B: Relative proportion (percent stacked)
pB <- ggplot(summary_df, aes(x = RichnessClass, y = Count, fill = Taxon, pattern = pattern)) +
  geom_bar_pattern(
    stat = "identity",
    position = "fill",
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.3,
    pattern_spacing = 0.05
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Set2") +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe")) +
  labs(y = "Proportion of Detections", fill = "Taxon") +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank()) +
  guides(pattern = "none") 

# Combine plots
combined_plot <- ggarrange(
  pA, pB,
  ncol = 2, nrow = 1,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "right"
)

combined_plot <- annotate_figure(
  combined_plot,
  bottom = text_grob("Species Richness Class", size = 12)
)

# Display plot
dev.new(); print(combined_plot)

outFile <- file.path(inDir, "figures/site_richness_taxon.png")
ggsave(outFile, combined_plot, width = 7, height = 3, units = "in", dpi = 600, bg = "white")

# ========================================================================
# Box plot of AGI by site richness class
# ========================================================================

# Read the input CSVs
abgi <- read_csv(file.path(inDir, "data/wildmon_hour-level_species_taxon_abgi.csv"))
richness <- read.csv(file.path(inDir,"data/wildmon_site-level_species_250907.csv")) %>% select(siteid, wm_richness_ge3)
sites <- read_csv(file.path(inDir,"data/biosoundscape_sites_daac_250507.csv")) %>% rename(siteid = SiteID)

# Filter minutes to dawn and dusk choruses
abgi <- abgi %>% filter((hour >= 5 & hour < 9) | (hour >= 16 & hour < 21))

# Merge datasets
merged <- abgi %>%
  inner_join(richness, by = "siteid")
merged <- merged %>%
  inner_join(sites, by = "siteid")

# Create RichnessClass factor using tertiles of wm_richness_ge3
merged <- merged %>%
  mutate(RichnessClass = cut(
    wm_richness_ge3,
    breaks = quantile(wm_richness_ge3, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("Low", "Moderate", "High"),
    include.lowest = TRUE
  ))

# Summarize by RichnessClass and Campaign using ABGI fields
summary_df <- merged %>%
  group_by(RichnessClass, Campaign) %>%
  summarise(
    Anthropophony = sum(Anthropophony, na.rm = TRUE),
    Geophony = sum(Geophony, na.rm = TRUE),
    Interference = sum(Interference, na.rm = TRUE),
    Biophony = sum(Biophony, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Anthropophony, Geophony, Interference, Biophony),
    names_to = "AcousticBin",
    values_to = "Count"
  ) %>%
  mutate(
    pattern = ifelse(Campaign == "Dry season", "stripe", "none")
  )

# Plot A: Detection counts with hatching for dry season
pA <- ggplot(summary_df, aes(x = RichnessClass, y = Count, fill = AcousticBin, pattern = pattern)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_stack(),
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.3,
    pattern_spacing = 0.05
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::label_comma()) + 
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe")) +
  labs(y = "Detection Count", fill = "Soundscape\nComponent") +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank()) + 
  guides(pattern = "none") 

# Plot B: Relative proportion (percent stacked)
pB <- ggplot(summary_df, aes(x = RichnessClass, y = Count, fill = AcousticBin, pattern = pattern)) +
  geom_bar_pattern(
    stat = "identity",
    position = "fill",
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.3,
    pattern_spacing = 0.05
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Set2") +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe")) +
  labs(y = "Proportion of Detections", fill = "Soundscape\nComponent") +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank()) +
  guides(pattern = "none") 

# Combine plots
combined_plot <- ggarrange(
  pA, pB,
  ncol = 2, nrow = 1,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "right"
)

combined_plot <- annotate_figure(
  combined_plot,
  bottom = text_grob("Species Richness Class", size = 12)
)

# Display plot
dev.new(); print(combined_plot)

# Save figure
outFile <- file.path(inDir, "figures/site_richness_acoustic_abgi.png")
ggsave(outFile, combined_plot, width = 7, height = 3, units = "in", dpi = 600, bg = "white")

# ========================================================================
# Box plot of pattern matching spectral-temporal statistics
# ========================================================================

# Read in data
pattern_df <- read_csv(file.path(inDir, "data/biosoundscape_pattern_matching_240809.csv"))
guild_df <- read_csv(file.path(inDir, "data/wildmon_scientific_names_with_guild_filled.csv"))

# Get unique species and calls
pattern_df$species_songtype <- as.factor(paste0(pattern_df$scientific_name,"_",pattern_df$y1))
pattern_df <- pattern_df %>% 
  group_by(species_songtype) %>%
  slice(1) %>%
  ungroup()

guild_df <- guild_df %>%
  mutate(Guild = case_when(
    Guild == "Frog" ~ "Amphibians",
    Guild == "Insect" ~ "Insects",
    TRUE ~ "Birds"
  )) %>%
  distinct()

# Merge on genus
merged_df <- pattern_df %>%
  left_join(guild_df, by = "scientific_name") %>%
  mutate(duration = x2 - x1,
         frequency_range = (y2 - y1)/1000,
         min_frequency = pmin(y1, y2, na.rm = TRUE)/1000,
         max_frequency = pmax(y1, y2, na.rm = TRUE)/1000)

# Reshape to long format
plot_df <- merged_df %>%
  select(Guild, duration, frequency_range, min_frequency, max_frequency) %>%
  pivot_longer(cols = c(duration, frequency_range, min_frequency, max_frequency),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         duration = "Duration (s)",
                         frequency_range = "Frequency Range (kHz)",
                         min_frequency = "Minimum Frequency (kHz)",
                         max_frequency = "Maximum Frequency (kHz)"))

# Remove statistical outliers (1.5x IQR rule) by Metric and Guild
plot_df_clean <- plot_df %>%
  group_by(Guild, Metric) %>%
  mutate(Q1 = quantile(Value, 0.25, na.rm = TRUE),
         Q3 = quantile(Value, 0.75, na.rm = TRUE),
         IQR = Q3 - Q1,
         lower = Q1 - 1.5 * IQR,
         upper = Q3 + 1.5 * IQR) %>%
  filter(Value >= lower & Value <= upper) %>%
  ungroup()

# Plot without outliers in data
p <- ggplot(plot_df_clean, aes(x = Guild, y = Value, fill = Guild)) +
  geom_boxplot() +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(
       x = "Taxonomic Group",
       y = "Value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        axis.title.x = element_blank())

# Display plot
dev.new(); print(p)

# Save figure
outFile <- file.path(inDir, "figures/boxplot_pattern_matching_spectral-temporal.png")
ggsave(outFile, p, width = 6, height = 6, units = "in", dpi = 600, bg = "white")

# Function to run two-way ANOVA controlling for site
run_anova_tests <- function(df, metric) {
  formula <- as.formula(paste(metric, "~ Guild"))
  model <- aov(formula, data = df)
  cat("\nANOVA for", metric, "\n")
  print(summary(model))
  
  cat("Tukey HSD post-hoc for Guild\n")
  print(TukeyHSD(model, "Guild"))
}

# Filter only needed columns
stats_df <- merged_df %>%
  filter(Guild %in% c("Birds", "Amphibians", "Insects")) %>%
  select(Guild, duration, frequency_range, min_frequency, max_frequency,species_songtype )

# Run tests
run_anova_tests(stats_df, "duration")
run_anova_tests(stats_df, "frequency_range")
run_anova_tests(stats_df, "min_frequency")
run_anova_tests(stats_df, "max_frequency")

# ========================================================================
# GCFR CNN precision threshold plot
# ========================================================================

# We load the final results for the optimization of the penalization for all taxa:
load(file=file.path(inDir, "data/allSpp_basePenalization_optimals_250413.RData"))

# The file contains the base table: all 51 species evaluated without any penalizations
head(base)

# The list object optresPars has three tables. Each is the evaluation of all parameters under the optimal penalization of one of them.
# For example, optresPars$Fbeta has the evaluation of all parameters (Fbeta, precision, sensitivity) under the optimal penalization for Fbeta for each species separately
# To make the plot for Figure 3, we need:
pdf<-optresPars$Prec[,c("GVspeciesDef","Prec","pen","numMatches")]

# We want to color by taxa type , so we need to merge with the lookup table
nlut<-read.csv(file.path(inDir, "data/WildMon_BirdNET_lookup.csv"), stringsAsFactors=FALSE)
names(pdf)[1]<-"WMshort"
pdf<-merge(pdf,nlut[,c("WMshort","AMsciName","Class")], by="WMshort", all.x=TRUE)
pdf$Class<-ifelse(pdf$Class=="Frog","Amphibian",pdf$Class)

## Let's make grobs and use gridextra to arrange them in a grid, so that we can have the same legend for all panels
precp <- ggplot(pdf, aes(x = AMsciName, y = Prec)) + 
  geom_bar(stat = "identity", aes(fill = Class)) + 
  coord_flip() +
  theme_bw() + 
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8)  # adjust size here
  ) +
  labs(x = "", y = "Precision", fill = "Taxa") +
  scale_fill_brewer(palette = "Set2")

penp <- ggplot(pdf, aes(x = AMsciName, y = pen, fill = Class)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "Threshold", fill = "Taxa") +
  scale_fill_brewer(palette = "Set2")

numMatchp <- ggplot(pdf, aes(x = AMsciName, y = numMatches, fill = Class)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "", y = "Num. Matches", fill = "Taxa") +
  scale_fill_brewer(palette = "Set2")

# Combine plots
combined_plot <- ggarrange(
  precp, penp, numMatchp,
  ncol = 3, nrow = 1,
  widths = c(1, 0.75, 0.75),
  labels = c("A", "B", "C"),
  common.legend = TRUE,
  legend = "bottom"
)

# Display plot
dev.new(); print(combined_plot)

# Save figure
outFile <- file.path(inDir, "figures/precision_threshold_species.png")
ggsave(outFile, combined_plot, width = 10, height = 6, units = "in", dpi = 600, bg = "white")

# Statistics
summary(pdf$Prec); sd(pdf$Prec)
summary(pdf$pen); sd(pdf$pen)
pdf %>%
  group_by(Class) %>%
  summarise(
    n = n(),
    mean_prec = mean(Prec, na.rm = TRUE),
    sd_prec = sd(Prec, na.rm = TRUE),
    mean_conf = mean(pen, na.rm = TRUE),
    sd_conf = sd(pen, na.rm = TRUE)
  )
summary(lm(Prec ~ numMatches, data = pdf))


# ========================================================================
# ABGI precision threshold plot
# ========================================================================

# We load and plot the results of the ABGI optimization
load(file=file.path(inDir, "data/evalBeta_Threshold_GVresults_4BAGI_BirdNET_0Conf_noAbs_250110.RData"))

# The data.frame bagi4bnmax has the optimals, column Precopt
optdf<- bagi4bnmax[,c("GVspeciesDef","Fbetaopt","Precopt","Sensopt")]
# Rename the first column to BAGI, so that we can merge with pendf later
names(optdf)[1]<-"BAGI"

# The data.frame bagi4bnpen has the penalizations. 
pendf<- bagi4bnpen[,c("BAGI","Pen","Fbeta","Prec","Sens")]

# Must reshape long, by BAGI and Pen, so that we have a column with the index value and a column with the index name
pendfL <- reshape(pendf, varying = c("Fbeta", "Prec", "Sens"), v.names = "Value", timevar = "Parameter", times = c("F0.5", "Precision", "Sensitivity"), direction = "long")
row.names(pendfL)<-NULL
# merging the optdf with pendfL by BAGI, so that we have the optimal values in the same data.frame as the penalizations
pendfL<-merge(pendfL, optdf, by="BAGI")

# rename the BAGI values to be more descriptive
pendfL$BAGI<-ifelse(pendfL$BAGI=="B","Biophony",ifelse(pendfL$BAGI=="A","Anthropophony",ifelse(pendfL$BAGI=="G","Geophony",ifelse(pendfL$BAGI=="I","Interference",NA))))

# plot it, but NOTE which optimal parameter we use for the vertical dashed lines
p <- ggplot(pendfL,aes(x=Pen,y=Value)) + geom_line(aes(color=Parameter),linewidth=1.2) +
  facet_wrap(~BAGI) + theme_bw() + theme(legend.position = "bottom") + 
  labs(x="Confidence",y="Value",color="Parameter") +
  scale_y_continuous(limits=c(0,1)) + 
  scale_color_brewer(palette = "Set2") +
  geom_vline(aes(xintercept=Fbetaopt),color="black",linetype="dashed", linewidth=1.1) # Using the optimal Fbeta
# Could use Precopt (optimal precision) or Sensopt (optimal sensitivity)

# Display plot
dev.new(); print(p)

# Save figure
outFile <- file.path(inDir, "figures/precision_threshold_abgi.png")
ggsave(outFile, p, width = 7, height = 5, units = "in", dpi = 600, bg = "white")

# Statistics
pendfL %>% group_by(BAGI) %>% filter(Pen == Fbetaopt)
