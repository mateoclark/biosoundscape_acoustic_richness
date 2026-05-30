library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(tidyr)
library(RColorBrewer)
library(tidyverse)
library(mgcv)

inDir <- "G:/Shared drives/BioSoundSCape/Paper development/RQ1.1 Animal-acoustic diversity relationships"


# Read in the data
richness <- read.csv(paste0(inDir,"/data/bySeason_RichnessEsts_F05_Q0.1_correctedFiltered_v20250602.csv")) 
embeddings <- read.csv(paste0(inDir,"/umap_tpd_modeling/BirdNET_pca_umap_tpd_stats/",
                              "birdnet_embeddings_minute_nonzero_mean_dawn-dusk_dry-wet_pca_umap_tpd_nn20_mdist0.99_UMAP5-UMAP1-UMAP3-UMAP2.csv"))
sites <- read_csv(paste0(inDir,"/data/biosoundscape_sites_daac_250507.csv"))

# Merge and filter
merged <- embeddings %>%
  inner_join(richness, by = "SiteID") %>%
  filter(!is.na(QuantRichness),
         QuantRichness >= 0)


df <- left_join(merged, sites, by = "SiteID")
df <- df %>%
  mutate(LogRecordingNum = log(RecordingNum)) %>%
  drop_na(QuantRichness, FRichness, FEvenness, FDivergence)

# GAM modeling with smooth spatial structure
gam_model_full <- gam(
  QuantRichness ~ FRichness + FEvenness + FDivergence + 
    s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)
summary(gam_model_full)

gam_model_fric <- gam(
  QuantRichness ~ FRichness + 
    s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)
summary(gam_model_fric)

gam_model_feve <- gam(
  QuantRichness ~ FEvenness + 
    s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)
summary(gam_model_feve)

gam_model_fdiv <- gam(
  QuantRichness ~ FDivergence + 
    s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)
summary(gam_model_fdiv)

# Predict from GAM models
df <- df %>%
  mutate(
    GAM_Full = predict(gam_model_full, type = "response"),
    GAM_FRic = predict(gam_model_fric, type = "response"),
    GAM_FEve = predict(gam_model_feve, type = "response"),
    GAM_FDiv = predict(gam_model_fdiv, type = "response")
  )

# Convert to long format
df_long <- df %>%
  select(QuantRichness, GAM_Full, GAM_FRic, GAM_FEve, GAM_FDiv) %>%
  pivot_longer(cols = starts_with("GAM_"), names_to = "Model", values_to = "Predicted") %>%
  mutate(
    Model = recode(Model,
                   GAM_Full = "Full",
                   GAM_FRic = "Fric",
                   GAM_FEve = "Feve",
                   GAM_FDiv = "Fdiv"),
    Model = fct_drop(as.factor(Model))
  )

# Calculate R² values per model
r2_values <- df_long %>%
  group_by(Model) %>%
  summarise(R2 = round(cor(QuantRichness, Predicted)^2, 2)) %>%
  mutate(
    label = paste0("R² = ", R2),
    Model = fct_drop(as.factor(Model))
  )

# Plot with facets
p <- ggplot(df_long, aes(x = QuantRichness, y = Predicted)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  facet_wrap(~Model) +
  geom_text(data = r2_values, aes(x = Inf, y = -Inf, label = label),
            hjust = 1.1, vjust = -1.2, inherit.aes = FALSE) +
  labs(
    x = "Observed Richness",
    y = "Predicted Richness",
  ) +
  theme_minimal(base_size = 12)

# Show and save plot
dev.new(); print(p)

outFile <- paste0(inDir, "/figures/AIrichness_tpd_gam_regression_scatterplot.png")
ggsave(outFile, p, width = 6, height = 4, units = "in", dpi = 600, bg = "white")



#################################
# Extra exploration
#################################

## More robust regression models
library(MASS)
library(rcompanion)
library(randomForest)
library(spaMM)

sites <- read_csv(paste0(inDir,"/data/biosoundscape_sites_daac_250507.csv"))
df <- left_join(merged, sites, by = "SiteID")
df <- df %>%
  mutate(LogRecordingNum = log(RecordingNum)) %>%
  drop_na(QuantRichness, FRichness, FEvenness, FDivergence)

glm_richness <- glm(QuantRichness ~ FRichness + FEvenness + FDivergence,
                    #glm_richness <- glm(QuantRichness ~ FRichness,
                    offset = LogRecordingNum,
                    family = "poisson",
                    data = df)
glm_null <- glm(QuantRichness ~ 1,
             offset = LogRecordingNum,
             family = "poisson",
             data = df)

R2 <- nagelkerke(glm_richness, null = glm_null)
R2$Pseudo.R.squared.for.model.vs.null

dispersion <- sum(residuals(glm_richness, type = "pearson")^2) / df.residual(glm_richness)
print(dispersion)  # >1 suggests overdispersion

nb_model <- glm.nb(QuantRichness ~ FRichness + FEvenness + FDivergence + offset(LogRecordingNum),
                   #nb_model <- glm.nb(QuantRichness ~ FRichness + offset(LogRecordingNum),
                   data = df)

summary(nb_model)

model_null <- glm.nb(QuantRichness ~ 1 + offset(LogRecordingNum), data = df)

# R²
model_null <- glm.nb(QuantRichness ~ 1 + offset(LogRecordingNum), data = df)
R2 <- nagelkerke(nb_model, null = model_null)
R2$Pseudo.R.squared.for.model.vs.null

avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spatial_model <- fitme(QuantRichness ~ FRichness + FEvenness + FDivergence +
                         Matern(1|Longitude + Latitude) + 
                         offset(LogRecordingNum),
                       family = negbin(), data = df, control.HLfit=list(NbThreads=max(avail_thr, 1L)))

# Spatial model R²
spatial_model_null <- fitme(QuantRichness ~ 1 +
                              Matern(1|Longitude + Latitude) + 
                              offset(LogRecordingNum),
                            family = negbin(), data = df, control.HLfit=list(NbThreads=max(avail_thr, 1L)))

# Extract log-likelihoods and sample size
loglik_null <- as.numeric(logLik(spatial_model_null))
loglik_full <- as.numeric(logLik(spatial_model))
n <- nobs(spatial_model)

# McFadden's R²
R2_mcfadden <- 1 - (loglik_full / loglik_null)
R2_mcfadden

# Compute Nagelkerke R²
R2_nagelkerke <- (1 - exp((2 / n) * (loglik_null - loglik_full))) /
  (1 - exp((2 / n) * loglik_null))
R2_nagelkerke


# GAM modeling with smooth spatial structure
gam_model <- gam(
  QuantRichness ~ FRichness + FEvenness + FDivergence + 
    s(Latitude, Longitude) + offset(LogRecordingNum),
  family = nb(),
  data = df
)

summary(gam_model)

# Fit a random forest regression model
rf_model <- randomForest(QuantRichness ~ FRichness + FEvenness + FDivergence,
                         data = df,
                         importance = TRUE,
                         ntree = 500)

# Summary of the model
print(rf_model)

# Variable importance
importance(rf_model)
varImpPlot(rf_model)