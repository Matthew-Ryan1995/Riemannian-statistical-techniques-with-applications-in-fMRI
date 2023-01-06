
# Load packages -----------------------------------------------------------

pacman::p_load(tidyverse, phdWork, MDMR)


# Set up parameters -------------------------------------------------------

set.seed(2022 + 1668286)
num_perm <- 999
df <- cobre

# Get distance matrices ---------------------------------------------------

if(file.exists("results/cobre_geometric_D.Rds")){
  D_geo <- readRDS("results/cobre_geometric_D.Rds")
}else{
  D_geo <- affine_dist(df, "cors", regularise = FALSE, cores = 4)
  saveRDS(D_geo, "results/cobre_geometric_D.Rds" )
}
if(file.exists("results/cobre_euclidean_D.Rds")){
  D_euc <- readRDS("results/cobre_euclidean_D.Rds")
}else{
  D_euc <- euclid_dist(df, "cors")
  saveRDS(D_euc, "results/cobre_euclidean_D.Rds" )
}
if(file.exists("results/cobre_correlation_D.Rds")){
  D_corr <- readRDS("results/cobre_correlation_D.Rds")
}else{
  D_corr <- corr_dist(df, "cors")
  saveRDS(D_corr, "results/cobre_correlation_D.Rds" )
}



# Model reduction on Geo --------------------------------------------------
# Consider an interaction term in gender and group there are less females in the patient group
# May or may not be interesting
# We will perform backwards selection based on p-values using alpha = 0.05

X_geo <- model.matrix(~age + handedness + group * gender, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 2))

# Least significant is the interaction term
X_geo <- model.matrix(~age + handedness + group + gender, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 2)) 

# Least significant is the handedness term
X_geo <- model.matrix(~age  + group + gender, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 2)) 

# Least significant is the gender term
X_geo <- model.matrix(~age  + group, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 2)) 

# Everything signficant
mdmr_geo <- summary(mdmr_geo) %>% 
  as_tibble(rownames = "Predictor") %>%
  mutate(method = "Geometric") %>% 
  janitor::clean_names()

# Model reduction on Euc --------------------------------------------------

X_euc <- model.matrix(~age + handedness + group * gender, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 2))

# Least significant is the gender term
# Due to principal of marginality, we keep this and instead remove the interaction
X_euc <- model.matrix(~age + handedness + group + gender, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 2))

# Least significant is the handedness term
X_euc <- model.matrix(~age  + group + gender, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 2))

# Least significant is the gender term
X_euc <- model.matrix(~age  + group, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 2))

# Everything signficant
mdmr_euc <- summary(mdmr_euc) %>% 
  as_tibble(rownames = "Predictor") %>% 
  mutate(method = "Euclidean") %>% 
  janitor::clean_names()

# Model reduction on Corr --------------------------------------------------

X_corr <- model.matrix(~age + handedness + group * gender, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 2))

# Least significant is the gender term
# Due to principal of marginality, we keep this and instead remove the interaction
X_corr <- model.matrix(~age + handedness + group + gender, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 2))

# Least significant is the handedness term
X_corr <- model.matrix(~age  + group + gender, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 2))

# Least significant is the gender term
X_corr <- model.matrix(~age  + group, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 2))

# Everything signficant
mdmr_corr <- summary(mdmr_corr) %>% 
  as_tibble(rownames = "Predictor") %>% 
  mutate(method = "Correlation") %>% 
  janitor::clean_names()


# Final results table -----------------------------------------------------

full_results <- mdmr_geo %>% 
  bind_rows(
    mdmr_euc,
    mdmr_corr
  )

saveRDS(full_results, "results/cobre_mdmr_results.Rds")
