
# Load packages -----------------------------------------------------------

pacman::p_load(tidyverse, phdWork, MDMR)


# Set up parameters -------------------------------------------------------

set.seed(1995 + 1668286)
num_perm <- 999
df <- abide_aal

# Get distance matrices ---------------------------------------------------

if(file.exists("results/abide_geometric_D.Rds")){
  D_geo <- readRDS("results/abide_geometric_D.Rds")
}else{
  D_geo <- affine_dist(df, "cors", regularise = TRUE, lambda = 1, cores = 4)
  saveRDS(D_geo, "results/abide_geometric_D.Rds" )
}
if(file.exists("results/abide_euclidean_D.Rds")){
  D_euc <- readRDS("results/abide_euclidean_D.Rds")
}else{
  D_euc <- euclid_dist(df, "cors")
  saveRDS(D_euc, "results/abide_euclidean_D.Rds" )
}
if(file.exists("results/abide_correlation_D.Rds")){
  D_corr <- readRDS("results/abide_correlation_D.Rds")
}else{
  D_corr <- corr_dist(df, "cors")
  saveRDS(D_corr, "results/abide_correlation_D.Rds" )
}


X_geo <- model.matrix(~age + group + sex + eye, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 4)) 

# Everything signficant
mdmr_geo <- summary(mdmr_geo) %>% 
  as_tibble(rownames = "Predictor") %>%
  mutate(method = "Geometric") %>% 
  janitor::clean_names()
(mdmr_euc <- mdmr(X = X_geo, D = D_euc, nperm = num_perm, ncores = 4))

# Everything signficant
mdmr_euc <- summary(mdmr_euc) %>% 
  as_tibble(rownames = "Predictor") %>% 
  mutate(method = "Euclidean") %>% 
  janitor::clean_names()
(mdmr_corr <- mdmr(X = X_geo, D = D_corr, nperm = num_perm, ncores = 4))

# Everything significant, sex is borderline, so we will keep it
mdmr_corr <- summary(mdmr_corr) %>% 
  as_tibble(rownames = "Predictor") %>% 
  mutate(method = "Correlation") %>% 
  janitor::clean_names()
full_results <- mdmr_geo %>% 
  bind_rows(
    mdmr_euc,
    mdmr_corr
  )
# Final results table -----------------------------------------------------

full_results <- mdmr_geo %>% 
  bind_rows(
    mdmr_euc,
    mdmr_corr
  )

saveRDS(full_results, "results/abide_mdmr_results_compare.Rds")
