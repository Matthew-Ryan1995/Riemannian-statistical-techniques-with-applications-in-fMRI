
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



# Model reduction on Geo --------------------------------------------------
# Consider interaction between group and age, viq, and fiq based on data exploration
# May or may not be interesting
# We will perform backwards selection based on p-values using alpha = 0.05

X_geo <- model.matrix(~age + group + sex + fiq + viq + piq + eye +
                        group:age + group:viq + group:fiq, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 4))

# Least significant is piq
X_geo <- model.matrix(~age + group + sex + fiq + viq  + eye +
                        group:age + group:viq + group:fiq, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 4)) 

# Least significant is group, but based on the principal of marginality we keep and remove age:group
X_geo <- model.matrix(~age + group + sex + fiq + viq  + eye +
                        group:viq + group:fiq, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 4)) 

# Least significant is group, but based on the principal of marginality we keep and remove fiq:group
X_geo <- model.matrix(~age + group + sex + fiq + viq  + eye +
                        group:viq, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 4)) 

# Least significant is group, but based on the principal of marginality we keep and remove viq:group
X_geo <- model.matrix(~age + group + sex + fiq + viq  + eye, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 4)) 

# Least significant is fiq
X_geo <- model.matrix(~age + group + sex + viq  + eye, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 4)) 

# Least significant is viq
X_geo <- model.matrix(~age + group + sex + eye, data = df)
(mdmr_geo <- mdmr(X = X_geo, D = D_geo, nperm = num_perm, ncores = 4)) 

# Everything signficant
mdmr_geo <- summary(mdmr_geo) %>% 
  as_tibble(rownames = "Predictor") %>%
  mutate(method = "Geometric") %>% 
  janitor::clean_names()

# Model reduction on Euc --------------------------------------------------

X_euc <- model.matrix(~age + group + sex + fiq + viq + piq + eye +
                        group:age + group:viq + group:fiq, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 4))

# Least significant is the gender term
# Due to principal of marginality, we keep this and instead remove the group:fiq
X_euc <- model.matrix(~age + group + sex + fiq + viq + piq + eye +
                        group:age + group:viq, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 4))

# Least significant is group:viq
X_euc <- model.matrix(~age + group + sex + fiq + viq + piq + eye +
                        group:age, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 4))

# Least significant is piq
X_euc <- model.matrix(~age + group + sex + fiq + viq + eye +
                        group:age, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 4))

# Least significant is fiq
X_euc <- model.matrix(~age + group + sex + viq + eye +
                        group:age, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 4))

# Least significant is viq
X_euc <- model.matrix(~age + group + sex  + eye +
                        group:age, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 4))

# Least significant is age:group
X_euc <- model.matrix(~age + group + sex  + eye, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 4))

# Least significant is sex
X_euc <- model.matrix(~age + group + eye, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 4))

# Least significant is group
X_euc <- model.matrix(~age + eye, data = df)
(mdmr_euc <- mdmr(X = X_euc, D = D_euc, nperm = num_perm, ncores = 4))

# Everything signficant
mdmr_euc <- summary(mdmr_euc) %>% 
  as_tibble(rownames = "Predictor") %>% 
  mutate(method = "Euclidean") %>% 
  janitor::clean_names()

# Model reduction on Corr --------------------------------------------------

X_corr <- model.matrix(~age + group + sex + fiq + viq + piq + eye +
                        group:age + group:viq + group:fiq, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 4))

# Least significant is piq
X_corr <- model.matrix(~age + group + sex + fiq + viq + eye +
                         group:age + group:viq + group:fiq, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 4))

# Least significant is group
# Due to the principal of marginality, we remove group:fiq
X_corr <- model.matrix(~age + group + sex + fiq + viq + eye +
                         group:age + group:viq, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 4))

# Least significant is the group:viq
X_corr <- model.matrix(~age + group + sex + fiq + viq + eye +
                         group:age, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 4))

# Least significant is the group:age
X_corr <- model.matrix(~age + group + sex + fiq + viq + eye, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 4))

# Least significant is the fiq
X_corr <- model.matrix(~age + group + sex + viq + eye, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 4))

# Least significant is the fiq
X_corr <- model.matrix(~age + group + sex + eye, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 4))

# Least significant is the group
X_corr <- model.matrix(~age + sex + eye, data = df)
(mdmr_corr <- mdmr(X = X_corr, D = D_corr, nperm = num_perm, ncores = 4))

# Everything significant, sex is borderline, so we will keep it
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

saveRDS(full_results, "results/abide_mdmr_results.Rds")
