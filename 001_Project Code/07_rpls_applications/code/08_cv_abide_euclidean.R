
# libraries ---------------------------------------------------------------
pacman::p_load(tidyverse, phdWork, rsample, yardstick)

walk(c("code/cross_validate_pls2.R"), source)

print(date())
## Detect cores from phoenix
slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (!is.na(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
}else {
  cores = detectCores() - 1 # Figure out how many cores there are
}

# job_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print(str_c("Doing job ABIDE"))


# Data --------------------------------------------------------------------
df <- abide_aal %>%
  select(age, group, sex, eye, cors) %>%
  mutate(group = factor(group, levels = c("Control", "Autism")),
         gender = as.factor(sex),
         eye = factor(eye, levels = c("Open", "Closed")),
         cors = map(cors, function(c){
           c <- c[upper.tri(c)]
           names(c) <- str_c("X", 1:length(c))
           c <- as_tibble(t(c))
           return(c)
         })
  ) %>%
  select(-sex) %>%
  unnest(cors)

model_formula <- ~ age + group + gender + eye


# seed and CV -------------------------------------------------------------

set.seed(1985)

folds <- rsample::vfold_cv(df, v = 10, strata = group)


# run cv ------------------------------------------------------------------

results <- cv_rpls_fit(folds, model_formula, mc.cores = cores)

write_rds(results, "results/abide_cv_results_euclidean.Rds")

# Fisher --------------------------------------------------------------------
df <- abide_aal %>%
  select(age, group, sex, eye, cors) %>%
  mutate(group = factor(group, levels = c("Control", "Autism")),
         gender = as.factor(sex),
         eye = factor(eye, levels = c("Open", "Closed")),
         cors = map(cors, function(c){
           c <- atanh(c[upper.tri(c)])
           names(c) <- str_c("X", 1:length(c))
           c <- as_tibble(t(c))
           return(c)
         })
  ) %>%
  select(-sex) %>%
  unnest(cors)

model_formula <- ~ age + group + gender + eye


# seed and CV -------------------------------------------------------------

set.seed(1985)

folds <- rsample::vfold_cv(df, v = 10, strata = group)


# run cv ------------------------------------------------------------------

results <- cv_rpls_fit(folds, model_formula, mc.cores = cores)

write_rds(results, "results/abide_cv_results_fisher.Rds")

print(str_c("Job complete"))

print(date())
