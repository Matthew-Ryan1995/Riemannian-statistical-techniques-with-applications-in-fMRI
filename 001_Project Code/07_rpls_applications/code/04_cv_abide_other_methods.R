
# libraries ---------------------------------------------------------------
pacman::p_load(tidyverse, phdWork, tidymodels, plsmod)

walk(c("code/cross_validate_other.R"), source)

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
  select(cors, age, group, sex, eye) %>%
  mutate(cors = map(cors, ~(.x + diag(1, nrow(.x)))),
         group = factor(group, levels = c("Control", "Autism")),
         sex = as.factor(sex),
         eye = factor(eye, levels = c("Open", "Closed")))



# seed and CV -------------------------------------------------------------

set.seed(1985)

folds <- rsample::vfold_cv(df, v = 10, strata = group)


# run cv ------------------------------------------------------------------

results <- cv_other_fit(folds, mc.cores = cores)

write_rds(results, "results/abide_cv_results_other.Rds")



print(str_c("Job complete"))

print(date())
