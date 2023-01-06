
# libraries ---------------------------------------------------------------
pacman::p_load(tidyverse, phdWork, rsample, yardstick)

walk(c("code/cross_validate_pls.R"), source)



print(date())
## Detect cores from phoenix
slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (!is.na(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
}else {
  cores = detectCores() - 1 # Figure out how many cores there are
}

# job_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print(str_c("Doing job COBRE"))


# Data --------------------------------------------------------------------
df <- cobre %>%
  select(cors, age, group) %>%
  mutate(
    group = fct_relevel(group, "Control")
  )

model_formula <- ~ age + group


# seed and CV -------------------------------------------------------------

set.seed(1668286)

folds <- rsample::vfold_cv(df, v = 10, strata = group)


# run cv ------------------------------------------------------------------

results <- cv_rpls_fit(folds, model_formula, mc.cores = cores)

write_rds(results, "results/cobre_cv_results.Rds")


print(str_c("Job complete"))

print(date())
