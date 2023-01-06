# tictoc::tic()
#### PREAMBLE ####
# Sys.sleep(sample(1:10, 1)) # Sleep the system so that no 2 jobs run at exactly the same time
print(date())
## Packages ##
pacman::p_load(tidyverse, MDMR, parallel, phdWork)

## Detect cores from phoenix
slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (!is.na(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
}else {
  cores = detectCores() - 1 # Figure out how many cores there are
}

job_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print(str_c("Doing job ", job_num))

# registerDoParallel(cores)


#### DATA PREPARATION ###
## Functions ##
# Load the necessary functions
# source("code/fc_simulation_functions.R")
# source("code/get_geo_dist6.R")
# source("code/other_distances.R")

## Defining Parameters ##
# Define all parameters needed for code
num_sims <- 999
num_perms <- 999

# sim_params <- expand_grid(block = 2:4,
#                           correlation = atanh(seq(-0.95, 0.95, by = 0.05)),
#                           se = c(0.01, 0.05, 0.1, 0.2))
sim_params <- expand_grid(block = 2:4,
                          correlation = atanh(c(-0.95, -0.5, 0, 0.5, 0.95)),
                          se = 0.267,
                          r = seq(0, 1, by = 0.1))

## Data ##
# Load the necessary data
df <- cobre
# df <- df %>% 
#   mutate(cors = map(msdl_roi, cor))

# X <- select(df, group)
X <- tibble(group = rep(c("Patient", "Control"), each = 75)) %>% 
  mutate(group = factor(group, c("Control", "Patient")))

#### SIMULATION CODE ####

dir_name <- str_c("power_results/sim_params_", job_num)

dir.create(dir_name) # Create directory

i <- job_num

correlation <- sim_params$correlation[i]
block <- sim_params$block[i]
se <- sim_params$se[i]
r <- sim_params$r[i]
# Set a seed for reproducibility
RNGkind("L'Ecuyer-CMRG")
# set.seed(correlation + block + r + i) # This gives a unique seed for each set of parameters
set.seed(2022 + i) # This gives a unique seed for each set of parameters
mc.reset.stream()
sims <- mclapply(1:num_sims, # Run the simulations
                 function(sim){
                   sub_sims <- map(X$group, # Generate simulation group
                                   ~simulate_subject(rho = correlation,
                                                     block_size = block,
                                                     group = .x,
                                                     method = "ar1",
                                                     vary = FALSE,
                                                     se = se,
                                                     r = r))
                   
                   # Calculate distance matrix
                   
                   D_geo <- affine_dist(sub_sims, regularise = FALSE, cores = 1)
                   
                   D_euc <- euclid_dist(sub_sims)
                   
                   D_cor <- corr_dist(sub_sims)
                   
                   
                   # Calculate MDMR results on 999 permutations
                   # Save results as a tibble
                   mdmr_output <- mdmr(X = X, D = D_geo, perm.p = TRUE,
                                       nperm = num_perms, ncores = 1) %>%
                     summary() %>%
                     as_tibble(rownames = "predictor") %>%
                     mutate(block = block,
                            rho = correlation,
                            se = se,
                            method = "geodesic",
                            sim = sim,
                            r = r)
                   
                   mdmr_output <- mdmr_output %>% 
                     bind_rows(
                       mdmr(X = X, D = D_euc, perm.p = TRUE,
                            nperm = num_perms, ncores = 1) %>%
                         summary() %>%
                         as_tibble(rownames = "predictor") %>%
                         mutate(block = block,
                                rho = correlation,
                                se = se,
                                method = "euclidean",
                                sim = sim,
                                r = r)
                     )
                   
                   mdmr_output <- mdmr_output %>% 
                     bind_rows(
                       mdmr(X = X, D = D_cor, perm.p = TRUE,
                            nperm = num_perms, ncores = 1) %>%
                         summary() %>%
                         as_tibble(rownames = "predictor") %>%
                         mutate(block = block,
                                rho = correlation,
                                se = se,
                                method = "correlation",
                                sim = sim,
                                r = r)
                     )
                   
                   return(mdmr_output)
                 },
                 mc.cores = cores)



sims <- map_dfr(sims, ~.x) # Convert all results into dataframe

saveRDS(sims, str_c(dir_name, "/sim_", i)) # Save dataframe

print(str_c("Job complete, power, sim number ", i))

print(date())