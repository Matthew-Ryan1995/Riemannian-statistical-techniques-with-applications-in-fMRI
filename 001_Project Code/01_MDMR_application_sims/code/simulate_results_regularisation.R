# tictoc::tic()
#### PREAMBLE ####
# Sys.sleep(sample(1:10, 1)) # Sleep the system so that no 2 jobs run at exactly the same time
print(date())
## Packages ##
pacman::p_load(tidyverse, phdWork, MDMR, parallel)

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
# # Load the necessary functions
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
# sim_params <- expand_grid(block = 2:4,
#                           correlation = atanh(c(-0.95, -0.5, 0, 0.5, 0.95)),
#                           se = 0.267,
#                           r = seq(0, 1, by = 0.1))
#                           

lambda <- seq(1e-5, 5, length.out = 100)

## Data ##
# Load the necessary data
df <- cobre #readRDS("data/df_msdl")
# df <- df %>% 
#   mutate(cors = map(msdl_roi, cor))

X <- select(df, group)

#### SIMULATION CODE ####

dir_name <- str_c("regularisation_results/sim_params_", job_num)

dir.create(dir_name) # Create directory
getwd()
i <- job_num
# 
# correlation <- sim_params$correlation[i]
# block <- sim_params$block[i]
# se <- sim_params$se[i]
# r <- sim_params$r[i]
lambda_i <- lambda[i]
# Set a seed for reproducibility
RNGkind("L'Ecuyer-CMRG")
# set.seed(correlation + block + r + i) # This gives a unique seed for each set of parameters
set.seed(1995 + i) # This gives a unique seed for each set of parameters
mc.reset.stream()
sims <- mclapply(1:num_sims, # Run the simulations
                 function(sim){
                   sub_sims <- map(1:nrow(X), # Generate simulation group
                                   ~simulate_subject(rho = 0,
                                                     block_size = 4,
                                                     group = "Control",
                                                     method = "ar1",
                                                     low_rank = TRUE,
                                                     rank_remove_max = 23,
                                                     rank_prop = 0.15))
                   
                   # Calculate distance matrix
                   
                   D_geo <- affine_dist(sub_sims, regularise = TRUE, lambda = lambda_i, cores = 1)
                   
                   
                   
                   # Calculate MDMR results on 999 permutations
                   # Save results as a tibble
                   mdmr_output <- mdmr(X = X, D = D_geo, perm.p = TRUE,
                                       nperm = num_perms, ncores = 1) %>%
                     summary() %>%
                     as_tibble(rownames = "predictor") %>%
                     mutate(sim = sim,
                            lam = lambda_i)
                   
                   return(mdmr_output)
                 },
                 mc.cores = cores)



sims <- map_dfr(sims, ~.x) # Convert all results into dataframe

saveRDS(sims, str_c(dir_name, "/sim_", i, ".Rds")) # Save dataframe

print(str_c("Job complete, regularisation, sim number ", i))

print(date())