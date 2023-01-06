#### PREAMBLE ####
Sys.sleep(sample(1:5, 1)) # Sleep the system so that no 2 jobs run at exactly the same time

## Packages ##
pacman::p_load(tidyverse, MDMR, parallel)

## Detect cores from phoenix
slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (!is.na(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
}else {
  cores = detectCores() - 1 # Figure out how many cores there are
}

#### DATA PREPARATION ###
## Functions ##
# Load the necessary functions
source("code/fc_simulation_functions.R")
source("code/get_geo_dist6.R")

## Defining Parameters ##
# Define all parameters needed for code
num_sims <- 100

sim_params <- expand_grid(block = 2:4,
                          correlation = atanh(seq(-0.95, 0.95, by = 0.05)),
                          se = c(0.01, 0.05, 0.1, 0.2))
sim_params <- expand_grid(block = 3:4,
                          correlation = -1.8,
                          se = 0.2)

## Data ##
# Load the necessary data
df <- readRDS("data/df_msdl")
df <- df %>% 
  mutate(cors = map(msdl_roi, cor))

X <- select(df, age:group)

#### SIMULATION CODE ####

# For loop to allow one job per array
for(i in 1:dim(sim_params)[1]){
  Sys.sleep(sample(5:10, 1)) # Rest system so no 2 jobs are running the same time
  
  dir_name <- str_c("results/sim_params_", i) # Name new directory
  if(!dir.exists(dir_name)){ # If job has not run, run it
    dir.create(dir_name) # Create directory
    
    correlation <- sim_params$correlation[i]
    block <- sim_params$block[i]
    se <- sim_params$se[i] # pull out necessary parameters
    
    sims <- mclapply(1:num_sims, # Run the simulations
                     function(sim){
                       sub_sims <- map(X$group, # Generate simulation group
                                       ~simulate_subject(rho = correlation, 
                                                         block_size = block,
                                                         group = .x,
                                                         method = "toeplitz",
                                                         vary = TRUE,
                                                         se = se))
                       
                       # Calculate distance matrix
                       D <- get_geo_dist(sub_sims, regularise = FALSE, cores = 1)
                       
                       # Calculate MDMR results on 5000 permutations
                       # Save results as a tibble
                       mdmr_output <- mdmr(X = X, D = D, perm.p = TRUE, 
                                           nperm = 5000, ncores = 1) %>% 
                         summary() %>% 
                         as_tibble(rownames = "predictor") %>% 
                         mutate(block = block, 
                                rho = correlation,
                                se = se,
                                sim = sim)
                       
                       return(mdmr_output)
                     },
                     mc.cores = cores)
    
    sims <- map_dfr(sims, ~.x) # Convert all results into dataframe
    
    saveRDS(sims, str_c(dir_name, "/sim_", i)) # Save dataframe
    break
  }
}
