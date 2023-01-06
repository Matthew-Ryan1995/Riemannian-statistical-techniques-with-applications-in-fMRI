## Matt Ryan
## 29/08/2022
## EDIT: 16/11/2022
# tictoc::tic()
#
# Packages ----------------------------------------------------------------
pacman::p_load(tidyverse, phdWork)

print(date())
## Detect cores from phoenix
slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (!is.na(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
}else {
  cores = detectCores() - 1 # Figure out how many cores there are
}

job_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print(str_c("Doing job ", job_num))




# Load functions ----------------------------------------------------------
walk(c("code/gen_data.R","code/perform_pls.R", "code/summary_stats.R", "code/fit_simulations.R"),
     source)

# Set up parameters -------------------------------------------------------

# Dimensions
# p <- q <- c(10, 50, 100, 200)
#
# # Sample size
# n <- c(10, 50, 100, 150)
#
# # Noise
# sig <- c(0, 0.1, 0.5, 1, 2)

n <- round(seq(10, 150, length.out = 5))

# Dimensions
# p <- q <- round(seq(5, 15, length.out = 5))
p <- q <- choose(round(seq(5, 15, length.out = 5)) + 1, 2)


# Noise
sig <- c(0, 0.1, 0.5, 1, 0.05)



# Number of sims
N_sims <- 100

predict_response_grid <- expand_grid(p, q)

i <- job_num

p_sim <- predict_response_grid$p[i]
q_sim <- predict_response_grid$q[i]

fit_simulation(sims = N_sims, n = n,
               p = p_sim, q = q_sim,
               sig = sig, L = 5, 10295674 * i)
# tictoc::toc()
#
print(str_c("Job complete, power, sim number ", i))

print(date())

