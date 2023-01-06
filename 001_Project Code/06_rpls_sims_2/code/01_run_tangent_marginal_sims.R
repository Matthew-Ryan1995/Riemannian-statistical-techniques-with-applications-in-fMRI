## Matt Ryan
## 05/09/2022
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
walk(c("code/collect_model_information.R","code/gen_data_and_fit_model.R",
       "code/measure_model_performance.R", "code/run_simulations.R", "code/prediction_functions.R"),
     source)

# Set up parameters -------------------------------------------------------

set.seed(5005) # This is to generate the list of seed values for my function

# Number of sims
N_sims <- 100
sim_type <- "tangent"

# Sample size
n <- round(seq(10, 150, length.out = 5))

# Dimensions
p <- q <- round(seq(5, 15, length.out = 5))


# Noise
sig <- c(0, 0.1, 0.5, 1, 2)


# Latent vars
L <- 5
K <- 1:(2*L)

# response and predictor types
type_x <- type_y <-  c("affine", "euclidean")

# Tolerance/iteration
tol <- 1e-4#c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
max.iter <- 20#seq(10, 50, by = 10)

# centre distance

r <- seq(0, 8, by = 0.5)

# Mariginalise
n <- n[3]
p <- p[3]
q <- q[3]
sig <- 0.05#sig[3]

predict_response_grid <- expand_grid(
  r = r,
  K = K,
  type_x = type_x,
  type_y = type_y,
) %>%
  mutate(seed_list = sample(1:100000, size = n(), replace = FALSE)) %>%
  filter(!(type_x == "euclidean" & type_y == "euclidean"), K<=5)

i <- job_num

r_sim <- predict_response_grid$r[i]
K_sim <- predict_response_grid$K[i]
type_x_sim <- predict_response_grid$type_x[i]
type_y_sim <- predict_response_grid$type_y[i]
seed_sim <- predict_response_grid$seed_list[i]


run_simulations(num_sims = N_sims, n = n,
                x_dim = p, y_dim = q, r = r_sim,
                sig = sig, L = L, seed = seed_sim, K = K_sim,
                type_x = type_x_sim, type_y = type_y_sim, sim_type = sim_type,
                mc.cores = cores, tol = tol, max.iter = max.iter)
# tictoc::toc()
#
print(str_c("Job complete, array number ", i))

print(date())

