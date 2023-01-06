
# libraries ---------------------------------------------------------------
pacman::p_load(tidyverse, phdWork)

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
M <- read_rds("results/cobre_rpls_model_K_3.Rds")

Y <- M$Y[[1]]
X <- M$X[[1]]
K <- M$L

set.seed(1668286)

j <- ncol(X)

res <- parallel::mclapply(
  1:j,
  function(i){
    X_tmp <- X
    beta_list <- list()
    v <- 1

    for(k in 1:v){
      s <- sample(1:nrow(X_tmp))
      X_tmp[, i] <- X_tmp[s, i]
      res <- riemannian_pls(X = X_tmp, Y = Y, L = K, tol = 1e-4, max.iter = 20,
                            type_x = "euclidean", type_y = "euclidean", method = "tangent",
                            scale = TRUE)

      P.hat <- do.call(cbind, res$loadingsX)
      W.hat <- do.call(cbind, res$weightsX)
      C.hat <- do.call(cbind, res$weightsY)
      B <- map_dbl(res$reg_steps, ~.x$b1) %>%
        diag()

      beta_list[[k]] <- (W.hat %*% solve(t(P.hat) %*% W.hat) %*% B %*% t(C.hat))[i,]
      beta_list[[k]] <- as.matrix(beta_list[[k]])
      colnames(beta_list[[k]]) <- colnames(X_tmp)[i]
    }
    return(beta_list)
  },
  mc.cores = cores
)

map(res, ~do.call(cbind, .x))

write_rds(res, "results/cobre_permute.Rds")

