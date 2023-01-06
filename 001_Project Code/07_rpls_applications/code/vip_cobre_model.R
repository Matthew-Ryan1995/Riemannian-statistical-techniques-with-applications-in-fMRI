# tictoc::tic()
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

job_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print(str_c("Doing job COBRE"))

source("code/calculate_vip.R")

# Data --------------------------------------------------------------------
M <- read_rds("results/cobre_rpls_model_K_3.Rds")

Y <- M$Y[[1]]
X <- M$X[[1]]
K <- M$L


# True VIP ----------------------------------------------------------------
Tmat <- do.call(cbind, M$scoresX)
W <- do.call(cbind, M$weightsX)
C <- do.call(cbind, M$weightsY)

true_vip <- calculate_vip(Tmat = Tmat, W = W, Y = Y)
true_vip_ij <- calculate_vip_per_y(Tmat = Tmat, W = W, Y = Y, C = C)

set.seed(1668286 + job_num)

j <- ncol(X)

res <- parallel::mclapply(
  1:j,
  function(i){
    X_tmp <- X

    s <- sample(1:nrow(X_tmp))
    X_tmp[, i] <- X_tmp[s, i]
    mod <- riemannian_pls(X = X_tmp, Y = Y, L = K, tol = 1e-4, max.iter = 20,
                          type_x = "euclidean", type_y = "euclidean", method = "tangent",
                          scale = TRUE)

    Tmat <- do.call(cbind, mod$scoresX)
    W <- do.call(cbind, mod$weightsX)
    C <- do.call(cbind, mod$weightsY)

    vip_score <- calculate_vip(Tmat = Tmat, W = W, Y = Y)
    vip_score_ij <- calculate_vip_per_y(Tmat = Tmat, W = W, Y = Y, C = C)

    ans <- cbind(vip_score[i] > true_vip[i], t(as.matrix(vip_score_ij[i, ] > true_vip_ij[i, ])))

    colnames(ans) <- c("full_vip", colnames(Y))

    rownames(ans)[1] <- colnames(X)[i]

    return(ans)
  },
  mc.cores = cores
)

res <- do.call(rbind, res)
# tictoc::toc()
write_rds(res, glue::glue("results/cobre_vip/cobre_vip_sim_{job_num}.Rds"))

