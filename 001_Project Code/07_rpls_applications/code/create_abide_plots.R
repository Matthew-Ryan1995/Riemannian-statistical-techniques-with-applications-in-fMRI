## Matthew Ryan
## 18/11/2022
# packages ----------------------------------------------------------------

pacman::p_load(tidyverse, phdWork)

source("code/plotting_functions.R")
source("code/calculate_vip.R")

# Data --------------------------------------------------------------------
dataset <- "abide"
height <- 15
M <- read_rds("results/abide_rpls_model_K_6.Rds")
regions <- read_csv("data/aal_to_rsn.csv",
                    col_types = cols())

Tmat <- do.call(cbind, M$scoresX)
P.hat <- do.call(cbind, M$loadingsX)
W.hat <- do.call(cbind, M$weightsX)
C.hat <- do.call(cbind, M$weightsY)
B <- map_dbl(M$reg_steps, ~.x$b1) %>%
  diag()
Y <- M$Y[[1]]

beta <- W.hat %*% solve(t(P.hat) %*% W.hat) %*% B %*% t(C.hat)

true_vip <- calculate_vip_per_y(Tmat = Tmat, W = W.hat, Y = Y, C = C.hat)



fl <- list.files("results/abide_vip/", full.names = T)
sig <- map(fl, read_rds)
sig <- do.call(rbind, sig)
sig <- sig/200
sig[1:116, ] <- 1
sig <- sig %>%
  apply(2, p.adjust, method = "fdr")

colnames(true_vip) <- colnames(sig)[-1]


# Create plots ------------------------------------------------------------


goodcopy <- TRUE
if(goodcopy){
  save_folder <- "../../002_Thesis/01_initial_draft/img/"

  walk(1:length(M$loadingsX),
       create_loading_plots,
       dataset = dataset,
       height = height, filter = TRUE,
       save_folder = save_folder)

  walk(1:ncol(beta),
       create_regression_coeff_plots,
       B = beta, dataset = dataset,
       height = height, filter = TRUE,
       save_folder = save_folder,
       vip = true_vip,
       sig = sig)


}else{
  ## Loading plots
  walk(1:length(M$loadingsX),
       create_loading_plots,
       dataset = dataset, height = height, filter = FALSE)
  walk(1:length(M$loadingsX),
       create_loading_plots,
       dataset = dataset, height = height, filter = TRUE)

  ## Regression coefficient plots
  walk(1:ncol(beta),
       create_regression_coeff_plots,
       B = beta, dataset = dataset, height = height, filter = FALSE)
  walk(1:ncol(beta),
       create_regression_coeff_plots,
       B = beta, dataset = dataset, height = height, filter = TRUE)
}
