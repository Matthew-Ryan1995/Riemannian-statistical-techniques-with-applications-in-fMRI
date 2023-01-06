## Matthew Ryan
## 18/11/2022
# packages ----------------------------------------------------------------

pacman::p_load(tidyverse, phdWork)

source("code/plotting_functions.R")
source("code/calculate_vip.R")

# Data --------------------------------------------------------------------
dataset <- "cobre"
height <- 7
M <- read_rds("results/cobre_rpls_model_K_3.Rds")
regions <- read_csv("data/msdl_rois_labels.csv",
                    col_types = cols())

Tmat <- do.call(cbind, M$scoresX)
P.hat <- do.call(cbind, M$loadingsX)
W.hat <- do.call(cbind, M$weightsX)
C.hat <- do.call(cbind, M$weightsY)
B <- map_dbl(M$reg_steps, ~.x$b1) %>%
  diag()
Y <- M$Y[[1]]

beta <- W.hat %*% solve(t(P.hat) %*% W.hat) %*% B %*% t(C.hat)

true_vip <- calculate_vip_per_y(Tmat = Tmat, W = W.hat,  Y = Y, C = C.hat)

# Significant values ------------------------------------------------------
fl <- list.files("results/cobre_vip/", full.names = T)
sig <- map(fl, read_rds)
sig <- (Reduce("+", sig)/200)
sig[1:39, ] <- 1
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
       vip = true_vip,
       save_folder = save_folder,
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
       B = beta, dataset = dataset, height = height, filter = FALSE, sig = sig)
  walk(1:ncol(beta),
       create_regression_coeff_plots,
       B = beta, dataset = dataset, height = height, filter = TRUE, sig = sig)

}
