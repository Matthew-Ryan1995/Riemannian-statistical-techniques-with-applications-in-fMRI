## Matthew Ryan
## 01/12/2022
# packages ----------------------------------------------------------------

pacman::p_load(tidyverse)

source("code/create_vip_plot.R")
source("code/calculate_vip.R")

# COBRE --------------------------------------------------------------------
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

true_vip <- calculate_vip_per_y(Tmat = Tmat, W = W.hat, Y = Y, C = C.hat)
single_vip <- calculate_vip(Tmat = Tmat, W = W.hat, Y = Y)


# Significant values ------------------------------------------------------
fl <- list.files("results/cobre_vip/", full.names = T)
sig <- map(fl, read_rds)
sig <- (Reduce("+", sig)/200)
sig[1:39, ] <- 1
sig <- sig %>%
  apply(2, p.adjust, method = "fdr")

p <- create_vip_plot(vip = true_vip, sig = sig[, -1], beta = beta)

ggsave("../../002_Thesis/01_initial_draft/img/cobre_vip_plot.png", plot = p,
       height = height, width = 1.618*height)

p2 <- create_single_vip_plot(vip = single_vip, sig = sig[, 1], beta = beta)

ggsave("../../002_Thesis/01_initial_draft/img/cobre_normal_vip_plot.png", plot = p2,
       height = height, width = 1.618*height)

# ABIDE -------------------------------------------------------------------
dataset <- "abide"
height <- 7
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

p <- create_vip_plot(vip = true_vip, sig = sig[, -1], beta = beta)

ggsave("../../002_Thesis/01_initial_draft/img/abide_vip_plot.png", plot = p,
       height = height, width = 1.618*height)

p2 <- create_single_vip_plot(vip = single_vip, sig = sig[, 1], beta = beta)

ggsave("../../002_Thesis/01_initial_draft/img/abide_normal_vip_plot.png", plot = p2,
       height = height, width = 1.618*height)

