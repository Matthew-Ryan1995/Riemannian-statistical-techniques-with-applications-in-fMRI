## Matt Ryan
## 31/10/2022
# packages ----------------------------------------------------------------
pacman::p_load(tidyverse, patchwork)
source("code/create_generic_sim_plot.R")

# data and parameters -----------------------------------------------------
files <- list.files(here::here("results"), full.names = T, recursive = T)
files <- files[str_detect(files, "summary")]
df <- map_dfr(files, read_rds)

x_dim <- y_dim <- round(seq(5, 15, length.out = 5))
sig <- unique(df$sig)

params <- expand_grid(
  p = x_dim,
  q = y_dim,
  sig = sig
)

walk(
  1:nrow(params),
  function(i){
    pars <- params %>%
      slice(i)
    create_generic_sim_results_euclidean(df = df,
                                         x_dim = pars$p,
                                         y_dim = pars$q,
                                         sigma = pars$sig,
                                         save = TRUE)
  }
)
