## Matt Ryan
## 30/08/2022
# packages ----------------------------------------------------------------
pacman::p_load(tidyverse)
source("code/create_facet_plots.R")

# data and parameters -----------------------------------------------------
files <- list.files(here::here("results"), full.names = T, recursive = T)
files <- files[str_detect(files, "summary")]
df <- map_dfr(files, read_rds)

N <- c(10, 50, 100, 150)

r2_stats <- c("P", "Q", "X", "Y")

r2_grid <- expand_grid(N, r2_stats)

rmse_stats <- c("X", "Y")

rmse_grid <- expand_grid(N, rmse_stats)


# Make r2 plots -----------------------------------------------------------

walk2(r2_grid$N, r2_grid$r2_stats, create_facet_plot_r2)

# Make rmse plots -----------------------------------------------------------

walk2(rmse_grid$N, rmse_grid$rmse_stats, create_facet_plot_rmse)

