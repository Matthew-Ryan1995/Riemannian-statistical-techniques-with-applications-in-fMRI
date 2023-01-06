## Matt Ryan
## 01/11/2022
# packages ----------------------------------------------------------------
pacman::p_load(tidyverse, patchwork)
source("code/create_key_finding_4.R")
source("code/unchoose.R")

# data and parameters -----------------------------------------------------
files <- list.files(here::here("results"), full.names = T, recursive = T)
files <- files[str_detect(files, "summary")]
df <- map_dfr(files, read_rds)

height <- 5

y_int <- expand_grid(
  val = c(0.2, 0.4, 0.6, 0.8, 1),
  metric = factor(c("P: R2", "Q: R2", "X: R2", "Y: R2",
                    "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE"),
                  levels = c("P: R2", "Q: R2", "X: R2", "Y: R2",
                             "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE"))
) %>%
  filter(!(str_detect(metric, "RMSE")))

vars <- c("P_r2", "Q_r2", "X_r2", "Y_r2", "X_rmse", "Y_rmse", "TP_rmse", "UQ_rmse")

walk(vars, create_finding_4_plots, df = df, height = height)
walk(vars, create_finding_4_plots_sub, df = df, height = height)
