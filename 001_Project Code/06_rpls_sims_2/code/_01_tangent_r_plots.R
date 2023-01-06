## Matt Ryan
## 03/11/2022
# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)

source("code/tangent_r_plot.R")
# load functions ----------------------------------------------------------

files <- list.files(here::here("results/tangent/"), full.names = T, recursive = T)
files <- files[!str_detect(files, "presim")]
df <- map(files, safely(read_rds))
df <- map_dfr(df, ~.x$result) %>%
  filter(K <= 5) %>%
  mutate(method = ifelse(method == "full", "Riemannian", method),
         method = str_to_title(method))

vars <- c("P_r2", "Q_r2", "X_r2", "Y_r2", "X_rmse", "Y_rmse", "uq_rmse", "tp_rmse")

walk(
  vars,
  function(v){
    tangent_r_plot(df = df,
                   metric = v,
                   x_type = "affine",
                   y_type = "affine",
                   save = TRUE)
    tangent_r_plot(df = df,
                   metric = v,
                   x_type = "euclidean",
                   y_type = "affine",
                   save = FALSE, save_git = TRUE)
    tangent_r_plot(df = df,
                   metric = v,
                   x_type = "affine",
                   y_type = "euclidean",
                   save = FALSE, save_git = TRUE)
  }
)
