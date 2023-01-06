## Matt Ryan
## 03/11/2022
# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)


# load functions ----------------------------------------------------------

source("code/get_data.R")
source("code/create_tangent_vs_riemann.R")
source("code/create_tangent_n.R")

tangent_df <- get_data("tangent", 0.05)
n_df <- get_data("sample_size", 0.05)
n_df <- n_df %>%
  mutate(method = "Riemannian")
df <- bind_rows(tangent_df, n_df) %>%
  mutate(method = str_to_title(method))


# parameters ----------------------------------------------------------------

vars <- c("P_r2", "Q_r2", "X_r2", "Y_r2", "X_rmse", "Y_rmse", "tp_rmse", "uq_rmse")

walk(
  vars,
  function(v){
    create_tangent_vs_riemann(df = df,
                              metric = v,
                              x_type = "affine",
                              y_type = "affine",
                              save = TRUE)
    create_tangent_vs_riemann(df = df,
                              metric = v,
                              x_type = "euclidean",
                              y_type = "affine",
                              save = FALSE, save_git = TRUE)
    create_tangent_vs_riemann(df = df,
                              metric = v,
                              x_type = "affine",
                              y_type = "euclidean",
                              save = FALSE, save_git = TRUE)
  }
)

walk(
  vars,
  function(v){
    create_tangent_n(df = df,
                     metric = v,
                     x_type = "affine",
                     y_type = "affine",
                     save = TRUE)
    create_tangent_n(df = df,
                     metric = v,
                     x_type = "euclidean",
                     y_type = "affine",
                     save = FALSE, save_git = TRUE)
    create_tangent_n(df = df,
                     metric = v,
                     x_type = "affine",
                     y_type = "euclidean",
                     save = FALSE, save_git = TRUE)
  }
)


