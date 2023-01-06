## Matt Ryan
## 03/11/2022
# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)


# load functions ----------------------------------------------------------

source("code/get_data.R")
source("code/create_noise_finding_plots.R")


# parameters ----------------------------------------------------------------

df <- get_data("noise", sig.filter = 1)

vars <- c("P_r2", "Q_r2", "X_r2", "Y_r2", "X_rmse", "Y_rmse", "tp_rmse", "uq_rmse")

walk(
  vars,
  function(v){
    create_noise_finding_plots(df = df,
                               vars = v,
                               x_type = "affine",
                               y_type = "affine",
                               save = TRUE)
    # create_noise_finding_plots(df = df,
    #                            vars = v,
    #                            x_type = "euclidean",
    #                            y_type = "affine",
    #                            save = TRUE)
    # create_noise_finding_plots(df = df,
    #                            vars = v,
    #                            x_type = "affine",
    #                            y_type = "euclidean",
    #                            save = TRUE)
  }
)
