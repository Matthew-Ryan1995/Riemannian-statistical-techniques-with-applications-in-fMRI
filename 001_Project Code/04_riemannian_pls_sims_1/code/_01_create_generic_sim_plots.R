## Matt Ryan
## 03/11/2022
# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)


# load functions ----------------------------------------------------------

source("code/get_data.R")
source("code/create_generic_sim_plot.R")


# parameters ----------------------------------------------------------------

folder_and_var_names <- tibble(
  folder = c("noise", "p", "q", "sample_size"),
  var = c("sig", "x_dim", "y_dim", "n"),
  filllab = c("sigma", "p", "q", "n"),
  sig.filter = c(1, rep(0.05, 3))
)

goodcopy <- FALSE

walk(
  1:nrow(folder_and_var_names),
  function(i){
    params <- folder_and_var_names %>%
      slice(i)
    df <- get_data(params$folder, params$sig.filter)

    create_generic_sim_plot(df = df,
                            vars = params$var,
                            x_type = "affine",
                            y_type = "affine",
                            save = TRUE)
    create_generic_sim_plot(df = df,
                            vars = params$var,
                            x_type = "euclidean",
                            y_type = "affine",
                            save = TRUE)
    create_generic_sim_plot(df = df,
                            vars = params$var,
                            x_type = "affine",
                            y_type = "euclidean",
                            save = TRUE)
  }
)
