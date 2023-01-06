## Matt Ryan
## 03/11/2022
# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)


# load functions ----------------------------------------------------------

source("code/get_data.R")
source("code/create_keyfinding3.R")


# parameters ----------------------------------------------------------------

df <- get_data("sample_size", sig.filter = 0.05) %>%
  filter(type_x == type_y, K == 5)

vars <- c("P_r2", "Q_r2", "X_r2", "Y_r2", "X_rmse", "Y_rmse", "tp_rmse", "uq_rmse")

walk(
  vars,
  function(v){
    create_keyfinding3(df = df, metric = v, save = TRUE)
  }
)
