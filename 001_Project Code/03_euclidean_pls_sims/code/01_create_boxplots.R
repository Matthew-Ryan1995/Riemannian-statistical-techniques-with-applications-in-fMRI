## Matt Ryan
## 30/08/2022
# packages ----------------------------------------------------------------
pacman::p_load(tidyverse, patchwork)
source("code/create_boxplots.R")

# data and parameters -----------------------------------------------------
files <- list.files(here::here("results"), full.names = T, recursive = T)
files <- files[str_detect(files, "summary")]
df <- map_dfr(files, read_rds)


vars <- c(str_c(c("P", "Q", "X", "Y"), "_r2"), str_c(c("X", "Y"), "_rmse"))
sig <- c(0, 0.05, 0.1, 0.5, 1)

params <- expand_grid(vars = vars,
                      s = sig) %>%
  mutate(title = map2(vars, s, function(vars, s){
    case_when(
      str_detect(vars, "P") ~ bquote(atop("P loading space recovery:"~R^2,
                                          sigma ==.(s))),
      str_detect(vars, "Q") ~ bquote(atop("Q loading space recovery:"~R^2,
                                          sigma ==.(s))),
      str_detect(vars, "X_r2") ~ bquote(atop("X recovery:"~R^2,
                                             sigma ==.(s))),
      str_detect(vars, "Y_r2") ~ bquote(atop("Y recovery:"~R^2,
                                             sigma ==.(s))),
      str_detect(vars, "X") ~ bquote(atop("X recovery: RMSE",
                                          sigma ==.(s))),
      str_detect(vars, "Y") ~ bquote(atop("Y recovery: RMSE",
                                          sigma ==.(s)))
    )
  })
  )

goodcopy <- FALSE

walk(1:nrow(params),
     function(i){
       s <- params$s[i]
       var <- params$vars[i]
       title <- params$title[i][[1]]

       create_boxplots(df = df, s = s, vars = var, title = title, height = 5,
                       savename = str_c(var, s, sep = "_"), goodcopy = goodcopy)
     })


walk(vars,
     function(v){
       title <- case_when(
         str_detect(v, "P") ~ expression("P loading space recovery:"~R^2),
         str_detect(v, "Q") ~ expression("Q loading space recovery:"~R^2),
         str_detect(v, "X_r2") ~ expression("X recovery:"~R^2),
         str_detect(v, "Y_r2") ~ expression("Y recovery:"~R^2),
         str_detect(v, "X") ~ expression("X recovery: RMSE"),
         str_detect(v, "Y") ~ expression("Y recovery: RMSE")
       )
       create_n_sig(df = df, savename = str_c(v, "n_sig", sep = "_"), goodcopy = goodcopy,
                    title = title, height = 5, vars = v)
     })

make_plots(df = df, var = "n", filllab = "Sample size", sig_filter = 0.5,
           savename = "n_sig05", goodcopy = goodcopy)
