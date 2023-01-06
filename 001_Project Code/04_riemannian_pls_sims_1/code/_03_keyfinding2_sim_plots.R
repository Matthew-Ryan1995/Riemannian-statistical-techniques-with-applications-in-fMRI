## Matt Ryan
## 03/11/2022
# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)


# load functions ----------------------------------------------------------

source("code/get_data.R")
source("code/create_keyfinding2a.R")


# parameters ----------------------------------------------------------------


folder_and_var_names <- tibble(
  folder = c( "p", "q", "sample_size"),
  var = c("x_dim", "y_dim", "n"),
  sig.filter = c(rep(0.05, 3))
)


vars <- c("P_r2", "Q_r2", "X_r2", "Y_r2", "X_rmse", "Y_rmse", "tp_rmse", "uq_rmse")

walk(
  1:nrow(folder_and_var_names),
  function(i){
    params <- folder_and_var_names %>%
      slice(i)
    df <- get_data(params$folder, params$sig.filter)

    for(v in vars){
      create_keyfinding2a(df = df,
                          vars = params$var,
                          metric = v,
                          x_type = "affine",
                          y_type = "affine",
                          save = TRUE)
    }
  }
)


# Comparing X and Y -------------------------------------------------------

df <- get_data("p", sig.filter = 0.05) %>%
  filter(x_dim == 10, type_x == type_y)

title <- str_c(
  "R-NIPALS\n",
  "Predictor: SPD; Response: SPD\n",
  "| p = 10 | q = 10 | sigma = 0.05 | n = 80 | "
)

### R2
###
height <- 5

p <- df %>%
  filter(K == 5) %>%
  select(X_r2, Y_r2) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x=  NULL, y = latex2exp::TeX("$R^2$"), fill = NULL,
       title = str_c(title, "K = 5 |")
  ) +
  scale_x_discrete(labels = c("X", "Y")) +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/riemannian_keyfinding2b_predictor_affine_response_affine_XY_r2_K5.png"),
       plot = p,
       height = height, width = 1.618 * height)

p <- df %>%
  filter(K <= 5) %>%
  mutate(K = as.factor(K)) %>%
  select(K, X_r2, Y_r2) %>%
  pivot_longer(-K) %>%
  ggplot(aes(x = K, y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x=  NULL, y = latex2exp::TeX("$R^2$"), fill = NULL,
       title = title
  ) +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/riemannian_keyfinding2b_predictor_affine_response_affine_XY_r2_Kall.png"),
       plot = p,
       height = height, width = 1.618 * height)

### RMSE
###

p <- df %>%
  filter(K == 5) %>%
  select(X_rmse, Y_rmse) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x=  NULL, y = latex2exp::TeX("RMSE"), fill = NULL,
       title = str_c(title, "K = 5 |")
  ) +
  scale_x_discrete(labels = c("X", "Y")) +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/riemannian_keyfinding2b_predictor_affine_response_affine_XY_rmse_K5.png"),
       plot = p,
       height = height, width = 1.618 * height)

p <- df %>%
  filter(K <= 5) %>%
  mutate(K = as.factor(K)) %>%
  select(K, X_rmse, Y_rmse) %>%
  pivot_longer(-K) %>%
  ggplot(aes(x = K, y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x=  NULL, y = latex2exp::TeX("RMSE"), fill = NULL,
       title = title
  ) +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/riemannian_keyfinding2b_predictor_affine_response_affine_XY_rmse_Kall.png"),
       plot = p,
       height = height, width = 1.618 * height)

### RMSE TP and UQ
###

p <- df %>%
  filter(K == 5) %>%
  select(tp_rmse, uq_rmse) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x=  NULL, y = latex2exp::TeX("RMSE"), fill = NULL,
       title = str_c(title, "K = 5 |")
  ) +
  scale_x_discrete(labels = c("TP", "UQ")) +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/riemannian_keyfinding2b_predictor_affine_response_affine_tpuq_rmse_K5.png"),
       plot = p,
       height = height, width = 1.618 * height)

p <- df %>%
  filter(K <= 5) %>%
  mutate(K = as.factor(K)) %>%
  select(K, tp_rmse, uq_rmse) %>%
  pivot_longer(-K) %>%
  mutate(name = ifelse(str_detect(name, "tp"), "TP", "UQ")) %>%
  ggplot(aes(x = K, y = value, fill = name)) +
  geom_boxplot(show.legend = TRUE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x=  NULL, y = latex2exp::TeX("RMSE"), fill = "Data type",
       title = title
  ) +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/riemannian_keyfinding2b_predictor_affine_response_affine_tpuq_rmse_Kall.png"),
       plot = p,
       height = height, width = 1.618 * height)
