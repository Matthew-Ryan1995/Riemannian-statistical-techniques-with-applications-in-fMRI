## Matt Ryan
## 01/11/2022
# packages ----------------------------------------------------------------
pacman::p_load(tidyverse, patchwork)

# data and parameters -----------------------------------------------------
files <- list.files(here::here("results"), full.names = T, recursive = T)
files <- files[str_detect(files, "summary")]
df <- map_dfr(files, read_rds)

height <- 5
x_dim <- 10
y_dim <- 10
sigma <- 0.5
N <- 80
tmp_df <- df %>%
  filter(sig == sigma,
         p == choose(x_dim + 1, 2),
         q == choose(y_dim + 1, 2),
         n == N
  )

y_int <- expand_grid(
  val = c(0.2, 0.4, 0.6, 0.8, 1),
  metric = factor(c("P: R2", "Q: R2", "X: R2", "Y: R2",
                    "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE"),
                  levels = c("P: R2", "Q: R2", "X: R2", "Y: R2",
                             "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE"))
) %>%
  filter(!(str_detect(metric, "RMSE")))




# R2 ----------------------------------------------------------------------

p <- tmp_df %>%
  filter(K == 5) %>%
  select(X_r2, Y_r2) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x=  NULL, y = latex2exp::TeX("$R^2$"), fill = NULL,
       title = str_c("Euclidean NIPALS\n",
                     "| p = ", x_dim,
                     " | q = ", y_dim,
                     " | sigma = ", sigma,
                     " | n = ", N,
                     " | K = 5 |"
       )
  ) +
  scale_x_discrete(labels = c("X", "Y")) +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

ggsave("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding5_r2_K5.png",
       plot = p, height = height, width = 1.618 * height)

p <- tmp_df %>%
  select(K, X_r2, Y_r2) %>%
  pivot_longer(-K) %>%
  ggplot(aes(x = as.factor(K), y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x = "K", y = latex2exp::TeX("$R^2$"), fill = NULL,
       title = str_c("Euclidean NIPALS\n",
                     "| p = ", x_dim,
                     " | q = ", y_dim,
                     " | sigma = ", sigma,
                     " | n = ", N, " |"
       )
  ) +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

ggsave("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding5_r2_all.png",
       plot = p, height = height, width = 1.618 * height)

# RMSE --------------------------------------------------------------------


p <- tmp_df %>%
  filter(K == 5) %>%
  select(X_rmse, Y_rmse) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x=  NULL, y = "RMSE", fill = NULL,
       title = str_c("Euclidean NIPALS\n",
                     "| p = ", x_dim,
                     " | q = ", y_dim,
                     " | sigma = ", sigma,
                     " | n = ", N,
                     " | K = 5 |"
       )
  ) +
  scale_x_discrete(labels = c("X", "Y")) +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

ggsave("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding5_rmse_K5.png",
       plot = p, height = height, width = 1.618 * height)

p <- tmp_df %>%
  select(K, X_rmse, Y_rmse) %>%
  pivot_longer(-K) %>%
  ggplot(aes(x = as.factor(K), y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x = "K", y = "RMSE", fill = NULL,
       title = str_c("Euclidean NIPALS\n",
                     "| p = ", x_dim,
                     " | q = ", y_dim,
                     " | sigma = ", sigma,
                     " | n = ", N, " |"
       )
  ) +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

ggsave("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding5_rmse_all.png",
       plot = p, height = height, width = 1.618 * height)

# RMSE 2 --------------------------------------------------------------------


p <- tmp_df %>%
  filter(K == 5) %>%
  select(TP_rmse, UQ_rmse) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x=  NULL, y = "RMSE", fill = NULL,
       title = str_c("Euclidean NIPALS\n",
                     "| p = ", x_dim,
                     " | q = ", y_dim,
                     " | sigma = ", sigma,
                     " | n = ", N,
                     " | K = 5 |"
       )
  ) +
  scale_x_discrete(labels = c("TP", "UQ")) +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

ggsave("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding5_tpuq_rmse_K5.png",
       plot = p, height = height, width = 1.618 * height)

p <- tmp_df %>%
  select(K, TP_rmse, UQ_rmse) %>%
  pivot_longer(-K) %>%
  mutate(name = ifelse(str_detect(name, "TP"), "TP", "UQ")) %>%
  ggplot(aes(x = as.factor(K), y = value, fill = name)) +
  geom_boxplot(show.legend = T) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  labs(x = "K", y = "RMSE", fill = "Data type",
       title = str_c("Euclidean NIPALS\n",
                     "| p = ", x_dim,
                     " | q = ", y_dim,
                     " | sigma = ", sigma,
                     " | n = ", N, " |"
       )
  ) +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

ggsave("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding5_tpuq_rmse_all.png",
       plot = p, height = height, width = 1.618 * height)
