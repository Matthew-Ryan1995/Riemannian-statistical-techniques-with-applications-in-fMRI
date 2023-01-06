## Matt Ryan
## 01/11/2022
# packages ----------------------------------------------------------------
pacman::p_load(tidyverse, patchwork)
source("code/unchoose.R")

# data and parameters -----------------------------------------------------
files <- list.files(here::here("results"), full.names = T, recursive = T)
files <- files[str_detect(files, "summary")]
df <- map_dfr(files, read_rds)

height <- 5
x_dim <- 15
y_dim <- 15
sigma <- 0
N <- 150
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
  filter(!(metric %in% c("X: RMSE", "Y: RMSE")))



# create R2 plots ---------------------------------------------------------

plot_df <- tmp_df %>%
  pivot_longer(-(sig:L), names_to = "metric") %>%
  mutate(metric = case_when(
    str_detect(metric, "P_r2") ~ "P: R2",
    str_detect(metric, "Q_r2") ~ "Q: R2",
    str_detect(metric, "X_r2") ~ "X: R2",
    str_detect(metric, "Y_r2") ~ "Y: R2",
    str_detect(metric, "X_rmse") ~ "X: RMSE",
    str_detect(metric, "Y_rmse") ~ "Y: RMSE",
    str_detect(metric, "TP_rmse") ~ "TP: RMSE",
    str_detect(metric, "UQ_rmse") ~ "UQ: RMSE"
  ),
  metric = factor(metric,
                  levels = c("P: R2", "Q: R2", "X: R2", "Y: R2",
                             "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE")),
  K = as.factor(K),
  n = as.factor(n)
  ) %>%
  group_by(metric, n, K) %>%
  summarise(
    median = median(value),
    Q1 = quantile(value, 1/4),
    Q3 = quantile(value, 3/4),
    .groups = "drop"
  ) %>%
  filter(!(str_detect(metric, "RMSE"))) %>%
  filter(!str_detect(metric, "X|Y")) %>%
  mutate(
    gp2 = ifelse(str_detect(metric, "P|Q"),
                 "Loading",
                 "Data"),
    gp = ifelse(str_detect(metric, "P|X"),
                "P",
                "Q"),
    label = case_when(
      str_detect(metric, "P") ~ "P",
      str_detect(metric, "Q") ~ "Q",
      str_detect(metric, "X") ~ "X",
      str_detect(metric, "Y") ~ "Y",
    ),
    gp = as.factor(gp))

p <- plot_df %>%
  ggplot(aes(x = K, y = median,
             colour = gp,
             group = metric,
             # shape = gp2,
             label = label
  )) +
  geom_hline(data = y_int,
             aes(yintercept = val),
             lty = 2, colour = "gray90") +
  # geom_linerange(aes(ymin = Q1, ymax = Q3),
  #                position = position_dodge(width = 0.5),
  #                show.legend = FALSE
  # ) +
  geom_point(position = position_dodge(width = 0.5), size = 4, shape = 17) +
  harrypotter::scale_color_hp_d("ravenclaw") +
  scale_shape(guide = "none") +
  labs(y = latex2exp::TeX("$R^2_S$"),
       colour = NULL,
       shape = NULL,
       title = str_c("Euclidean NIPALS\n",
                     "| p = ", x_dim,
                     " | q = ", y_dim,
                     " | sigma = ", sigma,
                     " | n = ", N, " |"
       )) +
  theme_classic() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "gray70"),
        plot.title = element_text(hjust = 0.5))

ggsave("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding1_r2.png",
       plot = p, height = height, width = 1.618 * height)

# create RMSE plots ---------------------------------------------------------

plot_df <- tmp_df %>%
  pivot_longer(-(sig:L), names_to = "metric") %>%
  mutate(metric = case_when(
    str_detect(metric, "P_r2") ~ "P: R2",
    str_detect(metric, "Q_r2") ~ "Q: R2",
    str_detect(metric, "X_r2") ~ "X: R2",
    str_detect(metric, "Y_r2") ~ "Y: R2",
    str_detect(metric, "X_rmse") ~ "X: RMSE",
    str_detect(metric, "Y_rmse") ~ "Y: RMSE",
    str_detect(metric, "TP_rmse") ~ "TP: RMSE",
    str_detect(metric, "UQ_rmse") ~ "UQ: RMSE"
  ),
  metric = factor(metric,
                  levels = c("P: R2", "Q: R2", "X: R2", "Y: R2",
                             "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE")),
  K = as.factor(K),
  n = as.factor(n)
  ) %>%
  group_by(metric, n, K) %>%
  summarise(
    median = median(value),
    Q1 = quantile(value, 1/4),
    Q3 = quantile(value, 3/4),
    .groups = "drop"
  ) %>%
  filter(str_detect(metric, "RMSE")) %>%
  filter(!str_detect(metric, "X|Y")) %>%
  mutate(
    gp2 = ifelse(str_detect(metric, "P|Q"),
                 "Loading",
                 "Data"),
    gp = ifelse(str_detect(metric, "P|X"),
                "TP",
                "UQ"),
    label = case_when(
      str_detect(metric, "P") ~ "P",
      str_detect(metric, "Q") ~ "Q",
      str_detect(metric, "X") ~ "X",
      str_detect(metric, "Y") ~ "Y",
      str_detect(metric, "TP") ~ "TP",
      str_detect(metric, "UQ") ~ "UQ",
    ),
    gp = as.factor(gp))

p <- plot_df %>%
  ggplot(aes(x = K, y = median,
             colour = gp,
             group = metric,
             shape = gp2,
             label = label
  )) +
  # geom_linerange(aes(ymin = Q1, ymax = Q3),
  #                position = position_dodge(width = 0.5),
  #                show.legend = FALSE
  # ) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  harrypotter::scale_color_hp_d("ravenclaw") +
  scale_shape(guide = "none") +
  labs(y = "RMSE",
       colour = NULL,
       shape = NULL,
       title = str_c("Euclidean NIPALS\n",
                     "| p = ", x_dim,
                     " | q = ", y_dim,
                     " | sigma = ", sigma,
                     " | n = ", N, " |"
       )
  ) +
  theme_classic() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "gray70"),
        plot.title = element_text(hjust = 0.5))

ggsave("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding1_rmse.png",
       plot = p, height = height, width = 1.618 * height)
