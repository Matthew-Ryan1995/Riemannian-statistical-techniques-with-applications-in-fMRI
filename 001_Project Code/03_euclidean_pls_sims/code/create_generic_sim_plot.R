create_generic_sim_results_euclidean <- function(df, x_dim,
                                                 y_dim, sigma, height = 7,
                                                 save = FALSE){
  tmp_df <- df %>%
    filter(sig == sigma,
           p == choose(x_dim + 1, 2),
           q == choose(y_dim + 1, 2))

  # return(tmp_df)
  y_int <- expand_grid(
    val = c(0.2, 0.4, 0.6, 0.8, 1),
    metric = factor(c("P: R2", "Q: R2", "X: R2", "Y: R2",
                      "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE"),
                    levels = c("P: R2", "Q: R2", "X: R2", "Y: R2",
                               "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE"))
  ) %>%
    filter(!(str_detect(metric, "RMSE")))

  p <- tmp_df %>%
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
    ggplot(aes(x = K, y = median, colour = n)) +
    geom_hline(data = y_int,
               aes(yintercept = val),
               lty = 2, colour = "gray90") +
    geom_linerange(aes(ymin = Q1, ymax = Q3),
                   position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5)) +
    facet_wrap(~metric, scale = "free_y") +
    harrypotter::scale_color_hp_d("ravenclaw") +
    labs(y = NULL, colour = "Sample size:",
         title = str_c("Euclidean PLS\n",
                       "| p = ", choose(x_dim + 1 , 2),
                       " | q = ", choose(y_dim + 1, 2),
                       " | sigma = ", sigma, " |"
         )) +
    theme_classic() +
    theme(legend.position = "bottom",
          text = element_text(size = 14),
          strip.text = element_text(size = 20),
          strip.background = element_rect(fill = "gray70"),
          plot.title = element_text(hjust = 0.5))

  if(save){
    ggsave(glue::glue("~/Desktop/Riemannian-statistical-techniques-with-applications-in-fMRI/Model Recovery/Euclidean PLS/fig/euclidean_p_{x_dim}_q_{y_dim}_sig_{sigma}.png"),
           height = height, width = 1.618 * height, plot = p)
  }else{
    return(p)
  }
}
