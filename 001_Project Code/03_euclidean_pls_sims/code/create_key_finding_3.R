
create_finding_3_plots <- function(df, var, height = 5){
  plot_df <- df %>%
    filter(n == 80, sig == 0.5) %>%
    select(K, p, q, vars = all_of(var))

  y_label <- case_when(
    str_detect(var, "P_r2") ~ "$R^2_S$ (P)",
    str_detect(var, "Q_r2") ~ "$R^2_S$ (Q)",
    str_detect(var, "X_r2") ~ "$R^2$ (X)",
    str_detect(var, "Y_r2") ~ "$R^2$ (Y)",
    str_detect(var, "X_rmse") ~ "RMSE (X)",
    str_detect(var, "Y_rmse") ~ "RMSE (Y)",
    str_detect(var, "TP_rmse") ~ "RMSE (TP)",
    str_detect(var, "UQ_rmse") ~ "RMSE (UQ)"
  ) %>%
    latex2exp::TeX()

  if(str_detect(var, "r2")){
    p1 <- plot_df %>%
      filter(q == 55) %>%
      mutate(p = map_dbl(p, unchoose)) %>%
      group_by(K, p) %>%
      summarise(median = median(vars),
                Q1 = quantile(vars, 1/4),
                Q3 = quantile(vars, 3/4),
                .groups = "drop") %>%
      ggplot(aes(x = K, y = median, col = as.factor(p), shape = as.factor(p))) +
      geom_hline(data = y_int,
                 aes(yintercept = val),
                 lty = 2, colour = "gray90") +
      # geom_linerange(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.5)) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }else{
    p1 <- plot_df %>%
      filter(q == 55) %>%
      mutate(p = map_dbl(p, unchoose)) %>%
      group_by(K, p) %>%
      summarise(median = median(vars),
                Q1 = quantile(vars, 1/4),
                Q3 = quantile(vars, 3/4),
                .groups = "drop") %>%
      ggplot(aes(x = K, y = median, col = as.factor(p), shape = as.factor(p))) +
      geom_linerange(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.5)) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }

  p1 <- p1 +
    harrypotter::scale_color_hp_d("ravenclaw") +
    theme_classic() +
    labs(x = "K", y = y_label, colour = "p:", shape = "p:",
         title = str_c("Euclidean PLS\n",
                       "| q = ", unchoose(55),
                       " | sigma = ", 0.5,
                       " | n = ", 80, " |"
         )
    ) +
    theme(legend.position = "bottom",
          text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE))



  ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding3_{var}_pvary.png"),
         plot = p1,
         height = height, width = 1.618 * height)

  if(str_detect(var, "r2")){
    p2 <- plot_df %>%
      filter(p == 55) %>%
      mutate(q = map_dbl(q, unchoose)) %>%
      group_by(K, q) %>%
      summarise(median = median(vars),
                Q1 = quantile(vars, 1/4),
                Q3 = quantile(vars, 3/4),
                .groups = "drop") %>%
      ggplot(aes(x = K, y = median, col = as.factor(q), shape = as.factor(q))) +
      geom_hline(data = y_int,
                 aes(yintercept = val),
                 lty = 2, colour = "gray90") +
      # geom_linerange(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.5)) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }else{
    p2 <- plot_df %>%
      filter(p == 55) %>%
      mutate(q = map_dbl(q, unchoose)) %>%
      group_by(K, q) %>%
      summarise(median = median(vars),
                Q1 = quantile(vars, 1/4),
                Q3 = quantile(vars, 3/4),
                .groups = "drop") %>%
      ggplot(aes(x = K, y = median, col = as.factor(q), shape = as.factor(q))) +
      geom_linerange(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.5)) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }

  p2 <- p2 +
    harrypotter::scale_color_hp_d("ravenclaw") +
    theme_classic() +
    labs(x = "K", y = y_label, colour = "q:", shape = "q:",
         title = str_c("Euclidean NIPALS\n",
                       "| p = ", unchoose(55),
                       " | sigma = ", 0.5,
                       " | n = ", 80, " |"
         )
    ) +
    theme(legend.position = "bottom",
          text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE))


  ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding3_{var}_qvary.png"),
         plot = p2,
         height = height, width = 1.618 * height)
}
