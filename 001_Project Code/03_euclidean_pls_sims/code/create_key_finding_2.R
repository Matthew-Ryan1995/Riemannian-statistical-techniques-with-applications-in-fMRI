create_finding_2_plots <- function(df, var, height = 5){
  plot_df <- df %>%
    filter(n == 80, K == 5) %>%
    select(sig, p, q, vars = all_of(var))

  if(str_detect(var, "P|X")){
    plot_df <- plot_df %>%
      filter(q == 55) %>%
      group_by(sig, p) %>%
      summarise(median = median(vars),
                .groups = "drop") %>%
      mutate(dim = map_dbl(p, unchoose), dim = as.factor(dim))
    col_label <- "p:"
    title_lab <- "| q = "
  }else{
    plot_df <- plot_df %>%
      filter(p == 55) %>%
      group_by(sig, q) %>%
      summarise(median = median(vars),
                .groups = "drop") %>%
      mutate(dim = map_dbl(q, unchoose), dim = as.factor(dim))
    col_label <- "q:"
    title_lab <- "| p = "
  }

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
    p <- plot_df %>%
      ggplot(aes(x = sig, y = median, col = dim, shape = dim)) +
      geom_hline(data = y_int,
                 aes(yintercept = val),
                 lty = 2, colour = "gray90") +
      geom_point(size = 3)
  }else{
    p <- plot_df %>%
      ggplot(aes(x = sig, y = median, col = dim, shape = dim)) +
      geom_point(size = 3)
  }

  p <- p +
    harrypotter::scale_color_hp_d("ravenclaw") +
    theme_classic() +
    labs(x= latex2exp::TeX("$\\sigma$"),
         y = y_label,
         colour = col_label,
         shape = col_label,
         title = str_c("Euclidean NIPALS\n",
                       title_lab, unchoose(55),
                       " | n = ", 80,
                       " | K = ", 5, " |"
         )
         ) +
    theme(legend.position = "bottom",
          text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE))

  ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding2_{var}.png"),
         plot = p,
         height = height, width = 1.618 * height)
}
