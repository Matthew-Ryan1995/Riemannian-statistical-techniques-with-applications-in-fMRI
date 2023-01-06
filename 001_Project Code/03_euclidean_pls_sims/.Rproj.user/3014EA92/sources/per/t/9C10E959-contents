create_finding_4_plots <- function(df, var, height = 5){
  plot_df <- df %>%
    filter(p == 55, q == 55, sig == 0.5) %>%
    select(K, n, vars = all_of(var)) %>%
    group_by(K, n) %>%
    summarise(median = median(vars),
              Q1 = quantile(vars, 1/4),
              Q3 = quantile(vars, 3/4),
              .groups = "drop")

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
    y_int <- tibble(val = c(0.2, 0.4, 0.6, 0.8, 1))
    p <- plot_df %>%
      ggplot(aes(x = K, y = median, col = as.factor(n), shape = as.factor(n))) +
      geom_hline(data = y_int,
                 aes(yintercept = val),
                 lty = 2, colour = "gray90") +
      geom_linerange(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.5),
                     show.legend = FALSE) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }else{
    p <- plot_df %>%
      ggplot(aes(x = K, y = median, col = as.factor(n), shape = as.factor(n))) +
      geom_linerange(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.5),
                     show.legend = FALSE) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }

  p <- p  +
    harrypotter::scale_color_hp_d("ravenclaw") +
    theme_classic() +
    labs(x = "K", y = y_label, colour = "Sample size:", shape = "Sample size:",
         title = str_c("Euclidean NIPALS\n",
                       "| p = ", unchoose(55),
                       " | q = ", unchoose(55),
                       " | sigma = ", 0.5, " |"
         )) +
    theme(legend.position = "bottom",
          text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE))


  ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding4_{var}.png"),
         plot = p,
         height = height, width = 1.618 * height)

}
create_finding_4_plots_sub <- function(df, var, height = 5){
  plot_df <- df %>%
    filter(p == 55, q == 55, sig == 0.05) %>%
    select(K, n, vars = all_of(var)) %>%
    group_by(K, n) %>%
    summarise(median = median(vars),
              Q1 = quantile(vars, 1/4),
              Q3 = quantile(vars, 3/4),
              .groups = "drop")

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
    y_int <- tibble(val = c(0.2, 0.4, 0.6, 0.8, 1))
    p <- plot_df %>%
      ggplot(aes(x = K, y = median, col = as.factor(n), shape = as.factor(n))) +
      geom_hline(data = y_int,
                 aes(yintercept = val),
                 lty = 2, colour = "gray90") +
      geom_linerange(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.5),
                     show.legend = FALSE) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }else{
    p <- plot_df %>%
      ggplot(aes(x = K, y = median, col = as.factor(n), shape = as.factor(n))) +
      geom_linerange(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.5),
                     show.legend = FALSE) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }

  p <- p  +
    harrypotter::scale_color_hp_d("ravenclaw") +
    theme_classic() +
    labs(x = "K", y = y_label, colour = "Sample size:", shape = "Sample size:",
         title = str_c("Euclidean NIPALS\n",
                       "| p = ", unchoose(55),
                       " | q = ", unchoose(55),
                       " | sigma = ", 0.05, " |"
         )) +
    theme(legend.position = "bottom",
          text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE))


  ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/euclidean_keyfinding4_{var}_appendix_version.png"),
         plot = p,
         height = height, width = 1.618 * height)

}
