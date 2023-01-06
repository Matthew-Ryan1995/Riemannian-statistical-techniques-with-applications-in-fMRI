create_keyfinding3 <- function(df, metric, height = 5, save = FALSE){

  y_label <- case_when(
    str_detect(metric, "P") ~ "$R^2_S$ (P)",
    str_detect(metric, "Q") ~ "$R^2_S$ (Q)",
    str_detect(metric, "X_r2") ~ "$R^2$ (X)",
    str_detect(metric, "Y_r2") ~ "$R^2$ (Y)",
    str_detect(metric, "X_rmse") ~ "RMSE (X)",
    str_detect(metric, "Y_rmse") ~ "RMSE (Y)",
    str_detect(metric, "tp_rmse") ~ "RMSE (TP)",
    str_detect(metric, "uq_rmse") ~ "RMSE (UQ)"
  ) %>%
    latex2exp::TeX()


  title <- str_c(
    "R-NIPALS\n",
    "Predictor: SPD; Response: SPD\n",
    "| p = 10 | q = 10 | sigma = 0.05 | K = 5 |"
  )

  y_int <- expand_grid(
    val = c(0.2, 0.4, 0.6, 0.8, 1),
  )

  plot_data <- df %>%
    select(gp = n, K, var = all_of(metric)) %>%
    mutate(K = as.factor(K),
           gp = as.factor(gp)) #%>%
    # group_by(gp, K) %>%
    # summarise(
    #   median = median(var),
    #   Q1 = quantile(var, 1/4),
    #   Q3 = quantile(var, 3/4),
    #   .groups = "drop"
    # )

  # if(str_detect(metric, "r2")){
  #   p <- plot_data %>%
  #     ggplot(aes(x = gp, y = median, colour = gp)) +
  #     geom_hline(data = y_int,
  #                aes(yintercept = val),
  #                lty = 2, colour = "gray90") +
  #     geom_linerange(aes(ymin = Q1, ymax = Q3),
  #                    position = position_dodge(width = 0.5),
  #                    show.legend = FALSE) +
  #     geom_point(position = position_dodge(width = 0.5))
  # }else{
    # p <- plot_data %>%
    #   ggplot(aes(x = gp, y = median, colour = gp)) +
    #   geom_linerange(aes(ymin = Q1, ymax = Q3),
    #                  position = position_dodge(width = 0.5),
    #                  show.legend = FALSE) +
    #   geom_point(position = position_dodge(width = 0.5))
    p <- plot_data %>%
      ggplot(aes(x = gp, y = var, fill = gp)) +
      geom_boxplot()
  # }
  p <- p +
    harrypotter::scale_fill_hp_d("ravenclaw") +
    labs(y = y_label,
         x = "Sample size",
         title = title) +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(size = 14),
          strip.text = element_text(size = 20),
          strip.background = element_rect(fill = "gray70"),
          plot.title = element_text(hjust = 0.5))



  if(save){
    ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/riemannian_keyfinding3_{metric}.png"),
           plot = p,
           height = height, width = 1.618 * height)
  }else{
    return(p)
  }

}
