create_noise_finding_plots <- function(df, x_type, y_type, vars, height = 5, save = FALSE){

  y_label <- case_when(
    str_detect(vars, "P_r2") ~ "$R^2_S$ (P)",
    str_detect(vars, "Q_r2") ~ "$R^2_S$ (Q)",
    str_detect(vars, "X_r2") ~ "$R^2$ (X)",
    str_detect(vars, "Y_r2") ~ "$R^2$ (Y)",
    str_detect(vars, "X_rmse") ~ "RMSE (X)",
    str_detect(vars, "Y_rmse") ~ "RMSE (Y)",
    str_detect(vars, "tp_rmse") ~ "RMSE (TP)",
    str_detect(vars, "uq_rmse") ~ "RMSE (UQ)"
  ) %>%
    latex2exp::TeX()


  col_lab <- latex2exp::TeX("$\\sigma$:")

  title <- str_c(
    "R-NIPALS\n",
    "Predictor: ", ifelse(x_type == "affine", "SPD", str_to_upper(x_type)),
    "; Response: ", ifelse(y_type == "affine", "SPD", str_to_upper(y_type)), "\n",
    "| p = 10 | q = 10 | n = 80 |"
  )

  y_int <- expand_grid(
    val = c(0.2, 0.4, 0.6, 0.8, 1),
    metric = factor(c("P: R2", "Q: R2", "X: R2", "Y: R2",
                      "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE"),
                    levels = c("P: R2", "Q: R2", "X: R2", "Y: R2",
                               "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE"))
  ) %>%
    filter(!(str_detect(metric, "RMSE")))

  plot_data <- df %>%
    filter(type_x == x_type,
           type_y == y_type,
           K <= 5) %>%
    select(sig, K, var = all_of(vars)) %>%
    mutate(K = as.factor(K),
           sig = as.factor(sig)) %>%
    group_by(sig, K) %>%
    summarise(
      median = median(var),
      Q1 = quantile(var, 1/4),
      Q3 = quantile(var, 3/4),
      .groups = "drop"
    )

  if(str_detect(vars, "r2")){
    p <- plot_data %>%
      ggplot(aes(x = K, y = median, colour = sig, shape = sig)) +
      geom_hline(data = y_int,
                 aes(yintercept = val),
                 lty = 2, colour = "gray90") +
      geom_linerange(aes(ymin = Q1, ymax = Q3),
                     position = position_dodge(width = 0.5),
                     show.legend = FALSE) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }else{
    p <- plot_data %>%
      ggplot(aes(x = K, y = median, colour = sig, shape = sig)) +
      geom_linerange(aes(ymin = Q1, ymax = Q3),
                     position = position_dodge(width = 0.5),
                     show.legend = FALSE) +
      geom_point(size = 3,position = position_dodge(width = 0.5))
  }
  p <- p +
    harrypotter::scale_color_hp_d("ravenclaw") +
    labs(y = y_label,
         colour = col_lab,
         shape = col_lab,
         title = title) +
    theme_classic() +
    theme(legend.position = "bottom",
          text = element_text(size = 14),
          strip.text = element_text(size = 20),
          strip.background = element_rect(fill = "gray70"),
          plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE))



  if(save){
    ggsave(glue::glue("../../002_Thesis/01_initial_draft/img/riemannian_keyfinding1_predictor_{x_type}_response_{y_type}_{vars}.png"),
           plot = p,
           height = height, width = 1.618 * height)
  }else{
    return(p)
  }

}
