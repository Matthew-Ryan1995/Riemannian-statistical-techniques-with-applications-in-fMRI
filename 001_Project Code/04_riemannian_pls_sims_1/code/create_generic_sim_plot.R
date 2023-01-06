create_generic_sim_plot <- function(df, x_type, y_type, vars,
                                    height = 5, save = FALSE){
  col_lab <- case_when(
    vars == "sig" ~ "$\\sigma$",
    vars == "x_dim" ~ "$p$",
    vars == "y_dim" ~ "$q$",
    vars == "n" ~ "$n$"
  ) %>%
  latex2exp::TeX()

  title <- str_c(
    "Riemannian PLS\n",
    "Predictor: ", ifelse(x_type == "affine", "SPD", str_to_upper(x_type)),
    "; Response: ", ifelse(y_type == "affine", "SPD", str_to_upper(y_type)), "\n",
    "| ",
    ifelse(vars == "x_dim", "", "p = 10 | "),
    ifelse(vars == "y_dim", "", "q = 10 | "),
    ifelse(vars == "sig", "", "sigma = 0.05 | "),
    ifelse(vars == "n", "", "n = 80 | ")
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
    select(var = all_of(vars), K, P_r2:uq_rmse) %>%
    pivot_longer(-(var:K), names_to = "metric") %>%
    mutate(metric = case_when(
      str_detect(metric, "P_r2") ~ "P: R2",
      str_detect(metric, "Q_r2") ~ "Q: R2",
      str_detect(metric, "X_r2") ~ "X: R2",
      str_detect(metric, "Y_r2") ~ "Y: R2",
      str_detect(metric, "X_rmse") ~ "X: RMSE",
      str_detect(metric, "Y_rmse") ~ "Y: RMSE",
      str_detect(metric, "tp_rmse") ~ "TP: RMSE",
      str_detect(metric, "uq_rmse") ~ "UQ: RMSE"
    ),
    metric = factor(metric,
                    levels = c("P: R2", "Q: R2", "X: R2", "Y: R2",
                               "X: RMSE", "Y: RMSE", "TP: RMSE", "UQ: RMSE")),
    K = as.factor(K),
    var = as.factor(var)
    ) %>%
    group_by(metric, var, K) %>%
    summarise(
      median = median(value),
      Q1 = quantile(value, 1/4),
      Q3 = quantile(value, 3/4),
      .groups = "drop"
    )

  p <- plot_data %>%
    ggplot(aes(x = K, y = median, colour = var)) +
    geom_hline(data = y_int,
               aes(yintercept = val),
               lty = 2, colour = "gray90") +
    geom_linerange(aes(ymin = Q1, ymax = Q3),
                   position = position_dodge(width = 0.5),
                   show.legend = FALSE) +
    geom_point(position = position_dodge(width = 0.5)) +
    facet_wrap(~metric, scale = "free_y") +
    harrypotter::scale_color_hp_d("ravenclaw") +
    labs(y = NULL,
         colour = col_lab,
         title = title) +
    theme_classic() +
    theme(legend.position = "bottom",
          text = element_text(size = 14),
          strip.text = element_text(size = 20),
          strip.background = element_rect(fill = "gray70"),
          plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE))

  if(save){
    ggsave(glue::glue("~/Desktop/Riemannian-statistical-techniques-with-applications-in-fMRI/Model Recovery/Riemannian PLS/fig/predictor_{x_type}_response_{y_type}_vary_{vars}.png"),
           height = height, width = 1.618 * height, plot = p)
  }else{
    return(p)
  }
}
