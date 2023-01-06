
# create_plot <- function(df, filllab = NULL, title = NULL){
#   if(isTRUE(filllab == "sigma")){
#     filllab <- latex2exp::TeX("$\\sigma$")
#   }
#   df %>%
#     mutate(name = str_remove(name, "_r2|_rmse")) %>%
#     group_by(K, var, name) %>%
#     summarise(m = median(value),
#               ymin = quantile(value, 1/4),
#               ymax = quantile(value, 3/4),
#               .groups = "drop") %>%
#     ggplot(aes(x = K, y = m)) +
#     # geom_vline(xintercept = 5, lty = 2, col = "gray70") +
#     geom_pointrange(aes(ymin = ymin, ymax = ymax,
#                         colour = as.factor(var),
#                         group = as.factor(var)),
#                     position = position_dodge(width = 0.5)) +
#     facet_wrap(~name, ncol = 2) +
#     theme_classic() +
#     harrypotter::scale_colour_hp_d("ravenclaw") +
#     labs(x = "K", y = NULL, colour = filllab, title = title) +
#     theme(
#       text = element_text(size = 14),
#       legend.position = "bottom",
#       strip.text = element_text(size = 20),
#       strip.background = element_rect(fill = "gray70"),
#       plot.title = element_text(hjust = 0.5)
#     )
# }

create_boxplots <- function(df, s, savename = NULL, goodcopy = FALSE,
                            title = "", height = 5, vars = "P_r2"){

  if(!is.null(savename)){
    savename <- str_c("euclidean", savename, "sim_results", sep = "_")
    if(goodcopy){
      folder <- "../../002_Thesis/01_initial_draft/img/"
    }else{
      folder <- "figs/"
    }
    savename <- str_c(folder, savename, sep = "/")
  }

  p <- df %>%
    filter(sig == s) %>%
    rename(var = all_of(vars)) %>%
    mutate(p = str_c("p: ", p),
           p = factor(p, levels = str_c("p: ", choose(c(5, 8, 10, 12, 15) + 1, 2))),
           q = str_c("q: ", q),
           q = factor(q, levels = str_c("q: ", choose(c(5, 8, 10, 12, 15) + 1, 2))),
           K = as.factor(K)) %>%
    ggplot(aes(x = K, y = var)) +
    geom_boxplot(aes(fill = as.factor(n))) +
    facet_grid(q~p) +
    theme_classic() +
    harrypotter::scale_fill_hp_d("ravenclaw") +
    labs(x = "K", y = NULL, fill = "Sample size", title = title) +
    theme(
      text = element_text(size = 14),
      legend.position = "bottom",
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "gray70"),
      plot.title = element_text(hjust = 0.5)
    )

  # p <- df %>%
  #   filter(sig == s) %>%
  #   rename(var = all_of(vars)) %>%
  #   mutate(p = str_c("p: ", p),
  #          p = factor(p, levels = str_c("p: ", choose(c(5, 8, 10, 12, 15) + 1, 2))),
  #          q = str_c("q: ", q),
  #          q = factor(q, levels = str_c("q: ", choose(c(5, 8, 10, 12, 15) + 1, 2))),
  #          K = as.factor(K)) %>%
  #   group_by(K, n, p, q) %>%
  #   summarise(m = median(var),
  #             ymin = quantile(var, 1/4),
  #             ymax = quantile(var, 3/4),
  #             .groups = "drop") %>%
  #   ggplot(aes(x = K, y = m)) +
  #   # geom_vline(xintercept = 5, lty = 2, col = "gray70") +
  #   geom_pointrange(aes(ymin = ymin, ymax = ymax,
  #                       colour = as.factor(n),
  #                       group = as.factor(n)),
  #                   position = position_dodge(width = 1)) +
  #   # ggplot(aes(x = K, y = var)) +
  #   # geom_boxplot(aes(fill = as.factor(n))) +
  #   facet_grid(q~p) +
  #   theme_classic() +
  #   harrypotter::scale_colour_hp_d("ravenclaw") +
  #   labs(x = "K", y = NULL, colour = "Sample size", title = title) +
  #   theme(
  #     text = element_text(size = 14),
  #     legend.position = "bottom",
  #     strip.text = element_text(size = 20),
  #     strip.background = element_rect(fill = "gray70"),
  #     plot.title = element_text(hjust = 0.5)
  #   )

  if(str_detect(vars, "P|Q")){
    best_subspace_recovery <- tibble(
      K = 1:5,
      var = map_dbl(K, ~min(.x/5, 1))
    ) %>%
      mutate(K = as.factor(K))

    p <- p +
      geom_point(data = best_subspace_recovery, shape = 4, size = 4, colour = "red")
  }

  if(!is.null(savename)){
    width <- 1.618 * height
    ggsave(str_c(savename, ".png"), plot = p, width = width, height = 10.5)
  }else{
    return(p)
  }
}

create_n_sig <- function(df, savename = NULL, goodcopy = FALSE,
                         title = "", height = 5, vars = "P_r2"){

  if(!is.null(savename)){
    savename <- str_c("euclidean", savename, "sim_results", sep = "_")
    if(goodcopy){
      folder <- "../../002_Thesis/01_initial_draft/img/"
    }else{
      folder <- "figs/"
    }
    savename <- str_c(folder, savename, sep = "/")
  }

  p <- df %>%
    filter(p == 15, q == 15) %>%
    rename(var = all_of(vars)) %>%
    mutate(K = as.factor(K)) %>%
    ggplot(aes(x = K, y = var)) +
    geom_boxplot(aes(fill = as.factor(n))) +
    facet_wrap(~sig) +
    theme_classic() +
    harrypotter::scale_fill_hp_d("ravenclaw") +
    labs(x = "K", y = NULL, fill = "Sample size", title = title) +
    theme(
      text = element_text(size = 14),
      legend.position = c(0.85, 0.25),
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "gray70"),
      plot.title = element_text(hjust = 0.5)
    )

  if(str_detect(vars, "P|Q")){
    best_subspace_recovery <- tibble(
      K = 1:5,
      var = map_dbl(K, ~min(.x/5, 1))
    ) %>%
      mutate(K = as.factor(K))

    p <- p +
      geom_point(data = best_subspace_recovery, shape = 4, size = 4, colour = "red")
  }

  if(!is.null(savename)){
    width <- 1.618 * height
    ggsave(str_c(savename, ".png"), plot = p, width = width, height = 10.5)
  }else{
    return(p)
  }
}


create_plot <- function(df, filllab = NULL, title = NULL){
  if(isTRUE(filllab == "sigma")){
    filllab <- latex2exp::TeX("$\\sigma$")
  }
  df %>%
    mutate(name = str_remove(name, "_r2|_rmse")) %>%
    ggplot(aes(x = K, y = value)) +
    geom_boxplot(aes(fill = var)) +
    # geom_vline(xintercept = 5, lty = 2, col = "gray70") +
    facet_wrap(~name, ncol = 2) +
    theme_classic() +
    harrypotter::scale_fill_hp_d("ravenclaw") +
    labs(x = "K", y = NULL, fill = filllab, title = title) +
    theme(
      text = element_text(size = 14),
      legend.position = "bottom",
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "gray70"),
      plot.title = element_text(hjust = 0.5)
    )
}

make_plots <- function(df, var, sig_filter = 0.5, filllab = NULL,
                       savename = NULL, goodcopy = FALSE, height = 5){
  if(!is.null(savename)){
    savename <- str_c(savename,"sim_results", sep = "_")
    if(goodcopy){
      folder <- "../../002_Thesis/01_initial_draft/img/"
    }else{
      folder <- "misc_analysis/fig/"
    }
    savename <- str_c(folder, savename, sep = "/")
  }

  df <- df %>%
    filter(sig == sig_filter,
           p == 120, q == 120) %>%
    rename(var = all_of(var)) %>%
    mutate(K = as.factor(K), var = as.factor(var))


  best_subspace_recovery <- tibble(
    K = 1:5,
    value = map_dbl(K, ~min(.x/5, 1))
  ) %>%
    mutate(K = as.factor(K))

  p1 <- df %>%
    select(K, var, P_r2, Q_r2) %>%
    pivot_longer(-c(K, var)) %>%
    create_plot(filllab = filllab,
                # title = latex2exp::TeX(
                #   str_c(
                #     "Loading space recovery: $R^2$, ",
                #     str_to_title(x_type), " predictor, ",
                #     str_to_title(y_type), " response."
                #   )
                # )
                # title = bquote(atop("Loading space recovery:"~R^2,
                #                     .(str_c(x_type, "predictor,",
                #                             y_type, "response.", sep = " "))))
                title = expression("Loading space recovery:"~R^2)
    ) +
    geom_point(data = best_subspace_recovery, shape = 4, size = 4, colour = "red")
  p2 <- df %>%
    select(K, var, X_r2, Y_r2) %>%
    pivot_longer(-c(K, var)) %>%
    create_plot(filllab = filllab,
                # title = latex2exp::TeX(
                #   str_c(
                #     "Response and predictor recovery: $R^2$, ",
                #     str_to_title(x_type), " predictor, ",
                #     str_to_title(y_type), " response."
                #   )
                # )
                # title = bquote(atop("Response and predictor recovery:"~R^2,
                #                     .(str_c(x_type, "predictor,",
                #                             y_type, "response.", sep = " "))))
                title = expression("Response and predictor recovery:"~R^2)
    )
  p3 <- df %>%
    select(K, var, X_rmse, Y_rmse) %>%
    pivot_longer(-c(K, var)) %>%
    create_plot(filllab = filllab,
                # title = latex2exp::TeX(
                #   str_c(
                #     "Response and predictor recovery: RMSE, ",
                #     str_to_title(x_type), " predictor, ",
                #     str_to_title(y_type), " response."
                #   )
                # )
                # title = bquote(atop("Response and predictor recovery: RMSE",
                #                     .(str_c(x_type, "predictor,",
                #                             y_type, "response.", sep = " "))))
                title = expression("Response and predictor recovery: RMSE")
    )

  if(!is.null(savename)){
    width <- 1.618 * height
    ggsave(str_c(savename, "_pq.png"), plot = p1, width = width, height = height)
    ggsave(str_c(savename, "_xy_r2.png"), plot = p2, width = width, height = height)
    ggsave(str_c(savename, "_xy_rmse.png"), plot = p3, width = width, height = height)
    p <- ((p1 + theme(legend.position = "none"))/(p2+ theme(legend.position = "none"))/p3) +
      plot_annotation(
        title = str_c(str_to_upper("euclidean"), "predictor,",
                                    str_to_upper("euclidean"), "response.", sep = " "),
                      theme = theme(plot.title = element_text(hjust = 0.55, size = 18)))
    ggsave(str_c(savename, ".png"), plot = p, width = width, height = 10.5)
  }else{
    return(p1/p2/p3)
  }
}
