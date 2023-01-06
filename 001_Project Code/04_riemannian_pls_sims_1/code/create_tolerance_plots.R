# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)

# functions ---------------------------------------------------------------
get_data <- function(folder){
  fl <- list.files(here::here(glue::glue("results/{folder}")), full.names = T, recursive = T)
  fl <- fl[!str_detect(fl, "presim")]
  fl <- fl[!str_detect(fl, "sig_0.5")]

  df <- map_dfr(fl, read_rds)
  return(df)
}

create_tolerance_plot <- function(df, title = NULL){
  df %>%
    mutate(name = str_remove(name, "_r2|_rmse")) %>%
    ggplot(aes(x = tol, y = value)) +
    geom_boxplot(aes(fill = max.iter)) +
    facet_wrap(~name, ncol = 2) +
    theme_classic() +
    harrypotter::scale_fill_hp_d("ravenclaw") +
    labs(x = "Tolerance", y = NULL, fill = "Maximum iterations", title = title) +
    theme(
      text = element_text(size = 16),
      legend.position = "bottom",
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "gray70"),
      plot.title = element_text(hjust = 0.5)
    )
}



make_tolerance_plots <- function(df, var, x_type, y_type, filllab = NULL,
                                 savename = NULL, goodcopy = FALSE, height = 4){
  if(!is.null(savename)){
    savename <- str_c(savename, x_type, y_type, "sim_results", sep = "_")
    if(goodcopy){
      folder <- "../../002_Thesis/01_initial_draft/img/"
    }else{
      folder <- "misc_analysis/fig/"
    }
    savename <- str_c(folder, savename, sep = "/")
  }

  df <- df %>%
    filter(type_x == x_type,
           type_y == y_type) %>%
    mutate(tol = as.factor(tol),
           tol = fct_rev(tol),
           max.iter = as.factor(max.iter))

  if(x_type == "affine"){
    x_type <- "SPD"
  }
  if(y_type == "affine"){
    y_type <- "SPD"
  }
  p1 <- df %>%
    select(tol, max.iter, P_r2, Q_r2) %>%
    pivot_longer(-c(tol, max.iter)) %>%
    create_tolerance_plot(
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
    )
  # p2 <- df %>%
  #   select(tol, max.iter, X_r2, Y_r2) %>%
  #   pivot_longer(-c(tol, max.iter)) %>%
  #   create_tolerance_plot(
  #     # title = latex2exp::TeX(
  #     #   str_c(
  #     #     "Response and predictor recovery: $R^2$, ",
  #     #     str_to_title(x_type), " predictor, ",
  #     #     str_to_title(y_type), " response."
  #     #   )
  #     # )
  #     # title = bquote(atop("Response and predictor recovery:"~R^2,
  #     #                     .(str_c(x_type, "predictor,",
  #     #                             y_type, "response.", sep = " "))))
  #     title = expression("Response and predictor recovery:"~R^2)
  #     #   expression(
  #     #   atop(
  #     #     "Response and predictor recovery: "^2 ", ",
  #     #     str_to_title(x_type) " predictor, "
  #     #     str_to_title(y_type) " response."
  #     #   )
      # )
    # )
  p3 <- df %>%
    select(tol, max.iter, tp_rmse, uq_rmse) %>%
    pivot_longer(-c(tol, max.iter)) %>%
    create_tolerance_plot(
      # title = latex2exp::TeX(
      #   str_c(
      #     "Response and predictor recovery: RMSE, ",
      #     str_to_title(x_type), " predictor, ",
      #     str_to_title(y_type), " response."
      #   )
      # title = bquote(atop("Response and predictor recovery: RMSE",
      #                     .(str_c(x_type, "predictor,",
      #                             y_type, "response.", sep = " "))))
      title = expression("Response and predictor signal: RMSE")

    )

  # return(
  #   ((p1 + theme(legend.position = "none"))/p3) +
  #     plot_annotation(title = str_c(x_type, "predictor,",
  #                                   y_type, "response.", sep = " "),
  #                     theme = theme(plot.title = element_text(hjust = 0.55, size = 18)))
  # )

  if(!is.null(savename)){
    width <- 1.618 * height
    ggsave(str_c(savename, "_pq.png"), plot = p1, width = width, height = height)
    # ggsave(str_c(savename, "_xy_r2.png"), plot = p2, width = width, height = height)
    ggsave(str_c(savename, "_tpuq_rmse.png"), plot = p3, width = width, height = height)
    p <- ((p1 + theme(legend.position = "none"))/p3) +
      plot_annotation(title = str_c(x_type, "predictor,",
                                    y_type, "response.", sep = " "),
                      theme = theme(plot.title = element_text(hjust = 0.55, size = 18)))
    ggsave(str_c(savename, ".png"), plot = , width = width, height = 11)
  }
}


# create plots ------------------------------------------------------------

df <- get_data("tolerance")
goodcopy <- TRUE

make_tolerance_plots(df = df, var = "tolerance",
                     x_type = "affine", y_type = "affine",
                     savename = "tolerance", goodcopy = goodcopy, height = 5)
