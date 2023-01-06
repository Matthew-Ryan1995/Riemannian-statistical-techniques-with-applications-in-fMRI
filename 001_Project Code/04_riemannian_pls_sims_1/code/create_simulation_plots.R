
# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)


# functions ---------------------------------------------------------------
get_data <- function(folder, sig.filter = 0.05){
  if(folder == "noise"){ # I accidentally have some double up in my data
    set.seed(round(2022 * 2019  - 1995)) # This is to generate the list of seed values for my function

    # Number of sims
    N_sims <- 100
    sim_type <- "noise"

    # Sample size
    n <- round(seq(10, 150, length.out = 5))

    # Dimensions
    p <- q <- round(seq(5, 15, length.out = 5))


    # Noise
    sig <- c(0, 0.1, 0.5, 1, 0.05)


    # Latent vars
    L <- 5
    K <- 1:(2*L)

    # response and predictor types
    type_x <- type_y <-  c("affine", "euclidean")

    # Tolerance/iteration
    tol <- 1e-4
    max.iter <- 20

    # Mariginalise
    n <- n[3]
    p <- p[3]
    q <- q[3]
    # sig <- sig[3]

    predict_response_grid <- expand_grid(
      n = n,
      x_dim = p,
      y_dim = q,
      sig = sig,
      K = K,
      type_x = type_x,
      type_y = type_y,
      tol = tol,
      max.iter = max.iter
    ) %>%
      mutate(seed_list = sample(1:100000, size = n(), replace = FALSE)) %>%
      filter(!(type_x == "euclidean" & type_y == "euclidean"))
    seed_list <- str_c(predict_response_grid$seed_list, collapse = "|")
  }

  fl <- list.files(here::here(glue::glue("results/{folder}")), full.names = T, recursive = T)
  fl <- fl[!str_detect(fl, "presim")]

  if(folder == "noise"){
    fl <-fl[str_detect(fl, seed_list)]
  }

  df <- map_dfr(fl, read_rds)
  df <- df %>%
    filter(sig <= sig.filter, K <=5)
  return(df)
}

create_plot2 <- function(df, filllab = NULL, title = NULL){
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
create_plot <- function(df, filllab = NULL, title = NULL){
  if(isTRUE(filllab == "sigma")){
    filllab <- latex2exp::TeX("$\\sigma$")
  }
  df %>%
    mutate(name = str_remove(name, "_r2|_rmse")) %>%
    group_by(K, var, name) %>%
    summarise(m = median(value),
              ymin = quantile(value, 1/4),
              ymax = quantile(value, 3/4),
              .groups = "drop") %>%
    ggplot(aes(x = K, y = m)) +
    # geom_vline(xintercept = 5, lty = 2, col = "gray70") +
      geom_pointrange(aes(ymin = ymin, ymax = ymax,
                          colour = as.factor(var),
                          group = as.factor(var)),
                      position = position_dodge(width = 0.5)) +
    facet_wrap(~name, ncol = 2) +
    theme_classic() +
    harrypotter::scale_colour_hp_d("ravenclaw") +
    labs(x = "K", y = NULL, colour = filllab, title = title) +
    theme(
      text = element_text(size = 14),
      legend.position = "bottom",
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "gray70"),
      plot.title = element_text(hjust = 0.5)
    )
}

make_plots <- function(df, var, x_type, y_type, filllab = NULL,
                       savename = NULL, goodcopy = FALSE, height = 5){
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
    rename(var = all_of(var)) %>%
    mutate(K = as.factor(K), var = as.factor(var))

  if(x_type == "affine"){
     x_type <- "SPD"
   }
  if(y_type == "affine"){
    y_type <- "SPD"
  }

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
    ) #+
    # geom_point(data = best_subspace_recovery, shape = 4, size = 4, colour = "red")
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
      plot_annotation(title = str_c(str_to_upper(x_type), "predictor,",
                                    str_to_upper(y_type), "response.", sep = " "),
                      theme = theme(plot.title = element_text(hjust = 0.55, size = 18)))
    ggsave(str_c(savename, ".png"), plot = p, width = width, height = 10.5)
  }


  return(p1/p2/p3)
}


# parameters ----------------------------------------------------------------

folder_and_var_names <- tibble(
  folder = c("noise", "p", "q", "sample_size"),
  var = c("sig", "x_dim", "y_dim", "n"),
  filllab = c("sigma", "p", "q", "n"),
  sig.filter = c(1, rep(0.05, 3))
)

goodcopy <- FALSE

walk(
  1:nrow(folder_and_var_names),
  function(i){
    params <- folder_and_var_names %>%
      slice(i)
    df <- get_data(params$folder, params$sig.filter)

    make_plots(df = df, var = params$var, x_type = "affine", y_type = "affine",
               filllab = params$filllab, savename = params$var, goodcopy = goodcopy)
    make_plots(df = df, var = params$var, x_type = "euclidean", y_type = "affine",
               filllab = params$filllab, savename = params$var, goodcopy = goodcopy)
    make_plots(df = df, var = params$var, x_type = "affine", y_type = "euclidean",
               filllab = params$filllab, savename = params$var, goodcopy = goodcopy)
  }
)

create_noise_by_dim <- function(){
  folder <- "p"
  fl <- list.files(here::here(glue::glue("results/{folder}")), full.names = T, recursive = T)
  fl <- fl[!str_detect(fl, "presim")]

  df <- map_dfr(fl, read_rds)
  df <- df %>%
    filter(sig == 0.5, K <=5)

  make_plots(df = df, var = "x_dim", x_type = "affine", y_type = "affine",
             filllab = "p", savename = "sig_noise", goodcopy = TRUE)
}

create_noise_by_dim()


# archive -----------------------------------------------------------------

# make_plots_orig <- function(df, var, x_type, y_type, filllab = NULL,
#                             savename = NULL, goodcopy = FALSE, height = 5){
#   if(!is.null(savename)){
#     savename <- str_c(savename, x_type, y_type, "sim_results", sep = "_")
#     if(goodcopy){
#       folder <- "../../002_Thesis/01_initial_draft/img/"
#     }else{
#       folder <- "misc_analysis/fig/"
#     }
#     savename <- str_c(folder, savename, sep = "/")
#   }
#
#   if(x_type == "affine"){
#     x_type <- "SPD"
#   }
#   if(y_type == "affine"){
#     y_type <- "SPD"
#   }
#
#   df <- df %>%
#     filter(type_x == x_type,
#            type_y == y_type) %>%
#     rename(var = all_of(var)) %>%
#     mutate(K = as.factor(K), var = as.factor(var))
#
#   best_subspace_recovery <- tibble(
#     K = 1:5,
#     value = map_dbl(K, ~min(.x/5, 1))
#   ) %>%
#     mutate(K = as.factor(K))
#
#   p1 <- df %>%
#     select(K, var, P_r2, Q_r2) %>%
#     pivot_longer(-c(K, var)) %>%
#     create_plot(filllab = filllab,
#                 # title = latex2exp::TeX(
#                 #   str_c(
#                 #     "Loading space recovery: $R^2$, ",
#                 #     str_to_title(x_type), " predictor, ",
#                 #     str_to_title(y_type), " response."
#                 #   )
#                 # )
#                 title = bquote(atop("Loading space recovery:"~R^2,
#                                     .(str_c(x_type, "predictor,",
#                                             y_type, "response.", sep = " "))))
#     ) +
#     geom_point(data = best_subspace_recovery, shape = 4, size = 4, colour = "red")
#   p2 <- df %>%
#     select(K, var, X_r2, Y_r2) %>%
#     pivot_longer(-c(K, var)) %>%
#     create_plot(filllab = filllab,
#                 # title = latex2exp::TeX(
#                 #   str_c(
#                 #     "Response and predictor recovery: $R^2$, ",
#                 #     str_to_title(x_type), " predictor, ",
#                 #     str_to_title(y_type), " response."
#                 #   )
#                 # )
#                 title = bquote(atop("Response and predictor recovery:"~R^2,
#                                     .(str_c(x_type, "predictor,",
#                                             y_type, "response.", sep = " "))))
#     )
#   p3 <- df %>%
#     select(K, var, X_rmse, Y_rmse) %>%
#     pivot_longer(-c(K, var)) %>%
#     create_plot(filllab = filllab,
#                 # title = latex2exp::TeX(
#                 #   str_c(
#                 #     "Response and predictor recovery: RMSE, ",
#                 #     str_to_title(x_type), " predictor, ",
#                 #     str_to_title(y_type), " response."
#                 #   )
#                 # )
#                 title = bquote(atop("Response and predictor recovery: RMSE",
#                                     .(str_c(x_type, "predictor,",
#                                             y_type, "response.", sep = " "))))
#     )
#
#   if(!is.null(savename)){
#     width <- 1.618 * height
#     ggsave(str_c(savename, "_pq.png"), plot = p1, width = width, height = height)
#     ggsave(str_c(savename, "_xy_r2.png"), plot = p2, width = width, height = height)
#     ggsave(str_c(savename, "_xy_rmse.png"), plot = p3, width = width, height = height)
#     p <- (p1 + theme(legend.position = "none"))/(p2+ theme(legend.position = "none"))/p3
#     ggsave(str_c(savename, ".png"), plot = , width = width, height = 11)
#   }
#
#
#   # return(p1/p2/p3)
# }
