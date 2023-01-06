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
