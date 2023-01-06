#' Fit and run the simulations for a given set of parameters
#'
#' @param sims - How many sims to run
#' @param n - Vector of sample sizes
#' @param p - Predictor space dimension
#' @param q - Response space dimension
#' @param sig - Vector of errors
#' @param L - True latent space size
#' @param seed - Seed
#'
#' @return
#' @export
#'
#' @examples
fit_simulation <- function(sims = 1, n, p, q, sig, L = 5, seed = 1234){

  file_name <- glue::glue("results/p_{p}_q_{q}_seed_{seed}/")

  if(!dir.exists(file_name)){
    dir.create(file_name)
  }

  if(file.exists(str_c(file_name, glue::glue("summaryStats_p_{p}_q_{q}_seed_{seed}.Rds")))){
    return(NULL)
  }

  if(p < 200){
    df <- gen_data(sims, n, p, q, sig, L, seed)

    # write_rds(df, str_c(file_name, glue::glue("simulations_p_{p}_q_{q}_seed_{seed}.Rds")))

    P <- df$P[[1]]
    Q <- df$Q[[1]]

    write_rds(P, str_c(file_name, glue::glue("P_p_{p}_q_{q}_seed_{seed}.Rds")))
    write_rds(Q, str_c(file_name, glue::glue("Q_p_{p}_q_{q}_seed_{seed}.Rds")))

    rm("P", "Q")

    df <- df %>%
      mutate(model = pmap(list(X, Y, K), perform_pls)) %>%
      unnest(model)

    # write_rds(df, str_c(file_name, glue::glue("model_p_{p}_q_{q}_seed_{seed}.Rds")))

    df <- df %>%
      mutate(P_r2 = pmap_dbl(list(P_fit, P, K), calculate_subspace_r2),
             Q_r2 = pmap_dbl(list(Q_fit, Q, K), calculate_subspace_r2),
             X_rmse = map2_dbl(X_fit, X, calculate_xy_rmse),
             X_r2 = map2_dbl(X_fit, X, calculate_xy_r2),
             Y_rmse = map2_dbl(Y_fit, Y, calculate_xy_rmse),
             Y_r2 = map2_dbl(Y_fit, Y, calculate_xy_r2),
             TP_rmse = map2_dbl(X_fit, TP, calculate_xy_rmse),
             UQ_rmse = map2_dbl(Y_fit, UQ, calculate_xy_rmse),
             ) %>%
      select(-X, -Y, -P, -Q, -TP, -UQ, -contains("fit"))

    write_rds(df, str_c(file_name, glue::glue("summaryStats_p_{p}_q_{q}_seed_{seed}.Rds")))
  }else{
    df_list <- list()
    k <- 0
    for(s in (sig)){
      for(ni in n){
        k <- k + 1
        df <- gen_data(sims, ni, p, q, s, L, seed + k)


        P <- df$P[[1]]
        Q <- df$Q[[1]]

        # if(k == 1){
          write_rds(P, str_c(file_name, glue::glue("P_p_{p}_q_{q}_seed_{seed + k}.Rds")))
          write_rds(Q, str_c(file_name, glue::glue("Q_p_{p}_q_{q}_seed_{seed + k}.Rds")))
        # }

        rm("P", "Q")

        df <- df %>%
          mutate(model = pmap(list(X, Y, K), perform_pls)) %>%
          unnest(model)


        df <- df %>%
          mutate(P_r2 = pmap_dbl(list(P_fit, P, K), calculate_subspace_r2),
                 Q_r2 = pmap_dbl(list(Q_fit, Q, K), calculate_subspace_r2),
                 X_rmse = map2_dbl(X_fit, X, calculate_xy_rmse),
                 X_r2 = map2_dbl(X_fit, X, calculate_xy_r2),
                 Y_rmse = map2_dbl(Y_fit, Y, calculate_xy_rmse),
                 Y_r2 = map2_dbl(Y_fit, Y, calculate_xy_r2),
                 TP_rmse = map2_dbl(X_fit, TP, calculate_xy_rmse),
                 UQ_rmse = map2_dbl(Y_fit, UQ, calculate_xy_rmse),) %>%
          select(-X, -Y, -P, -Q, -TP, -UQ, -contains("fit"))

        df_list[[k]] <- df
      }
    }
    df_list <- do.call(bind_rows, df_list)
    write_rds(df_list, str_c(file_name, glue::glue("summaryStats_p_{p}_q_{q}_seed_{seed}.Rds")))
  }
}
