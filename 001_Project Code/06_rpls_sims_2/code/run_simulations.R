run_simulations <- function(num_sims, n, x_dim, y_dim, sig, K, r, seed = NULL, L = 5,
                            muX = diag(1, x_dim), muY = diag(1, y_dim),
                            type_x = "affine", type_y = "affine", sim_type = "full_factorial",
                            tol = 1e-4, max.iter = 10, mc.cores = 1){

  if(is.null(seed)){
    seed <- 1994
  }

  if(!dir.exists(glue::glue("results/{sim_type}"))){
    dir.create(glue::glue("results/{sim_type}"))
  }

  dir_name <- glue::glue("results/{sim_type}/r_{r}_K_{K}_xtype_{type_x}_ytype_{type_y}_seedbase_{seed}")
  if(!dir.exists(dir_name)){
    dir.create(dir_name)
  }

  obj_name <- glue::glue("r_{r}_K_{K}_xtype_{type_x}_ytype_{type_y}_seedbase_{seed}")

  base_object <- tibble(
    sim = NULL,
    r = NULL,
    dist_error_X = NULL,
    dist_error_Y = NULL,
    K = NULL,
    type_x = NULL,
    type_y = NULL,
    L = NULL,
    method = NULL,
    P_r2 = NULL,
    Q_r2 = NULL,
    X_rmse = NULL,
    X_r2 = NULL,
    Y_rmse = NULL,
    Y_r2 = NULL,
    tp_rmse = NULL,
    uq_rmse = NULL
  )
  if(!file.exists(glue::glue("{dir_name}/{obj_name}.Rds"))){
    write_rds(base_object, glue::glue("{dir_name}/{obj_name}.Rds"))
  }

  res <- walk(
    1:num_sims,
    function(i){
      print(glue::glue("Running simulation {i} of {num_sims}"))
      if(file.exists(glue::glue("{dir_name}/{obj_name}_presim_{i}.Rds"))){
        print(glue::glue("Simulation {i} already complete"))
        return(NULL)
      }
      if(seed == 46666 & i == 52){
        return(NULL)
      }
      sim_params <- tibble(
        sim = i,
        r = r,
        K = K,
        type_x = type_x,
        type_y = type_y,
        L = L
      )
      # tryCatch(
      #   {
          data <- generate_rpls_data(n = n, x_dim = x_dim, y_dim = y_dim, sig = sig, r = r,
                                     seed = seed * i, # hopefully this seed is unique
                                     L = L, muX = muX, muY = muY,
                                     type_x = type_x, type_y = type_y)

          ## Full RPLS model
          model_fit <- fit_riemannian_pls(data, K = K,
                                          tol = tol, max.iter = max.iter,
                                          mc.cores = mc.cores)

          model_info <- collect_full_measures(model_fit)

          model_metrics <- get_metrics(data = data, measures = model_info, L = L)

          model_metrics <- model_metrics %>%
            mutate(method = "full")

          sim_run_metrics <- bind_cols(sim_params, model_metrics)

          # tangent RPLS model
          model_fit_tangent <- fit_riemannian_pls(data, K = K, method = "tangent",
                                                  tol = tol, max.iter = max.iter,
                                                  mc.cores = mc.cores)

          data_tangent <- list(
            X = model_fit_tangent$X[[1]],
            Y = model_fit_tangent$Y[[1]],
            TP = linearise_data(X = data$TP, mu = model_fit$muX),
            UQ = linearise_data(X = data$UQ, mu = model_fit$muY),
            type_x = model_fit_tangent$type_x,
            type_y = model_fit_tangent$type_y,
            trueP = data$trueP,
            trueQ = data$trueQ
          )



          model_info_tangent <- collect_full_measures(model_fit_tangent)



          model_metrics_tangent <- get_metrics(data = data_tangent, measures = model_info_tangent, L = L)

          model_metrics_tangent <- model_metrics_tangent %>%
            mutate(method = "tangent")

          sim_run_metrics_tangent <- bind_cols(sim_params, model_metrics_tangent)

          sim_run_metrics <- bind_rows(sim_run_metrics, sim_run_metrics_tangent)

          if(type_x == "affine"){
            dist_error_X <- mean(affine_dist(data$X)^2 - dist(data_tangent$X)^2)
          }else{
            dist_error_X <- NULL
          }
          if(type_y == "affine"){
            dist_error_Y <-  mean(affine_dist(data$Y)^2 - dist(data_tangent$Y)^2)
          }else{
            dist_error_Y <- NULL
          }

          sim_run_metrics <- sim_run_metrics %>%
            mutate(
              dist_error_X = dist_error_X,
              dist_error_Y = dist_error_Y
            )

          tmp <- read_rds(glue::glue("{dir_name}/{obj_name}.Rds"))
          write_rds(tmp, glue::glue("{dir_name}/{obj_name}_presim_{i}.Rds"))

          tmp <- bind_rows(tmp, sim_run_metrics)
          write_rds(tmp, glue::glue("{dir_name}/{obj_name}.Rds"))
      #   },
      #   error = function(e) "Simulation failed"
      # )
    }
  )

  # return(res)
}
