run_simulations <- function(num_sims, n, x_dim, y_dim, sig, K, seed = NULL, L = 5,
                            muX = diag(1, x_dim), muY = diag(1, y_dim),
                            type_x = "affine", type_y = "affine", sim_type = "full_factorial",
                            tol = 1e-4, max.iter = 10, mc.cores = 1, method = NULL){

  if(is.null(seed)){
    seed <- 1994
  }

  if(!dir.exists(glue::glue("results/{sim_type}"))){
    dir.create(glue::glue("results/{sim_type}"))
  }

  dir_name <- glue::glue("results/{sim_type}/n_{n}_p_{x_dim}_q_{y_dim}_sig_{sig}_K_{K}_xtype_{type_x}_ytype_{type_y}_seedbase_{seed}")
  if(!dir.exists(dir_name)){
    dir.create(dir_name)
  }

  obj_name <- glue::glue("n_{n}_p_{x_dim}_q_{y_dim}_sig_{sig}_K_{K}_xtype_{type_x}_ytype_{type_y}_seedbase_{seed}")

  base_object <- tibble(
    sim = NULL,
    n = NULL,
    x_dim = NULL,
    y_dim = NULL,
    sig = NULL,
    K = NULL,
    type_x = NULL,
    type_y = NULL,
    L = NULL,
    tol = NULL,
    max.iter = NULL,
    P_r2 = NULL,
    Q_r2 = NULL,
    X_rmse = NULL,
    X_r2 = NULL,
    Y_rmse = NULL,
    Y_r2 = NULL,
    tp_rmse = NULL,
    uq_rmse = NULL
  )


  if(!is.null(method)){
    base_object <- base_object %>%
      mutate(method = method)
  }
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
      sim_params <- tibble(
        sim = i,
        n = n,
        x_dim = x_dim,
        y_dim = y_dim,
        sig = sig,
        K = K,
        type_x = type_x,
        type_y = type_y,
        L = L,
        tol = tol,
        max.iter = max.iter
      )

      if(!is.null(method)){
        sim_params <- sim_params %>%
          mutate(method = method)
      }

      tryCatch(
        {
          data <- generate_rpls_data(n = n, x_dim = x_dim, y_dim = y_dim, sig = sig,
                                     seed = seed * i, # hopefully this seed is unique
                                     L = L, muX = muX, muY = muY,
                                     type_x = type_x, type_y = type_y)

          model_fit <- fit_riemannian_pls(data, K = K,
                                          tol = tol, max.iter = max.iter,
                                          mc.cores = mc.cores, method = method)

          if(isTRUE(method == "tangent")){
            data <- list(
                X = model_fit$X[[1]],
                Y = model_fit$Y[[1]],
                TP = linearise_data(X = data$TP, mu = model_fit$muX),
                UQ = linearise_data(X = data$UQ, mu = model_fit$muY),
                trueP = data$trueP,
                trueQ = data$trueQ,
                type_x = "euclidean",
                type_y = "euclidean"
              )
          }

          model_info <- collect_full_measures(model_fit)


          model_metrics <- get_metrics(data = data, measures = model_info, L = L)


          sim_run_metrics <- bind_cols(sim_params, model_metrics)

          tmp <- read_rds(glue::glue("{dir_name}/{obj_name}.Rds"))
          write_rds(tmp, glue::glue("{dir_name}/{obj_name}_presim_{i}.Rds"))

          tmp <- bind_rows(tmp, sim_run_metrics)
          write_rds(tmp, glue::glue("{dir_name}/{obj_name}.Rds"))
        },
        error = function(e) print("Simulation failed")
      )
    }
  )
}
