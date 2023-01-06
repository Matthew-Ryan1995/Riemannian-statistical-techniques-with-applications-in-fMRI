#' get_model_metrics
#'
#' @param M model
#' @param test_data test data
#' @param L number of components
#' @param formula model being fit for Y
#'
#' @return
#' @export
#'
#' @examples
get_model_metrics <- function(M, train_data, test_data, param, method){

  print(param)

  if(method == "pls_da"){
    M <- tryCatch(M %>%
      set_args(num_comp = param) %>%
      fit(group ~ ., data = train_data),
      error = function(e) "NO")

  }else if(method == "lr"){
    M <- tryCatch(M %>%
      set_args(penalty = param) %>%
      fit(group ~ ., data = train_data),
      error = function(e) "NO")
  }

  if(isTRUE(M == "NO")){
    return(
      tibble(
        param = param,
        method = method,
        accuracy = NA,
        spec = NA,
        sens = NA
      )
    )
  }

  # Get predicted values
  preds <- tryCatch(predict(M, new_data = test_data) %>%
    add_column(truth = test_data$group),
  error = function(e) "NO")

  if(isTRUE(preds == "NO")){
    return(
      tibble(
        param = param,
        method = method,
        accuracy = NA,
        spec = NA,
        sens = NA
      )
    )
  }

  met_set <- metric_set(accuracy, sens, spec)

  group_classification <- preds %>%
    met_set(truth = truth, estimate = .pred_class) %>%
    select(-.estimator) %>%
    pivot_wider(names_from = .metric, values_from = .estimate)

  res <- tibble(
    param = param,
    method = method
  ) %>%
    bind_cols(group_classification)

  return(res)
}

#' fit_and_test_single_fold
#'
#' @param cv_fold A single row of cv_folds
#' @param formula Model for Y
#'
#' @return
#' @export
#'
#' @examples
fit_and_test_single_fold <- function(cv_fold, mc.cores = 1){

  # Get split and ID from CV fold
  fold <- cv_fold$splits[[1]]
  id <- cv_fold$id

  # Extract training and testing data
  train_data <- analysis(fold)
  test_data <- assessment(fold)

  # Largest K we will consider
  # May need to change this to n - max(seg_length) - 1
  K <- min(nrow(train_data) - 1, 10)
           # choose(nrow(train_data$cors[[1]]) + 1, 2))

  # Get my X and Y
  X <- train_data$cors

  frechet_mean <- get_frechet_mean(X, cores = mc.cores)

  df_train <- linearise_data(X = X, mu = frechet_mean, cores = mc.cores) %>%
    as_tibble() %>%
    add_column(group = train_data$group)

  df_test <- linearise_data(X = test_data$cors, mu = frechet_mean, cores = mc.cores) %>%
    as_tibble() %>%
    add_column(group = test_data$group)

  # Fit the model on largest number of components considered
  pls_da_mod <- pls(mode = "classification")
  lr_mod <- logistic_reg()

  num_comp <- 1:K
  lambda <- c(0.01, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 3.0, 5.0)



  print(id)
  pls_fit <- parallel::mclapply(num_comp,
                                 function(l){
                                   get_model_metrics(M = pls_da_mod,
                                                     train_data = df_train,
                                                     test_data = df_test,
                                                     param = l,
                                                     method = "pls_da")
                                 },
                                 mc.cores = mc.cores)
  pls_fit <- do.call(bind_rows, pls_fit)
  lr_fit <- parallel::mclapply(lambda,
                                 function(l){
                                   get_model_metrics(M = lr_mod,
                                                     train_data = df_train,
                                                     test_data = df_test,
                                                     param = l,
                                                     method = "lr")
                                 },
                                 mc.cores = mc.cores)
  lr_fit <- do.call(bind_rows, lr_fit)

  fold_fit <- bind_rows(pls_fit, lr_fit)
  # fold_fit <- pls_fit

  fold_fit <- fold_fit %>%
    mutate(id = id) %>%
    select(id, everything())

  return(fold_fit)
}

#' cv_rpls_fit
#'
#' @param full_folds A v_folds from rsample
#' @param formula Model for Y
#'
#' @return
#' @export
#'
#' @examples
cv_other_fit <- function(full_folds, mc.cores = 1){
  V <- nrow(full_folds) # Total number of folds
  # CV ofr each fold
  res <- map_dfr(1:V, function(v){
    folds <- slice(full_folds, v)
    fit_and_test_single_fold(folds, mc.cores = mc.cores)
  }
  )

  return(res)
}
