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
get_model_metrics <- function(M, test_data, L, formula = ~ 1){

  print(L)
  # Get our X
  X_test <- test_data$cors
  # Get our response data
  Y_test <- model.matrix(formula, data = test_data) %>%
    as.matrix() %>%
    .[, -1]

  # Centre our data
  if(isTRUE(M$type_y == "euclidean")){
    Y_test <- t((t(Y_test) - M$muY)/M$Y_sd)
  }else{
    Y_test <- linearise_data(Y_test, M$muY)
  }
  # Get predicted values
  preds <- predict.riemannian_pls(object = M, newdata = X_test, tol = 1e-4, max.iter = 20,
                                  num_comp = L, method = "tangent")


  # Get R2 and RMSE
  r2_1 <- rsq.riemannian_pls(truth = Y_test, est = preds, muY = M$muY, m = M$muY,
                           type_y = "euclidean")
  r2_2 <- rsq.riemannian_pls(truth = Y_test, est = preds, muY = M$muY, m = M$muY,
                           type_y = "euclidean", test = TRUE)
  rmse <- rmse.riemannian_pls(truth = Y_test, est = preds,
                              type_y = "euclidean")
  met_set <- metric_set(accuracy, sens, spec)

  group_classification <- preds %>%
    as_tibble() %>%
    select(pred = contains("group")) %>%
    mutate(pred = ifelse(pred < 0, levels(test_data$group)[1], levels(test_data$group)[2]),
           pred = factor(pred, levels(test_data$group))) %>%
    add_column(truth = test_data$group) %>%
    met_set(truth = truth, estimate = pred) %>%
    select(-.estimator) %>%
    pivot_wider(names_from = .metric, values_from = .estimate)

  res <-
    tibble(
      K = L,
      r2_1 = r2_1,
      r2_2 = r2_2,
      rmse = rmse
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
fit_and_test_single_fold <- function(cv_fold, formula = ~ 1, mc.cores = 1){

  # Get split and ID from CV fold
  fold <- cv_fold$splits[[1]]
  id <- cv_fold$id

  # Extract training and testing data
  train_data <- analysis(fold)
  test_data <- assessment(fold)

  # Largest K we will consider
  # May need to change this to n - max(seg_length) - 1
  K <- 50#min(nrow(train_data) - 1,
           #choose(nrow(train_data$cors[[1]]) + 1, 2))

  # Get my X and Y
  X <- train_data$cors
  Y <- model.matrix(formula, data = train_data) %>%
    as.matrix() %>%
    .[, -1]

  print(id)
  print("Fitting model")
  # Fit the model on largest number of components considered
  M <- riemannian_pls(X = X, Y = Y, L = K, tol = 1e-4, max.iter = 20,
                      type_x = "affine", type_y = "euclidean", method = "tangent",
                      scale = TRUE)
  L <- 1:K

  # Get metrics under each sub-component selection
  # fold_fit <- map_dfr(L,
  #                     function(l){
  #                       get_model_metrics(M = M, test_data = test_data, L = l, formula = formula)
  #                     })
  fold_fit <- parallel::mclapply(1:K,
                                 function(l){
                                   get_model_metrics(M = M, test_data = test_data, L = l, formula = formula)
                                 },
                                 mc.cores = mc.cores)
  fold_fit <- do.call(bind_rows, fold_fit)

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
cv_rpls_fit <- function(full_folds, formula = ~ 1, mc.cores = 1){
  V <- nrow(full_folds) # Total number of folds
  # CV ofr each fold
  res <- map_dfr(1:V, function(v){
    folds <- slice(full_folds, v)
    fit_and_test_single_fold(folds, formula = formula, mc.cores = mc.cores)
  }
  )

  return(res)
}
