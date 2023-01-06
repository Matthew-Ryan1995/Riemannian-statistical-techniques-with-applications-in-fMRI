#' Perform PLS
#'
#' @param X - Matrix of predictors
#' @param Y  - Matrix of response
#' @param K - Number of variables to extract
#'
#' @return
#' @export
#'
#' @examples
#' NA
perform_pls <- function(X, Y, K){

  # Fit PLS model with my function
  pls_model <- riemannian_pls(X = X, Y = Y, L = K, tol = 1e-07, max.iter = 200,
                              type_x = "euclidean", type_y = "euclidean")

  # Extract relevant matrices
  P_fit <- do.call(cbind, pls_model$loadingsX)
  Q_fit <- do.call(cbind, pls_model$loadingsY)
  T_fit <- do.call(cbind, pls_model$scoresX)
  C_fit <- do.call(cbind, pls_model$weightsY)
  B_fit <- map_dbl(pls_model$reg_steps, ~.x$b1)

  # This stops tantrums
  if(K == 1){
    B_fit = matrix(B_fit)
  }else{
    B_fit <- diag(B_fit)
  }


  # Recover X and Y
  X_fit <- tcrossprod(T_fit, P_fit)
  Y_fit <- tcrossprod(tcrossprod(T_fit, B_fit), C_fit)

  # Return relevant matrices
  return(tibble(
    P_fit = list(P_fit),
    Q_fit = list(Q_fit),
    X_fit = list(X_fit),
    Y_fit = list(Y_fit),
  ))
}
