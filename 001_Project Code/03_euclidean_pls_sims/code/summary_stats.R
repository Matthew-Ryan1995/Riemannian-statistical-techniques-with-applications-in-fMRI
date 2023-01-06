#' calculate_subspace_distance
#'
#' @param U - Vectors spanning a subspace
#' @param V - Vectors spanning a subspace
#' @param K - Dimension of smaller subspace
#' @param L - Dimension of true space
#'
#' @return
#' @export
#'
#' @examples
calculate_subspace_distance <- function(U, V, K, L = 5){

  # Get orthonormal bases
  Uhat <- qr.Q(qr(U))
  Vhat <- qr.Q(qr(V))

  # Calculate angles between
  uv_svd <- svd(crossprod(Uhat, Vhat))
  uv_svd$d[uv_svd$d > 1 - 1e-8] <- 1 # Working at tolerance

  # Geo distance of the grassmannian
  dist <- sqrt(sum(acos(uv_svd$d)^2) + max(L - K, 0) * pi^2/4)

  return(dist)
}

#' R2 for subspaces
#'
#' @param U - Vectors spanning a subspace
#' @param V - Vectors spanning a subspace
#' @param K - Dimension of smaller subspace
#' @param L - Dimension of true space
#'
#' @return
#' @export
#'
#' @examples
#' NA
calculate_subspace_r2 <- function(U, V, K, L = 5){
  max_dist <- L * pi^2/4 # Largest possible distance we can see

  square_dist <- calculate_subspace_distance(U, V, K, L)^2 # Get distance between subspaces

  r2 <- 1 - square_dist/max_dist # R2 value

  return(r2)
}


#' Get RMSE for data X and Y
#'
#' @param est - The estimate matrix
#' @param truth - The true matrix
#'
#' @return
#' @export
#'
#' @details Gets RMSE defined by
#'
#' MSE = \frac{1}{n} \sum_{i = 1}^n \| y_i - \hat{y}_i\|^2
#'
#' @examples
#' NA
calculate_xy_rmse <- function(est, truth){

  # Number of subjects
  subs <- nrow(est)

  # Calc MSE
  mse <- sum((est - truth)^2)/subs

  return(sqrt(mse))
}

#' Get data R2
#'
#' @param est - The estimate matrix
#' @param truth - The true matrix
#'
#' @return
#' @export
#'
#' @details
#'
#' Defined by
#'
#' R^2 = 1 - ModelError/DataError
#'
#' @examples
#' NA
calculate_xy_r2 <- function(est, truth){

  # Frobenius from estimate to truth
  model_dist <- sum((est - truth)^2)

  # Column centre the truth
  centred_truth <- apply(truth, 2, scale, scale = FALSE)

  # Frobenius distance from truth to centre
  data_dist <- sum(centred_truth^2)

  # R2 value
  r2 <- 1 - model_dist/data_dist

  return(r2)
}
