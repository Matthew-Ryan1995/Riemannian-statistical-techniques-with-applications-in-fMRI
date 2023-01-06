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


#' get_metrics
#'
#' @param data - The data I have fit one
#' @param measures - The measures returned
#' @param L - The true number of latent variables
#'
#' @return
#' @export
#'
#' @examples
get_metrics <- function(data, measures, L = 5){
  X.rmse <- rmse.riemannian_pls(truth = data$X, est = measures$X.hat,
                                type_y = data$type_x, muY = measures$muX, ignore = TRUE)
  X.rsq <- rsq.riemannian_pls(truth = data$X, est = measures$X.hat,
                              type_y = data$type_x, muY = measures$muX, ignore = TRUE)

  tp.rmse <- rmse.riemannian_pls(truth = data$TP, est = measures$X.hat,
                                 type_y = data$type_x, muY = measures$muX, ignore = TRUE)

  Y.rmse <- rmse.riemannian_pls(truth = data$Y, est = measures$Y.hat,
                                type_y = data$type_y, muY = measures$muY, ignore = TRUE)
  Y.rsq <- rsq.riemannian_pls(truth = data$Y, est = measures$Y.hat,
                              type_y = data$type_y, muY = measures$muY, ignore = TRUE)
  uq.rmse <- rmse.riemannian_pls(truth = data$UQ, est = measures$Y.hat,
                                 type_y = data$type_y, muY = measures$muY, ignore = TRUE)

  K <- measures$K

  P.rsq <- calculate_subspace_r2(data$trueP, measures$P.hat, K = K, L = L)
  Q.rsq <- calculate_subspace_r2(data$trueQ, measures$Q.hat, K = K, L = L)

  return(
    tibble(
      P_r2 = P.rsq,
      Q_r2 = Q.rsq,
      X_rmse = X.rmse,
      X_r2 = X.rsq,
      Y_rmse = Y.rmse,
      Y_r2 = Y.rsq,
      tp_rmse = tp.rmse,
      uq_rmse = uq.rmse
    )
  )
}
