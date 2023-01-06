## NOTE TO SELF
## Clean up the Chol argument follow through.
## Clean up consistencies in my matrix inputs
## Maybe I set P, Q \in M, U, V \in T_I M?
# Make affine_exp and affine_log functions, just push through to affine_map


#' translate_to_identity
#'
#' A function to translate a point in the tangent space of P to the tangent space of the identity
#' by the affine group action.
#'
#' @param P Matrix: The point on the manifold
#' @param W Matrix: The point in TpM
#' @param chol Logical: Do you want ot use the Cholesky square root?
#' @param reverse Logical: Do you want to translate from the identity to P
#' @param P.sq Matrix: Optional. The square root of P, to help speed calculation.
#'
#' @return  The translated matrix at the identity.
#'
#' @details Performs the translation
#'
#'     P^{-1/2} W P^{-1/2}
#'
#' if reverse is FALSE to move W to the tangent space of the identity.  If reverse is TRUE, then
#' applies
#'
#'     P^{1/2} W P^{1/2}
#'
#' To move W from the tangent space at the identity to the tangent space at P.  Although this is
#' intended for use on tangent vectors W, this can be applied to any symmetric matrix.
#'
#' @export translate_to_identity
#'
#' @examples
#'
#' ## Get some matrices
#' ## We will look at the first two subjects in the COBRE dataset
#' P <- cobre$cors[[1]]
#' W <- cobre$cors[[2]]
#'
#' ## Translate W to a tangent vector at the identity:
#'
#' W_trans <- translate_to_identity(P, W)
#'
#' ## Translate back to W
#'
#' W_back <- translate_to_identity(P, W_trans, reverse = TRUE)
#'
#' ## Check it is the reverse (up to numerical error)
#' identical(zapsmall(W_back), zapsmall(W))
#'
translate_to_identity <- function(P, W, chol = FALSE,
                                  reverse = FALSE, P.sq = NULL){
  # The square root can be supplied to improve calculation speed
  # If not supplied, calculate either the symmetric square root of the Cholesky square root,
  # base on user input
  if(is.null(P.sq)){
    if(chol){
      P.sq <- get_sqrt(P, FALSE)
    }else{
      P.sq <- get_sqrt(P)
    }
  }

  # The reverse argument will translate from the tangent space at the identity
  # Back to the tangent space at P.
  # Otherwise, we translate from P to Identity
  if(reverse){
    inner.prod <- crossprod(P.sq, crossprod(W, P.sq))
  }else{
    P.sq.inv <- solve(P.sq)
    inner.prod <- crossprod(P.sq.inv, crossprod(W, P.sq.inv))
  }

  return(inner.prod)
}

#' affine_map
#'
#' This will apply either the Riemannian Log or Exp maps in the affine invariant geometry.
#'
#' @param P Matrix: The base point on the manifold
#' @param W Matrix: If applying the Exp map, W is a tangent vector at P.  If applying the Log map, W is another point on the manifold.
#' @param chol Logical: Should Cholesky square root be used?
#' @param FUN Character: Either exp or log, defining which map to apply.
#' @param P.sq Matrix: Optional. The square root of P, to help speed calculation.
#' @param ... Optional list:
#'    - method : The method to pass to expm or logm.  Only accepts Eigen or Higham08
#'
#' @return  A symmetric or positive definite symmetric matrix.
#'
#' @details When applying Log, you get that
#'
#'      Log_P(W) = P^{1/2} Log(P^{-1/2} W P^{-1/2}) P^{1/2}
#'
#' Where Log is the matrix logarithm.  This will take a point on the manifold W and return a
#' tangent vector in TpM.
#'
#' When applying Exp, you get
#'
#'      Exp_P(W) = P^{1/2} Exp(P^{-1/2} W P^{-1/2}) P^{1/2}
#'
#' Where Exp is the matrix logarithm.  This will take tangent vector W in TpM, and return a point in
#' the manifold.
#'
#' Both of these maps have the affect for translating to the identity, applying the matrix map, and
#' translating pack to the tangent space at P
#'
#' @export affine_map
#'
#' @examples
#' ## Get some matrices
#' ## We will look at the first two subjects in the COBRE dataset
#' P <- cobre$cors[[1]]
#' W <- cobre$cors[[2]]
#'
#' ## Push W to a tangent vector using the Logarithm
#'
#' W_vector <- affine_map(P, W, method = "log")
#'
#' ## Push W back down to the manifold using the exponential.
#'
#' W_back <- affine_map(P, W_vector, method = "exp")
#'
#' ## Check it is the reverse (up to numerical error)
#' identical(zapsmall(W_back), zapsmall(W))
#'
affine_map <- function(P, W, chol = FALSE, FUN = "exp", P.sq = NULL, ...){

  ARGS <- list(...)

  if(is.null(ARGS$method)){
    method <- "Higham08"
  }else{
    method <- ARGS$method
    if(!(method %in% c("Higham08", "Eigen"))){
      stop("Method must be Higham08 or Eigen")
    }
  }

  # Check to see they are performing a valid map
  FUN <- stringr::str_to_lower(FUN)
  if(!(FUN %in% c("exp", "log"))){
    stop("Please enter a valid method: either exp or log.")
  }

  # Calculate the square root either symmetric or Cholesky
  if(is.null(P.sq)){
    if(chol){
      P.sq <- get_sqrt(P, FALSE)
    }else{
      P.sq <- get_sqrt(P)
    }
  }

  # Translate W to the tangent space at the identity
  inner.prod <- translate_to_identity(P, W, chol = FALSE, P.sq = P.sq)

  # Apply the Riemannian map, which is matrix Exp/Log at the identity.
  if(FUN == "exp"){
    if(method == "Eigen"){
      method <- str_c("R_",method)
    }
    mat_trans <- expm::expm(inner.prod, method = method)
  }else if(FUN == "log"){
    mat_trans <- expm::logm(inner.prod, method = method)
  }else{
    stop("Enter a valid function")
  }

  # Translate the result back to P
  outer.prod <- translate_to_identity(P, mat_trans, reverse = TRUE, P.sq = P.sq)

  if(any(is.na(outer.prod))){
    affine_map(P, W, chol, FUN, P.sq, method = "Eigen")
  }

  return(outer.prod)
}

#' affine_exp
#'
#' Calculates the Riemannian exponential map under the affine-invariant geometry
#'
#' @param P Matrix: The base point on the manifold
#' @param W Matrix: If applying the Exp map, W is a tangent vector at P.  If applying the Log map, W is another point on the manifold.
#' @param chol Logical: Should Cholesky square root be used?
#' @param P.sq Matrix: Optional. The square root of P, to help speed calculation.
#' @param ... Optional arguments to be passed on to affine_map
#'
#' @return the Riemannian exponential map under the affine-invariant geometry
#' @export affine_exp
#'
#' @examples
#'
#' TBA
affine_exp <- function(P, W, chol = FALSE, P.sq = NULL, ...){
  affine_map(P, W, chol = FALSE, FUN = "exp", P.sq = NULL, ...)
}

#' affine_log
#'
#' Calculates the Riemannian log map under the affine-invariant geometry
#'
#' @param P Matrix: The base point on the manifold
#' @param W Matrix: If applying the Exp map, W is a tangent vector at P.  If applying the Log map, W is another point on the manifold.
#' @param chol Logical: Should Cholesky square root be used?
#' @param P.sq Matrix: Optional. The square root of P, to help speed calculation.
#' @param ... Optional arguments to be passed on to affine_map
#'
#' @return the Riemannian log map under the affine-invariant geometry
#' @export affine_log
#'
#' @examples
#'
#' TBA
affine_log <- function(P, W, chol = FALSE, P.sq = NULL, ...){
  affine_map(P, W, chol = FALSE, FUN = "log", P.sq = NULL, ...)
}

#' get_affine_m1
#'
#' This will calculate the Riemannian approximation of the first-order moment of the sample
#' distribution at P.  This function is used directly to calculate the Frechet mean.
#'
#' @param P Matrix: The point to find the moment around.
#' @param data list, array or data frame: The sample data
#' @param var.name Character: If you are using a dataframe, then the name of the variable in the
#' data frame.
#' @param P.sq Matrix: The square root of P, to speed up calculations. Optional
#' @param cores Integer > 0: This is set up for parallelisation on Mac through mclapply.
#' This will determine the number of cores to use. Defaults to 1.
#' @param ... Optional arguments to be passed on to affine_map
#'
#' @return A matrix of the same dimension as P that is the M1 approximation at P.
#'
#' @details Performs the calculation of
#'
#'     M1(P) = \frac{1}{m} \sum P^{1/2} Log(P^{-1/2} P_i P^{-1/2}) P^{1/2}
#'
#' where P_i is each observed point in data, and m is the number of observations in data.
#'
#' @export get_affine_m1
#'
#' @examples
#'
#' # We will perform the first step in calculating the Frechet mean for the Cobre data:
#'
#' P <- cobre$cors[[1]]
#' M1 <- get_affine_m1(P, data = cobre, var.name = "cors")
#'
#'
get_affine_m1 <- function(P, data, var.name = NULL, P.sq = NULL, cores = 1, ...){
  # Translate the data to the tangent space at P
  # This effectively removed P from the data
  # Think of it as centering on a manifold

  ## Make the data into a list
  tmp <- check_right_input(data, var.name)

  data <- tmp$x
  rm(tmp)

  if(is.null(P.sq)){
    P.sq <- get_sqrt(P)
  }
  logged_data <- parallel::mclapply(data,
                                    function(X){
                                      affine_map(P = P, W = X, FUN = "log", P.sq = P.sq, ...)
                                    },
                                    mc.cores = cores)

  # Calculate the first moment estimator, which is the Euclidean mean in this
  # tangent space.
  m1 <- 1/length(data) * Reduce("+", logged_data)

  return(m1)
}

#' affine_norm
#'
#' This will calculate the norm of a vector W in the tangent space of P.  This is used directly to
#' calculate the Frechet mean.
#'
#' @param P Matrix: the base point in the manifold
#' @param W Matrix: The tangent vector at P you wish to calculate the norm of.
#' @param P.sq Matrix: the square root of P, to speed calculations
#'
#' @return A number value of the norm of W at P.
#'
#' @details This performs the calculation
#'
#' ||W||_P = ||P^{-1/2} W P^{-1/2}||_I
#'
#' Which is calculated as the Frobenius norm of P^{-1/2} W P^{-1/2}
#'
#' @export affine_norm
#'
#' @examples
#'
#' ## We will push a vector up to the tangent space of P with the logarithm map
#' ## and then calculate the norm.
#'
#' P <- cobre$cors[[1]]
#' W <- cobre$cors[[2]]
#'
#' ## Push W to a tangent vector using the Logarithm
#'
#' W_vector <- affine_map(P, W, method = "log")
#'
#' ## Now get the norm of W_vector
#'
#' affine_norm(P, W_vector)
#'
#' ## Note, this method of pushing W up, and then getting the norm, is equivalent to calculating
#' ## the affine invariant distance
#'
#' affine_dist(slice(cobre, 1:2), "cors")
#'
affine_norm <- function(P, W, P.sq = NULL){

  # Norm given by ||W||_P = ||P^{-1/2}W P^{-1/2}||_I, which is the Frobenius norm.
  if(is.null(P.sq)){
    P.sq <- get_sqrt(P)
  }

  # Translate W to the tangent space at the identity
  prod <- translate_to_identity(P, W, P.sq = P.sq)

  # Calculate the Frobenius norm at the identity
  aff_norm <- sqrt(sum((prod)^2))

  return(aff_norm)
}

#' get_frechet_mean
#'
#' This will apply a Gauss-Newton gradient decent method to calculate the Frechet mean
#' on the cone of positive definite matrices with the affine-invariant metric for a given set
#' of data.
#'
#' @param data List, array, or dataframe: The data you want to find the mean of.
#' @param var.name Character: The name of the variable containing the correlation matrices. Defaults to NULL.
#' @param tol Numeric: The tolerance you want to work at for the mean calculation. Defaults to 1e-8.
#' @param nrep Integer: The total number of iterations before you stop searching for the mean. Defaults to 1000.
#' @param cores Integer > 0: This is set up for parallelisation on Mac through mclapply.
#' This will determine the number of cores to use. Defaults to 1.
#' @param symmetric Logical: Should the symmetric square root be used, or the Cholesky.
#' @param ... Optional arguments to be passed onto affine_map
#'     - convergence_properties : return convergence properties for algorithm
#'
#' @return
#' - Fmean: A matrix that is an estimation of the Frechet mean of the PSD matrices.
#' - error: The error achieved in the algorithm
#' - num_iter: The numer of iterations to convergence
#'
#' @details This function implements the gradient decent method on page 97 of
#' "Riemannian Geometric Statistics in Medical Image Analysis."  It aims to find
#'
#'     F = arg min_{F} \sum d(F, P_i).
#'
#' @export get_frechet_mean
#'
#' @examples
#'
#' ## Calculate the Frechet mean of the COBRE data.  This will be slow (take about 2 - 3 minutes)
#'
#' F <- get_frechet_mean(cobre, "cors")
#'
#' corrplot::corrplot(F, is.corr = FALSE)
#'
get_frechet_mean <- function(data, var.name = NULL, tol = 1e-8, nrep = 1000, cores = 1, symmetric = TRUE, ...){

  ## Make the data into a list
  tmp <- check_right_input(data, var.name)

  data <- tmp$x
  rm(tmp)

  # Initiate the algorithm with the first data point as the mean estimate
  # Calculate the first moment estimate
  # Start the error at the norm of the first moment in the tangent space
  # Start the timer at 0
  # Start our step size at 1/2
  # P.bar <- data[[1]]
  P.bar <- Reduce("+", data)/length(data)

  P.sq <- get_sqrt(P.bar, symmetric = symmetric)

  M1 <- get_affine_m1(P.bar, data, P.sq = P.sq, cores = cores, ...)
  e <- affine_norm(P.bar, M1, P.sq = P.sq)
  t <- 0

  tau <- 1/2

  # Perform Gauss-Newton search for mean
  while(isTRUE(e > tol & t < nrep)){ # terminate if error gets small or too many iterations
    # increase the timer
    t <- t + 1

    # Calculate a new mean by pushing along the geodesic
    tmp.P.bar <- affine_map(P.bar, 2*tau*M1, FUN = "exp", P.sq = P.sq, ...)
    P.sq <- get_sqrt(tmp.P.bar, symmetric = symmetric)
    # Recalculate the First moment
    M1 <- get_affine_m1(tmp.P.bar, data, P.sq = P.sq, cores = cores, ...)

    # Recalculate the error
    tmp.e <- affine_norm(tmp.P.bar, M1, P.sq = P.sq)

    if(isTRUE(tmp.e > e)){ # If the error has increased, then make your step size smaller
      tau <- tau/2
      P.sq <- get_sqrt(P.bar)
    }else{ # Otherwise, update error and mean estimate
      e <- tmp.e
      P.bar <- tmp.P.bar
    }
  }
  # Return mean, error, and timer.
  ARGS <- list(...)
  if(isTRUE(ARGS$convergence_properties)){
    return(
      list(
        Fmean = P.bar,
        error = e,
        num_iter = t
      )
    )
  }else{
    return(
      P.bar
    )
  }

}

#' vec
#'
#' This will perform the Vec function that is an explicit isometry between the tangent space
#' at P and Euclidean space.
#'
#' @param P Matrix: The base point in the manifold.
#' @param W Matrix: The tangent vector to vectorise.
#' @param chol Logical: Should the Cholesky square root be used?
#' @param P.sq Matrix: The square root of P to speed calculation. Optional.
#'
#' @return A vector of length n(n+1)/2
#'
#' @details First suppose P is the identity matrix.  Then the Vec operation is:
#'
#'      Vec_I(W) = c(diag(W), \sqrt{2} * upper.tri(W)).
#'
#' Now if P is not the identity, we translate to the identity such that
#'
#'      Vec_P(W) = Vec_I(P^{-1/2} W P^{-1/2}).
#'
#' @export vec
#'
#' @examples
#'
#' ## First grab some matrices
#'
#' P <- cobre$cors[[1]]
#' W <- cobre$cors[[2]]
#'
#' ## Push W to a tangent vector using the Logarithm
#'
#' W_v <- affine_map(P, W, method = "log")
#'
#' ## Now vectorise W_v
#'
#' (W_vec <- vec(P, W_v))
#'
#' ## Wow, what a large vector.
#' ## Now use unvec to undo this an make sure it worked:
#'
#' W_back <- unvec(P, W_vec)
#'
#' ## Check the same, up to rounding
#'
#' identical(zapsmall(W_v), zapsmall(W_back))
#'
vec <- function(P, W, chol = FALSE, P.sq = NULL){

  # If the matrices are names, then preserve this naming.
  if(!is.null(colnames(W))){
    name_matrix <- purrr::map2(colnames(W), purrr::map(1:ncol(W), ~colnames(W)),
                               str_c, sep = "-") %>%
      simplify2array() %>%
      t()
  }

  if(is.null(P.sq)){ # If square root supplied, use it, otherwise just translate normally
    W_mod <- translate_to_identity(P, W, chol = chol)
  }else{
    W_mod <- translate_to_identity(P, W, chol = chol, P.sq = P.sq)
  }

  # Apply the Vec transformation
  # Diagonal elements go in first slots, then upper triangle time sqrt(2)
  w1 <- diag(W_mod)
  w2 <- sqrt(2) * W_mod[upper.tri(W_mod)]

  # Preserve naming if possible
  if(!is.null(colnames(W))){
    names(w1) <- diag(name_matrix)
    names(w2) <- name_matrix[upper.tri(name_matrix)]
  }
  return(c(w1, w2))
}

#' unvec
#'
#' This is the reverse transformation of the vec function.
#'
#' @param P Matrix: The base point in the manifold.
#' @param vec Vector: THe vecotr to translate back to the tangent space at P.  Nust be length n(n+1)/2.
#' @param chol Logical: Should the Cholesky square root be used?
#' @param P.sq Matrix: The square root of P to speed calculation. Optional.
#'
#' @return A Matrix of size n x n.
#'
#' @details This will explicitly undo the vec function.  It will construct a symmetric matrix W such that
#'
#' diag(W) = vec[1:n]
#'
#' and
#'
#' upper.tri(W) = 1/\sqrt{2} vec[n:n(n+1)/2].
#'
#' @export unvec
#'
#' @examples
#'
#' ## First grab some matrices
#'
#' P <- cobre$cors[[1]]
#' W <- cobre$cors[[2]]
#'
#' ## Push W to a tangent vector using the Logarithm
#'
#' W_v <- affine_map(P, W, method = "log")
#'
#' ## Now vectorise W_v
#'
#' W_vec <- vec(P, W_v)
#'
#' ## Now use unvec to undo this an make sure it worked:
#' ## Definitely back to a matrix!
#'
#' (W_back <- unvec(P, W_vec))
#'
#' ## Check the same, up to rounding
#'
#' identical(zapsmall(W_v), zapsmall(W_back))
#'
unvec <- function(P, vec, chol = FALSE, P.sq = NULL){

  # Find how long the diagonal is
  n <- ncol(P)

  # If names, preserve the naming
  if(!is.null(names(vec))){
    col_names <- names(vec)[1:n]
    col_names <- stringr::str_remove(col_names,  "^[^-]*-")
  }

  # Start with a null matrix
  W_id <- matrix(0, n, n)
  # Undo the vec operation
  W_id[upper.tri(W_id)] <- 1/sqrt(2) * vec[(n+1):length(vec)]
  # Make symmetric
  W_id <- W_id + t(W_id)
  diag(W_id) <- vec[1:n]

  # Translate back to tangent space at P
  if(is.null(P.sq)){
    W <- translate_to_identity(P, W_id, chol = chol, reverse = TRUE)
  }else{
    W <- translate_to_identity(P, W_id, chol = chol, reverse = TRUE, P.sq = P.sq)
  }

  if(!is.null(names(vec))){
    colnames(W) <- rownames(W) <- col_names
  }


  return(W)
}

#' affine_inner_product
#'
#' @param P Matrix - the base point on the manifold.  Must be symmetric positive definite.
#' @param U Matrix - A symmetric matrix the same dimension as P
#' @param V Matrix - A symmetric matrix the same dimension as P
#' @param P.sq Matrix - If the square root of P has been calculated, you can input it for use here.
#'
#' @return A scalar inner product value
#'
#' @details This function calculates the affine invariant inner product value at P as
#'
#'      <U, V>_P = tr(U P^{-1} V P^{-1})
#'
#' @export affine_inner_product
#'
#' @examples
#'
#' set.seed(2021)
#'
#' P <- diag(1, 3)
#' A <- matrix(rnorm(9), 3, 3)
#' A <- (A + t(A))/2
#' B <- matrix(rnorm(9), 3, 3)
#' B <- (B + t(B))/2
#'
#' affine_inner_product(P = P, U = A, V = B)
#'
affine_inner_product <- function(P, U, V, P.sq = NULL){

  if(!isSymmetric(P, tol = sqrt(.Machine$double.eps))){
    stop("P must be a symmetric matrix")
  }

  if(!is_positive_definite(P)){
    stop("P must be positive definite, at least one eigenvalue is negative.")
  }

  if(is.null(P.sq)){
    P.sq <- get_sqrt(P)
  }

  # Translate W to the tangent space at the identity
  prod_U <- translate_to_identity(P, U, P.sq = P.sq)
  prod_V <- translate_to_identity(P, V, P.sq = P.sq)

  # Calculate the Frobenius norm at the identity
  aff_inner_prod <- sum(prod_U * prod_V)

  return(aff_inner_prod)

}

