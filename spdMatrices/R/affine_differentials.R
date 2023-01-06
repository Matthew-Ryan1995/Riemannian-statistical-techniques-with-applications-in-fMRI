
#' get_theta_matrix
#'
#' Calcualte the Theta matrix used to get the derivative of the matrix exponential.
#'
#' @param D Vector - the vector containing your eigenvalues.
#' @param t Numeric - The time component for the matrix derivative.
#' @param tol Numeric - The tolerance to decide whether eigenvalues are "close"
#'
#' @return The Theta matrix used to calculate the matrix derivative.
#'
#' @details Returns the matrix \Theta such that
#'
#' \theta_{ij} = \sinh(t(\lambda_i - \lambda_j)/2)/(t(\lambda_i - \lambda_j)/2)
#'
#' if \lambda_i \neq \lambda_j and 1 else.
#'
#' @export get_theta_matrix
#'
#' @examples
get_theta_matrix <- function(D, t = 1, tol = 1e-8){

  if(!is.numeric(D)){
    stop("D should be a numeric vector of eigenvalues.")
  }

  ## Convert D into a matrix
  D_mat <- matrix(D, length(D), length(D), byrow = TRUE)

  ## Use the vectorised form definition of the theta matrix
  Theta <- sinh(t * (D_mat - t(D_mat))/2)/(t * (D_mat - t(D_mat))/2)

  ## Change close eigen values value to be 1
  Theta[abs(D_mat - t(D_mat)) < tol] <- 1

  return(Theta)
}

#' get_omega_matrix
#'
#' Calculate the Omega matrix used to get the derivative of the matrix logarithm
#'
#' @param D - Vector - the vector containing your eigenvalues.
#' @param tol Numeric - The tolerance to decide whether eigenvalues are "close."
#'
#' @return A square matrix Omega used to calculate the derivative of the logarithm
#'
#' @details Calculates the square matrix Omega such that
#'
#' omega_{ij} = omega_{ji} = (log(lambda_i) - log(lambda_j))/(lambda_i - lambda_j) if lambda_i \neq lambda_j and 1/lambda_i else.
#'
#' @export get_omega_matrix
#'
#' @examples
get_omega_matrix <- function(D, tol = 1e-8){

  if(!is.numeric(D)){
    stop("D should be a numeric vector of eigenvalues.")
  }

  ## Convert D into a matrix
  D_mat <- matrix(D, length(D), length(D), byrow = TRUE)

  # Find which eigen values are "too close"
  idx <- apply(which(abs(D_mat - t(D_mat)) < tol, arr.ind = T),
               1,
               min)

  ## Use the vectorised form definition of the theta matrix
  Omega <- (log(D_mat) - log(t(D_mat)))/((D_mat - t(D_mat)))

  ## Change close eigen values value to be 1
  Omega[abs(D_mat - t(D_mat)) < tol] <- 1/D[idx]

  return(Omega)
}


#' exponential_derivative
#'
#' Calculate the exponential derivative for a semi-simple matrix as in Najfeld and Havel (1995)
#'
#' @param A Matrix: The evaluation point for the matrix derivative
#' @param V Matrix: The direction in which to take the directional derivative
#' @param t Numeric: The time at which to take the derivative
#' @param symm Logical: Is A a symmetric matrix or not
#' @param tol Numeric: Tolerance at which we determine if eigenvalues are "close"
#'
#' @return A matrix that is the directional derivative at A in the direction V
#'
#' @details Calculates the directional derivative of Exp(tA) in the direction V.  Specifically, it calculates
#'
#' D_v(t, A) = U (\bar{V} \odot \Theta(t)) U^{-1}
#'
#' where \bar{V} = U^{-1} V U, \odot represents the Hadamard product, and \Theta(t) is the output of get_theta_matrix.
#'
#' @export exponential_derivative
#'
#' @examples
#'
#' set.seed(2021)
#' A <- matrix(rnorm(3 * 3), 3, 3)
#' A <- (A + t(A))/2
#'
#' V <- matrix(rnorm(3 * 3), 3, 3)
#' V <- (V + t(V))/2
#'
#' # Make sure exponential derivative matches finite difference
#'
#' h <- 1e-8
#' (expm::expm(A + h*V) - expm::expm(A))/h
#' exponential_derivative(A, V)
#'
exponential_derivative <- function(A, V, t = 1, symm = TRUE, tol = 1e-8){

  # If performing symmetric operation, A and V must both be symmetric
  if(symm){
    if(!isSymmetric(A, tol = sqrt(.Machine$double.eps))){
      stop("A must be symmetric")
    }
    # if(!isSymmetric(V)){
    #   stop("V must be symmetric")
    # }
  }

  # Check that t is numeric

  if(!is.numeric(t)){
    stop("t must be a real number.")
  }

  # Find the Eigen-decomposition of the evaluation point P
  A_eig <- eigen(A, symmetric = symm)
  U <- A_eig$vectors

  ## Check that A is diagonalisable
  tmp <- tryCatch(
    expr = solve(U),
    error = function(e) "ERROR"
  )
  if("character" %in% class(tmp)){
    stop("A must be a diagonalisable matrix.")
  }

  lam <- A_eig$values
  # Remove waste variables
  rm(A_eig, tmp)

  if(symm){
    # If we are dealing with symmetric matrices, use the fact that U is orthogonal
    Vbar <- crossprod(U, crossprod(V, U))
  }else{
    Vbar <- solve(U) %*% V %*% (U)
  }

  # Calculate the theta matrix
  Theta <- get_theta_matrix(lam, t, tol = tol)

  # Calculate exp(t\Lambda/2)
  Lamda_exp <- diag(exp(t * lam/2))

  # Multiply on the left and the right by Lambda_exp
  Inner_prod <- translate_to_identity(Lamda_exp, Vbar * Theta, P.sq = Lamda_exp, reverse = TRUE)

  if(symm){
    # If we are dealing with symmetric matrices, use the fact that U is orthogonal
    Outer_prod <- tcrossprod(tcrossprod(U, Inner_prod), U)
  }else{
    U %*% Inner_prod %*% solve(U)
  }

  return(t * Outer_prod)
}


#' exponential_derivative_general
#'
#' Apply the chain rule to calculate the derivative of the affine exponential map
#'
#' @param P Matrix: a symmetric positive definite matrix
#' @param A Matrix: The evaluation point in the tangent space at P
#' @param V Matrix: The direction in the tangent space at P to take the derivative in
#' @param t Numeric: The time point at which to take the derivative.
#' @param symm Logical: Whether the matrix A is symmetric or not
#' @param P.sq Matrix: Optional square root of P to speed up calculation time
#' @param tol Numeric: Tolerance at which we determine if eigenvalues are "close"
#'
#' @return The directional derivative of Exp_P(tA) in the direction V
#' @export exponential_derivative_general
#'
#' @examples
#'
#' set.seed(2021)
#' P <- matrix(rnorm(3 * 3), 3, 3)
#' P <- (P + t(P))/2
#' P <- expm::expm(P) # Make sure P is positive definite
#'
#' A <- matrix(rnorm(3 * 3), 3, 3)
#' A <- (A + t(A))/2
#'
#' V <- matrix(rnorm(3 * 3), 3, 3)
#' V <- (V + t(V))/2
#'
#' # Make sure exponential derivative matches finite difference
#'
#' h <- 1e-8
#' (affine_map(P, A + h*V) - affine_map(P, A))/h
#' exponential_derivative_general(P, A, V)
#'
exponential_derivative_general <- function(P, A, V, t = 1, symm = TRUE, P.sq = NULL, tol = 1e-8){

  # P must be symmetric and positive definite
  if(!isSymmetric(P, tol = sqrt(.Machine$double.eps))){
    stop("P must be a symmetric matrix")
  }

  if(!is_positive_definite(P)){
    stop("P must be positive definite, at least one eigenvalue is negative.")
  }

  # If performing symmetric operation, A and V must both be symmetric
  if(symm){
    if(!isSymmetric(A, tol = sqrt(.Machine$double.eps))){
      stop("A must be symmetric")
    }
    # if(!isSymmetric(V)){
    #   stop("V must be symmetric")
    # }
  }

  # Check that t is numeric

  if(!is.numeric(t)){
    stop("t must be a real number.")
  }

  # Calculate the square root of our base point
  if(is.null(P.sq)){
    P.sq <- get_sqrt(P)
  }

  # Translate to the identity from P
  A_trans <- translate_to_identity(P, A, P.sq = P.sq)
  V_trans <- translate_to_identity(P, V, P.sq = P.sq)

  # Calculate the exponential derivative
  D_Id <- exponential_derivative(A_trans, V_trans, t, symm = symm, tol = tol)

  # Translate back to P.
  D_general <- translate_to_identity(P, D_Id, P.sq = P.sq, reverse = TRUE)

  return(D_general)
}

#' exponential_adjoint
#'
#' Calculate the adjoint of the derivative of the matrix exponential map
#'
#' @param A Matrix: Evaluation point.  Must be symmetric.
#' @param V Matrix: Directional point.  Must be symmetric
#' @param t Numeric: Time to evaluate at, i.e. exp(tA)
#' @param tol Numeric: Tolerance at which we determine if eigenvalues are "close"
#' @param ... Optional arguments
#'     -  method : The method to pass to expm or logm.  Only accepts Eigen or Higham08
#'
#'
#' @return A matrix that is the adjoint of d_A Exp applied to V.
#'
#' @details Calculates the adjoint of d_A Exp as
#'
#' (d_A Exp)^\dagger (V) = Exp(A)^{-1/2} (d_A Exp(Exp(A^{-1/2}) V Exp(A)^{-1/2})) Exp(A^{1/2}).
#'
#' This has the property that
#'
#' < d_AExp(W), V >_{Exp(A)} = < W, (d_A Exp)^\dagger(V) >_I
#'
#' For all symmetric matrices W and V.
#'
#' @export exponential_adjoint
#'
#' @examples
exponential_adjoint <- function(A, V, t = 1, tol = 1e-8, ...){
  # A and V must both be symmetric
  if(!isSymmetric(A, tol = sqrt(.Machine$double.eps))){
    stop("A must be symmetric")
  }
  # if(!isSymmetric(V)){
  #   stop("V must be symmetric")
  # }

  # Check that t is numeric

  if(!is.numeric(t)){
    stop("t must be a real number.")
  }

  ARGS <- list(...)

  if(is.null(ARGS$method)){
    method <- "Higham08"
  }else{
    method <- ARGS$method
    if(!(method %in% c("Higham08", "Eigen"))){
      stop("Method must be Higham08 or Eigen")
    }
  }

  if(method == "Eigen"){
    method <- str_c("R_", method)
  }

  E <- expm::expm(t * A, method = method)

  # Move A and V to the identity from P
  V_trans <- translate_to_identity(E, V)

  # Calculate the derivative at the identity, ending in A
  D_V <- exponential_derivative(A, V_trans, t = t, tol = tol)

  # Translate back to the identity
  D_V_adjoint <- translate_to_identity(E, D_V,  reverse = FALSE)

  if(any(is.na(D_V_adjoint)) & method != "R_Eigen"){
    return(exponential_adjoint(A, V, t, tol, method = "Eigen"))
  }

  return(D_V_adjoint)

}

#' exponential_adjoint_general
#'
#' Applies the chain rule to calculate the adjoint of the affine Riemannian exponential Exp_P(A) = P^{1/2} Exp(P^{-1/2} A P^{-1/2}) P^{1/2}.
#'
#' @param P Matrix: Base point.  Must be symmetric positive definite.
#' @param A Matrix: Evaluation point.  Must be symmetric.
#' @param V Matrix: Directional point.  Must be symmetric.
#' @param t Numeric: Time to evaluate at, i.e. exp(tA)
#' @param P.sq Matrix: Optional input of swaure root of P to speed up calculations.
#' @param tol Numeric: Tolerance at which we determine if eigenvalues are "close"
#' @param ... Optional argument to be passed to exponential_adjoint
#'
#' @return A matrix that is the adjoint of d_A Exp_P applied to V.
#'
#' @details Calculates the adjoint of d_A Exp_P as
#'
#' (d_A Exp_P)^\dagger (V) = P^{1/2} ((d_{P^{-1/2} A P^{-1/2}} Exp)^\dagger(P^{-1/2} V P^{-1/2}) P^{1/2}
#'
#' This has the property that
#'
#' < d_AExp_P(W), V >_{Exp_P(A)} = < W, (d_A Exp_P)^\dagger(V) >_{P}
#'
#' For all symmetric matrices W and V.
#'
#' @export exponential_adjoint_general
#'
#' @examples
exponential_adjoint_general <- function(P, A, V, t = 1, P.sq = NULL, tol = 1e-8, ...){

  # P must be symmetric and positive definite
  if(!isSymmetric(P, tol = sqrt(.Machine$double.eps))){
    stop("P must be a symmetric matrix")
  }

  if(!is_positive_definite(P)){
    stop("P must be positive definite, at least one eigenvalue is negative.")
  }

  # A and V must both be symmetric
    if(!isSymmetric(A, tol = sqrt(.Machine$double.eps))){
      stop("A must be symmetric")
    }
    # if(!isSymmetric(V)){
    #   stop("V must be symmetric")
    # }

  # Check that t is numeric

  if(!is.numeric(t)){
    stop("t must be a real number.")
  }

  # If the square root of P is not given, calculate it
  if(is.null(P.sq)){
    P.sq <- get_sqrt(P)
  }

  # Move A and V to the identity from P
  A_trans <- translate_to_identity(P, A, P.sq = P.sq)
  V_trans <- translate_to_identity(P, V, P.sq = P.sq)

  # Calculate the derivative at the identity, ending in A
  D_V <- exponential_adjoint(A_trans, V_trans, t = t, tol = tol, ...)

  # Translate back to the identity
  D_V_adjoint <- translate_to_identity(P, D_V, P.sq = P.sq, reverse = TRUE)

  return(D_V_adjoint)
}

#' log_derivative
#'
#' Calculates the differential of the matrix logarithm at A in the direction V.
#'
#' @param A  Matrix: Evaluation point.  Must be positive definite.
#' @param V  Matrix: Directional point.  Should be symmetric.
#' @param symm Logical: Is V a symmetric matrix.
#' @param tol Numeric: Tolerance at which we determine if eigenvalues are "close"
#'
#' @return A matrix of the same dimension as A that is the directional derivative of Log(A), the matrix logarithm, in the direction V.
#'
#' @details The differential d_A Log is found by using the fact that d_{Log(A)} Exp \circ d_A Log (V) = V.  This leverages the fact that
#' Log is the inverse of Exp, hence the differential of the composition is the identity mapping.
#'
#' @export log_derivative
#'
#' @examples
log_derivative <- function(A, V, symm = TRUE, tol = 1e-8){

  # A must be symmetric and positive definite
  # if(!isSymmetric(A)){
  #   stop("A must be a symmetric matrix")
  # }

  if(!isSymmetric(A, tol = sqrt(.Machine$double.eps))){
    stop("A must be a symmetric matrix")
    # warning("A wasn't symmetric, it has been adjusted")
    # A <- (A + t(A))/2
  }


  if(!is_positive_definite(A)){
    stop("A must be positive definite, at least one eigenvalue is negative.")
  }

  # # If performing symmetric operation, A and V must both be symmetric
  # if(symm){
  #   if(!isSymmetric(V)){
  #     stop("V must be symmetric")
  #   }
  # }

  # Find the Eigen-decomposition of the evaluation point P
  A_eig <- eigen(A, symmetric = symm)
  U <- A_eig$vectors

  ## Check that A is diagonalisable
  tmp <- tryCatch(
    expr = solve(U),
    error = function(e) "ERROR"
  )
  if("character" %in% class(tmp)){
    stop("A must be a diagonalisable matrix.")
  }

  lam <- A_eig$values
  # Remove waste variables
  rm(A_eig, tmp)

  #
  if(symm){# If we are dealing with symmetric matrices, use the fact that U is orthogonal
    Vbar <- crossprod(U, crossprod(V, U))
  }else{
    Vbar <- solve(U) %*% V %*% (U)
  }

  Omega <- get_omega_matrix(lam, tol = tol)

  Inner_prod <- Vbar * Omega

  if(symm){# If we are dealing with symmetric matrices, use the fact that U is orthogonal
    Outer_prod <- tcrossprod(tcrossprod(U, Inner_prod), U)
  }else{
    Outer_prod <- U %*% Inner_prod %*% solve(U)
  }

  return(Outer_prod)
}

#' log_derivative_general
#'
#' Applies the chain rule to calculate the derivative of Log_P(A) in the direction V.
#'
#' @param P  Matrix: Base point.  Must be symmetric positive definite.
#' @param A  Matrix: Evaluation point.  Must be symmetric positive definite.
#' @param V  Matrix: Directional point.  Should be symmetric.
#' @param symm Logical: Is V a symmetric matrix.
#' @param tol Numeric: Tolerance at which we determine if eigenvalues are "close"
#'
#' @return A matrix of the same dimension as A that is the directional derivative of Log_P(A), the matrix logarithm, in the direction V.
#'
#' @details Uses the chain rule to calculate that
#'
#' d_ALog_P(V) = P^{1/2} (d_{P^{-1/2} A P^{-1/2}} Log (P^{-1/2} V P^{-1/2})) P^{1/2}.
#'
#' @export log_derivative_general
#'
#' @examples
log_derivative_general <- function(P, A, V, symm = TRUE, tol = 1e-8){

  # P must be symmetric and positive definite
  if(!isSymmetric(P, tol = sqrt(.Machine$double.eps))){
    stop("P must be a symmetric matrix")
  }

  if(!is_positive_definite(P)){
    stop("P must be positive definite, at least one eigenvalue is negative.")
  }
  # A must be symmetric and positive definite
  if(!isSymmetric(A, tol = sqrt(.Machine$double.eps))){
    stop("A must be a symmetric matrix")
  }

  if(!is_positive_definite(A)){
    stop("A must be positive definite, at least one eigenvalue is negative.")
  }

  # # If performing symmetric operation, A and V must both be symmetric
  # if(symm){
  #   if(!isSymmetric(V)){
  #     stop("V must be symmetric")
  #   }
  # }

  # Translate Our evaluation and direction points back to the identity by P
  V_trans <- translate_to_identity(P, V)
  A_trans <- translate_to_identity(P, A)

  # Calculate the log derivative at the translated points
  log_der <- log_derivative(A_trans, V_trans, symm = symm, tol = tol)

  # Translate back to the tangent space of P
  final_dir <- translate_to_identity(P, log_der, reverse = TRUE)

  return(final_dir)
}

#' log_adjoint
#'
#' Calculates the adjoint of the linear map d_ALog, where Log is the matrix logarithm.
#'
#' @param A Matrix: Evaluation point.  Must be symmetric positive definite.
#' @param V Matrix: Directional point.  Should be symmetric.
#' @param symm Logical: Is V symmetric.
#' @param A.sq Matrix: Optional square root of A to speed up calculations.
#' @param tol Numeric: Tolerance to which eigenvalues are "the same."
#'
#' @return A matrix of the same dimension as A which is the adjoint of the linear map d_ALog.
#'
#' @details Calculates
#'
#' (d_A Log)^\dagger(V) = A^{1/2} (d_{A}Log(A^{1/2} V A^{1/2}))A^{1/2}
#'
#' Such that
#'
#' < d_A Log(V), W >_I = < V, (d_A Log)^\dagger(W) >_A
#'
#' @export log_adjoint
#'
#' @examples
log_adjoint <- function(A, V, symm = TRUE, A.sq = NULL, tol = 1e-8){

  # A must be symmetric and positive definite
  if(!isSymmetric(A, tol = sqrt(.Machine$double.eps))){
    stop("A must be a symmetric matrix")
    # warning("A wasn't symmetric, it has been adjusted")
    # A <- (A + t(A))/2
  }

  if(!is_positive_definite(A)){
    stop("A must be positive definite, at least one eigenvalue is negative.")
  }

  # # If performing symmetric operation, A and V must both be symmetric
  # if(symm){
  #   if(!isSymmetric(V)){
  #     stop("V must be symmetric")
  #   }
  # }

  if(is.null(A.sq)){
    A.sq <- get_sqrt(A, symmetric = symm)
  }

  # Move V into the tangent space at A
  V_trans <- translate_to_identity(A, V, P.sq = A.sq, reverse = TRUE)

  # Take the log derivative
  inner_derivative <- log_derivative(A, V_trans, tol = tol)

  # Translate back to the tangent space at A
  outer_adjoint <-  translate_to_identity(A, inner_derivative, P.sq = A.sq, reverse = TRUE)

  return(outer_adjoint)
}

#' log_adjoint_general
#'
#' @param P Matrix: The base point.  Must be symmetric positive definite.
#' @param A Matrix: The evaluation point.  Must be symmetric positive definite.
#' @param V Matrix: The directional point.  Should be symmetric.
#' @param P.sq Matrix: Optional square root of P to speed up calculation.
#' @param symm Logical: Is V symmetric?
#' @param tol Numeric: Tolerance to which eigenvalues are "the same."
#'
#' @return A matrix of the same dimension as A which is the adjoint of the linear map d_ALog_P.
#'
#' @details Calculates the adjoint of d_A Log_P as
#'
#' (d_A Log_P)^\dagger (V) = P^{1/2} ((d_{P^{-1/2} A P^{-1/2}} Log)^\dagger(P^{-1/2} V P^{-1/2}) P^{1/2}
#'
#' This has the property that
#'
#' < d_A Log_P(W), V >_{P} = < W, (d_A Log_P)^\dagger(V) >_{A}
#'
#' For all symmetric matrices W and V.
#'
#' @export log_adjoint_general
#'
#' @examples
log_adjoint_general <- function(P, A, V, P.sq = NULL, symm = TRUE, tol = 1e-8){

  # P must be symmetric and positive definite
  if(!isSymmetric(P, tol = sqrt(.Machine$double.eps))){
    stop("P must be a symmetric matrix")
  }

  if(!is_positive_definite(P)){
    stop("P must be positive definite, at least one eigenvalue is negative.")
  }

  # A must be symmetric and positive definite
  if(!isSymmetric(A, tol = sqrt(.Machine$double.eps))){
    stop("A must be a symmetric matrix")
  }

  if(!is_positive_definite(A)){
    stop("A must be positive definite, at least one eigenvalue is negative.")
  }

  # If performing symmetric operation, A and V must both be symmetric
  # if(symm){
  #   if(!isSymmetric(V)){
  #     stop("V must be symmetric")
  #   }
  # }

  if(is.null(P.sq)){
    P.sq <- get_sqrt(P, symmetric = symm)
  }

  # Translate to the identity via P
  V_trans <- translate_to_identity(P, V, P.sq = P.sq, reverse = FALSE)
  A_trans <- translate_to_identity(P, A, P.sq = P.sq, reverse = FALSE)

  # Calculate adjoint at the identity
  inner_adjoint <- log_adjoint(A_trans, V_trans, symm = symm, tol = tol)

  # Translate back to P
  outer_adjoint <- translate_to_identity(P, inner_adjoint, P.sq = P.sq, reverse = TRUE)

  return(outer_adjoint)
}
