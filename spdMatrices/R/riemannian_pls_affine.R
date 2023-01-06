#' Gradient of the score, affine variant
#'
#' Calculates the gradient for calculating the score when dealing with the affine-invariant geometry
#'
#' @param mu Matrix: Positive definite symmetric matrix, mean of the X.
#' @param W Matrix: Current estimate of weights matrix W.
#' @param X Matrix: Data matrix, i.e. correlation matrix.
#' @param t Numeric: The current estimate of the projection coefficient t.
#' @param ...
#'     - check: logical, should inputs be checked.
#'
#' @return A value that is the gradient for t at X and W.
#' @export gradient_score_affine
#'
#' @examples
#'
#' TBA
gradient_score_affine <- function(mu, W, X, t = 1, ...){

  ARGS <- list(...)
  if(isTRUE(ARGS$check)){
    # X and mu must be symmetric and positive definite
    if(!isSymmetric(mu, tol = sqrt(.Machine$double.eps))){
      stop("mu must be a symmetric matrix")
    }
    if(!isSymmetric(W, tol = sqrt(.Machine$double.eps))){
      stop("W must be symmetric")
    }
    if(!isSymmetric(X, tol = sqrt(.Machine$double.eps))){
      stop("X must be symmetric")
    }

    if(!is_positive_definite(X)){
      stop("X must be positive definite, .at least one eigenvalue is negative.")
    }
    if(!is_positive_definite(mu)){
      stop("mu must be positive definite, at least one eigenvalue is negative.")
    }

    # Check that t is numeric

    if(!is.numeric(t)){
      stop("t must be a real number.")
    }
  }

  # Calculate our S matrix, the exponential at t*W
  S <- affine_map(mu, t * W)
  # Calculate our log matrix, Log_X(S)
  L1 <- affine_map(X, S, FUN = "log")

  # Return our gradient
  return(
    2 * get_trace(solve(X) %*% L1 %*% solve(mu) %*% W)
  )
}

#' Gradient of direction, affine variant
#' Calculates a component of the gradient of the weights/loading vectors for the affine-invariant geometry.
#'
#' @param mu Matrix: Positive definite symmetric matrix, mean of the X.
#' @param P Matrix: Current estimate of direction vector.
#' @param X Matrix:  Data matrix, i.e correlation matrix
#' @param t Numeric: Current estimate of projection coefficient
#' @param ...
#'     - check: logical, should inputs be checked.
#'
#' @return A matrix that is the portion of gradient of P attributed to X.
#' @export gradient_direction_affine
#'
#' @examples
#'
#' TBA
gradient_direction_affine <- function(mu, P, X, t = 1, ...){

  ARGS <- list(...)
  if(isTRUE(ARGS$check)){
    # X and mu must be symmetric and positive definite
    if(!isSymmetric(mu, tol = sqrt(.Machine$double.eps))){
      stop("mu must be a symmetric matrix")
    }
    if(!isSymmetric(P, tol = sqrt(.Machine$double.eps))){
      stop("P must be symmetric")
    }
    if(!isSymmetric(X, tol = sqrt(.Machine$double.eps))){
      stop("X must be symmetric")
    }

    if(!is_positive_definite(X)){
      stop("X must be positive definite, at least one eigenvalue is negative.")
    }
    if(!is_positive_definite(mu)){
      stop("mu must be positive definite, at least one eigenvalue is negative.")
    }

    # Check that t is numeric

    if(!is.numeric(t)){
      stop("t must be a real number.")
    }
  }

  # get square root
  if(is.null(ARGS$mu.sq)){
    mu.sq <- get_sqrt(mu)
  }else{
    mu.sq <- ARGS$mu.sq
  }

  # Calculate the exponential Exp_mu(t * P)
  E1 <- affine_map(mu, t*P)

  # Calculate the "difference" between X and E1.
  L1 <- affine_map(E1, X, FUN = "log")

  # Calculate the exponential adjoint, based at mu, evaluates at p, in direction L1 at time t.
  # Similar to geodesic regression.
  adjoint_map <- exponential_adjoint_general(P = mu, A = P, V = L1, t = t, P.sq = mu.sq)

  return(
    -2 * adjoint_map
  )
}


#' Calculate the gradient of the direction, affine variant.
#'
#' Calculates the gradient for the direction (weight or loading) in the affine-invariant geometry.
#'
#' @param mu Matrix: Mean of the X.
#' @param P Matrix: Current estimate of direction matrix.
#' @param X List: List of the X matrices.
#' @param t Numeric: Numeric vector of the current t values.
#' @param ... Optional:
#'     - mc.cores: Integer, to run in parallel, how many cores.
#'     - mu.sq: Matrix, the square root of mu to speed up computation.
#'     - check: Logical, should we check the variable?
#'
#' @return The full gradient of P at X and t.
#' @export gradient_direction_affine_full
#'
#' @examples
#'
#' TBA
gradient_direction_affine_full <- function(mu, P, X, t, ...){

  ARGS <- list(...)

  # Get number of cores
  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  # get square root
  if(is.null(ARGS$mu.sq)){
    mu.sq <- get_sqrt(mu)
  }else{
    mu.sq <- ARGS$mu.sq
  }

  # Do we check our paramters?
  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }

  # Convert X and t into a paired list.
  params <- map2(X, t, ~list(X = .x, t = .y))

  # Calculate gradient value for each X and t
  direction_gradients <- parallel::mclapply(params,
                                            function(l){
                                              X <- l$X
                                              t <- l$t

                                              return(
                                                gradient_direction_affine(mu = mu, P = P, X = X, t = t,
                                                                          mu.sq = mu.sq, check = check)
                                              )
                                            },
                                            mc.cores = mc.cores)

  # Add them all up
  direction_gradient <- Reduce("+", direction_gradients)

  return(
    (t(direction_gradient) + direction_gradient)/2 # Make sure symmetric.
  )
}



#' Calculate the score for a subject, affine variant
#'
#' Using gradient descent, calculate the score for a single given subject.
#'
#' @param t Numeric: Initial estimate of t.
#' @param W Matrix: Current estimate of weight direction W.
#' @param X Matrix: Data point X.
#' @param muX Matrix: The mean of the X.
#' @param tau Numeric: Initial step size.
#' @param max.iter Numeric: Maximum number of iterations for while loop.
#' @param tol Numeric: Desired tolerance.
#' @param ...
#'     - check: Logical, should we check our variables.
#'
#' @return A new estimate of the score.
#' @export calculate_score_affine
#'
#' @examples
#' TBA
calculate_score_affine <- function(t, W, X, muX, tau = 1/2, max.iter = 50, tol = 1e-6, ...){

  # t individual,
  # X individual

  # Convergence info
  #
  ARGS <- list(...)


  # Do we check our parameters?
  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }


  # initialize values
  k <- 0
  e <- 1
  t_old <- t

  while(k < max.iter){
    # Increase counter
    k <- k + 1

    # Calculate gradient, stop if you can't
    t_grad <- tryCatch(gradient_score_affine(mu = muX, W = W, X = X, t = t_old, check = check),
                       error = function(e) "stop")
    if(identical(t_grad, "stop")){
      t_new <- t_old
      break
    }

    # updat t and error
    t_new <- t_old - tau * t_grad
    e_new <- abs(t_new - t_old)

    # If error has increased, keep decreasing tau until it decreases.
    while(e_new > e & k < max.iter){
      k <- k + 1
      tau <- tau/2

      t_new <- t_old - tau * t_grad
      e_new <- abs(t_new - t_old)
    }

    # If tau got too small, stop
    if(tau < 1e-12){
      break
    }

    if(k%%5 == 0){
      check <- TRUE
      if(isTRUE(ARGS$verbose)){
        print(str_c(
          "Current Iteration: ", k,
          "; Current error: ", round(e, 8),
          "; Current step size: ", round(tau, 8)
        ))
      }
    }else{
      check <- FALSE
    }

    # Update t and e
    t_old <- t_new
    e <- e_new

    # Stop if error small
    if(isTRUE(e < tol)){
      break
    }
  }

  if(isTRUE(ARGS$verbose)){
    print(str_c(
      "Maximum iteration: ", max.iter,
      "; Iterations completed: ", k,
      "; Final error: ", round(e, 8),
      "; Final step size: ", round(tau, 8)
    ))
  }
  return(t_new)

}

#' Calculate all scores, affine variant
#'
#' @param mu Matrix: Mean of the X.
#' @param W Matrix: Current estimate of direction W.
#' @param X List: List of the data matrices X.
#' @param t Numeric: Vector of current estimates of t.
#' @param tau Numeric: Initial Step size
#' @param max.iter Numeric: Maximum number of iterations
#' @param tol Numeric: Tolerance to work at.
#' @param ... Optional arguments
#'    - check: Logical, should we check our estimates
#'    - mc.cores: Numeric, if parallelising.
#'
#' @return An updated vector of scores.
#' @export calculate_score_affine_full
#'
#' @examples
#' TBA
calculate_score_affine_full <- function(mu, W, X, t, tau = 1/2, max.iter = 50, tol = 1e-6, ...){

  ARGS <- list(...)

  # Do we parallelise?
  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  # Do we check our parameters?
  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }

  # Set up parameters as a paired list
  params <- map2(X, t, ~list(X = .x, t = .y))

  # Parallelise over all scores
  score_gradients <- parallel::mclapply(params,
                                        function(l){
                                          X <- l$X
                                          t <- l$t

                                          return(
                                            calculate_score_affine(mu = mu, W = W, X = X, t = t,
                                                                   tau = tau, max.iter = max.iter, tol = tol,
                                                                   check = check, ...)
                                          )
                                        },
                                        mc.cores = mc.cores) %>%
    unlist()

  return(score_gradients)
}


#' Calculate the direction vector, affine invariant
#'
#' Using gradient descent, calculate the direction vector in the affine invariant geometry.
#'
#' @param mu Matrix: Mean of the X
#' @param P Matrix: Current estimate of direction vector
#' @param X List: List of data matrices X. If as data frame, provide optional argument var.name.X
#' @param t Numeric: Vector of current scores.
#' @param max.iter Numeric: Maximum number of iterations.
#' @param tol Numeric: Tolerance to work at.
#' @param ... Optional arguments
#'    - var.name.X: Character, if X is data frame, name of column with data.
#'    - mu.sq: Matrix: Square root of mu to speed up computation.
#'    - tau: Numeric, initial step size.
#'    - mc.cores: Numeric, if using parallel, number of cores.
#'    - check: Logical, should we check our inputs.
#'    - verbose: Logical, should we display progress output.
#'
#' @return A new, updated direction
#' @export calculate_direction_affine
#'
#' @examples
#' TBA
#'
calculate_direction_affine <- function(mu, P, X, t, max.iter = 50, tol = 1e-6, ...){

  ## THIS HAS BEEN CHANGED TO TEST GRADIENT 0
  ARGS <- list(...)

  # Make sure X is a lists.
  if(is.null(ARGS$var.name.X)){
    var.name.X <- NULL
  }else{
    var.name.X <- ARGS$var.name.X
  }

  X <- check_right_input(X, var.name = var.name.X)$x


  if(is.null(ARGS$mu.sq)){
    mu.sq <- get_sqrt(mu)
  }else{
    mu.sq <- ARGS$mu.sq
  }


  if(is.null(ARGS$tau)){
    tau <- 1/2
  }else{
    tau <- ARGS$tau
  }


  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }

  # Initialise counters and error
  P_old <- P
  P_grad_old <- tryCatch(gradient_direction_affine_full(mu = mu, P = P_old, X = X, t = t, mc.cores = mc.cores,
                                                        mu.sq = mu.sq, check = check),
                         error = function(e) "stop")

  # If gradient throws error, stop
  if(identical(P_grad_old, "stop")){
    tmp <- (P_old + t(P_old))/2
    warning("error in gradient")
    return(tmp)
  }

  k <- 0

  error <- affine_norm(P = mu, W = P_grad_old, P.sq = mu.sq)

  if(isTRUE(error < tol)){
    tmp <- (P_old + t(P_old))/2
    print("No iteration perfromed, solution found at initial value.")
    if(isTRUE(ARGS$convergence_results)){
      return(list(
        ANS = tmp,
        tau = tau,
        k = k,
        e = error
      ))
    }
    return(tmp)
  }



  while(k < max.iter){
    # Increase counter
    k <- k + 1

    # Update P
    P_new <- P_old - tau * P_grad_old

    # Calculate gradient
    P_grad <- tryCatch(gradient_direction_affine_full(mu = mu, P = P_new, X = X, t = t, mc.cores = mc.cores,
                                                      mu.sq = mu.sq, check = check),
                       error = function(e) "stop")

    # If gradient throws error, stop
    if(identical(P_grad, "stop")){
      warning("error in gradient")
      P_new <- P_old
      break
    }

    if(isTRUE(error < tol)){
      break
    }

    error_new <- affine_norm(P = mu, W = P_grad, P.sq = mu.sq)

    if(isTRUE(error_new > error)){
      k <- k + 1
      # P_old <- P_old + tau * P_grad_old - tau * P_grad_old / 2
      tau <- tau/2
    }else{
      P_old <- P_new
      error <- error_new
      P_grad_old <- P_grad
    }

    if(k%%5 == 0){
      check <- TRUE
      if(isTRUE(ARGS$verbose)){
        print(str_c(
          "Current Iteration: ", k,
          "; Current error: ", round(error, 8),
          "; Current step size: ", round(tau, 8)
        ))
      }
    }else{
      check <- FALSE
    }

    # # Update P
    # P_new <- P_old - tau * P_grad
    #
    #
    # # Update error
    # error_new <- affine_norm(P = mu, W = P_new - P_old, P.sq = mu.sq)
    #
    # while(error_new > error & k < max.iter){ # If error increases, decrease step size until decrease.
    #   k <- k + 1
    #   tau <- tau/2
    #
    #   P_new <- P_old - tau * P_grad
    #   error_new <- affine_norm(P = mu, W = P_new, P.sq = mu.sq)
    # }
    #
    #
    # if(tau < 1e-12){
    #   break
    # }
    #
    # P_old <- P_new
    # error <- error_new
    #
    #
    # if(k%%5 == 0){
    #   check <- TRUE
    #   if(isTRUE(ARGS$verbose)){
    #     print(str_c(
    #       "Current Iteration: ", k,
    #       "; Current error: ", round(error, 8),
    #       "; Current step size: ", round(tau, 8)
    #     ))
    #   }
    # }else{
    #   check <- FALSE
    # }
  }


  if(isTRUE(ARGS$verbose)){
    print(str_c(
      "Maximum iteration: ", max.iter,
      "; Iterations completed: ", k,
      "; Final error: ", round(error, 8),
      "; Final step size: ", round(tau, 8)
    ))
  }
  # Fix numerical impurities
  tmp <- (P_new + t(P_new))/2

  if(isTRUE(ARGS$convergence_results)){
    return(list(
      ANS = tmp,
      tau = tau,
      k = k,
      e = error
    ))
  }

  return(tmp)
}



#' Calculate the score t, affine variant
#'
#' A wrapper for easy calculation of t.
#'
#' @param mu Matrix: mean of X
#' @param W Matrix: Current estimate of W
#' @param X List: List of data matrices X
#' @param t Numeric: Numeric vector of current t scores
#' @param tau Numeric: Initial step size.
#' @param max.iter Numeric: Maximum number of iterations
#' @param tol Numeric: Tolerance to work at.
#' @param ... Optional arguments
#'    - mc.cores: Numeric, number of cores to run on
#'    - check: Logical, Should we check our inputs
#'    - verbose: Logical, Do we want to print progress output?
#'
#' @return The new scores t
#' @export calculate_t_affine
#'
#' @examples
#' TBA
calculate_t_affine <- function(mu, W, X, t, tau = 1/2, max.iter = 50, tol = 1e-6, ...){
  ARGS <- list(...)

  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }

  if(isTRUE(ARGS$verbose)){
    print("Calculating scores t:")
  }

  return(
    calculate_score_affine_full(mu = mu, W = W, X = X, t = t, tau = tau,
                                max.iter = max.iter, tol = tol, mc.cores = mc.cores,
                                check = check, ...)
  )
}


#' Calculate the score u, affine variant
#'
#' A wrapper for easy calculation of u.
#'
#' @param mu Matrix: mean of X
#' @param C Matrix: Current estimate of C
#' @param Y List: List of data matrices Y
#' @param u Numeric: Numeric vector of current u scores
#' @param tau Numeric: Initial step size.
#' @param max.iter Numeric: Maximum number of iterations
#' @param tol Numeric: Tolerance to work at.
#' @param ... Optional arguments
#'    - mc.cores: Numeric, number of cores to run on
#'    - check: Logical, Should we check our inputs
#'    - verbose: Logical, Do we want to print progress output?
#'
#' @return The new scores u
#' @export calculate_u_affine
#'
#' @examples
#' TBA
calculate_u_affine <- function(mu, C, Y, u, tau, max.iter = 50, tol = 1e-6, ...){
  ARGS <- list(...)

  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }

  if(isTRUE(ARGS$verbose)){
    print("Calculating scores u:")
  }

  return(
    calculate_score_affine_full(mu = mu, W = C, X = Y, t = u, tau = tau,
                                max.iter = max.iter, tol = tol, mc.cores = mc.cores,
                                check = check, ...)
  )
}


#' Calculate the weights W, affine variant
#'
#' A wrapper for easy calculation of W.
#'
#' @param mu Matrix: mean of X
#' @param W Matrix: Current estimate of W
#' @param X List: List of data matrices X. If X is a data frame, please use var.name.X
#' @param u Numeric: Numeric vector of current u scores
#' @param max.iter Numeric: Maximum number of iterations
#' @param tol Numeric: Tolerance to work at.
#' @param ... Optional arguments
#'    - mc.cores: Numeric, number of cores to run on
#'    - check: Logical, Should we check our inputs
#'    - verbose: Logical, Do we want to print progress output?
#'    - var.name.X: Character, name of column of data if X is data frame.
#'    - mu.sq: Matrix, square root of mu to speed computation
#'    - tau: Numeric, initial step size.
#'
#' @return The new weights W.
#' @export calculate_w_affine
#'
#' @examples
#' TBA
calculate_w_affine <- function(mu, W, X, u, max.iter = 50, tol = 1e-6, ...){

  ARGS <- list(...)

  # Make sure X is a list.
  if(is.null(ARGS$var.name.X)){
    var.name.X <- NULL
  }else{
    var.name.X <- ARGS$var.name.X
  }

  X <- check_right_input(X, var.name = var.name.X)$x

  if(is.null(ARGS$mu.sq)){
    mu.sq <- get_sqrt(mu)
  }else{
    mu.sq <- ARGS$mu.sq
  }


  if(is.null(ARGS$tau)){
    tau <- 1/2
  }else{
    tau <- ARGS$tau
  }


  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }

  if(isTRUE(ARGS$verbose)){
    print("Calculating weights W:")
  }

  tmp <- calculate_direction_affine(mu = mu, P = W, X = X, t = u, max.iter = max.iter, tol = tol,
                                    mu.sq = mu.sq, tau = tau, mc.cores = mc.cores, check = check, ...)

  # Scale to norm 1
  tmp <- tmp/affine_norm(P = mu, W = tmp, P.sq = mu.sq)

  return(
    tmp
  )

}

#' Calculate the weights C, affine variant
#'
#' A wrapper for easy calculation of C.
#'
#' @param mu Matrix: mean of Y
#' @param C Matrix: Current estimate of C
#' @param Y List: List of data matrices Y. If Y is a data frame, please use var.name.Y
#' @param t Numeric: Numeric vector of current t scores
#' @param max.iter Numeric: Maximum number of iterations
#' @param tol Numeric: Tolerance to work at.
#' @param ... Optional arguments
#'    - mc.cores: Numeric, number of cores to run on
#'    - check: Logical, Should we check our inputs
#'    - verbose: Logical, Do we want to print progress output?
#'    - var.name.Y: Character, name of column of data if Y is data frame.
#'    - mu.sq: Matrix, square root of mu to speed computation
#'    - tau: Numeric, initial step size.
#'
#' @return The new Weights C
#' @export calculate_c_affine
#'
#' @examples
#' TBA
calculate_c_affine <- function(mu, C, Y, t, max.iter = 50, tol = 1e-6, scale_c = TRUE, ...){

  ARGS <- list(...)

  # Make sure X and Y are lists.
  if(is.null(ARGS$var.name.Y)){
    var.name.Y <- NULL
  }else{
    var.name.Y <- ARGS$var.name.Y
  }

  Y <- check_right_input(Y, var.name = var.name.Y)$x


  if(is.null(ARGS$mu.sq)){
    mu.sq <- get_sqrt(mu)
  }else{
    mu.sq <- ARGS$mu.sq
  }


  if(is.null(ARGS$tau)){
    tau <- 1/2
  }else{
    tau <- ARGS$tau
  }


  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }

  if(isTRUE(ARGS$verbose)){
    print("Calculating weights C:")
  }

  tmp <- calculate_direction_affine(mu = mu, P = C, X = Y, t = t, max.iter = max.iter, tol = tol,
                                    mu.sq = mu.sq, tau = tau, mc.cores = mc.cores, check = check, ...)

  if(scale_c){
    tmp <- tmp/affine_norm(P = mu, W = tmp, P.sq = mu.sq)
  }


  return(
    tmp
  )

}

#' Calculate the loadings P, affine variant
#'
#' A wrapper for easy calculation of P.
#'
#' @param mu Matrix: mean of X
#' @param P Matrix: Current estimate of P
#' @param X List: List of data matrices X. If X is a data frame, please use var.name.X
#' @param t Numeric: Numeric vector of current t scores
#' @param max.iter Numeric: Maximum number of iterations
#' @param tol Numeric: Tolerance to work at.
#' @param ... Optional arguments
#'    - mc.cores: Numeric, number of cores to run on
#'    - check: Logical, Should we check our inputs
#'    - verbose: Logical, Do we want to print progress output?
#'    - var.name.X: Character, name of column of data if X is data frame.
#'    - mu.sq: Matrix, square root of mu to speed computation
#'    - tau: Numeric, initial step size.
#'
#' @return The new Loadings P.
#' @export calculate_p_affine
#'
#' @examples
#' TBA
calculate_p_affine <- function(mu, P, X, t, max.iter = 50, tol = 1e-6, j = 0, ...){

  ARGS <- list(...)

  s <- round(seq(1, length(X), length.out = min(length(X), max.iter)))

  # Make sure X and Y are lists.
  if(is.null(ARGS$var.name.X)){
    var.name.X <- NULL
  }else{
    var.name.X <- ARGS$var.name.X
  }

  X <- check_right_input(X, var.name = var.name.X)$x


  if(is.null(ARGS$mu.sq)){
    mu.sq <- get_sqrt(mu)
  }else{
    mu.sq <- ARGS$mu.sq
  }


  if(is.null(ARGS$tau)){
    tau <- 1/2
  }else{
    tau <- ARGS$tau
  }


  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }

  if(is.null(ARGS$P0)){
    P0 <- P
  }else{
    P0 <- ARGS$P0
  }

  if(isTRUE(ARGS$verbose)){
    print("Calculating Loadings P:")
  }


  res <- tryCatch(calculate_direction_affine(mu = mu, P = P, X = X, t = t, max.iter = max.iter, tol = tol,
                                             mu.sq = mu.sq, tau = tau, mc.cores = mc.cores, check = check,
                                             convergence_results = TRUE, ...),
                  error = function(e){ "stop" },
                  warning = function(e){ "stop" })

  if(identical(res, "stop")){
    j <- j + 1
    if(j > length(X)){
      print("P did not converge.")
      return(P0)
    }

    P <- X[[s[j]]]
    return(
      calculate_p_affine(mu = mu, P = P, X = X, t = t, max.iter = max.iter, tol = tol, j = j, P0 = P0, ...)
    )
  }

  if((res$e > tol & res$tau < 1e-12)){
    j <- j + 1
    if(j > length(X)){
      print("P did not converge.")
      return(res$ANS)
    }

    P <- X[[s[j]]]
    return(
      calculate_p_affine(mu = mu, P = P, X = X, t = t, max.iter = max.iter, tol = tol, j = j, P0 = P0, ...)
    )
  }
  # Something for when J hits the end

  return(
    # calculate_direction_affine(mu = mu, P = P, X = X, t = t, max.iter = max.iter, tol = tol,
    #                            mu.sq = mu.sq, tau = tau, mc.cores = mc.cores, check = check, ...)
    return(res$ANS)
  )

}
#' Calculate the loadings Q, affine variant
#'
#' A wrapper for easy calculation of Q.
#'
#' @param mu Matrix: mean of Y
#' @param Q Matrix: Current estimate of Q
#' @param Y List: List of data matrices Y. If Y is a data frame, please use var.name.Y
#' @param u Numeric: Numeric vector of current u scores
#' @param max.iter Numeric: Maximum number of iterations
#' @param tol Numeric: Tolerance to work at.
#' @param ... Optional arguments
#'    - mc.cores: Numeric, number of cores to run on
#'    - check: Logical, Should we check our inputs
#'    - verbose: Logical, Do we want to print progress output?
#'    - var.name.Y: Character, name of column of data if Y is data frame.
#'    - mu.sq: Matrix, square root of mu to speed computation
#'    - tau: Numeric, initial step size.
#'
#' @return The new Loadings Q.
#' @export calculate_q_affine
#'
#' @examples
#' TBA
calculate_q_affine <- function(mu, Q, Y, u, max.iter = 50, tol = 1e-6, j = 0, ...){

  ARGS <- list(...)

  s <- round(seq(1, length(Y), length.out = min(length(Y), max.iter)))

  # Make sure X and Y are lists.
  if(is.null(ARGS$var.name.Y)){
    var.name.Y <- NULL
  }else{
    var.name.Y <- ARGS$var.name.Y
  }

  Y <- check_right_input(Y, var.name = var.name.Y)$x


  if(is.null(ARGS$mu.sq)){
    mu.sq <- get_sqrt(mu)
  }else{
    mu.sq <- ARGS$mu.sq
  }


  if(is.null(ARGS$tau)){
    tau <- 1/2
  }else{
    tau <- ARGS$tau
  }


  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }

  if(isTRUE(ARGS$verbose)){
    print("Calculating Loadings Q:")
  }

  if(is.null(ARGS$Q0)){
    Q0 <- Q
  }else{
    Q0 <- ARGS$Q0
  }

  res <- tryCatch(calculate_direction_affine(mu = mu, P = Q, X = Y, t = u, max.iter = max.iter, tol = tol,
                                             mu.sq = mu.sq, tau = tau, mc.cores = mc.cores, check = check,
                                             convergence_results = TRUE, ...),
                  error = function(e){ "stop" },
                  warning = function(e){ "stop" })

  if(identical(res, "stop")){
    j <- j + 1
    if(j > length(Y)){
      print("Q did not converge.")
      return(Q0)
    }

    Q <- Y[[s[j]]]
    return(
      calculate_q_affine(mu = mu, Q = Q, Y = Y, u = u, max.iter = max.iter, tol = tol, j = j, Q0 = Q0, ...)
    )
  }

  if((res$e > tol & res$tau < 1e-12)){
    j <- j + 1
    if(j > length(Y)){
      print("Y did not converge.")
      return(res$ANS)
    }

    Q <- Y[[s[j]]]
    return(
      calculate_q_affine(mu = mu, Q = Q, Y = Y, u = u, max.iter = max.iter, tol = tol, j = j, Q0 = Q0, ...)
    )
  }
  # Something for when J hits the end

  return(
    # calculate_direction_affine(mu = mu, P = P, X = X, t = t, max.iter = max.iter, tol = tol,
    #                            mu.sq = mu.sq, tau = tau, mc.cores = mc.cores, check = check, ...)
    return(res$ANS)
  )

}

#' Calculate the scores and weights for geodesic PLS, affine variant
#'
#' The title says it all.
#'
#' @param X List: Data predictor matrices X.  If dataframe, please use optional input var.name.X
#' @param Y List: Data predictor matrices Y.  If dataframe, please use optional input var.name.Y
#' @param muX Matrix: Mean of the X.
#' @param muY Matrix: Mean of Y.
#' @param max.iter Numeric: Total number of iterations to do.
#' @param tol Numeric: Tolerance to work at.
#' @param ... Optional arguments
#'     - var.name.X, var.name.Y: Character, the column containing the data matrices
#'     - muX.sq, muY.sq: Matrix, the square root of mu to speed computation
#'     - tau: Numeric, Initial step size
#'     - mc.cores: Numeric, if running parallel, how many cores.
#'     - check: Logical, should we check our parameters
#'     - t, u, W, C: Initial estimates of scores and weights.
#'     - verbose: Logical, should run time information be printed.
#'     - save: Logical, to be implemented.
#'
#' @return A list with W = weights of X, C = Weights of Y, t = scores of X, u = scores of Y.
#' @export scores_and_weights_affine
#'
#' @examples
#' TBA
scores_and_weights_affine <- function(X, Y, muX, muY, max.iter = 50, tol = 1e-6, tau = 1/2, ...){

  ARGS <- list(...)

  ## Add saving options

  if(is.null(ARGS$muX.sq)){
    muX.sq <- get_sqrt(muX)
  }else{
    muX.sq <- ARGS$muX.sq
  }
  if(is.null(ARGS$muY.sq)){
    muY.sq <- get_sqrt(muY)
  }else{
    muY.sq <- ARGS$muY.sq
  }


  if(is.null(ARGS$tau)){
    tau <- 1/2 # Experimentation has found this is an appropriate step size.
  }else{
    tau <- ARGS$tau
  }


  if(is.null(ARGS$mc.cores)){
    mc.cores <- 1
  }else{
    mc.cores <- ARGS$mc.cores
  }

  if(!isTRUE(ARGS$check)){
    check <- FALSE
  }else{
    check <- TRUE
  }


  # Intitialise inputs
  if(is.null(ARGS$u)){
    u_old <- initiate_u(Y, muY)
  }else{
    u_old <- ARGS$u
  }
  if(is.null(ARGS$t)){
    t <- initiate_u(X, muX)
  }else{
    t <- ARGS$t
  }
  if(is.null(ARGS$W)){
    W <- muX
  }else{
    W <- ARGS$W
  }
  if(is.null(ARGS$C)){
    C <- muY
  }else{
    C <- ARGS$C
  }



  e <- 1
  k <- 0

  if(isTRUE(ARGS$verbose)){
    print(str_c(
      "Begins calculating Scores/Weights:"
    ))
  }

  while(k < max.iter){

    k <- k + 1

    W <- calculate_w_affine(mu = muX, W = W, X = X, u = u_old, max.iter = max.iter, tol = tol,
                            tau = tau, mu.sq = muX.sq, mc.cores = mc.cores, check = check, ...)

    t <- calculate_t_affine(mu = muX, W = W, X = X, t = t, tau = tau, max.iter = max.iter, tol = tol,
                            mc.cores = mc.cores, check = check, ...)

    C <- calculate_c_affine(mu = muY, C = C, Y = Y, t = t, max.iter = max.iter, tol = tol,
                            tau = tau, mu.sq = muY.sq, mc.cores = mc.cores, check = check, ...)

    u_new <- calculate_u_affine(mu = muY, C = C, Y = Y, u = u_old, tau = tau, max.iter = max.iter, tol = tol,
                                mc.cores = mc.cores, check = check, ...)

    e <- sqrt(sum((u_new - u_old)^2))

    # Add in verbose
    if(isTRUE(ARGS$verbose)){
      print(str_c(
        "Scores/Weights convergence information:"
      ))
      print(str_c(
        "Iteration: ", k,
        "; Current change in u: ", round(e, 8),
        "; Current step size: ", round(tau, 8)
      ))
    }

    u_old <- u_new

    if(isTRUE(e < tol)){
      break
    }

    # Decrease step size to help convergence.
    tau <- tau/2
  }

  if(isTRUE(ARGS$verbose)){
    print(str_c(
      "End calculating Scores/Weights:"
    ))
    print(str_c(
      "Total iterations: ", k,
      "; Final change in u: ", round(e, 8),
      "; Final step size: ", round(tau, 8)
    ))
  }

  return(
    list(
      W = W,
      C = C,
      t = t,
      u = u_new
    )
  )
}


#' Parallel transport in the affine-invariant geometry
#'
#' Calculate the affine invariant parallel transport of W from T_PM to T_Q M
#'
#' @param P Matrix: base point matrix
#' @param Q Matrix: Endpoint matrix
#' @param W Matrix: Vector in T_P M
#' @param method Character: Either symmetric or non-symmetric.  Changes the matrix used to calculated the parallel transport.
#'
#' @return A vector in T_Q M
#' @export calculate_parallel_transport
#'
#' @examples
#' TBA
calculate_parallel_transport <- function(P, Q, W, method = c("symmetric", "non-symmetric")){
  method <- str_to_lower(method[1])
  if(method == "symmetric"){
    translate_W <- translate_to_identity(P = P, W = W)
    translate_Q <- translate_to_identity(P = P, W = Q)

    translate_W_by_Q <- translate_to_identity(P = translate_Q,
                                              W = translate_W,
                                              reverse = TRUE)

    final_translation <- translate_to_identity(P = P,
                                               W = translate_W_by_Q,
                                               reverse = TRUE)
  }else if(method == "non-symmetric"){
    P.inv.Q <-  solve(P) %*% Q
    P.inv.Q.sq <- get_sqrt(P.inv.Q, method = "expm")
    final_translation <- translate_to_identity(P = P.inv.Q,
                                               W = W,
                                               reverse = TRUE,
                                               P.sq = P.inv.Q.sq)

  }else{
    stop("Please enter a method as either symmetric or non-symmetric.")
  }


  return(final_translation)
}




#' Deflate X, affine-variant
#'
#' Deflate the X-data in the affine-PLS model.
#'
#' @param X List: List of PSD matrices, the X data.
#' @param P Matrix: The loadings from the PLS model.
#' @param t Numeric: Vector of estimated scores from the PLS model.
#' @param muX Matrix: Frechet mean of the X
#' @param par Logical: Use parallel transport or more direct group action?
#' @param ... Optional arguments
#'     - mu.sq: Square root of the matrix muX to speed up computation.
#'
#' @return A list the same size as X containing deflated data
#' @export deflate_X_affine
#'
#' @examples
#' TBA
deflate_X_affine <- function(X, P, t, muX, par = TRUE, ...){

  ARGS <- list(...)

  if(is.null(ARGS$mu.sq)){
    mu.sq <- get_sqrt(muX)
  }else{
    mu.sq <- ARGS$mu.sq
  }

  # Calculate the geodesic approximation
  E <- map(t,
           function(t){
             affine_exp(muX, t * P, P.sq = mu.sq)
           }
  )

  # Map to tangent space at geodesic approximation
  L <- map2(E, X, affine_log)

  if(par){
    eta <- map2(E, L, calculate_parallel_transport, Q = muX)
  }else{
    # Translate to the identity
    eta <- map2(E, L, translate_to_identity)

    # Translate back out to mean
    eta <- map(eta, translate_to_identity, P = muX, reverse = TRUE, P.sq = mu.sq)
  }


  # # Map L to the identity
  # eta <- map2(E, L, translate_to_identity)
  #
  # # Map back out to the mean muX
  # eta <- map(eta, translate_to_identity, P = muX, reverse = TRUE, P.sq = mu.sq)

  # Push this back down to the manifold with exponential
  deflated <- map(eta, affine_exp, P = muX, P.sq = mu.sq)

  # Return deflated
  return(deflated)
}

#' Deflate Y, affine variant
#'
#' Deflate the response data Y in the affine-PLS model
#'
#' @param Y
#' @param Q
#' @param B
#' @param muY
#' @param par
#' @param ...
#'
#' @return
#' @export deflate_Y_affine
#'
#' @examples
deflate_Y_affine <- function(Y, Q, B, muY, par = TRUE, ...){
  # add mc.cores? Depends how slow it is.

  ARGS <- list(...)

  if(is.null(ARGS$mu.sq)){
    mu.sq <- get_sqrt(muY)
  }else{
    mu.sq <- ARGS$mu.sq
  }

  # Calculate fitted values from LM
  fits <- fitted(B)

  # Get geodesic approximation
  E <- map(fits, function(t){
    affine_exp(muY, t * Q, P.sq = mu.sq)
  })

  # Find "difference"
  L <- map2(E, Y, affine_log)

  if(par){
    eta <- map2(E, L, calculate_parallel_transport, Q = muY)
  }else{
    # Translate to the identity
    eta <- map2(E, L, translate_to_identity)

    # Translate back out to mean
    eta <- map(eta, translate_to_identity, P = muY, reverse = TRUE, P.sq = mu.sq)
  }

  # Push back down to manifold
  deflated <- map(eta, affine_exp, P = muY, P.sq = mu.sq)

  # Return deflated values.
  return(deflated)
}

