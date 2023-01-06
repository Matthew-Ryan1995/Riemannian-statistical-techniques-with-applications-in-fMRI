#' Initiate u
#'
#' Initiate the score vector u in the NIPALS algorithm
#'
#' @param Y List: List of correlation matrices.  The data.
#' @param muY Matrix: Frechet mean of Y.
#' @param ... Optional arguments:
#'     - check: logical, should inputs be checked.
#'     - muY.sq: Matrix, square root of muY to speed up computation
#'
#' @return A vector u of initial estimates of the scores
#'
#' @details First, each Y_i is mapped to the tangent space of muY via Log_{muY}.  Then these are vectorised
#' using the vec operation.  Finally, the first column of the resulting matrix is returned.
#'
#' @export initiate_u
#'
#' @examples
#'
#' TBA
initiate_u <- function(Y, muY, ...){

  ARGS <- list(...)

  if(isTRUE(ARGS$check)){
    if(any(map_lgl(Y, ~!is.matrix(.x)), map_lgl(Y, ~!is.numeric(.x)))){
      stop("Y must be a list of numeric matrices")
    }
    if(any(map_lgl(Y, ~!isSymmetric(.x, tol = sqrt(.Machine$double.eps))))){
      stop("The Y matrices must be symmetric")
    }

    if(any(map_lgl(Y, ~{any(dim(muY) != dim(.x))}))){
      stop("muY and Y should be the same dimension.")
    }

    if(!isSymmetric(muY, tol = sqrt(.Machine$double.eps))){
      stop("muY must be a symmetric matrix")
    }

    if(!is_positive_definite(mu)){
      stop("mu must be positive definite, at least one eigenvalue is negative.")
    }
  }

  if(!is.null(ARGS$type)){
    ARGS$type <- str_to_lower(ARGS$type)
    if(ARGS$type == "euclidean"){
      return(Y[, 1])
    }
  }

  if(is.null(ARGS$muY.sq)){ # Get square root to speed up computation
    muY.sq <- get_sqrt(muY)
  }else{
    muY.sq <- ARGS$muY.sq
  }

  n <- nrow(Y[[1]])


  # Take the data Y to the tangent space at muY
  logged_Y <- map(Y, affine_log, P = muY, P.sq = muY.sq)

  # Vectorise
  vector_Y <- map(logged_Y, vec, P = muY, P.sq = muY.sq)

  # Return first column not related to the diagonal.
  u <- map_dbl(vector_Y, ~.x[n + 1])
  return(u)
}

#' Calculate the scores and weights for Riemannian PLS
#'
#' @param X List or matrix: Predictor data.  List for affine, matrix for euclidean
#' @param Y List or matrix: Response data.  List for affine, matrix for euclidean
#' @param muX Matrix:  If X is affine, the frechet mean of X
#' @param muY Matrix:  If Y is affine, the frechet mean of Y
#' @param max.iter Numeric: The maximum number of iterations in the while loops
#' @param tol Numeric: The tolerance/precision level to work at
#' @param tau Numeric: Step size in the gradient descent.
#' @param ... Optional list:
#'     - type_x, type_y: String. The data type.  Can be "euclidean" or "affine".  Defualts to affine.
#'     - muX.sq, muY.sq: Matrix.  The square roots of muX and muY.
#'     - mc.cores: Numeric.  The number of cores to run on.
#'     - check: Logical.  Should we check the input data types.
#'     - u, t, W, C: Initial estimates of u, t, W, and C.
#'     - verbose: Logical.  Should we print convergence output.
#'
#' @return  A list of values containing the estimates of u, t, W and C.
#' @export scores_and_weights
#'
#' @examples
#' TBA
scores_and_weights <- function(X, Y, muX = NULL, muY = NULL, max.iter = 50, tol = 1e-6, tau = 1/2, ...){

  ARGS <- list(...)

  ## Add saving options
  # Find X type
  if(!is.null(ARGS$type_x)){
    ARGS$type_x <- str_to_lower(ARGS$type_x)
    if(ARGS$type_x == "affine"){
      if(any(map_lgl(X, ~!is.matrix(.x)), map_lgl(X, ~!is.numeric(.x)))){
        stop("Since X is affine, it must be a list of numeric matrices")
      }
      if(is.null(muX)){
        stop("Enter the frechet mean of X")
      }
      if(is.null(ARGS$muX.sq)){
        muX.sq <- get_sqrt(muX)
      }else{
        muX.sq <- ARGS$muX.sq
      }
      calculate_w <- calculate_w_affine
      calculate_t <- calculate_t_affine
    }else if(ARGS$type_x == "euclidean"){
      if(!is.matrix(X)){
        stop("Since X is Euclidean, X should be a matrix.")
      }
      calculate_w <- calculate_w_euclidean
      calculate_t <- calculate_t_euclidean
    }else{
      stop("We can only currently consider predictors of either affine type or euclidean type.")
    }
  }else{
    if(any(map_lgl(X, ~!is.matrix(.x)), map_lgl(X, ~!is.numeric(.x)))){
      stop("Since X is affine, it must be a list of numeric matrices")
    }
    if(is.null(muX)){
      stop("Enter the frechet mean of X")
    }
    if(is.null(ARGS$muX.sq)){
      muX.sq <- get_sqrt(muX)
    }else{
      muX.sq <- ARGS$muX.sq
    }
    calculate_w <- calculate_w_affine
    calculate_t <- calculate_t_affine
  }

  # Find Y type
  if(!is.null(ARGS$type_y)){
    ARGS$type_y <- str_to_lower(ARGS$type_y)
    if(ARGS$type_y == "affine"){
      if(any(map_lgl(Y, ~!is.matrix(.x)), map_lgl(Y, ~!is.numeric(.x)))){
        stop("Since Y is affine, it must be a list of numeric matrices")
      }
      if(is.null(muY)){
        stop("Enter the frechet mean of Y")
      }

      if(is.null(ARGS$muY.sq)){
        muY.sq <- get_sqrt(muY)
      }else{
        muY.sq <- ARGS$muY.sq
      }
      calculate_c <- calculate_c_affine
      calculate_u <- calculate_u_affine
    }else if(ARGS$type_y == "euclidean"){
      if(!is.matrix(Y)){
        stop("Since Y is Euclidean, Y should be a matrix.")
      }
      calculate_c <- calculate_c_euclidean
      calculate_u <- calculate_u_euclidean
    }else{
      stop("We can only currently consider response variables of either affine type or euclidean type.")
    }
  }else{
    if(any(map_lgl(Y, ~!is.matrix(.x)), map_lgl(Y, ~!is.numeric(.x)))){
      stop("Since Y is affine, it must be a list of numeric matrices")
    }
    if(is.null(muY)){
      stop("Enter the frechet mean of Y")
    }

    if(is.null(ARGS$muY.sq)){
      muY.sq <- get_sqrt(muY)
    }else{
      muY.sq <- ARGS$muY.sq
    }
    calculate_c <- calculate_c_affine
    calculate_u <- calculate_u_affine
  }




  # if(is.null(ARGS$tau)){
  #   tau <- 1/2 # Experimentation has found this is an appropriate step size.
  # }else{
  #   tau <- ARGS$tau
  # }


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
  if(is.null(ARGS$u0)){
    u_old <- initiate_u(Y, muY, type = ARGS$type_y, ...)
  }else{
    u_old <- ARGS$u0
  }
  if(is.null(ARGS$t0)){
    t <- initiate_u(X, muX, type = ARGS$type_x, ...)
  }else{
    t <- ARGS$t0
  }
  if(is.null(ARGS$W0)){
    W <- muX
  }else{
    W <- ARGS$W0
  }
  if(is.null(ARGS$C0)){
    C <- muY
  }else{
    C <- ARGS$C0
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

    W <- calculate_w(mu = muX, W = W, X = X, u = u_old, max.iter = max.iter, tol = tol,
                     tau = tau, mu.sq = muX.sq, mc.cores = mc.cores, check = check, ...)
    t <- calculate_t(mu = muX, W = W, X = X, t = t, tau = tau, max.iter = max.iter, tol = tol,
                     mc.cores = mc.cores, check = check, ...)

    C <- calculate_c(mu = muY, C = C, Y = Y, t = t, max.iter = max.iter, tol = tol,
                     tau = tau, mu.sq = muY.sq, mc.cores = mc.cores, check = check, ...)

    u_new <- calculate_u(mu = muY, C = C, Y = Y, u = u_old, tau = tau, max.iter = max.iter, tol = tol,
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


#' Calculate the X-loadings in Riemannian PLS
#'
#' @param X List or matrix: Predictor data.  List for affine, matrix for euclidean
#' @param t Numeric: Vector of scores for the X-data.
#' @param muX Matrix:  If X is affine, the frechet mean of X
#' @param max.iter Numeric: The maximum number of iterations in the while loops
#' @param tol Numeric: The tolerance/precision level to work at
#' @param tau Numeric: Step size in the gradient descent.
#' @param P Matrix: Initial estimate of the loadings P.
#' @param ... Optional arguments.
#'     - type_x: String. The data type.  Can be "euclidean" or "affine".  Defualts to affine.
#'     - muX.sq: Matrix.  The square roots of muX.
#'     - other arguments passed to calculate_p.
#'
#' @return The loading vector P
#' @export calculate_loading_X
#'
#' @examples
#' TBA
calculate_loading_X <- function(X, t, muX = NULL, max.iter = 50, tol = 1e-6, tau = 1/2, P = NULL, ...){
  # type_x - important

  ARGS <- list(...)

  ## Add saving options
  # Find X type
  if(!is.null(ARGS$type_x)){
    ARGS$type_x <- str_to_lower(ARGS$type_x)
    if(ARGS$type_x == "affine"){
      if(any(map_lgl(X, ~!is.matrix(.x)), map_lgl(X, ~!is.numeric(.x)))){
        stop("Since X is affine, it must be a list of numeric matrices")
      }
      if(is.null(muX)){
        stop("Enter the frechet mean of X")
      }

      if(is.null(ARGS$muX.sq)){
        mu.sq <- get_sqrt(muX)
      }else{
        mu.sq <- ARGS$muX.sq
      }
      calculate_p <- calculate_p_affine
    }else if(ARGS$type_x == "euclidean"){
      if(!is.matrix(X)){
        stop("Since X is Euclidean, X should be a matrix.")
      }
      calculate_p <- calculate_p_euclidean
    }else{
      stop("We can only currently consider predictors of either affine type or euclidean type.")
    }
  }else{
    if(any(map_lgl(X, ~!is.matrix(.x)), map_lgl(X, ~!is.numeric(.x)))){
      stop("Since X is affine, it must be a list of numeric matrices")
    }
    if(is.null(muX)){
      stop("Enter the frechet mean of X")
    }
    if(is.null(ARGS$muX.sq)){
      mu.sq <- get_sqrt(muX)
    }else{
      mu.sq <- ARGS$muX.sq
    }
    calculate_p <- calculate_p_affine
  }

  return(
    calculate_p(X = X, mu = muX, P = P, t = t, max.iter = max.iter, tol = tol, tau = tau,
                mu.sq = mu.sq, ...)
  )
}

#' Calculate the X-loadings in Riemannian PLS
#'
#' @param X List or matrix: Predictor data.  List for affine, matrix for euclidean
#' @param t Numeric: Vector of scores for the X-data.
#' @param muX Matrix:  If X is affine, the frechet mean of X
#' @param max.iter Numeric: The maximum number of iterations in the while loops
#' @param tol Numeric: The tolerance/precision level to work at
#' @param tau Numeric: Step size in the gradient descent.
#' @param P Matrix: Initial estimate of the loadings P.
#' @param ... Optional arguments.
#'     - type_x: String. The data type.  Can be "euclidean" or "affine".  Defualts to affine.
#'     - muX.sq: Matrix.  The square roots of muX.
#'     - other arguments passed to calculate_p.
#'
#' @return The loading vector P
#' @export calculate_loading_Y
#'
#' @examples
#' TBA
calculate_loading_Y <- function(Y, u, muY = NULL, max.iter = 50, tol = 1e-6, tau = 1/2, Q = NULL, ...){
  # type_x - important

  ARGS <- list(...)

  # if(any(map_lgl(Y, ~!is.matrix(.x)), map_lgl(Y, ~!is.numeric(.x)))){
  #   stop("Since Y is affine, it must be a list of numeric matrices")
  # }
  # if(is.null(muY)){
  #   stop("Enter the frechet mean of Y")
  # }
  #
  # if(is.null(ARGS$muY.sq)){
  #   mu.sq <- get_sqrt(muY)
  # }else{
  #   mu.sq <- ARGS$muY.sq
  # }
  # calculate_q <- calculate_q_affine
  #
  # return(
  #   calculate_q(Y = Y, mu = muY, Q = Q, u = u, max.iter = max.iter, tol = tol, tau = tau,
  #               mu.sq = mu.sq, ...)
  # )

  # Add saving options
  # Find X type
  if(!is.null(ARGS$type_y)){
    ARGS$type_y <- str_to_lower(ARGS$type_y)
    if(ARGS$type_y == "affine"){
      if(any(map_lgl(Y, ~!is.matrix(.x)), map_lgl(Y, ~!is.numeric(.x)))){
        stop("Since Y is affine, it must be a list of numeric matrices")
      }
      if(is.null(muY)){
        stop("Enter the frechet mean of Y")
      }

      if(is.null(ARGS$muY.sq)){
        mu.sq <- get_sqrt(muY)
      }else{
        mu.sq <- ARGS$muY.sq
      }
      calculate_q <- calculate_q_affine
    }else if(ARGS$type_y == "euclidean"){
      if(!is.matrix(Y)){
        stop("Since Y is Euclidean, Y should be a matrix.")
      }
      calculate_q <- calculate_q_euclidean
    }else{
      stop("We can only currently consider predictors of either affine type or euclidean type.")
    }
  }else{
    if(any(map_lgl(Y, ~!is.matrix(.x)), map_lgl(Y, ~!is.numeric(.x)))){
      stop("Since Y is affine, it must be a list of numeric matrices")
    }
    if(is.null(muY)){
      stop("Enter the frechet mean of Y")
    }
    if(is.null(ARGS$muY.sq)){
      mu.sq <- get_sqrt(muY)
    }else{
      mu.sq <- ARGS$muY.sq
    }
    calculate_q <- calculate_q_affine
  }

  # return(
  #   calculate_p(X = X, mu = muX, P = P, t = t, max.iter = max.iter, tol = tol, tau = tau,
  #               mu.sq = mu.sq, ...)
  # )

  return(
    calculate_q(Y = Y, mu = muY, Q = Q, u = u, max.iter = max.iter, tol = tol, tau = tau,
                mu.sq = mu.sq, ...)
  )
}

#' Calculate the inner relationship for Riemannian PLS
#'
#' @param t Numeric: Vector of scores for X
#' @param u Numeric: Vector of scores for Y
#'
#' @return A list containing the linear model u ~ t, as well as the coefficients.
#' @export calculate_regression
#'
#' @examples
#' TBA
calculate_regression <- function(t, u){

  df <- tibble(t = as.numeric(t),
               u = as.numeric(u))

  M <- lm(u ~ t, data = df)
  return(
    list(
      M = M,
      b0 = coefficients(M)[1],
      b1 = coefficients(M)[2]
    )
  )
}

#' Deflate the X-data in Riemannian PLS
#'
#' @param X List or matrix: Predictor data.  List for affine, matrix for euclidean
#' @param t Numeric: Vector of scores for the X-data.
#' @param P Matrix: The loading matrix for the X-data and scores t.
#' @param muX Matrix:  If X is affine, the frechet mean of X
#' @param ... Optional list
#'     - type_x: String. The data type.  Can be "euclidean" or "affine".  Defualts to affine.
#'     - muX.sq: Matrix.  The square roots of muX.
#'     - other arguments passed to deflate
#'
#' @return A deflated version of the X-data, removing the information from the previous PLS step.
#' @export deflate_X
#'
#' @examples
#' TBA
deflate_X <- function(X, t, P, muX = NULL, ...){
  ARGS <- list(...)

  if(!is.null(ARGS$type_x)){
    ARGS$type_x <- str_to_lower(ARGS$type_x)
    if(ARGS$type_x == "affine"){
      if(any(map_lgl(X, ~!is.matrix(.x)), map_lgl(X, ~!is.numeric(.x)))){
        stop("Since X is affine, it must be a list of numeric matrices")
      }
      if(is.null(muX)){
        stop("Enter the frechet mean of X")
      }

      if(is.null(ARGS$muX.sq)){
        mu.sq <- get_sqrt(muX)
      }else{
        mu.sq <- ARGS$muX.sq
      }
      deflate <- deflate_X_affine
    }else if(ARGS$type_x == "euclidean"){
      if(!is.matrix(X)){
        stop("Since X is Euclidean, X should be a matrix.")
      }
      deflate <- deflate_X_euclidean
    }else{
      stop("We can only currently consider predictors of either affine type or euclidean type.")
    }
  }else{
    if(any(map_lgl(X, ~!is.matrix(.x)), map_lgl(X, ~!is.numeric(.x)))){
      stop("Since X is affine, it must be a list of numeric matrices")
    }
    if(is.null(muX)){
      stop("Enter the frechet mean of X")
    }
    if(is.null(ARGS$muX.sq)){
      mu.sq <- get_sqrt(muX)
    }else{
      mu.sq <- ARGS$muX.sq
    }

    deflate <- deflate_X_affine
  }

  return(deflate(X = X, t = t, P = P, muX = muX, mu.sq = mu.sq, ...))
}

#' Deflate the Y-data in Riemannian PLS
#'
#' @param Y List or matrix: Response data.  List for affine, matrix for euclidean
#' @param C Matrix: The weights matrix for the Y-data (may need to change to Y-loadings)
#' @param Q Matrix: The loadings matrix for the Y-data (may need to change to Y-loadings)
#' @param B Model: The inner relation model of u ~ t
#' @param muY Matrix:  If Y is affine, the frechet mean of Y
#' @param ... Optional list
#'     - type_y: String. The data type.  Can be "euclidean" or "affine".  Defualts to affine.
#'     - muY.sq: Matrix.  The square roots of muY.
#'     - other arguments passed to deflate
#'
#' @return A deflated version of the Y-data, removing the information from the previous PLS step.
#' @export deflate_Y
#'
#' @examples
#' TBA
deflate_Y <- function(Y, C, Q, B, muY = NULL, ...){
  # B is the linear model of u ~ t
  ARGS <- list(...)

  if(!is.null(ARGS$type_y)){
    ARGS$type_y <- str_to_lower(ARGS$type_y)
    if(ARGS$type_y == "affine"){
      if(any(map_lgl(Y, ~!is.matrix(.x)), map_lgl(Y, ~!is.numeric(.x)))){
        stop("Since Y is affine, it must be a list of numeric matrices")
      }
      if(is.null(muY)){
        stop("Enter the frechet mean of Y")
      }

      if(is.null(ARGS$muY.sq)){
        mu.sq <- get_sqrt(muY)
      }else{
        mu.sq <- ARGS$muY.sq
      }
      deflate <- deflate_Y_affine
    }else if(ARGS$type_y == "euclidean"){
      if(!is.matrix(Y)){
        stop("Since Y is Euclidean, Y should be a matrix.")
      }
      deflate <- deflate_Y_euclidean
    }else{
      stop("We can only currently consider predictors of either affine type or euclidean type.")
    }
  }else{
    if(any(map_lgl(Y, ~!is.matrix(.x)), map_lgl(Y, ~!is.numeric(.x)))){
      stop("Since Y is affine, it must be a list of numeric matrices")
    }
    if(is.null(muY)){
      stop("Enter the frechet mean of Y")
    }

    if(is.null(ARGS$muY.sq)){
      mu.sq <- get_sqrt(muY)
    }else{
      mu.sq <- ARGS$muY.sq
    }

    deflate <- deflate_Y_affine
  }

  return(deflate(Y = Y, C = C, Q = Q, B = B, muY = muY, mu.sq = mu.sq, ...))
}

#' Fit the Riemannian PLS model
#'
#' @param X List or matrix: Predictor data.  List for affine, matrix for euclidean
#' @param Y List or matrix: Response data.  List for affine, matrix for euclidean
#' @param L Numeric: Number of PLS components to calculate
#' @param tol Numeric: The tolerance/precision level to work at
#' @param max.iter Numeric: The maximum number of iterations in the while loops
#' @param muX Matrix:  If X is affine, the frechet mean of X
#' @param muY Matrix:  If Y is affine, the frechet mean of Y
#' @param tau Numeric: step size in gradient descent
#' @param ... Optional list
#'     - type_x, type_y: String. The data type.  Can be "euclidean" or "affine".  Defualts to affine.
#'     - method: String.  Are we performing the tangent-space method.  Only applicable if data is affine.  Only available choice is method = "tangent"
#'     - centre: Logical.  If X, Y are Euclidean, do we need to centre?  Suggested no if method = "tangent"
#'     - muX.sq, muY.sq: Matrix.  The square roots of muX and muY.
#'     - mc.cores: Numeric.  The number of cores to run on.
#'     - check: Logical.  Should we check the input data types.
#'     - u, t, W, C: Initial estimates of u, t, W, and C.
#'     - verbose: Logical.  Should we print convergence output.
#'     - Other options to be passed to scores, loadings, regression, and deflate functions.
#'
#' @return A list containing:
#'    - scoresX: List of scores for X at each PLS step
#'    - scoresY: List of scores for Y at each PLS step
#'    - weightsX: List of weights for X at each PLS step
#'    - weightsY: List of weights for Y at each PLS step
#'    - loadingsX: List of loadings for X at each PLS step
#'    - reg_steps: List containing the inner relationship data at each step
#'    - L: Numeric: How many PLS directions were calculated.
#'    - X: List: contains the deflated X-data, position l + 1 contains the data used to calculate PLS direction l
#'    - Y: List: contains the deflated Y-data, position l + 1 contains the data used to calculate PLS direction l
#' @export riemannian_pls
#'
#' @examples
#' TBA
riemannian_pls <- function(X, Y, L = 1, tol = 1e-7, max.iter = 50, muX = NULL, muY = NULL, tau = 1/2, ...){

  ARGS <- list(...)

  if(is.null(mc.cores <- ARGS$mc.cores)){
    mc.cores <- 1
  }
  # if(is.null(tau <- ARGS$tau)){
  #   tau <- 1/2
  # }
  #
  #

  if(isTRUE(ARGS$type_x == "affine") | is.null(ARGS$type_x)){
    if(is.null(muX)){
      muX <- get_frechet_mean(X, cores = mc.cores)
    }
    X_sd <- NULL
  }else{
    muX <- colMeans(X)
    # if(isTRUE(ARGS$centre)){
    if(is.null(scale <- ARGS$scale)){
      scale <- FALSE
      X_sd <- NULL
    }else{
      X_sd <- apply(X, 2, sd)
    }
    X <- apply(X, 2, scale, scale = scale)
  }

  if(isTRUE(ARGS$type_y == "affine") | is.null(ARGS$type_y)){
    if(is.null(muY)){
      muY <- get_frechet_mean(Y, cores = mc.cores)
    }
    Y_sd <- NULL
  }else{
    muY <- colMeans(Y)
    if(is.null(scale <- ARGS$scale)){
      scale <- FALSE
      Y_sd <- NULL
    }else{
      Y_sd <- apply(Y, 2, sd)
    }
    # if(isTRUE(ARGS$centre)){
    Y <- apply(Y, 2, scale, scale = scale)
    # }
  }

  if(isTRUE(ARGS$method == "tangent")){
    if(isTRUE(ARGS$type_x == "affine") | is.null(ARGS$type_x)){
      if(is.null(muX)){
        muX <- get_frechet_mean(X, cores = mc.cores)
      }
      muX.sq <- get_sqrt(muX)
      X <- linearise_data(X = X, mu = muX, mc.cores = mc.cores, mu.sq = muX.sq) ## WRITE THIS FUNCTION
      ARGS$type_x <- "euclidean"
      # muX_tangent <- muX
    }
    if(isTRUE(ARGS$type_y == "affine") | is.null(ARGS$type_y)){
      if(is.null(muY)){
        muY <- get_frechet_mean(Y, cores = mc.cores)
      }
      muY.sq <- get_sqrt(muY)
      Y <- linearise_data(X = Y, mu = muY, mc.cores = mc.cores, mu.sq = muY.sq)
      ARGS$type_y <- "euclidean"
      # muY_tangent <- muY
    }
  }


  scores_t <- list()
  scores_u <- list()
  weights_W <- list()
  weights_C <- list()
  loadings_P <- list()
  loadings_Q <- list()
  regression_B <- list()
  deflated_X <- list(X)
  deflated_Y <- list(Y)

  for(l in 1:L){
    tmp_scores_weights <- scores_and_weights(X = deflated_X[[l]], Y = deflated_Y[[l]], muX = muX, muY = muY,
                                             max.iter = max.iter, tol = tol, tau = tau, mc.cores = mc.cores,
                                             type_x = ARGS$type_x, type_y = ARGS$type_y, ...)

    scores_t[[l]] <- tmp_scores_weights$t
    scores_u[[l]] <- tmp_scores_weights$u
    weights_W[[l]] <- tmp_scores_weights$W
    weights_C[[l]] <- tmp_scores_weights$C

    loadings_P[[l]] <- calculate_loading_X(X = deflated_X[[l]], t = scores_t[[l]], P = weights_W[[l]], muX = muX,
                                           max.iter = max.iter, tau = tau, tol = tol,
                                           type_x = ARGS$type_x, type_y = ARGS$type_y, ...)

    loadings_Q[[l]] <- calculate_loading_Y(Y = deflated_Y[[l]], u = scores_u[[l]], Q = weights_C[[l]], muY = muY,
                                           max.iter = max.iter, tau = tau, tol = tol,
                                           type_x = ARGS$type_x, type_y = ARGS$type_y, ...)

    regression_B[[l]] <- calculate_regression(t = scores_t[[l]], u = scores_u[[l]])

    deflated_X[[l + 1]] <- deflate_X(X = deflated_X[[l]], t = scores_t[[l]], P = loadings_P[[l]], muX = muX,
                                     type_x = ARGS$type_x, type_y = ARGS$type_y, ...)
    ## Fix this up so that it does the right thing with C and Q
    deflated_Y[[l + 1]] <- deflate_Y(Y = deflated_Y[[l]], C = weights_C[[l]], Q = loadings_Q[[l]], B = regression_B[[l]]$M,  muY = muY,
                                     type_x = ARGS$type_x, type_y = ARGS$type_y, ...)
  }


  # if(isTRUE(ARGS$method == "tangent")){
  #   if(!is.null(muX_tangent)){
  #     muX <- muX_tangent
  #   }
  #   if(!is.null(muY_tangent)){
  #     muY <- muY_tangent
  #   }
  # }
  res <-
    list(
      scoresX = scores_t,
      scoresY = scores_u,
      weightsX = weights_W,
      weightsY = weights_C,
      loadingsX = loadings_P,
      loadingsY = loadings_Q,
      reg_steps = regression_B,
      L = L,
      X = deflated_X,
      Y = deflated_Y,
      muX = muX,
      muY = muY,
      X_sd = X_sd,
      Y_sd = Y_sd,
      type_x = ARGS$type_x,
      type_y = ARGS$type_y,
      methods = ARGS$method
    )

  return(res)
}


#' Linearise the data to apply the tangent space method.
#'
#' @param X List: List of SPD matrices
#' @param mu Matrix: Mean of X
#' @param ... Optional list
#'    - mc.cores: Numeric.  Number of cores to run on.
#'    - mu.sq: Matrix.  Square root of muX.
#'
#' @return A matrix that is the tangent-space vectorisation of the X-data
#' @export linearise_data
#'
#' @examples
#' TBA
linearise_data <- function(X, mu, ...){

  ARGS <- list(...)

  if(is.null(mc.cores <- ARGS$mc.cores)){
    mc.cores <- 1
  }

  if(is.null(mu.sq <- ARGS$mu.sq)){
    mu.sq <- get_sqrt(mu)
  }

  logged_data <- parallel::mclapply(X,
                                    function(x){
                                      affine_log(P = mu, W = x, P.sq = mu.sq)
                                    },
                                    mc.cores = mc.cores)
  vec_data <- parallel::mclapply(logged_data,
                                 function(x){
                                   vec(P = mu, W = x, P.sq = mu.sq)
                                 },
                                 mc.cores = mc.cores) %>%
    simplify2array() %>%
    t()

  return(vec_data)
}
