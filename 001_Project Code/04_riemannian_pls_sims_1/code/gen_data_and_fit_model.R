#' generate_rpls_data
#' Generate data from the Rpls model
#'
#' @param n - number of subjects
#' @param x_dim - Dimension of X-space, either matrix dim or euclidean dim
#' @param y_dim  - Dimension of Y-space, either matrix dim or euclidean dim
#' @param sig - sigma noise level
#' @param seed - seed to generate data from
#' @param L - True number of latent variables
#' @param muX - mean to generate x about
#' @param muY - mean to generate y about
#' @param type_x - either affine or euclidean
#' @param type_y - either affine or euclidean
#'
#' @return
#' @export
#'
#' @examples
generate_rpls_data <- function(n, x_dim, y_dim, sig, seed = NULL, L = 5,
                               muX = diag(1, x_dim), muY = diag(1, y_dim),
                               type_x = "affine", type_y = "affine"){

  # if(type_x == "affine"){
  p <- choose(x_dim + 1, 2)
  # }else{
  #   p <- x_dim
  # }
  # if(type_y == "affine"){
  q <- choose(y_dim + 1, 2)
  # q <- y_dim
  # }else{
  #   q <- y_dim
  # }

  ## Fix a seed to generate the loading matrices
  set.seed(2022^2)
  # Fixed P and Q across simulation parameters p, q and L
  P <- matrix(rnorm(L*p), nrow = p, ncol = L)
  Q <- matrix(rnorm(L*q), nrow = q, ncol = L)
  # Make P and Q orthogonal
  if(p < L){
    P <- t(svd(P)$v)
  }else{
    P <- svd(P)$u
  }
  if(q < L){
    Q <- t(svd(Q)$v)
  }else{
    Q <- svd(Q)$u
  }

  ## Change the seed for the data generation
  if(!is.null(seed)){
    set.seed(seed)
  }

  ## Create the scores
  Tmat <- matrix(rnorm(n * L), nrow = n)
  ## Scale the scores about 0
  Tmat <- apply(Tmat, 2, scale, scale = FALSE)

  if(type_x == "affine"){
    P_mat <- map(1:L,
                 function(l){
                   p_vec <- P[, l]
                   unvec(muX, vec = p_vec)
                 })

    TP <- map(1:n,
              function(i){
                t_l <- Tmat[i, ] # Scores for subject

                p_l <- map2(P_mat, t_l, ~.x * .y) # Multiply scores by loading and add
                p_l <- Reduce("+", p_l)

                return(
                  affine_exp(muX, p_l) # Push down to manifold
                )
              })

    X <- map(TP,
             function(tp){

               # Make symmetric noise matrix
               E <- matrix(rnorm(x_dim^2, sd = sig), ncol = x_dim)
               E <- (E + t(E))/2
               return(affine_exp(tp, E))
             })
#
#     X <- map(1:n,
#              function(i){
#                t_l <- Tmat[i, ] # Scores for subject
#
#                p_l <- map2(P_mat, t_l, ~.x * .y) # Multiply scores by loading and add
#                p_l <- Reduce("+", p_l)
#
#                # Make symmetric noise matrix
#                E <- matrix(rnorm(x_dim^2, sd = sig), ncol = x_dim)
#                E <- (E + t(E))/2
#
#                return(
#                  affine_exp(affine_exp(muX, p_l), E) # Push down to manifold
#                )
#              })
  }else{
    E <-  matrix(rnorm(p * n, sd = sig), ncol = p)
    TP <- Tmat %*% t(P)
    X <-  TP + E
    # Scale data
    X <- apply(X, 2, scale, scale = FALSE)
  }

  if(type_y == "affine"){
    Q_mat <- map(1:L,
                 function(l){
                   q_vec <- Q[, l]
                   unvec(muY, vec = q_vec)
                 })

    UQ <- map(1:n,
              function(i){
                t_l <- Tmat[i, ] # Scores for subject, assuming B = I

                p_l <- map2(Q_mat, t_l, ~.x * .y) # Multiply scores by loading and add
                p_l <- Reduce("+", p_l)

                return(
                  affine_exp(muY, p_l) # Push down to manifold
                )
              })

    Y <- map(UQ,
             function(uq){
               E <- matrix(rnorm(y_dim^2, sd = sig), ncol = y_dim)
               E <- (E + t(E))/2
               return(affine_exp(uq, E))
             })

    # Y <- map(1:n,
    #          function(i){
    #            t_l <- Tmat[i, ] # Scores for subject, assuming B = I
    #
    #            p_l <- map2(Q_mat, t_l, ~.x * .y) # Multiply scores by loading and add
    #            p_l <- Reduce("+", p_l)
    #
    #            # Make symmetric noise matrix
    #            E <- matrix(rnorm(y_dim^2, sd = sig), ncol = y_dim)
    #            E <- (E + t(E))/2
    #
    #            return(
    #              affine_exp(affine_exp(muY, p_l), E) # Push down to manifold
    #            )
    #          })
  }else{
    E <-  matrix(rnorm(q * n, sd = sig), ncol = q)
    UQ <- Tmat %*% t(Q)
    Y <- UQ + E
    # Scale data
    Y <- apply(Y, 2, scale, scale = FALSE)
  }


  res <- list(
    X = X,
    Y = Y,
    TP = TP,
    UQ = UQ,
    trueP = P,
    trueQ = Q,
    type_x = type_x,
    type_y = type_y
  )

  return(res)
}


#' fit_riemannian_pls
#' fit Rpls on data returned from generate_rpls_data
#'
#' @param data - object from generate_rpls_data
#' @param K - Number of latent variables to fit
#' @param tol - Tolerance to work at
#' @param max.iter - max iter to work at.
#'
#' @return
#' @export
#'
#' @examples
fit_riemannian_pls <- function(data, K = 1, tol = 1e-4, max.iter = 10,
                               mc.cores = 1, method = NULL){
  model <- suppressWarnings(riemannian_pls(X = data$X, Y = data$Y, L = K,
                                           tol = tol, max.iter = max.iter,
                                           type_x = data$type_x, type_y = data$type_y,
                                           mc.cores = mc.cores, verbose = F, method = method))

  return(model)
}
