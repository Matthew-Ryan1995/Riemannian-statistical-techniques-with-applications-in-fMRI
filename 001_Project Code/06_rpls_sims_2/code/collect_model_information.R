#' collect_X
#' Returns the fitted X-values from the model
#'
#' @param model - Fitted riemannian pls
#'
#' @return
#' @export
#'
#' @examples
collect_X <- function(model){
  T.hat <- do.call(cbind, model$scoresX)
  if(model$type_x == "affine"){
    P.hat <- model$loadingsX
    X.hat <- map(1:nrow(T.hat),
                 function(i){
                   t_i <- T.hat[i, ]
                   tmp <- map2(t_i, P.hat, ~.x*.y)
                   return(
                     Reduce("+", tmp)
                   )
                 })

    X.hat <- map(X.hat, affine_exp, P = model$muX)
  }else{
    P.hat <- do.call(cbind, model$loadingsX)
    X.hat  <- T.hat %*% t(P.hat)
  }

  return(X.hat)
}

#' collect_Y
#' Returns the fitted Y-values from the model
#'
#' @param model - Fitted riemannian pls
#'
#' @return
#' @export
#'
#' @examples
collect_Y <- function(model){
  Y.hat <- fitted.riemannian_pls(model)
  return(Y.hat)
}

#' collect_P
#' Returns the fitted P matrix
#'
#' @param model - Fitted riemannian pls
#'
#' @return
#' @export
#'
#' @examples
collect_P <- function(model){
  P.hat <- model$loadingsX
  if(model$type_x == "affine"){
    P.hat <- map(P.hat, vec, P = model$muX)
  }
  P.hat <- do.call(cbind, P.hat)

  return(P.hat)
}

#' collect_Q
#' Returns the fitted Q matrix
#'
#' @param model - Fitted riemannian pls
#'
#' @return
#' @export
#'
#' @examples
collect_Q <- function(model){
  Q.hat <- model$loadingsY
  if(model$type_y == "affine"){
    Q.hat <- map(Q.hat, vec, P = model$muY)
  }
  Q.hat <- do.call(cbind, Q.hat)

  return(Q.hat)
}

#' collect_full_measures
#' collects all relevant data from RPLS model
#'
#' @param model - Fitted Riemannian PLS
#'
#' @return
#' @export
#'
#' @examples
collect_full_measures <- function(model){
  X.hat <- collect_X(model)
  Y.hat <- collect_Y(model)
  P.hat <- collect_P(model)
  Q.hat <- collect_Q(model)

  # X.hat.tangent <- collect_X(model$model2)
  # Y.hat.tangent <- collect_Y(model$model2)
  # P.hat.tangent <- collect_P(model$model2)
  # Q.hat.tangent <- collect_Q(model$model2)
  #
  # X.affine.distance <- affine_dist(model$model1$X[[1]])
  # Y.affine.distance <- affine_dist(model$model1$Y[[1]])
  # X.euc.distance <- affine_dist(model$model2$X[[1]])
  # Y.euc.distance <- affine_dist(model$model2$Y[[1]])

  return(
    list(
      X.hat = X.hat,
      Y.hat = Y.hat,
      P.hat = P.hat,
      Q.hat = Q.hat,
      # X.hat.tangent = X.hat.tangent,
      # Y.hat.tangent = Y.hat.tangent,
      # P.hat.tangent = P.hat.tangent,
      # Q.hat.tangent = Q.hat.tangent,
      K = model$L,
      muX = model$muX,
      muY = model$muY
    )
  )
}
