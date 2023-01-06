#' Calculate scores in Euclidean NIPALS
#'
#' @param X Matrix: matrix of data values of size n x p
#' @param W Numeric: Vector of weights of size p x 1
#' @param ...
#'
#' @return Numeric vector of length n containing the scores for X in the direction W.
#' @export calculate_score_euclidean
#'
#' @examples
#' TBA
calculate_score_euclidean <- function(X, W, ...){
  return(
    X %*% W
  )
}

#' Calculate the direction of interest in Euclidean NIPALS
#'
#' @param X Matrix: matrix of data values of size n x p
#' @param u Numeric: Current estimate of scores for Y
#' @param ...
#'
#' @return Numeric vector of size p X 1, it is a direction in X-space
#' @export calculate_direction_euclidean
#'
#' @examples
#' TBA
calculate_direction_euclidean <- function(X, u, ...){
  w <- t(X) %*% u / sum(u^2)
  return(w)
}

#' Calculate X-scores, Euclidean
#'
#' @param X Matrix: matrix of data values of size n x p
#' @param W Numeric: Vector of weights of size p x 1
#' @param ...
#'
#' @return Numeric: Vector of weights of size p x 1
#' @export calculate_t_euclidean
#'
#' @examples
#' TBA
calculate_t_euclidean <- function(X, W, ...){
  t <- calculate_score_euclidean(X = X, W = W)
  return(t)
}

#' Calculate Y-scores, Euclidean
#'
#' @param Y Matrix: matrix of data values of size n x q
#' @param C Numeric: Vector of weights of size q x 1
#' @param ...
#'
#' @return Numeric: Vector of weights of size q x 1
#' @export calculate_u_euclidean
#'
#' @examples
#' TBA
calculate_u_euclidean <- function(Y, C, ...){
  u <- calculate_score_euclidean(X = Y, W = C)
  return(u)
}

#' Calculate X-weights, Euclidean
#'
#' @param X Matrix: matrix of data values of size n x p
#' @param u Numeric: Scores for Y of size n x 1
#' @param ...
#'
#' @return Numeric vector of size p X 1, it is a direction in X-space
#' @export calculate_w_euclidean
#'
#' @examples
#' TBA
calculate_w_euclidean <- function(X, u, ...){
  w <- calculate_direction_euclidean(X = X, u = u)
  # Scale to length 1
  w <- w/sqrt(sum(w^2))
  return(w)
}

#' Calculate Y-weights, Euclidean
#'
#' @param Y Matrix: matrix of data values of size n x q
#' @param t Numeric: Scores for X of size n x 1
#' @param ...
#'
#' @return Numeric vector of size q X 1, it is a direction in Y-space
#' @export calculate_c_euclidean
#'
#' @examples
#' TBA
calculate_c_euclidean <- function(Y, t, ...){
  c <- calculate_direction_euclidean(X = Y, u = t)
  # Scale to length 1
  c <- c/sqrt(sum(c^2))
  return(c)
}

#' Calculate X-loadings, Euclidean
#'
#' @param X Matrix: matrix of data values of size n x p
#' @param t Numeric: Scores for X of size n x 1
#' @param ...
#'
#' @return Numeric vector of size p X 1, it is a direction in X-space
#' @export calculate_p_euclidean
#'
#' @examples
#' TBA
calculate_p_euclidean <- function(X, t, ...){
  p <- calculate_direction_euclidean(X = X, u = t)
  return(p)
}

#' Calculate Y-loadings, Euclidean
#'
#' @param Y Matrix: matrix of data values of size n x q
#' @param u Numeric: Scores for U of size n x 1
#' @param ...
#'
#' @return Numeric vector of size q X 1, it is a direction in X-space
#' @export calculate_q_euclidean
#'
#' @examples
#' TBA
calculate_q_euclidean <- function(Y, u, ...){
  q <- calculate_direction_euclidean(X = Y, u = u)
  return(q)
}

#' Deflate the X-values in Euclidean space
#'
#' @param X Matrix: matrix of data values of size n x p
#' @param t Numeric: Scores for X of size n x 1
#' @param P Numeric: Vector of loadings of size p x 1
#' @param ...
#'
#' @return A new set of X-values that have removed the information from the previous direction
#' @export deflate_X_euclidean
#'
#' @examples
#' TBA
deflate_X_euclidean <- function(X, t, P, ...){
  X_deflate <- X - t %*% t(P)

  return(X_deflate)
}

#' Deflate the Y-values in Euclidean space
#'
#' @param Y Matrix: matrix of data values of size n x q
#' @param C Numeric: Vector of weights of size q x 1
#' @param B Linear model: A linear model fitting the Y-scores to the X-scores.  The inner relation
#' @param ...
#'
#' @return A new set of Y-values that have removed the information from the previous direction
#' @export deflate_Y_euclidean
#'
#' @examples
#' TBA
deflate_Y_euclidean <- function(Y, C, B, ...){

  preds <- fitted(B)

  Y_deflate <- Y - preds %*% t(C)

  return(Y_deflate)
}
