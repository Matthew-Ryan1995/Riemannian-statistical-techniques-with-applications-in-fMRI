#' fitted.riemannian_pls
#' Returns fitted values from riemannian pls model
#'
#' @param M
#' @param ...
#'
#' @return TBA
#' @export fitted.riemannian_pls
#'
#' @examples
#' TBA
fitted.riemannian_pls <- function(M, ...){

  # Get output from model
  ARGS <- list(...)

  if(!is.null(ARGS$num_comp)){
    L <- ARGS$num_comp
  }else{
    L <- M$L
  }


  # Get predicted u values
  regs <- map(1L:L, ~M$reg_steps[[.x]])
  # Get predictions from inner relation
  pred_scores <- map(regs,
                     function(m){
                       model = m$M
                       return(
                         fitted(model)
                       )
                     }
  )

  if(isTRUE(M$type_y == "euclidean")){
    C <- do.call(cbind, M$weightsY)
    C <- C[, 1L:L]


    u_vals <- do.call(cbind, pred_scores)

    preds <- u_vals %*% t(C)

    return(preds)
  }

  muY <- M$muY
  mu.sq <- get_sqrt(muY)

  n <- length(M$Y[[1]])

  Q <- map(1:L, ~M$loadingsY[[.x]])
  # regs <- map(1:L, ~M$reg_steps[[.x]])
  #
  #
  # # Get predictions from inner relation
  # pred_scores <- map(regs,
  #                    function(m){
  #                      model = m$M
  #                      return(
  #                        fitted(model)
  #                      )
  #                    }
  # )

  pred1 <- map2(Q, pred_scores,
                function(q, t){
                  map(t, ~.x * q) # Multiply each predicted score by corresponding loading
                })

  pred2 <- map(1:n, # Convert list to subject wise list, list element i contains all information for subject i.
               function(i){
                 map(pred1,
                     function(l){
                       l[[i]]
                     }
                 )
               })

  pred3 <- map(pred2,  # Add across individual subjects
               function(l){
                 Reduce("+", l)
               })

  # Map down with Exponential map
  final_pred <- map(pred3,
                    function(q){
                      affine_exp(P = muY, W = q, P.sq = mu.sq)
                    })

  return(final_pred)
}

#' rmse.riemannian_pls
#' Calculates RMSE from Riemannian PLS model
#'
#' @param truth - True data
#' @param est - Estimated data
#' @param ...
#'
#' @return TBA
#' @export rmse.riemannian_pls
#'
#' @examples
#' TBA
rmse.riemannian_pls <- function(truth, est, ...){

  ARGS <- list(...)

  if(is.null(ARGS$type_y) | isTRUE(ARGS$type_y == "affine") ){
    if(any(map_lgl(truth, ~!is.matrix(.x)), map_lgl(truth, ~!is.numeric(.x)))){
      stop("Since truth is affine, it must be a list of numeric matrices")
    }
    if(any(map_lgl(est, ~!is.matrix(.x)), map_lgl(est, ~!is.numeric(.x)))){
      stop("Since est is affine, it must be a list of numeric matrices")
    }
    squared_dist <- map2_dbl(truth, est,
                             function(t, e){
                               affine_dist(list(t, e))^2
                             })
    mse <- mean(squared_dist)
  }else if (isTRUE(ARGS$type_y == "euclidean")){
    if(!(is.matrix(truth) | is.matrix(est))){
      stop("Since truth and estimate are Euclidean, they should be matrices.")
    }
    squared_dist <- sum((truth - est)^2)
    mse <- squared_dist / nrow(truth)
  }else{
    warning("Unrecognised data type")
    return(NULL)
  }

  return(
    sqrt(mse)
  )
}

#' rsq.riemannian_pls
#' Returns R^2 value for Riemannian PLS model
#'
#' @param truth - True data
#' @param est - Estimated data
#' @param muY - True muY
#' @param m - column mean if euclidean
#' @param ignore - Logical, should I calculate means anyway
#' @param ...
#'
#' @return TBA
#' @export rsq.riemannian_pls
#'
#' @examples
#' TBA
rsq.riemannian_pls <- function(truth, est, muY = NULL, m = NULL, ignore = TRUE,
                               test = FALSE,
                               ...){

  ARGS <- list(...)

  if(is.null(ARGS$type_y) | isTRUE(ARGS$type_y == "affine") ){
    if(any(map_lgl(truth, ~!is.matrix(.x)), map_lgl(truth, ~!is.numeric(.x)))){
      stop("Since truth is affine, it must be a list of numeric matrices")
    }
    if(any(map_lgl(est, ~!is.matrix(.x)), map_lgl(est, ~!is.numeric(.x)))){
      stop("Since est is affine, it must be a list of numeric matrices")
    }
    if(is.null(muY)){
      if(isTRUE(ignore)){
        muY = get_frechet_mean(truth)
      }else{
        stop("Please enter the Frechet mean muY.")
      }
    }
    squared_dist <- map2_dbl(truth, est,
                             function(t, e){
                               affine_dist(list(t, e))^2
                             })
    ssd <- sum(squared_dist)

    squared_error <- map_dbl(truth,
                             function(t){
                               affine_dist(list(t, muY))^2
                             })
    sse <- sum(squared_error)

    squared_model_error <- map_dbl(est,
                                   function(t){
                                     affine_dist(list(t, muY))^2
                                   })
    ssr <- sum(squared_model_error)
  }else{
    if(!(is.matrix(truth) | is.matrix(est))){
      stop("Since truth and estimate are Euclidean, they should be a matrices.")
    }
    # if(is.null(m)){
    #   if(isTRUE(ignore)){
        m <- colMeans(truth)
    #   }else{
    #     stop("Please enter the mean of the response data")
    #   }
    # }
    ssd <- sum((truth - est)^2)

    M <- matrix(m, nrow = nrow(truth), ncol = ncol(truth), byrow = TRUE)

    sse <- sum((truth - M)^2)
    ssr <- sum((est - M)^2)
  }


  return(
    if(isTRUE(test)){
      ssr/sse
    }else{
      1 - ssd/sse
    }
  )
}


## Need to add in tolerance, iterations, parallel over newdata
#' predict.riemannian_pls
#'
#' @param object - the Riemannian PLS model
#' @param newdata - New X data to predict on
#' @param tol - Tolerance at which to calc t-values
#' @param max.iter - Max iter at which to calc t-values
#' @param tau - Step size
#' @param ...
#'
#' @return TBA
#' @export predict.riemannian_pls
#'
#' @examples
#' TBA
predict.riemannian_pls <- function(object, newdata, tol = 1e-5, max.iter = 50, tau = 1/2, ...){

  ARGS <- list(...)
  M <- object

  if(!is.null(ARGS$num_comp)){
    L <- ARGS$num_comp
  }else{
    L <- M$L
  }

  if(is.matrix(newdata)){
    n <- nrow(newdata)
  }else{
    n <- length(newdata)
  }

  if(is.null(mc.cores <- ARGS$mc.cores)){
    mc.cores <- 1
  }

  if(isTRUE(ARGS$method == "tangent")){
    # if(isTRUE(ARGS$type_x == "affine") | is.null(ARGS$type_x)){
      # if(is.null(muX)){
      #   muX <- get_frechet_mean(X, cores = mc.cores)
      # }
      muX <- object$muX
      muX.sq <- get_sqrt(muX)
      newdata <- linearise_data(X = newdata, mu = muX,
                                mc.cores = mc.cores, mu.sq = muX.sq) ## WRITE THIS FUNCTION
      # ARGS$type_x <- "euclidean"
    # }
  }


  if(!is.null(M$type_x)){
    M$type_x <- str_to_lower(M$type_x)
    if(M$type_x == "affine"){
      if(any(map_lgl(newdata, ~!is.matrix(.x)), map_lgl(newdata, ~!is.numeric(.x)))){
        stop("Since newdata is affine, it must be a list of numeric matrices")
      }
      calculate_t <- calculate_t_affine
    }else if(M$type_x == "euclidean"){
      if(!is.matrix(newdata)){
        stop("Since newdata is Euclidean, newdata should be a matrix.")
      }
      calculate_t <- calculate_t_euclidean
    }else{
      stop("We can only currently consider predictors of either affine type or euclidean type.")
    }
  }else{
    if(any(map_lgl(newdata, ~!is.matrix(.x)), map_lgl(newdata, ~!is.numeric(.x)))){
      stop("Since X is affine, it must be a list of numeric matrices")
    }
    calculate_t <- calculate_t_affine
  }

  # Review this
  # t_values <- map(newdata,
  #                 function(x){
  #                   t <- numeric()
  #                   for(i in 1:L){
  #                     t[i] <- calculate_t(mu = M$muX, X = list(x), W = M$weightsX[[i]], t = 1, tau = tau,
  #                                         max.iter = max.iter, tol = tol, ...)
  #                     x <- deflate_X(X = list(x), t = t[i], P = M$loadingsX[[i]], muX = M$muX,
  #                                    type_x = ARGS$type_x, ...)[[1]]
  #                   }
  #                   return(t)
  #                 })

  t_values <- matrix(0, ncol = L, nrow = n)

  x <- newdata
  for(i in 1L:L){
    t_values[, i] <- calculate_t(mu = M$muX, X = x, W = M$weightsX[[i]], t = rep(1, n), tau = tau,
                                 max.iter = max.iter, tol = tol, ...)
    x <- deflate_X(X = x, t = t_values[, i], P = M$loadingsX[[i]], muX = M$muX,
                   type_x = M$type_x, ...)
  }

  u_preds <- matrix(0, ncol = L, nrow = n)

  for(i in 1L:L){
    u_preds[, i] <- predict(M$reg_steps[[i]]$M, newdata = data.frame(t = t_values[, i]))
  }

  # u_preds <- map(t_values,
  #                function(t){
  #                  map_dbl(1:L, function(l){
  #                    predict(M$reg_steps[[l]]$M, newdata = data.frame(t = t[l]))
  #                  })
  #                })

  if(is.null(M$type_y) | isTRUE(M$type_y == "affine") ){
    u_preds <- map(1L:nrow(u_preds), ~u_preds[.x, ])

    Q <- map(1L:L, ~M$loadingsY[[.x]])
    Q_preds <- map(u_preds,
                   function(u){
                     q1 <- map2(Q, u,
                                function(q, u){
                                  q * u
                                }
                     )
                     return(
                       Reduce("+", q1)
                     )
                   })

    muY <- M$muY
    mu.sq <- get_sqrt(muY)

    pred_vals <- map(Q_preds,
                     function(q){
                       affine_exp(P = muY, W = q, P.sq = mu.sq)
                     })

  }else if (isTRUE(M$type_y == "euclidean")){
    C <- map(1L:L, ~M$weightsY[[.x]])
    C <- do.call(cbind, C)

    pred_vals <- u_preds %*% t(C)
  }else{
    return(NULL)
  }


  return(pred_vals)
}

## To review this function
# get_metrics <- function(object, truth = NULL, est = NULL, newdata = NULL, type = "training",
#                         tol = 1e-5, max.iter = 50, tau = 1/2, ...){
#
#   ARGS <- list(...)
#
#   if(!is.null(ARGS$num_comp)){
#     L <- ARGS$num_comp
#   }else{
#     L <- object$L
#   }
#
#   if(isTRUE(type == "training")){
#     if(is.null(truth)){
#       truth <- object$Y[[1]]
#     }
#     if(is.null(est)){
#       est <- fitted.riemannian_pls(object, ...)
#     }
#   }else if(isTRUE(type == "testing")){
#     if(is.null(truth)){
#       stop("Enter the true response values of the testing set")
#     }
#     if(is.null(est)){
#       if(is.null(newdata)){
#         stop("Either input the estimate values, or the new data to predict on.")
#       }
#       est <- predict.riemannian_pls(object, newdata, tol = tol, max.iter = max.iter, tau = tau, ...)
#     }
#   }else{
#     stop("type must be testing or training.")
#   }
#
#   rmse <- rmse.riemannian_pls(truth = truth, est = est,
#                               type_x = object$type_x,type_y = object$type_y, ...)
#   rsq <- rsq.riemannian_pls(truth = truth, est = est, muY = object$muY,
#                             type_x = object$type_x, type_y = object$type_y, ...)
#
#   return(
#     tibble(
#       type = type,
#       rmse = rmse,
#       r.sq = rsq,
#       L = L
#     )
#   )
# }
