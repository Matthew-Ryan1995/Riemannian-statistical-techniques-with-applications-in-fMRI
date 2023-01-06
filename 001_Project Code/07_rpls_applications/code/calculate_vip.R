# calculate_vip <- function(Tmat, W, C, b){
#   p <- nrow(W)
#   K <- ncol(W)
#   q <- ncol(C)
#
#   res <- map(1:ncol(C),
#              function(i){
#                c <- C[, i]
#                denom <- sum(diag(t(Tmat) %*% Tmat) * b^2)
#
#                numer <- W^2 %*% (diag(t(Tmat) %*% Tmat) * b^2 * c^2)
#
#                ans <- sqrt(p * q * numer/denom)
#                return(ans)
#              })
#   res <- do.call(cbind, res)
#   colnames(res) <- colnames(C)
#
#   return(res)
# }
calculate_vip <- function(Tmat, W, Y){

  p <- nrow(W)

  RD_yi_tk <-  cor(Y, Tmat)^2
  RD_tk <- colMeans(RD_yi_tk)

  RD <- sum(RD_tk)

  vip_scores <- map_dbl(1:p,
                        function(j){
                          sqrt(p * sum(RD_tk * W[j, ]^2) / RD)
                        }
  )
  return(vip_scores)
}

calculate_vip_per_y <- function(Tmat, W, Y, C){

  p <- nrow(W)
  q <- nrow(C)

  RD_yi_tk <-  cor(Y, Tmat)^2
  RD_tk <- colMeans(RD_yi_tk)

  RD <- sum(RD_tk)

  vip_scores <- map(1:q,
                        function(i){
                          top_sum <-t(t(W^2) * C[i, ]^2) %*% RD_tk
                          sqrt(p * q * top_sum / RD)
                        }
  )
  vip_scores <- do.call(cbind, vip_scores)


  return(vip_scores)
}


