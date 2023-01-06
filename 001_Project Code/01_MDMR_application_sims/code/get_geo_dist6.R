
## This is a function to calculate the geodesic distances 
## on the cone of positive definite symmetric matrices

## 

#' get_geo_dist
#'
#' This is a function to calculate the geodesic distances
#'   on the cone of positive definite symmetric matrices.  Dependent 
#'   packages are dplyr, purrr, parallel, geigen, Matrix.  See details for more information.
#'
#' @param x Subject array/list/df.
#' @param var.name Variable name containing correlations, df only, must be string.
#' @param regularise Do we want to regularise the matrices.  Default is TRUE.
#' @param lambda Regularisation value.
#' @param cores Can run in parallel using mclapply, number of cores to run on.  Defualts to 1.
#' @param ranks A vector of length n containing the ranks of the matrices. Optional.
#'
#' @return A dist object containing the geodesic distances.
#' 
#' @details For two PSD matrices Q1 and Q2, we calculate Q = Q1^{-1/2} Q2 Q1^{-1/2}.
#'   Then the distance d( Q1, Q2 ) = sum( log( lambda )^2 ), where lambda are the
#'   eigenvalues of Q.  This amounts to finding the generalised eigenvalues of 
#'   Q1 and Q2. 
#'   
#'   When regularise = TRUE, we translate all matrices Q to Q + lambda x I, 
#'   where I is the identity matrix.
#'   
#'   When regularise = FALSE, we only translate the "troublesome" matrices.  A matrix is troublesome
#'   if it is rank defficient or almost trank defficient (eigenvalues < 1e-7).
#'   Distances are then caluculated as d( Q1, Q2 ) when Q1 and Q2 are not troublesome and
#'   d( Q1 + lambda x I, Q2 + lambda x I ) when Q1 OR Q2 are troublesome.
#' @export
#'
#' @examples
#' 
#' require(tidyverse) 
#' 
#' #### Making some random normal data, calculating the correlations
#' 
#' sig <- matrix( c(1, 0.4, -0.1,
#'                  0.4, 1, 0.5,
#'                  -0.1, 0.5, 1),
#'                nrow = 3,
#'                ncol = 3 )
#' 
#' mu <- 1:3
#' 
#' df <- tibble( id = 1:171,
#'               obs = map(id, ~MASS::mvrnorm(n = 30, mu - mu, Sigma = sig ) ),
#'               cors = map( obs, cor ) )
#' 
#' #### Get the geodisic distance with my function
#' 
#' df.distances.geo <- get_geo_dist( x = df, var.name = "cors",
#'                                   regularise = FALSE, lambda = 1,  cores = 1 )
#' 
#' 
#' #### Transform the data into the upper triangle to let me use euclidean distance
#' df.Y <- df %>% 
#'   transmute( id,
#'             upptri = map(cors, ~.x[upper.tri(.x)]),
#'              axes = map(1:n(), ~c("x", "y", "z"))) %>% 
#'   tidyr::unnest(cols = everything()) %>% 
#'   tidyr::pivot_wider( names_from = "axes", values_from = "upptri" )
#' 
#' #### Get Euclidean distance to see how they are different.
#' 
#' df.distances.euc <- dist( df.Y %>% select( -id ) )
#' 
get_geo_dist <- function( x, var.name = NULL, regularise = TRUE, lambda = 1, cores = 1, ranks = NULL ){
  
  ##################  Required packages ###########################
  
  require(dplyr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(parallel, quietly = TRUE)
  require(geigen, quietly = TRUE)
  
  ##################  Convert data to desired format ###########################
  
  # Converts input from an array or df into a list.
  # Gets the number of subjects and the dimension of the correlation matrices
  # Assumes all matrices are the same dimensions
  if( is.array(x) ){
    n <- dim( x )[3]
    roi.num <- dim( x )[1]
    x <- map(1:n, ~x[,,.])
  }else if( is.data.frame( x ) ){
    if( is.null( var.name ) ){
      stop("When supplying a dataframe, please provide the variable name containing correlations")
    }
    n <- dim( x )[1]
    roi.num <- dim( x[ , var.name ][[1]][[1]] )[1]
    x <- x[ , var.name ] %>% 
      pull() 
  }else if( is.list( x ) ){
    n <- length( x )
    roi.num <- dim( x[[1]] )[1]
    x <- x 
  }else{
    stop("x must be an array, a list, or a dataframe containing the correlation values")
  }
  
  
  ##################  Regularisation ###########################
  
  # Regularises the matrices by the lambda amount.
  # Regularisation is given by A -> A + lambda * I, where I is the roi x roi identity matrix
  
  if( regularise ){
    x <- map( x, ~.x + diag( lambda, roi.num ) )
    ranks <- rep( roi.num, n )
  }
  
  
  
  ##################  Ranks ###########################
  
  # Calculate the ranks of the matrices
  # This is to deal with the "troublesome" matrices.
  # Matrices are troublesome if they are either rank defficient, 
  # or very close to rank defficient.
  if( is.null( ranks ) ){
    ranks <- mclapply(x, Matrix::rankMatrix, tol = 1e-7, mc.cores = cores ) %>% 
      unlist()
  }
  
  ##################  Get distance matrix ###########################
  
  # We only calculate the lower triangle to save computation time
  # the following calculates a list of distances from position i to position n
  
  # If the matrices are troublesome, we will automatically regularise them.
  # i.e. if A or B is troublesome, we will find the distance  
  # d( A + lambda*I, B + lambda*I ) instead of d( A, B )
  geo_dist <- mclapply( 1:(n-1), 
                        function( i ){ 
                          
                          if( ranks[i] == roi.num ){
                            dist_vec<-map_dbl( (i+1):n, 
                                               function(k){
                                                 
                                                 if(ranks[k] == roi.num){
                                                   eigs <- geigen( A = x[[i]], 
                                                                   B = x[[k]], 
                                                                   symmetric = T, 
                                                                   only.values = TRUE )$values
                                                   dist <- sqrt( sum( log( eigs )^2 ) )
                                                 }else{
                                                   eigs <- geigen( A = x[[i]] +  diag( lambda, roi.num ), 
                                                                   B = x[[k]] + diag( lambda, roi.num ), 
                                                                   symmetric = T, 
                                                                   only.values = TRUE )$values
                                                   dist <- sqrt( sum( log( eigs )^2 ) )
                                                 }
                                                 
                                                 return(dist)
                                               } )
                          }else{
                            dist_vec<-map_dbl( (i+1):n, 
                                               function(k){
                                                 
                                                 eigs <- geigen( A = x[[i]] +  diag( lambda, roi.num ), 
                                                                 B = x[[k]] + diag( lambda, roi.num ), 
                                                                 symmetric = T, 
                                                                 only.values = TRUE )$values
                                                 dist <- sqrt( sum( log( eigs )^2 ) )
                                                 
                                                 return(dist)
                                               } )
                          }
                          return(dist_vec)
                        }, 
                        mc.cores = cores )
  
  # We need to convert from a list to a matrix to be able to convert it into a dist object
  # We need to add a column of 0's to make it a square matrix
  geo_dist <- map(1:(n-1), ~c( rep(0, .x), geo_dist[[.x]] )) %>% 
    bind_cols() %>% 
    mutate( n = 0 ) %>% 
    as.dist()
  
  return(geo_dist)
}
