#' affine_dist
#'
#' This is a function to calculate the geodesic distances
#'   on the cone of positive definite symmetric matrices.  Dependent
#'   packages are dplyr, purrr, parallel, geigen, Matrix.  See details for more information.
#'
#' @param df Subject array/list/df.
#' @param var.name String. Variable name containing correlations, df only, must be string.
#' @param regularise Logical. Do we want to regularise the matrices?  Default is FALSE.
#' @param lambda Numeric > 0.  Regularisation value.
#' @param cores Numeric. Can run in parallel using mclapply, number of cores to run on.  Defaults to 1.
#' @param ranks Numeric vector. A vector of length n containing the ranks of the matrices. Optional.
#' @param sub_ids Character vector.  A list of subject labels. Optional.
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
#'   if it is rank deficient or almost rank deficient (eigenvalues < 1e-7).
#'   Distances are then calculated as d( Q1, Q2 ) when Q1 and Q2 are not troublesome and
#'   d( Q1 + lambda x I, Q2 + lambda x I ) when Q1 OR Q2 are troublesome.
#'
#' @export affine_dist
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
#' df.distances.geo <- affine_dist( x = df, var.name = "cors",
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
affine_dist <- function(df, var.name = NULL, regularise = FALSE, lambda = 1, cores = 1, ranks = NULL, sub_ids = NULL, tol = 1e-7){

  ##################  Check correct inputs ###########################
  if(lambda <= 0){
    stop("You cannot regularise by a negative value of lambda.")
  }
  if(cores < 1 | cores > parallel::detectCores()){
    stop("You are trying to run on an incorrect number of cores. This must be an integer between 1 and the maximum number of available cores on your system.")
  }
  if(!is.logical(regularise) | is.na(regularise)){
    stop("The regularise variable should either be a TRUE or FALSE Boolean.")
  }


  ##################  Required packages ###########################
  ##################  Need to remove this.

  require(dplyr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(parallel, quietly = TRUE)
  require(geigen, quietly = TRUE)

  # ##################  Convert data to desired format ###########################
  #
  # # Converts input from an array or df into a list.
  # # Gets the number of subjects and the dimension of the correlation matrices
  # # Assumes all matrices are the same dimensions
  # if(is.array(df)){
  #   n <- dim(df)[3]
  #   roi.num <- dim(df)[1]
  #   x <- map(1:n, ~df[, , .])
  # }else if(is.data.frame(df)){
  #   if(is.null(var.name)){
  #     stop("When supplying a dataframe, please provide the variable name containing correlations")
  #   }
  #   n <- dim(df)[1]
  #   if("id" %in% colnames(df)){
  #     sub_ids <- df %>%
  #       slice(-n()) %>%
  #       pull(id) %>%
  #       as.character()
  #   }else{
  #     sub_ids <- 1:(n-1) %>%
  #       as.character()
  #   }
  #   roi.num <- dim(df[ , var.name][[1]][[1]])[1]
  #   x <- df[ , var.name] %>%
  #     pull()
  # }else if(is.list(df)){
  #   n <- length(df)
  #   roi.num <- dim(df[[1]])[1]
  #   x <- df
  # }else{
  #   stop("x must be an array, a list, or a dataframe containing the correlation values")
  # }

  ##################  Convert data to desired format ###########################

  tmp <- check_right_input(df, var.name)

  x <- tmp$x
  roi.num <- tmp$roi.num
  n <- tmp$n
  rm(tmp)


  ##################  Regularisation ###########################

  # Regularises the matrices by the lambda amount.
  # Regularisation is given by A -> A + lambda * I, where I is the roi x roi identity matrix

  if(regularise){
    x <- map(x, ~.x + diag(lambda, roi.num))
    ranks <- rep(roi.num, n)
  }



  ##################  Ranks ###########################

  # Calculate the ranks of the matrices
  # This is to deal with the "troublesome" matrices.
  # Matrices are troublesome if they are either rank deficient,
  # or very close to rank deficient.
  if(is.null(ranks)){
    ranks <- mclapply(x, Matrix::rankMatrix, tol = tol, mc.cores = cores) %>%
      unlist()
    if(any(ranks < roi.num)){
      low_rank_subjects <- which(ranks < roi.num)
      rank_message <- str_c(low_rank_subjects, collapse = ", ")
      warning(glue::glue("There are subjects with low rank functional connectivity profiles.  Regularised distances will be calculated for distances involving subjects {rank_message}."))
    }
  }

  ##################  Get distance matrix ###########################

  # We only calculate the lower triangle to save computation time
  # the following calculates a list of distances from position i to position n

  # If the matrices are troublesome, we will automatically regularise them.
  # i.e. if A or B is troublesome, we will find the distance
  # d( A + lambda*I, B + lambda*I ) instead of d( A, B )
  geo_dist <- mclapply(1:(n-1),
                        function(i){

                          if(ranks[i] == roi.num){
                            dist_vec <- map_dbl((i+1):n,
                                               function(k){

                                                 if(ranks[k] == roi.num){
                                                   eigs <- geigen(A = x[[i]],
                                                                   B = x[[k]],
                                                                   symmetric = T,
                                                                   only.values = TRUE)$values
                                                   dist <- sqrt(sum( log( eigs )^2))
                                                 }else{
                                                   eigs <- geigen(A = x[[i]] +  diag(lambda, roi.num),
                                                                   B = x[[k]] + diag(lambda, roi.num),
                                                                   symmetric = T,
                                                                   only.values = TRUE)$values
                                                   dist <- sqrt(sum(log(eigs )^2))
                                                 }

                                                 return(dist)
                                               } )
                          }else{
                            dist_vec <- map_dbl((i+1):n,
                                               function(k){

                                                 eigs <- geigen(A = x[[i]] +  diag( lambda, roi.num ),
                                                                 B = x[[k]] + diag( lambda, roi.num ),
                                                                 symmetric = T,
                                                                 only.values = TRUE)$values
                                                 dist <- sqrt(sum(log( eigs )^2))

                                                 return(dist)
                                               } )
                          }
                          return(dist_vec)
                        },
                        mc.cores = cores )

  if(is.data.frame(df) & "id" %in% colnames(df)){ # Pull out subject ID's if input is dataframe.
    sub_ids <- df %>%
      pull(id) %>%
      as.character()
  }

  if(is.null(sub_ids)){
    sub_ids <- 1:(n-1) %>%
      as.character()
  }
  # We need to convert from a list to a matrix to be able to convert it into a dist object
  # We need to add a column of 0's to make it a square matrix
  geo_dist <- 1:(n-1) %>%
    set_names(sub_ids[1:n-1]) %>%
    map(~c(rep(0, .x), geo_dist[[.x]])) %>%
    bind_cols() %>%
    mutate(n = 0) %>%
    as.dist()

  return(geo_dist)
}
