#' euclid_dist
#'
#' This is a function to calculate the Euclidean distance between correlation matrices.  Dependent
#'   packages are dplyr, purrr.  See details for more information.
#'
#' @param df Subject array/list/df.
#' @param var.name Variable name containing correlations, df only, must be string.
#' @param regularise Do we want to regularise the matrices.  Default is FALSE.
#' @param lambda Regularisation value.
#'
#' @return A dist object containing the Euclidean distances.
#'
#' @details This will take a list of correlation matrices, take the upper triangle of each,
#' and calculates the Euclidean distance between subjects based on this vector.
#'
#' @export euclid_dist
#'
#' @examples
#'
#'
#'
euclid_dist <- function(df, var.name = NULL, sub_ids = NULL){

  ##################  Required packages ###########################

  require(dplyr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  # ##################  Convert data to desired format ###########################
  tmp <- check_right_input(df, var.name)

  x <- tmp$x
  roi.num <- tmp$roi.num
  n <- tmp$n
  rm(tmp)

  ##################  Get distance matrix ###########################

  # Much easier for Euclidean distance
  # convert to a matrix that is n x p
  # Use dist function
  x <- x %>%
    map(~.x[upper.tri(.x)]) %>%
    simplify2array() %>%
    t()

  if(nrow(x) == 1){
    x <- t(x)
  }

  if(is.data.frame(df) & "id" %in% colnames(df)){ # Pull out subject ID's if input is dataframe.
    sub_ids <- df %>%
      pull(id) %>%
      as.character()
  }

  if(is.null(sub_ids)){
    sub_ids <- 1:n %>%
      as.character()
  }

  euc_dist <- x %>%
    dist() %>%
    usedist::dist_setNames(sub_ids)

  return(euc_dist)
}

#' corr_dist
#'
#' This is a function to calculate the Euclidean distance between correlation matrices.  Dependent
#'   packages are dplyr, purrr.  See details for more information.
#'
#' @param df Subject array/list/df.
#' @param var.name Variable name containing correlations, df only, must be string.
#' @param sub_ids A character list of subject id's to keep track of.  Defaults to NULL.
#'
#' @return A dist object containing the Correlation distances.
#'
#' @details This will take a list of correlation matrices, take the upper triangle of each,
#' and calculates the correlation distance between subjects based on this vector as defined by
#'
#'    d(x, y) = \sqrt{2(1 - \rho(x, y))}
#'
#' where \rho(x, y) is the correlation between vectors x and y.
#'
#' @export corr_dist
#'
#' @examples
#'
#'
#'
corr_dist <- function(df, var.name = NULL, sub_ids = NULL){

  ##################  Required packages ###########################

  require(dplyr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  # ##################  Convert data to desired format ###########################
  tmp <- check_right_input(df, var.name)

  x <- tmp$x
  roi.num <- tmp$roi.num
  n <- tmp$n
  rm(tmp)



  ##################  Get distance matrix ###########################

  # Much easier for Euclidean distance
  # convert to a matrix that is n x p
  # Use dist function
  x <- x %>%
    map(~.x[upper.tri(.x)]) %>%
    simplify2array()

  if(is.null(ncol(x))){
    stop("You cannot calculate the correlation distance for 2 x 2 matrices.")
  }

  if(is.data.frame(df) & "id" %in% colnames(df)){ # Pull out subject ID's if input is dataframe.
    sub_ids <- df %>%
      pull(id) %>%
      as.character()
  }

  if(is.null(sub_ids)){
    sub_ids <- 1:n %>%
      as.character()
  }

  corr_dist <- x %>%
    cor()
  corr_dist <- sqrt(2 * (1 - corr_dist))
  corr_dist <- corr_dist %>%
    as.dist() %>%
    usedist::dist_setNames(sub_ids)

  return(corr_dist)
}
