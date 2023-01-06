#' check_right_input
#'
#' @param df - data the user input, should be a list, array, or dataframe
#' @param var.name - variable name when df is a dataframe
#'
#' @return the data as a list
#'
#' @details Not for external use.  This function will convert the user input into a list.
#' It will also extract the number of ROI and the number of subject
#'
#' @export check_right_input
#'
#' @examples
#' df <- check_right_input(cobre, "cors")
check_right_input <- function(df, var.name = NULL){
  ##################  Convert data to desired format ###########################

  # Converts input from an array or df into a list.
  # Gets the number of subjects and the dimension of the correlation matrices
  # Assumes all matrices are the same dimensions
  if(is.array(df)){
    n <- dim(df)[3]
    roi.num <- dim(df)[1]
    x <- map(1:n, ~df[, , .])
  }else if(is.data.frame(df)){
    if(is.null(var.name)){
      stop("When supplying a dataframe, please provide the variable name containing correlations")
    }
    n <- dim(df)[1]
    roi.num <- dim(df[ , var.name][[1]][[1]])[1]
    x <- df[ , var.name] %>%
      pull()
  }else if(is.list(df)){
    n <- length(df)
    roi.num <- dim(df[[1]])[1]
    x <- df
  }else{
    stop("x must be an array, a list, or a dataframe containing the correlation values")
  }
  return(list(
    x = x,
    n = n,
    roi.num = roi.num
  ))
}

#' check_input_values
#'
#' This function throws an error if any of the parameters in Riemannian PLS
#' are not defined as they should be.  This is not for external use.
#'
#' @param t
#' @param u
#' @param W
#' @param C
#' @param X
#' @param Y
#' @param muX
#' @param muY
#' @param lambda_t
#' @param lambda_u
#' @param lambda_1
#' @param lambda_2
#' @param nu
#' @param tau
#'
#' @return
#' @export check_input_values
#'
#' @examples
check_input_values <- function(t, u, W, C, X, Y, muX, muY, lambda_t, lambda_u, lambda_1, lambda_2, nu, tau){
  # Check for missing data
  if(any(is.na(t))){
    stop("There are missing values in t.")
  }
  if(any(is.na(u))){
    stop("There are missing values in u.")
  }
  if(any(is.na(lambda_t))){
    stop("There are missing values in lambda_t.")
  }
  if(any(is.na(lambda_u))){
    stop("There are missing values in lambda_u.")
  }
  if(any(is.na(W))){
    stop("There are missing values in W.")
  }
  if(any(is.na(C))){
    stop("There are missing values in C.")
  }
  if(any(is.na(muX))){
    stop("There are missing values in muX.")
  }
  if(any(is.na(muY))){
    stop("There are missing values in muY.")
  }
  if(any(map_lgl(X, ~any(is.na(.x))))){
    stop("There are missing values in X.")
  }
  if(any(map_lgl(Y, ~any(is.na(.x))))){
    stop("There are missing values in Y.")
  }
  if(any(is.na(lambda_1))){
    stop("lambda_1 is missing.")
  }
  if(any(is.na(lambda_2))){
    stop("lambda_2 is missing.")
  }
  if(is.na(nu)){
    stop("nu is missing.")
  }
  if(is.na(tau)){
    stop("tau is missing.")
  }

  # Check right type
  if(!is.numeric(t)){
    stop("t should be a numeric vector.")
  }
  if(!is.numeric(u)){
    stop("u should be a numeric vector.")
  }
  if(!is.numeric(lambda_t)){
    stop("lambda_t should be a numeric vector.")
  }
  if(!is.numeric(lambda_u)){
    stop("lambda_u should be a numeric vector.")
  }
  if(!is.numeric(lambda_1)){
    stop("lambda_1 should be a numeric value.")
  }
  if(!is.numeric(lambda_2)){
    stop("lambda_2 should be a numeric value.")
  }

  if(any(!is.matrix(W), !is.numeric(W))){
    stop("W should be a numeric matrix.")
  }
  if(any(!is.matrix(muX), !is.numeric(muX))){
    stop("muX should be a numeric matrix.")
  }
  if(any(!is.matrix(C), !is.numeric(C))){
    stop("C should be a numeric matrix.")
  }
  if(any(!is.matrix(muY), !is.numeric(muY))){
    stop("muY should be a numeric matrix.")
  }

  if(any(!is.numeric(nu), nu<=0)){
    stop("nu should be a positive number.")
  }
  if(any(!is.numeric(tau), tau<0)){
    stop("tau should be a non-negative number.")
  }

  if(any(map_lgl(X, ~!is.matrix(.x)), map_lgl(X, ~!is.numeric(.x)))){
    stop("X must be a list of numeric matrices")
  }
  if(any(map_lgl(Y, ~!is.matrix(.x)), map_lgl(Y, ~!is.numeric(.x)))){
    stop("Y must be a list of numeric matrices")
  }


  # Check right sizes.
  if(length(t) != length(u)){
    stop("t and u should be the same length.")
  }
  if(length(lambda_t) != length(t)){
    stop("The number of constraints on t does not much the number of t-values.")
  }
  if(length(lambda_u) != length(u)){
    stop("The number of constraints on u does not much the number of u-values.")
  }
  if(length(lambda_1) != 1){
    stop("lambda_1 should be a single value.")
  }
  if(length(lambda_2) != 1){
    stop("lambda_2 should be a single value.")
  }

  if(any(dim(W) != dim(muX))){
    stop("W and muX should be the same dimension.")
  }
  if(any(dim(C) != dim(muY))){
    stop("C and muY should be the same dimension.")
  }

  if(any(map_lgl(X, ~{any(dim(W) != dim(.x))}))){
    stop("W and X should be the same dimension.")
  }
  if(any(map_lgl(Y, ~{any(dim(C) != dim(.x))}))){
    stop("C and Y should be the same dimension.")
  }

}


#' check_right_input_loadings
#'
#' Check that the loading inputs are correct
#'
#' @param X
#' @param muX
#' @param P
#' @param t
#' @param ...
#'
#' @return
#' @export check_right_input_loadings
#'
#' @examples
#'
#' TBA
check_right_input_loadings <- function(X, muX, P, t, ...){

  ARGS <- list(...)

  if(is.null(ARGS$tau)){
    tau <- 1/2
  }else{
    tau <- ARGS$tau
  }

  # Check missing

  if(any(map_lgl(X, ~any(is.na(.x))))){
    stop("There are missing values in X.")
  }

  if(any(is.na(muX))){
    stop("There are missing values in muX.")
  }

  if(any(is.na(P))){
    stop("There are missing values in P.")
  }

  if(any(is.na(t))){
    stop("There are missing values in t.")
  }

  if(is.na(tau)){
    stop("tau is missing.")
  }

  # Check input type
  if(any(map_lgl(X, ~!is.matrix(.x)), map_lgl(X, ~!is.numeric(.x)))){
    stop("X must be a list of numeric matrices")
  }

  if(any(!is.matrix(muX), !is.numeric(muX))){
    stop("muX should be a numeric matrix.")
  }

  if(any(!is.matrix(P), !is.numeric(P))){
    stop("P should be a numeric matrix.")
  }

  if(!is.numeric(t)){
    stop("t should be a numeric vector.")
  }

  if(any(!is.numeric(tau), tau<0)){
    stop("tau should be a non-negative number.")
  }

  # Check right sizes
  if(length(t) != length(X)){
    stop("t and X should be the same length.")
  }
  if(any(dim(P) != dim(muX))){
    stop("P and muX should be the same dimension.")
  }

  if(any(map_lgl(X, ~{any(dim(P) != dim(.x))}))){
    stop("P and X should be the same dimension.")
  }

}

