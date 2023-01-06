#' is_positive_definite
#'
#' Check to see if a matrix is positive definite.
#'
#' @param P - A matrix
#' @param tol - A tolerance level
#'
#' @return Logical: TRUE if the matrix is positive definite, FALSE otherwise
#'
#' @details Uses the eigen-value decomposition to determine is a matrix is positive definite up to a certain tolerance.
#' @export is_positive_definite
#'
#' @examples
is_positive_definite <- function(P, tol = 1e-10){
  eigenvalues <- eigen(P, only.values = TRUE)$values
  pos_eig <- map_dbl(1:length(eigenvalues),
                     function(l){
                       e <- eigenvalues[l]
                       if(abs(e) < 1e-10){
                         return(0)
                       }
                       return(e)
                     }
  )

  if(any(pos_eig <= 0)){
    return(FALSE)
  }else{
    return(TRUE)
  }
}


#' get_sqrt
#'
#' This function will obtain the square root of a matrix for you.  Currently two choices are available: The symmetric
#' sqaure root or Cholesky square root.  This is not optomised for large or poorly behaved matrices.
#'
#' @param A - A square matrix to find the square root of.
#' @param symmetric - Logical: Is the matrix A symmetric or not.  Defaults to isSymmetric(A)
#' @param method - Character: Either custom or expm, to choose how you want to calculate the square root.
#'
#' @return A matrix B such that B^2 = A (or B^T B = A in the non-symmetric case).
#' @export get_sqrt
#'
#' @examples
#'
#' ## Symmetric square root
#' A <- matrix(
#' c(
#'   1, 0.5,
#'   0.5, 1)
#'  ),
#'  2, 2
#')
#'
#' get_sqrt(A)
#'
#' ## Non-symmetric square root
#'
#' get_sqrt(A, symmetric = FALSE)
#'
#'
get_sqrt <- function(A, symmetric = isSymmetric(A, tol = sqrt(.Machine$double.eps)), method = c("custom", "expm")){
  tmp.names <- colnames(A)
  method <- method[1]
  if(method == "custom"){
    if(symmetric){
      # tmp.names <- colnames(A)
      A_eig <- eigen(A, symmetric = TRUE)

      if(any(A_eig$values < 0)){
        warning("Numeric instabilities caused negative eigenvalues, using methods from expm")
        A.sqr <- expm::sqrtm(A) %>% Re() # Can get numerically complex square root, hopefully this doesn't bite me.
      }else{
        A.sqr <-  tcrossprod(tcrossprod(A_eig$vectors, diag(sqrt(A_eig$values))), A_eig$vectors)
      }
      # colnames(A.sqr) <- rownames(A.sqr) <- tmp.names
    }else{
      A.sqr <- chol(A)
    }
  }else if(method == "expm"){
    A.sqr <- expm::sqrtm(A) %>% Re()
  }else{
    stop("Please enter a method that is either custom or expm.")
  }

  colnames(A.sqr) <- rownames(A.sqr) <- tmp.names

  return(A.sqr)
}

#' get_trace
#'
#' Calculate the trace if A
#'
#' @param A - A square matrix
#'
#' @return  The trace of A
#' @export get_trace
#'
#' @examples
get_trace <- function(A){
  if(nrow(A) != ncol(A)){
    stop("A must be a square matrix.")
  }

  return(sum(diag(A)))
}


#' double_centre
#'
#' This function will double centre the input matrix, that is, it will make the row and column means 0.
#'
#' @param X Matrix: A square, numeric matrix.
#'
#' @return A double-centred version of X.
#'
#' @details Returns the matrix J X J where J is the centring matrix
#'
#' J = I - 1/n 1 t(1).
#'
#' @export double_centre
#'
#' @examples
#'
#' set.seed(2022)
#' X <- matrix(rnorm(9), 3, 3)
#'
#' centred_X <- double_centre(X)
#'
#' all.equal(apply(centred_X, 1, mean), rep(0, 3))
#' all.equal(apply(centred_X, 2, mean), rep(0, 3))
#'
double_centre <- function(X){
  if(!is.matrix(X)){
    stop("X must be a square numeric matrix")
  }

  if(!is.numeric(X)){
    stop("X must be a square numeric matrix")
  }

  if(nrow(X) != ncol(X)){
    stop("X must be a square numeric matrix")
  }

  if(any(is.na(X))){
    stop("I don't know how to deal with missing values, please forgive me.")
  }

  n <- nrow(X)

  X_row <- matrix(rowMeans(X), n, n)
  X_col <- t(matrix(colMeans(X), n, n))
  X_mean <- mean(X)

  return(
    X - X_row - X_col + X_mean
  )
}



#' multiply_and_reduce
#'
#' Take a set of numeric values and a list of matrices, multiply them, and add them up.
#'
#' @param t Numeric: a numeric vector of coefficient values
#' @param W List: A list of numeric matrices
#'
#' @return A matrix S the same size as the matrices in W.
#'
#' @details Performs the calculation
#'
#' S = sum( t_i * W_i )
#'
#' @export multiply_and_reduce
#'
#' @examples
#'
#' t <- 1:3
#'
#' W <- map(1:3, ~diag(1, 3))
#'
#' multiply_and_reduce(t, W)
#'
multiply_and_reduce <- function(t, W){

  if(!is.list(W)){
    stop("W must be a list of numeric matrices")
  }

  if(any(map_lgl(W, ~!is.matrix(.x)))){
    stop("W must be a list of numeric matrices")
  }

  if(any(map_lgl(W, ~!is.numeric(.x)))){
    stop("W must be a list of numeric matrices")
  }

  if(length(t) == 1){
    t <- rep(t, length(W))
  }

  if(length(t) != length(W)){
    stop("t and W must be the same length")
  }


  # Use map to create new list
  step_1_multiply <- map2(t, W,
                          function(t, W){
                            t * W
                          })

  # Add up the results
  step_2_reduce <- Reduce("+", step_1_multiply)

  return(
    step_2_reduce
  )
}

