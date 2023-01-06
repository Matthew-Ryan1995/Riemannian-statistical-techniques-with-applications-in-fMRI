#' Implant signal matrix
#'
#' @param A Matrix. The matrix we want to implant signal in.  The dimension of A must be larger than the dimension of B.
#' @param B Matrix. The signal we wish to implant
#' @param network Character.  Do we place the signal in a particular network.  Default is DMN, which
#' is the only implemented network.
#'
#' @return Matrix of the same dimension as A with B implanted into the top left block.
#' @export
#' 
#' @details Given
#' 
#' $$
#' A = \begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}
#' $$
#'
#' Let 
#' 
#' $$
#' C = \begin{bmatrix} B_c' ((A_{11})_c^{-1})' & 0 \\ 0 & I \end{bmatrix}
#' $$
#' 
#' where the subscript c represents the Cholesky square root.  The resulting matrix given by 
#' 
#' $$
#' \tilde{A} = C A C'
#' $$
#' 
#' where the prime indicates matrix transpose, is such that $$\tilde{A}_{11} = B\, .$$
#'
#' @examples
#' A <- diag(1, 10)
#' B <- matrix(5, nrow = 2, ncol = 2)
#' 
#' A.tilde <- implant_signal(A, B, network = NULL)
implant_signal <- function(A, B, network = "DMN", r = NULL){
  block_size <- dim(B)[1]
  
  if(network == "DMN"){
    if(block_size>4){
      stop("Your block is too big")
    }
    A_chol <- chol(A[4:(block_size+3), 4:(block_size+3)])
  }else{
    A_chol <- chol(A[1:block_size, 1:block_size])
  }
  
  if(is.numeric(r)){
    if(r < 0 | r > 1){
      stop("r must be between 0 and 1")
    }
    B_chol <- chol((1-r) * A[4:(block_size+3), 4:(block_size+3)] + r * B)
  }else{
    B_chol <- chol(B)
  }
  
  if(network == "DMN"){
    C <- Matrix::bdiag(diag(1,3),
                       crossprod(B_chol, t(solve(A_chol))),
                       diag(1, dim(A)[1] - block_size - 3)) %>% 
      as.matrix()
  }else{
    C <- Matrix::bdiag(crossprod(B_chol, t(solve(A_chol))),
                       diag(1, dim(A)[1] - block_size)) %>% 
      as.matrix()
  }
  
  return(tcrossprod(crossprod(t(C), A), C))
}


#' Generate a signal matrix
#'
#' @param rho Numeric. The average strength of the signal to implant.  Should be between -1 and 1.
#' @param block_size Numeric. The size of the matrix to implant.  An integer greater than 0.
#' @param method Character. What method of correlation structure should be used.  The current options are
#' "constant", "toeplitz", or "hub".
#' @param rho_max Numeric. Used only for method = "hub".  The largest signal strength we want to measure.  
#' Should be between -1 and 1. 
#' @param rho_min Numeric. Used only for method = "hub".  The smallest signal strength we want to measure. 
#' Should be between -1 and rho_max.
#' @param vary Logical. Do we add variation to the average signal strength.  
#' This is used for more subject-to-subject variability.
#' @param se Numeric.  The standard error to be used if vary = TRUE.
#' @param seed Numeric.  If a value is given, this will be used as the seed for the signal generation.
#'
#' @return A matrix of size block_size x block_size used to implant signal into an FC matrix
#' 
#' @details TO BE FILLED OUT 
#' @export
#'
#' @examples
generate_signal_matrix <- function(rho, block_size, 
                                   method = "constant",
                                   rho_max = NULL, rho_min = NULL, vary = FALSE,
                                   se = 0.05, seed = NULL){
  if(!is_null(seed)){
    set.seed(seed)
  }
  if(vary){
    signal <- tanh(rnorm(1, rho, se))
  }else{
    signal <- tanh(rho)
  }
  
  if(method == "constant"){
    A <- matrix(data = signal, 
                nrow = block_size,
                ncol = block_size)
  }else if(method == "toeplitz"){
    A <- diag(1, block_size)
    rho <- map_dbl(1:block_size, ~(signal^{.x}))
    
    for(i in 1:(block_size-1)){
      A[i, (i+1):block_size] <- rho[-((block_size-(i-1)):block_size)]
    }
    
    A <- A + t(A)
  }else if(method == "hub"){
    if(is_null(rho_max)){
      rho_max <- signal
    }
    if(is_null(rho_min)){
      rho_min <- 0
    }
    
    A <- diag(1, block_size)
    alpha <- map_dbl(2:block_size, 
                     function(b){
                       return(
                         rho_max - ((b - 2)/(block_size - 2)) * (rho_max - rho_min)
                       )
                     })
    
    for(i in 1:(block_size-1)){
      A[i, (i+1):block_size] <- alpha[1:(block_size - i)]
    }
    
    A <- A + t(A)
  }else{
    stop("Enter a valid method")
  }
  
  diag(A) <- 1
  
  eigs <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  if(sum(eigs <=0) != 0){
    stop("Matrix is not positive definite")
  }
  return(A)
}



#' Simulate a subject
#'
#' @param rho Numeric. The average signal strength to be implanted. Should be between -1 and 1.
#' @param block_size Numeric. The size of the matrix to implant.  An integer greater than 0.
#' @param group Character. Is the subject a "patient" or "control."  
#' Patients will have an implanted signal, controls will not.
#' @param method Character. What method of correlation structure should be used.  The current options are
#' "constant", "toeplitz", or "hub".
#' @param test Logical. This is for testing purposes to check simulations against truth.
#' @param ... Other inputs for generate_signal_matrix.
#'
#' @return A single simulated functional connectivity matrix, with signal implanted if they are a patient.
#' @export
#' 
#' @details As a base, the subject will be randomly assigned a functional connectivity matrix from
#' a set of true data (our method used the COBRE data and the MSDL atlas).  A signal matrix will then
#' be created using the generate_signal_matrix function using the variables specified in the function
#' input.  This signal will be implanted using the implant_signal function.  Finally, the resulting
#' matrix with the implanted signal will be used as the scale matrix in a Wishart distribution
#' to act as the simulated data.
#'
#' @examples
simulate_subject <- function(rho, block_size, group,
                             method = "constant", test = FALSE, ...){
  
  if(!exists("df", inherits = FALSE)){
    df <- readRDS("data/df_msdl")
    df <- df %>% 
      mutate(cors = map(msdl_roi, cor))
  }
  
  ARGS <- list(...)
  if(!("vary" %in% names(ARGS))){
    ARGS$vary <- FALSE
  }
  if(!("se" %in% names(ARGS))){
    ARGS$se <- 0.05
  }
  if(!("r" %in% names(ARGS))){
    ARGS$r <- NULL
  }
  
  s <- sample(1:dim(df)[1], 1)
  A <- df$cors[[s]]
  
  if(group == "Patient"){
    if(method == "hub"){
      if(block_size < 3){
        stop("This wil not work, block size must be >2 for the hub method.")
      }
      B <- generate_signal_matrix(rho = rho, block_size = block_size, method = method,
                                  rho_max = ARGS$rho_max, rho_min = ARGS$rho_min,
                                  vary = ARGS$vary, se = ARGS$se, seed = ARGS$seed)
    }else{
      B <- generate_signal_matrix(rho = rho, block_size = block_size, method = method,
                                  vary = ARGS$vary, se = ARGS$se, seed = ARGS$seed)
    }
    
    A <- implant_signal(A, B, r = ARGS$r)
  }
  
  sim <- rWishart(1, 
                  df = sample(dim(A)[1]:150, 1),
                  Sigma = A)[,,1]
  
  sim <- cov2cor(sim)
  
  if(test){
    return(
      list(simulation = sim,
           base = s,
           args = ARGS)
    )
  }else{
    return(sim)
  }
}
