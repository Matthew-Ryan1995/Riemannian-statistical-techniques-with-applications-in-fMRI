#' Generate data for a given set of parameters
#'
#' @param sims - Number of simulations to run per parameter set
#' @param n - vector of sample sizes to consider
#' @param p - Dimension of predictor space
#' @param q - Dimension of response space
#' @param sig - Vector of noise values to consider
#' @param L - True latent space size, defaults to 5
#' @param seed - Seed for reproducibi;ity
#'
#' @return
#' @export
#'
#' @examples
#' NA
gen_data <- function(sims = 1, n, p, q, sig, L = 5, seed = 1234){

  # Need sig ordered for my code to work
  sig <- sort(sig)

  if(!is.null(seed)){
    set.seed(seed)
  }

  # Fixed P and Q across simulation parameters p, q and L
  P <- matrix(rnorm(L*p), nrow = p, ncol = L)
  Q <- matrix(rnorm(L*q), nrow = q, ncol = L)


  # Make P and Q orthogonal
  if(p < L){
    P <- t(svd(P)$v)
  }else{
    P <- svd(P)$u
  }
  if(q < L){
    Q <- t(svd(Q)$v)
  }else{
    Q <- svd(Q)$u
  }

  # Set up simulation parameters, this repeats the parameters for each simulation per sample size
  params <- map_dfr(n,
                    function(n){
                      expand_grid(sig = sig, n = n, sim = rep(1:sims, n), K = 1:(L)) #Changed to L
                    }) %>%
    arrange(sig)

  N <- sum(sims*n*(L)*(length(sig))) # Total number of obs to construct #Changed to L

  ## Create the scores
  Tmat <- matrix(rnorm(N * L), nrow = N)
  ## Taking B = I so T = U

  ## Create error matrices
  E <-  do.call(rbind, map(sig,
                           ~matrix(rnorm(N * p/(length(sig)), sd = .x),
                                   nrow = N/(length(sig)))))
  Fmat <-  do.call(rbind, map(sig,
                              ~matrix(rnorm(N * q/(length(sig)), sd = .x),
                                      nrow = N/(length(sig)))))

  X <- tcrossprod(Tmat, P) + E
  Y <- tcrossprod(Tmat, Q) + Fmat
  TP <- tcrossprod(Tmat, P)
  UQ <- tcrossprod(Tmat, Q)

  colnames(X) <- str_c("X", 1:ncol(X))
  colnames(Y) <- str_c("Y", 1:ncol(Y))
  colnames(TP) <- str_c("TP", 1:ncol(TP))
  colnames(UQ) <- str_c("UQ", 1:ncol(UQ))

  dat <- cbind(X, Y, TP, UQ) %>%
    as_tibble() %>%
    bind_cols(params) %>% # Join on sim params
    group_by(sig, sim, n, K) %>%
    nest(X = starts_with("X"), Y = starts_with("Y"),
         TP = starts_with("TP"), UQ = starts_with("UQ")) %>% # nest my data
    mutate(X = map(X, as.matrix), # Convert to matrices
           X = map(X, ~apply(.x, 2, scale, scale = FALSE)), # I am applying centering here, this may be an issue.
           Y = map(Y, as.matrix),
           Y = map(Y, ~apply(.x, 2, scale, scale = FALSE)),
           TP = map(TP, as.matrix),
           TP = map(TP, ~apply(.x, 2, scale, scale = FALSE)),
           UQ = map(UQ, as.matrix),
           UQ = map(UQ, ~apply(.x, 2, scale, scale = FALSE)),
           ) %>%
    ungroup() %>%
    mutate(p = p, q = q, L = L, P = list(P), Q = list(Q)) # add extra params, plus P and Q

  return(dat)
}
