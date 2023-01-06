#' affine_pga
#'
#' Performs principle geodesic analysis on the cone of positive definite matrices using the affine
#' invariant metric.
#'
#'
#' @param data List, array, or dataframe: The data you want to find the mean of.
#' @param var.name Character: The name of the variable containing the correlation matrices. Defaults to NULL.
#' @param tol Numeric: The tolerance you want to work at for the mean calculation. Defaults to 1e-8.
#' @param nrep Integer: The total number of iterations before you stop searching for the mean. Defaults to 1000.
#' @param num_comp Integer: The number of geodesic directions you want to calculate.
#' @param symmetric: Logical: Should the symmetric square root be used, or the Cholesky? Defaults to TRUE.
#' @param cores Integer > 0: This is set up for parallelisation on Mac through mclapply.
#' This will determine the number of cores to use. Defaults to 1.
#' @param scale: Logical: Should scaling be performed before PGA components are calculated. Defaults to FALSE.
#'
#' @return
#' - score: the lower dimensional representation of the data as scores.
#' - loadings: The tangent vectors at the Frechet mean as matrices.  Used for interpretation.
#' - sdev: The varaince along each principle geodesic.
#'
#' @details This performs the tangent space approximation algorithm of PGA found in
#' "Riemannian Geometric Statistics in Medical Image Analysis."  This function makes use of the
#' prcomp function built into R.
#'
#' @export affine_pga
#'
#' @examples
#' ## Calculate the first tow PGA directions for the COBRE dataset.
#' ## This will take a while.
#'
#' data(cobre)
#'
#' cobre_pga <- affine_pga(cobre, "cors", num_comp = 2)
#'
#' # Look at the scores against each other:
#' plot(cobre_pga$scores, xlab = "PGA 1", ylab = "PGA 2")
#'
#' # Look at the loadings
#' corrplot::corrplot(cobre_pga$loadings[[1]], is.corr = FALSE)
#' corrplot::corrplot(cobre_pga$loadings[[2]], is.corr = FALSE)
affine_pga <- function(data, var.name = NULL, tol = 1e-4, nrep = 1000, num_comp = 5,
                       symmetric = TRUE, cores = 1, scale = FALSE){
  ## Warning for scaling
  if(scale){
    warning("Scaling can lead to an inflation of impact due to the diagonal elements.")
  }
  ## Make the data into a list
  tmp <- check_right_input(data, var.name)

  data <- tmp$x
  rm(tmp)

  # First calculate the Frechet mean
  Fmean <- get_frechet_mean(data, tol = tol, nrep = nrep, cores = cores, symmetric = symmetric)

  # Calculate the square root to speed up calculations
  F.sq <- get_sqrt(Fmean, symmetric = symmetric)

  # Push data into the tangent space at the Frechet mean
  # "Centering the data"
  logged_data <- map(data, ~affine_map(P = Fmean, W = .x, FUN = "log", P.sq = F.sq))

  # Vectorise, or make Euclidean, the tangent vectors
  vec_data <- map(logged_data, ~vec(P = Fmean, W = .x, P.sq = F.sq)) %>%
    simplify2array() %>%
    t()

  num.comp <- min(ncol(vec_data), nrow(vec_data) - 1, num_comp)

  if(num.comp < num_comp){
    warning(glue::glue("The maximum number of componets we can calculate is {num.comp}.  We will aim for this instead."))
  }

  # Perform PCA in the tangent space
  # No centering, the data has already been centered by the log map
  pca <- prcomp(vec_data, center = FALSE, scale = scale, rank = num.comp)

  # Move loadings back to tangent space of mean for interpretability.
  loadings <- map(1:num.comp, ~unvec(Fmean, pca$rotation[, .x], P.sq = F.sq))


  return(list(
    scores = pca$x[, 1:num.comp],
    loadings = loadings,
    sdev = pca$sdev[1:num.comp]
  ))
}

