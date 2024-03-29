% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/riemannian_pls_affine.R
\name{scores_and_weights_affine}
\alias{scores_and_weights_affine}
\title{Calculate the scores and weights for geodesic PLS, affine variant}
\usage{
scores_and_weights_affine(
  X,
  Y,
  muX,
  muY,
  max.iter = 50,
  tol = 1e-06,
  tau = 1/2,
  ...
)
}
\arguments{
\item{X}{List: Data predictor matrices X.  If dataframe, please use optional input var.name.X}

\item{Y}{List: Data predictor matrices Y.  If dataframe, please use optional input var.name.Y}

\item{muX}{Matrix: Mean of the X.}

\item{muY}{Matrix: Mean of Y.}

\item{max.iter}{Numeric: Total number of iterations to do.}

\item{tol}{Numeric: Tolerance to work at.}

\item{...}{Optional arguments
- var.name.X, var.name.Y: Character, the column containing the data matrices
- muX.sq, muY.sq: Matrix, the square root of mu to speed computation
- tau: Numeric, Initial step size
- mc.cores: Numeric, if running parallel, how many cores.
- check: Logical, should we check our parameters
- t, u, W, C: Initial estimates of scores and weights.
- verbose: Logical, should run time information be printed.
- save: Logical, to be implemented.}
}
\value{
A list with W = weights of X, C = Weights of Y, t = scores of X, u = scores of Y.
}
\description{
The title says it all.
}
\examples{
TBA
}
