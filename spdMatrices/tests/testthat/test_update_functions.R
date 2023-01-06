library(phdWork)


# h1 ----------------------------------------------------------------------

# test_that(
#   "Test that h1 gives right errors",
#   {
#     mu <- matrix(1, 2, 2)
#     W <- diag(1, 2)
#     mu[1, 2] <- 2
#     expect_error(h1(mu, W), "mu must be a symmetric matrix")
#
#     mu <- matrix(1, 2, 2)
#     W <- diag(1, 2)
#     W[1, 2] <- 2
#     expect_error(h1(mu, W), "W must be a symmetric matrix")
#
#     mu <- matrix(1, 2, 2)
#     W <- diag(1, 2)
#     expect_error(h1(mu, W), "mu must be positive definite, at least one eigenvalue is negative or zero.")
#   }
# )
#
# test_that(
#   "Test that h1 gives correct answer",
#   {
#     data(cobre)
#     mu <- cobre$cors[[1]]
#     W <- cobre$cors[[2]]
#     expect_equal(round(h1(mu, W), 4), 1502.2298)
#   }
# )
