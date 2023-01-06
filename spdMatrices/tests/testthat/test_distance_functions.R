library(testthat)
library(phdWork)

A <- diag(1, 2)
B <- matrix(
  c(
    1, 0.5,
    0.5, 1
  ),
  2, 2
)
C <- matrix(
  c(
    1, -0.5,
    -0.5, 1
  ),
  2, 2
)

matrix_list <- list(A, B, C)

## normal distance
D <- diag(0, 3, 3)

D[1, 2] <- D[2, 1] <- D[1, 3] <- D[3, 1] <- sqrt(log(1.5)^2 + log(0.5)^2)
D[2, 3] <- D[3, 2] <- sqrt(log(3)^2 + log(1/3)^2)

dimnames(D) <-  map(1:2, ~c("1", "2", "n"))

## Regularised distance
D_reg <- diag(0, 3, 3)

D_reg[1, 2] <- D_reg[2, 1] <-  sqrt(log(5/4)^2 + log(3/4)^2) -> D_reg[1, 3] -> D_reg[3, 1]
D_reg[2, 3] <- D_reg[3, 2] <- sqrt(log(5/3)^2 + log(3/5)^2)

dimnames(D_reg) <-  map(1:2, ~c("1", "2", "n"))

## Regularised distance l = 1/2
D_reg2 <- diag(0, 3, 3)

D_reg2[1, 2] <- D_reg2[2, 1] <-  sqrt(log(4/3)^2 + log(2/3)^2) -> D_reg2[1, 3] -> D_reg2[3, 1]
D_reg2[2, 3] <- D_reg2[3, 2] <- sqrt(log(2)^2 + log(1/2)^2)

dimnames(D_reg2) <-  map(1:2, ~c("1", "2", "n"))



# affine distance ---------------------------------------------------------


test_that("Test that distance is as expected.",
          {
            expect_equal(as.matrix(affine_dist(matrix_list)), D)
          }
)

test_that("Test that regulariesd distance is as expected.",
          {
            expect_equal(as.matrix(affine_dist(matrix_list, regularise = TRUE)), D_reg)
            expect_equal(as.matrix(affine_dist(matrix_list, regularise = TRUE, lambda = 1/2)), D_reg2)
          }
)

test_that("Test that we get warning from regularisation.",
          {
            expect_warning(
              abide_aal %>%
                slice(1:5) %>%
                affine_dist("cors"),
              "There are subjects with low rank functional connectivity profiles.  Regularised distances will be calculated for distances involving subjects 1, 5."
            )
          }
)

test_that("Test that function stops on incorrect entries",
          {
            # Regularise
            expect_error(affine_dist(matrix_list, regularise = "FALSE"),
                         "The regularise variable should either be a TRUE or FALSE Boolean.")
            expect_error(affine_dist(matrix_list, regularise = NA),
                         "The regularise variable should either be a TRUE or FALSE Boolean.")

            # Cores
            expect_error(affine_dist(matrix_list, cores = 9),
                         "You are trying to run on an incorrect number of cores. This must be an integer between 1 and the maximum number of available cores on your system.")
            expect_error(affine_dist(matrix_list, cores = -1),
                         "You are trying to run on an incorrect number of cores. This must be an integer between 1 and the maximum number of available cores on your system.")

            # Lambda
            expect_error(affine_dist(matrix_list, lambda = 0),
                         "You cannot regularise by a negative value of lambda.")
            expect_error(affine_dist(matrix_list, lambda = -5),
                         "You cannot regularise by a negative value of lambda.")
          }
)

# Euclidean and Corr ------------------------------------------------------

test_that("Test that Euclidean distance works",
          {
            D <- matrix(0, 3, 3)
            D[1, 2] <- D[1, 3] <- D[2, 1] <- D[3, 1] <-  1/2
            D[2, 3] <- D[3, 2] <- 1
            dimnames(D) <- map(1:2, ~c("1", "2", "3"))
            expect_equal(as.matrix(euclid_dist(matrix_list)), D)

            A <- diag(1, 3)
            A[1, 3] <- A[3, 1] <- 1
            B <- matrix(
              c(1, 0.5, -1,
                0.5, 1, 0.5,
                -1, 0.5, 1),
              3, 3, byrow = TRUE
            )
            C <- matrix(
              c(1, -0.5, 1,
                -0.5, 1, -0.5,
                1, -0.5, 1),
              3, 3, byrow = TRUE
            )
            matrix_list <- list(A, B, C)

            D <- diag(0, 3)
            D[1, 2] <- D[2, 1] <- sqrt(9/2)
            D[1, 3] <- D[3, 1] <-  sqrt(1/2)
            D[2, 3] <- D[3, 2] <- sqrt(6)
            dimnames(D) <- map(1:2, ~c("1", "2", "3"))

            expect_equal(as.matrix(euclid_dist(matrix_list)), D)
          }
)

test_that("Test that error in correlation distances works.",
          {
            expect_error(corr_dist(matrix_list),
                         "You cannot calculate the correlation distance for 2 x 2 matrices.")
          }
)

test_that("Test that correlation distance works.",
          {
            A <- diag(1, 3)
            A[1, 3] <- A[3, 1] <- 1
            B <- matrix(
              c(1, 0.5, -1,
                0.5, 1, 0.5,
                -1, 0.5, 1),
              3, 3, byrow = TRUE
            )
            C <- matrix(
              c(1, -0.5, 1,
                -0.5, 1, -0.5,
                1, -0.5, 1),
              3, 3, byrow = TRUE
            )
            matrix_list <- list(A, B, C)

            D <- diag(0, 3)
            D[1, 2] <- D[2, 1] <- 2
            D[1, 3] <- D[3, 1] <-  0
            D[2, 3] <- D[3, 2] <- 2
            dimnames(D) <- map(1:2, ~c("1", "2", "3"))

            expect_equal(as.matrix(corr_dist(matrix_list)), D)
          }
)

rm(list = ls())

