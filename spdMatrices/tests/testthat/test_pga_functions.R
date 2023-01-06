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
M1 <- 1/3 * diag(log(3/4), 2)
Pbar <- diag(3/4, 2)^(1/3)

logged_data <- map(matrix_list, ~affine_map(P = Pbar, W = .x, FUN = "log"))

# Vectorise, or make Euclidean, the tangent vectors
vec_data <- map(logged_data, ~vec(P = Pbar, W = .x)) %>%
  simplify2array() %>%
  t()


# Check Riemannian maps ---------------------------------------------------

test_that(
  "Test that Riemannian maps work.",
  {
    expect_equal(
      affine_map(B, C, FUN = "log"),
      matrix(
        c(
          -log(3)/2, -log(3),
          -log(3), -log(3)/2
        ),
        2, 2
      )
    )
    expect_equal(
      affine_map(B, C, FUN = "exp"),
      1/4 * matrix(
        c(
          exp(3) + 3*exp(1/3), 3 * exp(1/3) - exp(3),
          3 * exp(1/3) - exp(3), exp(3) + 3*exp(1/3)
        ),
        2, 2
      )
    )
    expect_equal(affine_map(A, B, FUN = "log"), expm::logm(B))
    expect_equal(affine_map(A, B, FUN = "exp"), expm::expm(B))
  }
)

test_that(
  "Test that affine maps errors work, and inputs work.",
  {
    expect_equal(affine_map(B, C, FUN = "log"), affine_map(B, C, FUN = "log"))
    expect_equal(affine_map(B, C, FUN = "exp"), affine_map(B, C, FUN = "exp"))
    expect_error(affine_map(B, C, FUN = "Gary"),
                 "Please enter a valid method: either exp or log.")
  }
)



# Check Frechet mean algorithm and function -------------------------------

test_that(
  "Test that M1 does as expected",
  {
    expect_equal(
      get_affine_m1(A, matrix_list),
      M1
    )
  }
)

test_that(
  "Test that error calculation is correct",
  {
    expect_equal(affine_norm(A, M1), sqrt(2 * log(3/4)^2/9))
    }
)

test_that(
  "Test that output is list in get_frechet_mean.",
  {
    expect_type(get_frechet_mean(matrix_list, convergence_properties = TRUE), "list")
  }
)

test_that(
  "Test that output is correct list length in get_frechet_mean.",
  {
    expect_length(get_frechet_mean(matrix_list, convergence_properties = TRUE), 3)
  }
)

test_that(
  "Test that Frechet Mean is correct.",
  {
    expect_equal(get_frechet_mean(matrix_list), Pbar)
  }
)



# Affine PGA --------------------------------------------------------------


test_that(
  "Test that vectorisation is as expected",
  {
    expect_equal(
      vec_data,
      matrix(
        c(
          -1/3 * log(3) + 2/3 * log(2), -1/3 * log(3) + 2/3 * log(2), 0,
          1/6 * log(3) - 1/3 * log(2), 1/6 * log(3) - 1/3 * log(2), sqrt(2)/2 * log(3),
          1/6 * log(3) - 1/3 * log(2), 1/6 * log(3) - 1/3 * log(2), -sqrt(2)/2 * log(3)
        ),
        nrow = 3, ncol = 3, byrow = TRUE
      )
    )
  }
)

test_that(
  "Test that PGA components match eigenvectors and variances match eigenvalues",
  {
    ## Faff and define things
    eig <- eigen(cov(vec_data))
    tmp_eig <- map(1:ncol(eig$vectors), ~unvec(Pbar, eig$vectors[, .x]))
    tmp <- affine_pga(matrix_list, num_comp = 2)
    expect_equal_loadings <- function(){
      test_equal <- logical()
      for(i in 1:length(tmp$loadings)){
        if(sum(zapsmall((tmp_eig[[i]] - tmp$loadings[[i]]))) | sum(zapsmall((tmp_eig[[i]] + tmp$loadings[[i]]))) ){
          test_equal[i] <- TRUE
        }else{
          test_equal[i] <- FALSE
        }
      }
      return(
        all(test_equal)
      )
    }

    # Actuall test
    expect_true(expect_equal_loadings())

    expect_equal(tmp$sdev, sqrt(eig$values[1:2]))
  }
)

test_that(
  "Test that we get sufficient errors in affine_pga.",
  {
    expect_warning(affine_pga(matrix_list),
                   "The maximum number of componets we can calculate is 2.  We will aim for this instead.")
    expect_warning(affine_pga(matrix_list, num_comp = 2, scale = TRUE),
                   "Scaling can lead to an inflation of impact due to the diagonal elements.")
  }
)

test_that(
  "Test that output is list in affine_pga.",
  {
    expect_type(affine_pga(matrix_list, num_comp = 2), "list")
  }
)

test_that(
  "Test that output is correct list length in affine_pga.",
  {
    expect_length(affine_pga(matrix_list, num_comp = 2), 3)
  }
)
