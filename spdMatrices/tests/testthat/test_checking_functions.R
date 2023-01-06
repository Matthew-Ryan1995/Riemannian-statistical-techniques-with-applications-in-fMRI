library(phdWork)

A <- matrix(c(1, 0.5,
              0.5, 1), 2, 2)
B <- matrix(c(1, 1.5,
              1.5, 1), 2, 2)


# Positive definite -------------------------------------------------------


test_that(
  "Test that is_positive_definite() works.",
  {
    expect_true(is_positive_definite(A))
    expect_true(!is_positive_definite(B))
  }
)


# Square root -------------------------------------------------------------


test_that(
  "Test that square root works",
  {
    tmp <- get_sqrt(A)
    expect_equal(tmp %*% tmp, A)

    tmp <- get_sqrt(A, symmetric = FALSE)
    expect_equal(t(tmp) %*% tmp, A)
  }
)


# Trace -------------------------------------------------------------------



test_that(
  "Test that get_trace returns the trace",
  {
    expect_equal(get_trace(A), 2)
    expect_equal(get_trace(B), 2)

    p <- 10
    X <- matrix(rnorm(100 * p), ncol = p)
    H <- X %*% solve(t(X) %*% X) %*% t(X)
    expect_equal(get_trace(H), p)
  }
)

test_that(
  "Test that trace does not work for non-square matrices.",
  {
    p <- 10
    X <- matrix(rnorm(100 * p), ncol = p)
    expect_error(get_trace(X), "A must be a square matrix.")
  }
)


# Double centre -----------------------------------------------------------

# Check all error messages
# Check double centres

test_that(
  "That that errors in double_centre works",
  {
    x <- 1:3
    expect_error(double_centre(x), "X must be a square numeric matrix")

    X <- matrix("A", 3, 3)
    expect_error(double_centre(X), "X must be a square numeric matrix")

    X <- matrix(1, 2, 3)
    expect_error(double_centre(X), "X must be a square numeric matrix")

    X <- diag(1, 4)
    X[1, 2] <- NA
    expect_error(double_centre(X), "I don't know how to deal with missing values, please forgive me.")
  }
)

test_that(
  "Test that double_centre actually works",
  {
    set.seed(2022)
    X <- matrix(rnorm(9), 3, 3)

    centred_X <- double_centre(X)
    expect_equal(apply(centred_X, 1, mean), rep(0, 3))
    expect_equal(apply(centred_X, 2, mean), rep(0, 3))
    expect_equal(mean(centred_X), 0)
  }
)

# Multiply and reduce -----------------------------------------------------

# Check error messages
# Check it works


test_that(
  "Test that errors in multiply_and_reduce work",
  {
    t <- 1:4
    W <- map(1:3, ~diag(1, 2))
    expect_error(multiply_and_reduce(t, W), "t and W must be the same length")

    t <- 1:3
    W <- list(1, 2, 3)
    expect_error(multiply_and_reduce(t, W), "W must be a list of numeric matrices")


    t <- 1:3
    W <- map(1:3, ~matrix("a", 2, 2))
    expect_error(multiply_and_reduce(t, W), "W must be a list of numeric matrices")

    t <- 1:3
    W <- matrix(1, 2, 2)
    expect_error(multiply_and_reduce(t, W), "W must be a list of numeric matrices")
  }
)

test_that(
  "Test that multiply_and_reduce works as expected",
  {
    t <- 1:3
    W <- map(1:3, ~diag(1,2))
    expect_equal(multiply_and_reduce(t, W), diag(6, 2))
  }
)
