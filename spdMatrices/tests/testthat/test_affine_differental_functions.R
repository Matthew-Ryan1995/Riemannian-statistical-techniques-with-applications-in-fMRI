library(phdWork)


set.seed(1234)

A <- matrix(
  rnorm(16),
  4, 4
)
A <- (A + t(A))/10
A <- expm::expm(A)

I <- diag(1, 4)

P <- matrix(
  rnorm(16),
  4, 4
)
P <- (P + t(P))/10
P <- expm::expm(P)

V <- matrix(
  rnorm(16),
  4, 4
)
V <- V + t(V)

W <- matrix(
  rnorm(16),
  4, 4
)
W <- W + t(W)

test_that(
  "Test that matrix exponential matches",
  {
    tol <- 1e-6
    h <- 1e-8
    expect_lt(
      sqrt(sum((expm::expm(A + h*V) - expm::expm(A) - h * exponential_derivative(A, V))^2))/h,
      tol
    )

    t <- 0.5
    expect_lt(
      sqrt(sum((expm::expm(t * (A + h*V)) - expm::expm(t * A) - h * exponential_derivative(A, V, t = t))^2))/h,
      tol
    )

    t <- -2
    expect_lt(
      sqrt(sum((expm::expm(t * (A + h*V)) - expm::expm(t * A) - h * exponential_derivative(A, V, t = t))^2))/h,
      tol
    )
  }
)

test_that(
  "Test that my exponential derivative errors work as expected",
  {
    # A non-diagonalisable matrix
    expect_error(exponential_derivative(A = matrix(c(1, 0, 1, 1), 2, 2, byrow = TRUE), V, symm = FALSE),
                 "A must be a diagonalisable matrix.")

    tmp_A <- A
    tmp_A[4,3] <- 0
    tmp_V <- V
    tmp_V[4,3] <- 0
    expect_error(exponential_derivative(A = tmp_A, V = V, symm = TRUE),
                 "A must be symmetric")
    # expect_error(exponential_derivative(A, tmp_V),
    #              "V must be symmetric")
    expect_error(exponential_derivative(A, V, t = "1"),
                 "t must be a real number.")
  }
)

test_that(
  "Test that the general matrix exponential derivative works as expected.",
  {
    h <- 1e-8
    expect_lt(
      (sqrt(sum((affine_map(P, A + h*V) - affine_map(P, A) - h * exponential_derivative_general(P, A, V))^2)))/h,
      1e-6
    )

    t <- 0.5
    expect_lt(
      sqrt(sum((affine_map(P, t*(A + h*V)) - affine_map(P, t*A) - h *exponential_derivative_general(P, A, V, t = t))^2))/h,
      1e-6
    )

    t <- -2
    expect_lt(
      sqrt((sum((affine_map(P, t*(A + h*V)) - affine_map(P, t*A) - h *exponential_derivative_general(P, A, V, t = t))^2)))/h,
      1e-6
    )
  }
)

test_that(
  "Test that adjoint works",
  {
    t <- 1
    expect_equal(
      affine_inner_product(
        expm::expm(t * A),
        exponential_derivative(A, V, t = t),
        W
      ),
      affine_inner_product(
        diag(1, 4),
        V,
        exponential_adjoint(A, W, t = t)
      )
    )
    t <- 0.5
    expect_equal(
      affine_inner_product(
        expm::expm(t * A),
        exponential_derivative(A, V, t = t),
        W
      ),
      affine_inner_product(
        diag(1, 4),
        V,
        exponential_adjoint(A, W, t = t)
      )
    )
    t <- -2
    expect_equal(
      affine_inner_product(
        expm::expm(t * A),
        exponential_derivative(A, V, t = t),
        W
      ),
      affine_inner_product(
        diag(1, 4),
        V,
        exponential_adjoint(A, W, t = t)
      )
    )
  }
)

test_that(
  "Test that general adjoint works",
  {
    t <- 1
    expect_equal(
      affine_inner_product(
        affine_map(P, t*A),
        exponential_derivative_general(P, A, V, t = t),
        W
      ),
      affine_inner_product(
        P,
        V,
        exponential_adjoint_general(P, A, W, t = t)
      )
    )
    t <- -0.5
    expect_equal(
      affine_inner_product(
        affine_map(P, t*A),
        exponential_derivative_general(P, A, V, t = t),
        W
      ),
      affine_inner_product(
        P,
        V,
        exponential_adjoint_general(P, A, W, t = t)
      )
    )
    t <- 2
    expect_equal(
      affine_inner_product(
        affine_map(P, t*A),
        exponential_derivative_general(P, A, V, t = t),
        W
      ),
      affine_inner_product(
        P,
        V,
        exponential_adjoint_general(P, A, W, t = t)
      )
    )
  }
)


test_that(
  "Test that the matrix logarithm differential works.",
  {
    expect_equal(
      exponential_derivative(
        expm::logm(A),
        log_derivative(A, V)
      ),
      V
    )
    expect_equal(
      log_derivative(
        expm::expm(A),
        exponential_derivative(A, V)
      ),
      V
    )
  }
)

test_that(
  "Test that the general matrix logarithm differential works.",
  {
    expect_equal(
      exponential_derivative_general(
        P,
        affine_map(P, A, FUN = "log"),
        log_derivative_general(P, A, V)
      ),
      V
    )
    expect_equal(
      log_derivative_general(
        P,
        affine_map(P, A, FUN = "exp"),
        exponential_derivative_general(P, A, V)
      ),
      V
    )
  }
)

test_that(
  "Test that matrix logarithm adjoint works.",
  {
    expect_equal(
    affine_inner_product(
      I,
      log_derivative(A, V),
      W
    ),
    affine_inner_product(
      A,
      V,
      log_adjoint(
        A,
        W
      )
    )
    )

    expect_equal(
    affine_inner_product(
      P,
      log_derivative_general(P, A, V),
      W
    ),
    affine_inner_product(
      A,
      V,
      log_adjoint_general(
        P,
        A,
        W
      )
    )
    )
  }
)
