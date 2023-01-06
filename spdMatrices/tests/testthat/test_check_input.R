library(phdWork)

df <- cobre %>%
  slice(1:5)
df_array <- map(df$cors, ~.x) %>%
  simplify2array()


# check_right_input -------------------------------------------------------


test_that(
  "Test that output is list.",
  {
    # Full output
    expect_type(check_right_input(df, "cors"), "list")
    expect_type(check_right_input(pull(df, "cors")), "list")
    expect_type(check_right_input(df_array), "list")
    # Data output
    expect_type(check_right_input(df, "cors")$x, "list")
    expect_type(check_right_input(pull(df, "cors"))$x, "list")
    expect_type(check_right_input(df_array)$x, "list")
  }
)

test_that(
  "Test that output names are correct.",
  {
    # DF entry
    df_convert_names <- names(check_right_input(df, "cors"))
    expect_match(df_convert_names[1], "x")
    expect_match(df_convert_names[2], "n")
    expect_match(df_convert_names[3], "roi.num")
    # DF entry
    list_convert_names <- names(check_right_input(pull(df, "cors")))
    expect_match(list_convert_names[1], "x")
    expect_match(list_convert_names[2], "n")
    expect_match(list_convert_names[3], "roi.num")
    # DF entry
    array_convert_names <- names(check_right_input(df_array))
    expect_match(array_convert_names[1], "x")
    expect_match(array_convert_names[2], "n")
    expect_match(array_convert_names[3], "roi.num")
    rm("df_convert_names", "list_convert_names", "array_convert_names")
  }
)

test_that(
  "Test that output right length.",
  {
    # Full length
    expect_length(check_right_input(df, "cors"), 3)
    expect_length(check_right_input(pull(df, "cors")), 3)
    expect_length(check_right_input(df_array), 3)
    # Data length
    expect_length(check_right_input(df, "cors")$x, 5)
    expect_length(check_right_input(pull(df, "cors"))$x, 5)
    expect_length(check_right_input(df_array)$x, 5)
  }
)

test_that(
  "Test for correct number of subjects.",
  {
    expect_equal(check_right_input(df, "cors")$n, 5)
    expect_equal(check_right_input(pull(df, "cors"))$n, 5)
    expect_equal(check_right_input(df_array)$n, 5)
  }
)

test_that(
  "Test for correct number of ROI.",
  {
    expect_equal(check_right_input(df, "cors")$roi.num, 39)
    expect_equal(check_right_input(pull(df, "cors"))$roi.num, 39)
    expect_equal(check_right_input(df_array)$roi.num, 39)
  }
)

test_that(
  "Test that you get errors",
  {
    expect_error(check_right_input(df$age),
                 "x must be an array, a list, or a dataframe containing the correlation values")
    expect_error(check_right_input(df),
                 "When supplying a dataframe, please provide the variable name containing correlations")
    # expect_error(check_right_input(df, cors),
    # "When supplying a dataframe, please provide the variable name containing correlations")
  }
)

rm("df", "df_array")


# check_input_values ------------------------------------------------------


test_that(
  "Test that you get correct errors for wrong inputs.",
  {
    t <- 1:3
    u <- 1:3
    W <- C <- matrix(1, 2, 2)
    X <- Y <- map(1:3, ~diag(1, 2))
    muX <- muY <- diag(1, 2)
    lambda_t <- lambda_u <- 1:3
    lambda_1 <- lambda_2 <- 1
    nu <- tau <- 1

    # test t
    t <- "a"
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "t should be a numeric vector."
    )
    t <- 1:3

    #test u
    u <- "a"
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "u should be a numeric vector."
    )
    u <- 1:3

    #test lambda_t
    lambda_t <- "a"
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "lambda_t should be a numeric vector."
    )
    lambda_t <- 1:3

    #test lambda_u
    lambda_u <- "a"
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "lambda_u should be a numeric vector."
    )
    lambda_u <- 1:3

    #test lambda_1
    lambda_1 <- "a"
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "lambda_1 should be a numeric value."
    )
    lambda_1 <- 1

    #test lambda_2
    lambda_2 <- "a"
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "lambda_2 should be a numeric value."
    )
    lambda_2 <- 1

    #test nu
    nu <- "a"
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "nu should be a positive number."
    )
    nu <- -1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "nu should be a positive number."
    )
    nu <- 1

    #test tau
    tau <- "a"
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "tau should be a non-negative number."
    )
    tau <- -1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "tau should be a non-negative number."
    )
    tau <- 1

    #test W
    W <- 1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "W should be a numeric matrix."
    )
    W <- matrix("a", 2, 2)
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "W should be a numeric matrix."
    )
    W <- matrix(1, 2, 2)

    #test C
    C <- 1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "C should be a numeric matrix."
    )
    C <- matrix("a", 2, 2)
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "C should be a numeric matrix."
    )
    C <- matrix(1, 2, 2)

    #test muX
    muX <- 1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "muX should be a numeric matrix."
    )
    muX <- matrix("a", 2, 2)
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "muX should be a numeric matrix."
    )
    muX <- diag(1, 2)

    #test muY
    muY <- 1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "muY should be a numeric matrix."
    )
    muY <- matrix("a", 2, 2)
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "muY should be a numeric matrix."
    )
    muY <- diag(1, 2)

    #test X
    X <- list(1, matrix(2, 2, 2))
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "X must be a list of numeric matrices"
    )
    X <- list(matrix("2", 2, 2), matrix(2, 2, 2))
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "X must be a list of numeric matrices"
    )
    X <-  map(1:3, ~diag(1, 2))

    #test Y
    Y <- list(1, matrix(2, 2, 2))
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "Y must be a list of numeric matrices"
    )
    Y <- list(matrix("2", 2, 2), matrix(2, 2, 2))
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "Y must be a list of numeric matrices"
    )
    Y <-  map(1:3, ~diag(1, 2))
  }
)

test_that(
  "Test that inputs match each other correctly",
  {
    t <- 1:3
    u <- 1:3
    W <- C <- matrix(1, 2, 2)
    X <- Y <- map(1:3, ~diag(1, 2))
    muX <- muY <- diag(1, 2)
    lambda_t <- lambda_u <- 1:3
    lambda_1 <- lambda_2 <- 1
    nu <- tau <- 1

    # test t and u same length
    t <- 1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "t and u should be the same length."
    )
    t <- 1:3

    u <- 1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "t and u should be the same length."
    )
    u <- 1:3

    # Test t and lambda_t same length
    lambda_t <- 1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "The number of constraints on t does not much the number of t-values."
    )
    lambda_t <- 1:3

    # Test u and lambda_u same length
    lambda_u <- 1
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "The number of constraints on u does not much the number of u-values."
    )
    lambda_u <- 1:3

    # Test lambda_1 and lambda_2 are single numbers
    lambda_1 <- 1:2
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "lambda_1 should be a single value."
    )
    lambda_1 <- 1

    lambda_2 <- 1:2
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "lambda_2 should be a single value."
    )
    lambda_2 <- 1

    # Test W and muX same dimension
    W <- matrix(1, 1, 1)
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "W and muX should be the same dimension."
    )
    W <- matrix(1, 2, 2)

    muX <- matrix(1, 1, 1)
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "W and muX should be the same dimension."
    )
    muX <- diag(1, 2)

    # Test C and muY same dimension
    C <- matrix(1, 1, 1)
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "C and muY should be the same dimension."
    )
    C <- matrix(1, 2, 2)

    muY <- matrix(1, 1, 1)
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "C and muY should be the same dimension."
    )
    muY <- diag(1, 2)

    # Check W and X same dimension
    X <- map(1:3, ~diag(1, 1))
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "W and X should be the same dimension."
    )

    X <- list(diag(1, 2), diag(1, 1))
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "W and X should be the same dimension."
    )
    X <- map(1:3, ~diag(1, 2))

    # Check C and Y same dimension
    Y <- map(1:3, ~diag(1, 1))
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "C and Y should be the same dimension."
    )

    Y <- list(diag(1, 2), diag(1, 1))
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "C and Y should be the same dimension."
    )
    Y <- map(1:3, ~diag(1, 2))
  }
)

test_that(
  "Test that check_input_values picks up missing values",
  {

    t <- 1:3
    u <- 1:3
    W <- C <- matrix(1, 2, 2)
    X <- Y <- map(1:3, ~diag(1, 2))
    muX <- muY <- diag(1, 2)
    lambda_t <- lambda_u <- 1:3
    lambda_1 <- lambda_2 <- 1
    nu <- tau <- 1

    # test t for missing
    t <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in t."
    )
    t <- 1:3

    # test u for missing
    u <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in u."
    )
    u <- 1:3

    # test lambda_t for missing
    lambda_t <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in lambda_t."
    )
    lambda_t <- 1:3

    # test lambda_u for missing
    lambda_u <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in lambda_u."
    )
    lambda_u <- 1:3

    # test W for missing
    W[1, 1] <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in W."
    )
    W <- matrix(1, 2, 2)

    # test C for missing
    C[1, 1] <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in C."
    )
    C <- matrix(1, 2, 2)

    # test muX for missing
    muX[1, 1] <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in muX."
    )
    muX <- diag(1, 2)

    # test muY for missing
    muY[1, 1] <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in muY."
    )
    muY <- diag(1, 2)

    # test X for missing
    X[[1]][1, 1] <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in X."
    )
    X <- map(1:3, ~diag(1, 2))

    # test Y for missing
    Y[[1]][1, 1] <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "There are missing values in Y."
    )
    Y <- map(1:3, ~diag(1, 2))

    # test lambda_1 for missing
    lambda_1 <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "lambda_1 is missing."
    )
    lambda_1 <- 1

    # test lambda_2 for missing
    lambda_2 <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "lambda_2 is missing."
    )
    lambda_2 <- 1

    # test nu for missing
    nu <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "nu is missing."
    )
    nu <- 1

    # test nu for missing
    tau <- NA
    expect_error(
      check_input_values(t = t, u = u, W = W, C = C,
                         X = X, Y = Y, muX = muX, muY = muY,
                         lambda_t = lambda_t, lambda_u = lambda_u, lambda_1 = lambda_1,
                         lambda_2 = lambda_2, nu = nu, tau = tau),
      "tau is missing."
    )
    tau <- 1

  }
)

