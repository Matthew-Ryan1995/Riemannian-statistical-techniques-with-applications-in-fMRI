
# libraries ---------------------------------------------------------------
pacman::p_load(tidyverse, phdWork, rsample, yardstick)

walk(c("code/cross_validate_pls.R"), source)

# job_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print(str_c("Doing job COBRE"))

set.seed(12345)

# Data --------------------------------------------------------------------
df <- cobre %>%
  select(cors, age, group) %>%
  mutate(
    group = fct_relevel(group, "Control")
  )

model_formula <- ~ age + group

lambda <- seq(0, 5, by = 0.5)

s <- sample(1:nrow(df), round(nrow(df))/10)

train_data <- df %>%
  slice(-s)
test_data <- df %>%
  slice(s)

res <- numeric()

for(i in 1:length(lambda)){
  X <- train_data %>%
    mutate(cors = map(cors, ~.x + diag(lambda[i], nrow(.x)))) %>%
    pull(cors)
  Y <- model.matrix(model_formula, data = train_data) %>%
    as.matrix() %>%
    .[, -1]

  # Fit the model on largest number of components considered
  M <- riemannian_pls(X = X, Y = Y, L = 2, tol = 1e-4, max.iter = 20,
                      type_x = "affine", type_y = "euclidean", method = "tangent",
                      scale = TRUE, mc.cores = 6)

  X_test <- test_data %>%
    mutate(cors = map(cors, ~.x + diag(lambda[i], nrow(.x)))) %>%
    pull(cors)

  Y_test <- model.matrix(model_formula, data = test_data) %>%
    as.matrix() %>%
    .[, -1]

  Y_test <- t((t(Y_test) - M$muY)/M$Y_sd)
  preds <- predict.riemannian_pls(object = M, newdata = X_test, tol = 1e-4, max.iter = 20,
                                  num_comp = 2, method = "tangent")

  res[i] <- rsq.riemannian_pls(truth = Y_test, est = preds, muY = M$muY, m = M$muY,
                               type_y = "euclidean")
}
