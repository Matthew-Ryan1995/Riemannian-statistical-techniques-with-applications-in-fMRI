
# libraries ---------------------------------------------------------------
pacman::p_load(tidyverse, phdWork)

# Data --------------------------------------------------------------------
df <- cobre %>%
  select(cors, age, group) %>%
  mutate(
    group = fct_relevel(group, "Control")
  )

model_formula <- ~ age + group

K <- 2 # Best from CV
cores <- 6

# Full data -------------------------------------------------------------

X <- df$cors
Y <- model.matrix(model_formula, data = df) %>%
  as.matrix() %>%
  .[, -1]


# Fit model ------------------------------------------------------------------

results <- riemannian_pls(X = X, Y = Y, L = K, tol = 1e-4, max.iter = 20,
                          type_x = "affine", type_y = "euclidean", method = "tangent",
                          scale = TRUE, mc.cores = cores)


write_rds(results, "results/cobre_rpls_model_K_3.Rds")

