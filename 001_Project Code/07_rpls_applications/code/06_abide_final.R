
# libraries ---------------------------------------------------------------
pacman::p_load(tidyverse, phdWork)

# Data --------------------------------------------------------------------
df <- abide_aal %>%
  select(cors, age, group, sex, eye) %>%
  mutate(cors = map(cors, ~(.x + diag(1, nrow(.x)))),
         group = factor(group, levels = c("Control", "Autism")),
         sex = as.factor(sex),
         eye = factor(eye, levels = c("Open", "Closed")))

model_formula <- ~ age + group + sex + eye



K <- 3 # Best from CV
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


write_rds(results, "results/abide_rpls_model_K_6.Rds")

