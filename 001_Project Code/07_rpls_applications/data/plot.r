pacman::p_load(tidyverse, targets)
theme_set(theme_minimal())
df <- crossing(
  x = 1:20, y = 1:20,
  x_grp = LETTERS[1:5],
  y_grp = LETTERS[1:5]
)
n <- nrow(df)
df$r <- runif(n, min = -1, max = 1)
df$r[sample(1:n, 0.8*n)] <- NA
df |>
  ggplot(aes(x, y, fill = r)) +
  geom_tile(show.legend = FALSE) +
  facet_grid(y_grp ~ x_grp, space = "free") +
  scale_fill_gradient2(na.value = "white") +
  labs( x = NULL, y = NULL) +
  theme(
    axis.text = element_blank(),
    panel.spacing = unit(0, "mm"),
    panel.border = element_rect(colour = "black", fill = NA)
  )
