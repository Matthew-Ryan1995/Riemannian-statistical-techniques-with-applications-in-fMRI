create_facet_plot_r2 <- function(N, interest = "P", dest_folder = NULL){

  max_possible_r2 <- expand_grid(tibble(K = 1:10,
                                        r2 = map_dbl(K, ~min(1, 1 - ((5 - .x)/5)))),
                                 p = unique(df$p),
                                 q = unique(df$q)
  )

  p <- df %>%
    filter(n == N) %>%
    group_by(sig, n, K, p, q) %>%
    summarise(across(c(contains("rmse"), contains("r2")), mean), .groups = "drop") %>%
    mutate(sig = as.factor(sig)) %>%
    rename(r2 = glue::glue("{interest}_r2")) %>%
    ggplot(aes(x = K, y = r2)) +
    geom_line(aes(colour = sig)) +
    geom_vline(xintercept = 5, colour = "red", lty = 2) +
    facet_grid(p ~ q, labeller = label_both) +
    harrypotter::scale_color_hp_d("ravenclaw") +
    scale_x_continuous(labels = 1:10, breaks = 1:10) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Number of PLS components",
         y = latex2exp::TeX("$R^2$"),
         colour = "Error",
         title = glue::glue("Recovery of {interest}, n = {N}")) +
    theme_classic() +
    theme(panel.spacing = unit(2, "lines"),
          # strip.text = element_text(size = 10),
          strip.background = element_rect(fill = "gray70"))

  if(interest %in% c("P", "Q")){
    p <- p +
      geom_point(data = max_possible_r2, shape = 4, colour = "black", size = 2)
  }

  if(is.null(dest_folder)){
    ggsave(glue::glue("figs/r2_{interest}_{N}.png"),
           plot = p, width = 10, height = 10)
  }else{
    ggsave(glue::glue("{dest_folder}/figs/r2_{interest}_{N}.png"),
           plot = p, width = 10, height = 10)
  }

  return(p)
}

create_facet_plot_rmse <- function(N, interest = "X", dest_folder = NULL){
  p <- df %>%
    filter(n == N) %>%
    group_by(sig, n, K, p, q) %>%
    summarise(across(c(contains("rmse"), contains("r2")), mean), .groups = "drop") %>%
    mutate(sig = as.factor(sig)) %>%
    rename(r2 = glue::glue("{interest}_rmse")) %>%
    ggplot(aes(x = K, y = r2)) +
    geom_line(aes(colour = sig)) +
    # geom_point(data = max_possible_r2, shape = 4, colour = "black", size = 2) +
    geom_vline(xintercept = 5, colour = "red", lty = 2) +
    facet_grid(p ~ q, labeller = label_both) +
    harrypotter::scale_color_hp_d("ravenclaw") +
    scale_x_continuous(labels = 1:10, breaks = 1:10) +
    labs(x = "Number of PLS components",
         y = "RMSE",
         colour = "Error",
         title = glue::glue("Recovery of {interest} (RMSE), n = {N}")) +
    theme_classic() +
    theme(panel.spacing = unit(2, "lines"),
          # strip.text = element_text(size = 10),
          strip.background = element_rect(fill = "gray70"))

  if(is.null(dest_folder)){
    ggsave(glue::glue("figs/rmse_{interest}_{N}.png"),
           plot = p, width = 10, height = 10)
  }else{
    ggsave(glue::glue("{dest_folder}/figs/rmse_{interest}_{N}.png"),
           plot = p, width = 10, height = 10)
  }

  return(p)
}
