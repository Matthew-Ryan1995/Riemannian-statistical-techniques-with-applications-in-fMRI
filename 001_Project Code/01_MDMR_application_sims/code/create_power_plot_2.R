## Matthew Ryan
## 15/08/2022
# packages ----------------------------------------------------------------
pacman::p_load(tidyverse)

# data --------------------------------------------------------------------
power_files <- list.files("power_results_2", full.names = T, recursive = T)
power_sims <- map_dfr(power_files, readRDS)


# create plot -------------------------------------------------------------

p <- power_sims %>% 
  filter(predictor == "groupControl") %>% # select the group variable
  rename(p_val = `Permutation p-value`) %>% # Simplify p_value name
  group_by(block, rho, method, r) %>% 
  summarise(power = sum(p_val < 0.05)/max(sim), .groups = "drop") %>% # Calculate power for simulation parameters
  mutate(
    rho = round(rho, 2),
    rho = as.factor(rho)
  ) %>% 
  ggplot(aes(x = r, y = power, colour = method)) +
  geom_line() +
  geom_hline(yintercept = 0.05, lty = 2) +
  theme_classic() +
  theme(panel.spacing=unit(1,"lines"),
        strip.text = element_text(size = 14),
        text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "bottom") +
  labs(
    colour = "Distance metric",
    x = "r",
    y = "Power"
  ) +
  harrypotter::scale_color_hp_d("ravenclaw") +
  facet_grid(block ~ rho)

p

ggsave("../../002_Thesis/01_initial_draft/img/power_sim_plot.png", plot = p, width = 10, height = 10)
