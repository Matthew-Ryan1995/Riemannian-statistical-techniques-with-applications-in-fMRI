## Matthew Ryan
## 15/08/2022
# packages ----------------------------------------------------------------
pacman::p_load(tidyverse)

# data --------------------------------------------------------------------
reg_files <- list.files("regularisation_results", full.names = T, recursive = T)
reg_sims <- map_dfr(reg_files, readRDS)


# create plot -------------------------------------------------------------

p <- reg_sims %>% 
  filter(predictor == "groupControl") %>% 
  rename(p_val = `Permutation p-value`) %>% 
  group_by(lam) %>% 
  summarise(power = sum(p_val < 0.05)/max(sim), .groups = "drop") %>% 
  ggplot(aes(x = lam, y = power)) +
  geom_line() +
  geom_hline(yintercept = 0.05, lty = 2) +
  theme_classic() +
  labs(
    x = latex2exp::TeX("$\\lambda$"),
    y = "Type I Error"
  ) 

p

ggsave("../../002_Thesis/01_initial_draft/img/regularisation_sim_plot.png", plot = p, width = 7, height = 5)

reg_sims %>% 
  filter(predictor == "groupControl") %>% 
  rename(p_val = `Permutation p-value`) %>% 
  group_by(lam) %>% 
  summarise(power = sum(p_val < 0.05)/max(sim), .groups = "drop") %>% 
  rstatix::get_summary_stats(power)
