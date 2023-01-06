pacman::p_load(tidyverse, phdWork, rstatix, ggpubr)

# Get data
df <- cobre %>%
  select(-cors, -msdl_roi, -id)

height <- 7

# Perfrom wilcox test on age
cobre_group_wilcox_test <- df %>%
  wilcox_test(age ~ group)

# Create boxplot with p-values
cobre_age_group_plot <- df %>%
  ggplot() +
  geom_boxplot(aes(x = group, y = age, fill = group),
               show.legend = FALSE) +
  stat_pvalue_manual(data = cobre_group_wilcox_test,
                     label = "W = {statistic}; p = {p}",
                     y.position = 68,
                     size = 6) +
  labs(
    x =  NULL,
    y = "Age"
  ) +
  ylim(17, 72) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )

ggsave("~/Desktop/PhD/PhD/002_Thesis/01_initial_draft/img/cobre_age_group_boxplot.png",
       plot = cobre_age_group_plot,
       height = .618 * height,
       width = height)
