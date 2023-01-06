pacman::p_load(tidyverse, phdWork, rstatix, ggpubr)

# Get data
df <- abide_aal %>%
  select(-sub_id, -roi, - cors)

height <- 7

# Perfrom wilcox test on age
abide_age_group_wilcox <- df %>%
  wilcox_test(age ~ group)

# Create boxplot with p-values
abide_age_group_plot <- df %>%
  ggplot() +
  geom_boxplot(aes(x = group, y = age, fill = group),
               show.legend = FALSE) +
  stat_pvalue_manual(data = abide_age_group_wilcox,
                     label = "W = {statistic}; p = {p}",
                     y.position = 35,
                     size = 6) +
  labs(
    x =  NULL,
    y = "Age"
  ) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

ggsave("~/Desktop/PhD/PhD/002_Thesis/01_initial_draft/img/abide_age_group_boxplot.png",
       plot = abide_age_group_plot,
       height = .618 * height,
       width = height)


height <- 9
# Perform wilcox test on IQ
abide_iq_group_wilcox <- df %>%
  select(group, contains("iq")) %>%
  pivot_longer(-group, values_to = "iq") %>%
  mutate(name = str_to_upper(name)) %>%
  group_by(name) %>%
  wilcox_test(iq ~ group)

# Create a boxplot for each type of IQ, add p-values
abide_iq_group_plot <- df %>%
  select(group, contains("iq")) %>%
  pivot_longer(-group, values_to = "iq") %>%
  mutate(name = str_to_upper(name)) %>%
  ggplot() +
  geom_boxplot(aes(x = group, y = iq, fill = group),
               show.legend = FALSE) +
  facet_wrap(~name) +
  stat_pvalue_manual(data = abide_iq_group_wilcox,
                     label = "W = {statistic};\n p = {p}",
                     y.position = 155, size = 6) +
  ylim(60, 170) +
  labs(
    x = NULL,
    y = "IQ"
  ) +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "gray80"),
    axis.text = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text =  element_text(size = 20)
  )

ggsave("~/Desktop/PhD/PhD/002_Thesis/01_initial_draft/img/abide_iq_group_boxplot.png",
       plot = abide_iq_group_plot,
       height =  .618 * height,
       width =  height)





