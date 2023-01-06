### Matt Ryan
### 20/12/2022
### 
### 


# Packages ----------------------------------------------------------------
pacman::p_load(tidyverse, patchwork, phdWork)


# Functions ---------------------------------------------------------------

create_corrplot <- function(P, title, height = 5, vip = NULL){
  data <- P %>%
    as_tibble(rownames = "source") %>%
    pivot_longer(-source, names_to = "target") %>%
    drop_na() %>%
    mutate(source = fct_inorder(source),
           target = factor(target, levels(source)),
           pred = str_c(target, source, sep = "-"))
  max_lim <- max(abs(data$value))
  
  p <- data %>%
    ggplot(aes(x = target, y = fct_rev(source))) +
    geom_tile(aes(fill = value), colour = "black") +
    geom_text(aes(label = round(value, 2)), colour = "black", size = 10) +
    scale_fill_distiller(palette = "RdBu", direction = 1,
                         na.value = "white",
                         limits = c(-1.001, 1.001),
                         guide = guide_colorbar(frame.colour = "black",
                                                frame.linewidth = 2,
                                                ticks.colour = "black",
                                                ticks.linewidth = 2,
                                                barheight = unit(1, "cm"),
                                                barwidth = unit(5*height, "cm"))) +
    labs(x = NULL, y = NULL, fill = "", title = title)
  
  
  p <- p +
    theme_void() +
    theme(title = element_text(size = 30),
          legend.text = element_text(size = 25),
          axis.text = element_blank(),
          # axis.text.y = element_text(angle = 0),
          plot.title = element_text(hjust = 0.5, size = 30),
          legend.position = "bottom", )
  
  return(p)
}


# Data --------------------------------------------------------------------

height <- 8
df <- cobre
s <- sample(1:nrow(df), 1)
C <- df$cors[[s]]

B <- generate_signal_matrix(rho = atanh(-0.5), block_size = 4, method = "ar1", vary = FALSE)
F0 <- implant_signal(C, B, r = 0)
F25 <- implant_signal(C, B, r = 0.25)
F5 <- implant_signal(C, B, r = 0.5)
F75 <- implant_signal(C, B, r = 0.75)
F1 <- implant_signal(C, B, r = 1)

C0 <- F0[4:7, 4:7] 
rownames(C0) <- colnames(C0) <- str_c("C", 1:4)
C1 <- F25[4:7, 4:7] 
rownames(C1) <- colnames(C1) <- str_c("C", 1:4)
C2 <- F75[4:7, 4:7] 
rownames(C2) <- colnames(C2) <- str_c("C", 1:4)
C3 <- F1[4:7, 4:7] 
rownames(C3) <- colnames(C3) <- str_c("C", 1:4)


# Make plot ---------------------------------------------------------------

p1 <- create_corrplot(C0, title = "")
p2 <- create_corrplot(C1, title = "")
p3 <- create_corrplot(C2, title = "")
p4 <- create_corrplot(C3, title = "")

p <- ((p1 | p2)/(p3 | p4) + 
  plot_layout(guides='collect') &
  theme(legend.position='bottom')) +
  plot_annotation(tag_levels = 'A', tag_suffix = ")", 
                  theme = theme(text = element_text(size = 6)))

ggsave("../../002_Thesis/01_initial_draft/img/implant_signal_example.png", 
       plot = p, height = 1.618 *height, width = 1.618 * height)
