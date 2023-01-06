## Hemodynamic Response Function from SPM
library(ggplot2)
a1 <- 6
a2 <- 12
  
b1 <- 1
b2 <- 1

A <- 1
c <- 1/6
  
t <- seq( 0, 20, length.out = 1000 )

G1 <- t ^ ( a1 - 1 ) * b1 ^ ( a1 ) * exp( -b1 * t ) / gamma( a1 ) 

G2 <- t ^ ( a2 - 1 ) * b2 ^ ( a2 ) * exp( -b2 * t ) / gamma( a2 ) 



hrf <- A * ( G1 - c * G2 )

hrf_plot <- ggplot( mapping = aes( x = t, y = hrf ) ) +
  geom_line() +
  geom_hline(yintercept = 0, colour = "red", lty = 2) +
  theme_minimal() +
  labs( x = "Time", y = "Hemodynamic Response Function" ) +
  theme( 
    axis.text = element_blank()
    )
  
hrf_plot

ggsave(
  "hrf_plot.png",
  plot = hrf_plot,
  device = "png"
)
