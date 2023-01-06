get_euc_dist <- function(cor_mats){
  require(tidyverse, quietly = TRUE)
  
  cor_mats <- map(cor_mats, ~.x[upper.tri(.x)]) %>% 
    simplify2array() %>% 
    t()
  
  D <- dist(cor_mats)
  return(D)
}

get_cor_dist <- function(cor_mats){
  require(tidyverse, quietly = TRUE)
  
  cor_mats <- map(cor_mats, ~.x[upper.tri(.x)]) %>% 
    simplify2array() 
  
  
  
  D <- sqrt(2-2*cor(cor_mats)) %>% 
    as.dist()
  
  return(D)
}
