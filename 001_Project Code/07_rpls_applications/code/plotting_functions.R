# Functions ---------------------------------------------------------------

find_network <- function(s){
  regions %>%
    filter(region == s) %>%
    pull(rsn)
}

create_subnetwork_fc <- function(P, sig = NULL){
  P_mod <- P %>%
    as_tibble(rownames = "source") %>%
    pivot_longer(-source, names_to = "target") %>%
    rowwise() %>%
    mutate(
      net1 = map_chr(source,
                     find_network),
      net2 = map_chr(target,
                     find_network)
    ) %>%
    ungroup() %>%
    mutate(net1 = fct_inorder(net1),
           net2 = factor(net2, levels = levels(net1)),
           pred = str_c(target, source, sep = "-"))

  if(!is.null(sig)){
    P_mod <- left_join(P_mod, sig, by = "pred") %>%
      mutate(value = ifelse(V1 < 0.05, value, NA))
  }

  P_mod <- P_mod %>%
    select(-pred) %>%
    group_by(net1, net2) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(value = ifelse(is.nan(value), 0, value)) %>%
    pivot_wider(names_from = net2, values_from = value) %>%
    select(-net1) %>%
    as.matrix()
  rownames(P_mod) <- colnames(P_mod)
  return(P_mod)
}

create_vip_subnetwork <- function(vip, sig){

  P_mod <- vip %>%
    as_tibble(rownames = "source") %>%
    pivot_longer(-source, names_to = "target") %>%
    rowwise() %>%
    mutate(
      net1 = map_chr(source,
                     find_network),
      net2 = map_chr(target,
                     find_network)
    ) %>%
    ungroup() %>%
    mutate(net1 = fct_inorder(net1),
           net2 = factor(net2, levels = levels(net1)),
           pred = str_c(target, source, sep = "-"))

  if(!is.null(sig)){
    P_mod <- left_join(P_mod, sig, by = "pred") %>%
      mutate(value = ifelse(V1 < 0.05, value, NA))
  }

  P_mod <- P_mod %>%
    select(-pred) %>%
    group_by(net1, net2) %>%
    summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(value = ifelse(is.nan(value), 0, value)) %>%
    pivot_wider(names_from = net2, values_from = value) %>%
    select(-net1) %>%
    as.matrix()
  rownames(P_mod) <- colnames(P_mod)
  return(P_mod)
}

create_corrplot <- function(P, title, filter = FALSE, sig = NULL,
                            plot_type = "", height = 5, vip = NULL){
  data <- P %>%
    as_tibble(rownames = "source") %>%
    pivot_longer(-source, names_to = "target") %>%
    drop_na() %>%
    mutate(source = fct_inorder(source),
           target = factor(target, levels(source)),
           pred = str_c(target, source, sep = "-"))

  if(!is.null(sig)){
    if(plot_type != "rsn"){
      data <- left_join(data, sig, by = "pred") %>%
        mutate(value = ifelse(V1 < 0.05, value, NA))
    }
  }
  if(filter){
    if(plot_type == "rsn" & !is.null(vip)){
      vip_dat <- vip %>%
        as_tibble(rownames = "source") %>%
        pivot_longer(-source, names_to = "target", values_to = "vip") %>%
        drop_na() %>%
        mutate(source = fct_inorder(source),
               target = factor(target, levels(source)),
               pred = str_c(target, source, sep = "-"),
               filter_val = ifelse(vip == 0, NA, vip))
      data_sig <- left_join(data, vip_dat, by = c("source", "target", "pred")) %>%
        mutate(value = ifelse(vip > quantile(filter_val, 0.75, na.rm = TRUE), value, NA)) %>%
        drop_na()

    }else{
      data <- data %>%
        mutate(value = ifelse(abs(value) < quantile(abs(value), 0.9, na.rm = TRUE),
                              NA,
                              value))
    }
  }
  max_lim <- max(abs(data$value))

  if(isFALSE(plot_type=="rsn")){
    p <- data %>%
      mutate(value = ifelse(value == 0, NA, value)) %>%
      ggplot(aes(x = target, y = fct_rev(source))) +
      geom_tile(fill = "white", colour = "black") +
      geom_point(aes(colour = value), size = 2) +
      scale_colour_distiller(palette = "RdBu", direction = 1,
                             na.value = "white",
                             limits = c(-max_lim, max_lim),
                             guide = guide_colorbar(frame.colour = "black",
                                                    frame.linewidth = 2,
                                                    ticks.colour = "black",
                                                    ticks.linewidth = 2,
                                                    barheight = unit(0.5, "cm"),
                                                    barwidth = unit(2*height, "cm"))) +
      labs(x = NULL, y = NULL, colour = NULL,  title = title)
  }else{
    p <- data %>%
      ggplot(aes(x = target, y = fct_rev(source))) +
      geom_tile(aes(fill = value), colour = "black") +
      scale_fill_distiller(palette = "RdBu", direction = 1,
                           na.value = "white",
                           limits = c(-max_lim, max_lim),
                           guide = guide_colorbar(frame.colour = "black",
                                                  frame.linewidth = 2,
                                                  ticks.colour = "black",
                                                  ticks.linewidth = 2,
                                                  barheight = unit(0.5, "cm"),
                                                  barwidth = unit(2*height, "cm"))) +
      labs(x = NULL, y = NULL, fill = "", title = title)
    if(exists("data_sig")){
      p <- p +
        geom_tile(data = data_sig, aes(fill = value), colour = "black", size = 2)
    }
  }


  p <- p +
    theme_void() +
    theme(title = element_text(size = 30),
          axis.text.x = element_text(angle = 90),
          axis.text.y = element_text(angle = 0),
          plot.title = element_text(hjust = 0.5, size = 30),
          legend.position = "bottom")

  return(p)
}

create_loading_plots <- function(L, dataset = "cobre",
                                 height = 5, filter = FALSE,
                                 save_folder = NULL){
  if(is_null(save_folder)){
    save_folder <- "img"
  }
  P_true <- unvec(P = M$muX, vec = M$loadingsX[[L]])
  Q <- M$loadingsY[[L]]
  colnames(Q) <- "value"
  Q <- Q %>%
    as_tibble(rownames = "predictor") %>%
    mutate(predictor = str_remove(predictor, "Patient|Autism|Male|Closed"),
           title = str_c(predictor, " = ", round(value, 3)))
  title <- str_c("Y Loadings:\n", str_c(Q$title[1:2], collapse = "; "))
  if(dataset == "abide"){
    title <- str_c(title, "\n", str_c(Q$title[3:4], collapse = "; "))
  }

  # p <- corrplot(P, is.corr = FALSE, type = "lower",
  #               col = COL2("RdBu", 200), title = title,
  #               mar=c(0,0,2,0))
  P <- P_true
  P[upper.tri(P)] <- NA
  p <- create_corrplot(P, title, filter = filter, height = height)
  ggsave(glue::glue("{save_folder}/{dataset}_loading_{filter}_{L}.png"),
         plot = p,
         height = height,
         width =  height)
  ## Subnetwork plots
  P <- create_subnetwork_fc(P_true)
  P[upper.tri(P)] <- NA
  p <- create_corrplot(P, title, filter = filter, plot_type = "rsn", height = height)
  # p <- corrplot(P, is.corr = FALSE, type = "lower",
  #               col = COL2("RdBu", 200), title = title,
  #               mar=c(0,0,2,0))
  ggsave(glue::glue("{save_folder}/{dataset}_loading_subnetworks_{filter}_{L}.png"),
         plot = p,
         height = height,
         width = height)
}

create_regression_coeff_plots <- function(B, L, dataset = "cobre",
                                          height = 5, filter = FALSE, sig = NULL,
                                          vip_vals = NULL,
                                          save_folder = NULL){
  if(is_null(save_folder)){
    save_folder <- "img"
  }
  P_true <- unvec(P = M$muX, vec = B[, L])

  if(!is.null(sig)){
    var_name <- colnames(B)[L]
    sig <- sig %>%
      as_tibble(rownames = "pred") %>%
      select(pred, V1 = all_of(var_name))
  }
  if(!is.null(vip_vals)){
    vip_mat <- vip_vals[, which(colnames(vip_vals) == var_name)]

    vip_mat <- vip_mat %>%
      as_tibble(rownames = "pred") %>%
      mutate(R = str_split(pred, pattern = "-"),
             source = map_chr(R, ~.x[1]),
             target = map_chr(R, ~.x[2])
      ) %>%
      select(-R, -pred) %>%
      pivot_wider(names_from = target, values_from = value, values_fill = 0) %>%
      select(-source) %>%
      as.matrix()
    vip_mat <- vip_mat + t(vip_mat)
    colnames(vip_mat) <- rownames(vip_mat) <- colnames(P_true)
  }
  Q <- colnames(B)[L]
  Q <- str_remove(Q, "Patient|Autism|Male|Closed")
  title <- str_to_title(Q)

  # p <- corrplot(P, is.corr = FALSE, type = "lower",
  #               col = COL2("RdBu", 200), title = title,
  #               mar=c(0,0,2,0))
  P <- P_true
  P[upper.tri(P)] <- NA
  p <- create_corrplot(P, title, filter = filter, height = height, sig = sig)
  ggsave(glue::glue("{save_folder}/{dataset}_coefficient_{filter}_{L}.png"),
         plot = p,
         height = height,
         width =  height)
  ## Subnetwork plots
  vip_sub <- create_vip_subnetwork(vip = vip_mat, sig = sig)
  P <- create_subnetwork_fc(P_true, sig = sig)
  P[upper.tri(P)] <- NA
  vip_sub[upper.tri(vip_sub)] <- NA
  p <- create_corrplot(P, title, filter = filter, plot_type = "rsn", height = height, vip = vip_sub)
  # p <- corrplot(P, is.corr = FALSE, type = "lower",
  #               col = COL2("RdBu", 200), title = title,
  #               mar=c(0,0,2,0))
  ggsave(glue::glue("{save_folder}/{dataset}_coefficient_subnetworks_{filter}_{L}.png"),
         plot = p,
         height = height,
         width = height)
}
