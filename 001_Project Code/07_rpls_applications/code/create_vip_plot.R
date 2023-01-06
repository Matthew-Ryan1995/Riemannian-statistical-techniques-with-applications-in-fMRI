find_na <- function(x){
  all(is.na(x))
}
create_vip_plot <- function(vip, sig, beta){
  vip[which(sig >= 0.05, arr.ind = T)] <- NA # Get rid of non-significant

  colnames(vip) <- colnames(sig)
  vip <- vip %>%
    as_tibble(rownames = "pred") %>%
    rowwise() %>%
    mutate(filt = find_na(c_across(-pred))) %>%
    filter(!filt) %>%
    select(-filt) %>%
    ungroup() %>%
    pivot_longer(-pred) %>%
    group_by(name) %>%
    drop_na() %>%
    slice_max(value, n = 10) %>%
    nest() %>%
    mutate(
      data = map(data, function(df) arrange(df, pred)),
      sign = map2(data,name,
                  function(df, n){
                    df <- df %>%
                      arrange(pred)
                    j <- which(colnames(beta) == n)
                    tmp <- beta[which(rownames(beta) %in% df$pred), j]
                    sign(tmp[sort(names(tmp))])
                  }),
      # strength = map2(data,name,
      #             function(df, n){
      #               j <- which(colnames(beta) == n)
      #               abs(beta[which(rownames(beta) %in% df$pred), j])
      #             }),
    ) %>%
    unnest(cols = c(data, sign)) %>%
    ungroup() %>%
    mutate(sign = ifelse(sign < 0, "negative", "positive"),
           pred = fct_reorder(pred, value),
           name = str_remove(name, "Patient|Autism|Male|Closed"),
           name = str_to_title(name))

  p <- vip %>%
    ggplot(aes(x = value, y = pred, fill = sign)) +
    geom_col() +
    # geom_label(aes(label = round(strength, 2))) +
    facet_wrap(~name, ncol = 2, scales = "free_y") +
    labs(x = NULL, y = NULL, fill = "Coefficient sign") +
    theme_classic() +
    harrypotter::scale_fill_hp_d("ravenclaw", direction = -1) +
    theme(legend.position = "bottom",
          text = element_text(size = 20),
          axis.text.y = element_text(size = 10),
          strip.background = element_rect(fill = "gray80"))

  return(p)
}
create_single_vip_plot <- function(vip, sig, beta){
  vip[which(sig >= 0.05, arr.ind = T)] <- NA # Get rid of non-significant

  # colnames(vip) <- "value"
  vip <- vip %>%
    as_tibble(rownames = "pred") %>%
    mutate(pred = names(sig)) %>%
    drop_na() %>%
    slice_max(value, n = 10)

  # return(vip)
  p <- vip %>%
    mutate(pred = fct_reorder(pred, value)) %>%
    ggplot(aes(x = value, y = pred)) +
    geom_col() +
    labs(x = "VIP", y = NULL) +
    theme_classic() +
    harrypotter::scale_fill_hp_d("ravenclaw", direction = -1) +
    theme(legend.position = "bottom",
          text = element_text(size = 20),
          axis.text.y = element_text(size = 10),
          strip.background = element_rect(fill = "gray80"))

  return(p)
}
