pacman::p_load(tidyverse, phdWork, igraph, ggraph)

set.seed(12345)

# Create adjacency matrix
adj <- cobre %>%
  filter(id == "40012") %>%
  pull(cors) %>%
  .[[1]]
diag(adj) <- 0
adj[abs(adj) < 0.28] <- 0 # Threshold so full connected

# Create graph, add network to attributes
g <- graph_from_adjacency_matrix(adjmatrix = adj,
                                 mode = "undirected",
                                 weighted = TRUE,
                                 diag = FALSE)
labels <- read_csv(
  "~/nilearn_data/msdl_atlas/MSDL_rois/msdl_rois_labels.csv", col_types = cols()
) %>%
  pull(`net name`)

g <- g %>%
  set_vertex_attr("network", value = labels)


## make plot
p <- ggraph((g), layout = "dh") +
  geom_edge_link(alpha = 0.5) +
  geom_node_label(aes(fill = network, label = name), size = 4) +
  theme_void() +
  labs(fill = "Network") +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme(legend.position = "bottom")

ggsave("~/Desktop/PhD/PhD/002_Thesis/01_initial_draft/img/network_graph.png", plot = p,
       width = 12, height = 12)

