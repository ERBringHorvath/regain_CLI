args <- commandArgs(trailingOnly = TRUE)

boot_RDS <- args[1]
data_matrix <- args[2]
metadata_file <- args[3]
stats_file <- args[4]

pkgs <- c("visNetwork", "igraph", "bnlearn", "RColorBrewer")
for (pkg in pkgs) {
  suppressMessages({
    if (!require(pkg, character.only=TRUE)) {
      install.packages(pkg)
      suppressMessages(library(pkg, character.only=TRUE))
    }
  })
}

boot = readRDS(boot_RDS)
averaged.network(boot)
avg_boot = averaged.network(boot, threshold = 0.5)
arcs(avg_boot) <- directed.arcs(avg_boot)

data <- read.csv(data_matrix, row.names = 1)
metadata <- read.csv(metadata_file)
gene_names <- metadata[,1]
lookup <- setNames(as.character(metadata[, 2]), metadata[, 1])
valid_genes <- intersect(gene_names, colnames(data))

stats <- read.csv(stats_file)

net <- igraph.from.graphNEL(as.graphNEL(avg_boot))

if(vcount(net) < 2 || ecount(net) == 0) {
  stop("Network is too small or has no edges for visualization.")
}

visZach <- toVisNetworkData(net)

nodes <- visZach$nodes
edges <- visZach$edges

# Apply function to each edge

get_width <- function(node_from, node_to, stats) {
  row <- stats %>% filter((stats$Gene_A == node_from & stats$Gene_B == node_to) |
                            (stats$Gene_A == node_to & stats$Gene_B == node_from))
  if(nrow(row) > 0) {
    width_value <- stats$Conditional_Probability_CI_high[1]
    # Check for negative or unexpected values
    if(width_value <= 0) {
      stop("Encountered negative or zero width value in get_width function.")
    }
    return(width_value)
  } else {
    return(0)  # Default value if no match is found
  }
}

# Apply the function to each pair of nodes in 'edges'
edges$width <- sapply(1:nrow(edges), function(i) get_width(edges$from[i], edges$to[i], stats))

# Rescale the weights
edges$width <- scales::rescale(edges$width, to = c(1, 5))

nodes <- nodes[nodes$id %in% valid_genes, ]

nodes$group <- lookup[nodes$id]

palette1 <- RColorBrewer::brewer.pal(8, "Set3")
palette2 <- RColorBrewer::brewer.pal(8, "Set2")
palette3 <- RColorBrewer::brewer.pal(12, "Paired")
palette4 <- RColorBrewer::brewer.pal(8, "Dark2")

color_palette <- c(palette1, palette2, palette3, palette4)

unique_groups <- unique(nodes$group)

color_palette <- color_palette[1:length(unique_groups)]

if (length(unique_groups) > length(color_palette)) {
  color_palette <- rep(color_palette, length.out = length(unique_groups))
}

group_color_lookup <- setNames(color_palette, unique_groups)
nodes$color <- group_color_lookup[nodes$group]

# Prepare the legend (lnodes) to match the unique groups and their respective colors
lnodes <- data.frame(color = color_palette,
                     label = unique_groups,
                     font.color = 'white')
# Create and save the network
network <- visNetwork(nodes = nodes, edges = edges, width = '100%', height = 900) %>%
  visNodes(size = 20, color = list(highlight = 'yellow'), font = list(size = 25)) %>%
  visEdges(smooth = list(enabled = T, type = 'diagonalCross', roundness = 0.1),
           physics = F, color = "black") %>%
  visIgraphLayout(layout = "layout_with_fr", type = 'full') %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T, labelOnly = T),
             nodesIdSelection = T) %>%
  visLegend(addNodes = lnodes, width = 0.1, position = 'left', main = 'Gene Class', ncol = 2, useGroups = F)
visSave(network, "Bayesian_Network.html", background = "#F5F4F4")

cat("\n \033[32mAnalysis complete.\033[39m\n \n")