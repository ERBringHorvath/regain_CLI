args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
metadata_file <- args[2]
output_boot <- args[3]
threads <- if (length(args) >= 4 && !is.na(as.integer(args[4]))) as.integer(args[4]) else 2
number_of_bootstraps <- as.integer(args[5])
resamples <- as.integer(args[6])

print("User Inputs:")
print(paste("Matrix:", args[1]))
print(paste("Metadata:", args[2]))
print(paste("Bootstrapped Network Name:", args[3]))
print(paste("Number of Threads:", threads))
print(paste("Number of Bootstraps:", args[5]))
print(paste("Number of Resamples:", args[6]))

# Install required packages if not already installed
pkgs <- c('dplyr', 'parallel', 'pbapply', 'BiocManager', 'RColorBrewer', 'visNetwork', 'igraph', 'reshape2')
for (pkg in pkgs) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

pkgs_b <- c('bnlearn', 'gRain')
for (pkgs in pkgs_b) {
  if(!require(pkgs, character.only = TRUE)) {
    BiocManager::install(pkgs)
  }
}

data <- read.csv(input_file, row.names = 1)
d_fact <- data %>% mutate_if(is.numeric, as.factor)

n_cores <- threads
cl = parallel::makeCluster(n_cores)
clusterSetRNGStream(cl, 13245)

required_bootstraps = number_of_bootstraps

cat("\n \033[32mBootstrapping started.\033[39m\n \n")

boot = boot.strength(data = d_fact, R = number_of_bootstraps, algorithm = "hc",
                     algorithm.args = list(score="bde", iss=10), cluster = cl)

output_boot <- ifelse(grepl("\\.rds$", args[3]), args[3], paste0(args[3], ".rds"))
saveRDS(boot, file = output_boot)

averaged.network(boot)
avg_boot = averaged.network(boot, threshold = 0.5)
arcs(avg_boot) <- directed.arcs(avg_boot)

##Prepare to query the network
metadata <- read.csv(metadata_file)
gene_names <- metadata[,1]
lookup <- setNames(as.character(metadata[, 2]), metadata[, 1])
valid_genes <- intersect(gene_names, colnames(data))

Nlists <- resamples

boosts = function(d_fact, Nlists, avg_boot) {
  sample_data_list <- lapply(1:Nlists, function(i) {
    sample_data <- dplyr::slice_sample(d_fact, n=nrow(d_fact), replace = T)
    bnlearn::bn.fit(avg_boot, sample_data, method = 'bayes')
  })
  return(sample_data_list)
}

boosts_list <- boosts(d_fact, Nlists, avg_boot)

N <- length(valid_genes)
epsilon <- ((N + 0.5) / (N + 1))

combinations <- expand.grid(Gene_1 = valid_genes, Gene_2 = valid_genes)
combinations <- subset(combinations, Gene_1 != Gene_2)
cat(paste("\n \033[32mNumber of genes in dataset:\033[39m", N, "\n"))

max_combinations <- nrow(combinations)

probs_list <- vector("list", max_combinations * length(boosts_list))
risk_list <- vector("list", max_combinations * length(boosts_list))

compute_gene_stats <- function(gene1, gene2, bn, epsilon, query_number = 1) {
  if (query_number %% 1000 == 0) {
    message(sprintf("Calculated %d", query_number))
  }
  ##Compute conditional probability
  P_ij <- cpquery(bn, 
                  event = eval(parse(text = paste0(gene1, ' == ', 1))),
                  evidence = eval(parse(text = paste0(gene2, ' == ', 1))),
                  n = 10000)
  
  ##Compute relative risk
  exposed = P_ij
  
  exp_adj <- (exposed * epsilon)
  
  unexposed = cpquery(bn, event = eval(parse(text = paste0(gene1, ' == ', 1))),
                      evidence = eval(parse(text = paste0(gene2, ' == ', 0))),
                      n = 10000)
  unexp_adj <- (unexposed * epsilon)
  
  relodds <- (exposed / unexposed)
  
  ##Return a list of dataframes
  list(probs_data = data.frame(Gene_1 = gene1, Gene_2 = gene2, Conditional_Probability = P_ij),
       risk_data = data.frame(Gene_1 = gene1, Gene_2 = gene2, Relative_Risk = relodds))
}

n_queries <- (N * (N - 1) * Nlists)

cat(paste("\n \033[32mCores registered:\033[39m", n_cores, "\n"))
cat(paste("\n \033[32mNumber of queries:\033[39m", n_queries, "\n"))
cat("\n \033[35mQuerying network. Please be patient.\033[39m\n \n")

counter <- 1
for (i in 1:max_combinations) {
  gene1 <- combinations$Gene_1[i]
  gene2 <- combinations$Gene_2[i]
  
  for (bn in boosts_list) {
    result <- compute_gene_stats(gene1, gene2, bn, epsilon, query_number = counter)
    probs_list[[counter]] <- result$probs_data
    risk_list[[counter]] <- result$risk_data
    counter <- counter + 1
  }
}

probs_data <- do.call(rbind, probs_list)
risk_data <- do.call(rbind, risk_list)

# Compute statistics for conditional probabilities
probs_stats <- probs_data %>%
  group_by(Gene_1, Gene_2) %>%
  summarise(Conditional_Probability_Mean = mean(Conditional_Probability),
            Conditional_Probability_SD = sd(Conditional_Probability),
            Conditional_Probability_CI_low = Conditional_Probability_Mean - qt(0.975, n() - 1) * Conditional_Probability_SD / sqrt(n()),
            Conditional_Probability_CI_high = Conditional_Probability_Mean + qt(0.975, n() - 1) * Conditional_Probability_SD / sqrt(n()))

# Compute statistics for absolute risks
risk_stats <- risk_data %>%
  group_by(Gene_1, Gene_2) %>%
  summarise(Relative_Risk_Mean = mean(Relative_Risk),
            Relative_Risk_SD = sd(Relative_Risk),
            Relative_Risk_CI_low = Relative_Risk_Mean - qt(0.975, n() - 1) * Relative_Risk_SD / sqrt(n()),
            Relative_Risk_CI_high = Relative_Risk_Mean + qt(0.975, n() - 1) * Relative_Risk_SD / sqrt(n()))

# Join the two sets of statistics into one data frame
stats <- full_join(probs_stats, risk_stats, by = c("Gene_1", "Gene_2"))
stats <- na.omit(stats)

write.csv(stats, "Results.csv", row.names = FALSE)

###BDPS
calculate_ratio <- function(stats, gene1, gene2) {
  prob1 <- stats %>%
    filter(Gene_1 == gene1, Gene_2 == gene2) %>%
    pull(Conditional_Probability_Mean)
  prob2 <- stats %>%
    filter(Gene_1 == gene2, Gene_2 == gene1) %>%
    pull(Conditional_Probability_Mean)
  
  if (length(prob1) > 0 && length(prob2) > 0) {
    return(prob1 / prob2) ###BDPS
  } else {
    return(NA)
  }
}

##Run BDPS
result <- stats %>%
  distinct(Gene_1, Gene_2) %>%
  rowwise() %>%
  mutate(BDPS = calculate_ratio(stats, Gene_1, Gene_2)) %>%
  select(Gene_A = Gene_1, Gene_B = Gene_2, BDPS)

result <- result[!is.na(result$BDPS),]

##Fold Change
calculate_fold_change <- function(stats, gene1, gene2) {
  fc1 <- stats %>%
    filter(Gene_1 == gene1, Gene_2 == gene2) %>%
    pull(Relative_Risk_Mean)
  fc2 <- stats %>%
    filter(Gene_1 == gene2, Gene_2 == gene1) %>%
    pull(Relative_Risk_Mean)
  
  if (length(fc1) > 0 && length(fc2) > 0) {
    return((fc1 / fc2) / 2) ###Fold Change
  } else {
    return(NA)
  }
}

fold_change_results <- stats %>%
  distinct(Gene_1, Gene_2) %>%
  rowwise() %>%
  mutate(Fold_Change = calculate_fold_change(stats, Gene_1, Gene_2)) %>%
  select(Gene_A = Gene_1, Gene_B = Gene_2, Fold_Change)

fold_change_results <- fold_change_results[!is.na(fold_change_results$Fold_Change),]

post_hoc <- full_join(result, fold_change_results, by = c("Gene_A", "Gene_B"))
post_hoc <- na.omit(post_hoc)

write.csv(post_hoc, "post_hoc_analysis.csv", row.names = FALSE)

##Stop cluster
parallel::stopCluster(cl)
cat("\n \033[032mStatistics calculated.\033[39m \n")

## Prepare data for the network visualization
net <- igraph.from.graphNEL(as.graphNEL(avg_boot))

# Check if the network has enough nodes and edges
if(vcount(net) < 2 || ecount(net) == 0) {
  stop("Network is too small or has no edges for visualization.")
}

visZach <- toVisNetworkData(net)

nodes <- visZach$nodes
edges <- visZach$edges

# Apply function to each edge

get_width <- function(node_from, node_to, stats) {
  row <- stats %>% filter((Gene_1 == node_from & Gene_2 == node_to) |
                            (Gene_1 == node_to & Gene_2 == node_from))
  if(nrow(row) > 0) {
    width_value <- row$Conditional_Probability_CI_high[1]
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
