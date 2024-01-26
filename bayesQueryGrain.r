args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
metadata_file <- args[2]
output_boot <- args[3]
threads <- as.integer(args[4])
number_of_bootstraps <- as.integer(args[5])
resamples <- as.integer(args[6])

print("User Inputs:")
print(paste("Matrix:", args[1]))
print(paste("Metadata:", args[2]))
print(paste("Bootstrapped Network Name:", args[3]))
print(paste("Number of Threads:", args[4]))
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
        BiocManager::instal(pkgs)
    }
}

data <- read.csv(input_file, row.names = 1)
d_fact <- data %>% mutate_if(is.numeric, as.factor)

n_cores <- threads
print(paste("Using", threads, "threads for parallel processing."))
cl = parallel::makeCluster(n_cores)
clusterSetRNGStream(cl, 13245)

required_bootstraps = number_of_bootstraps

boot = boot.strength(data = d_fact, R = number_of_bootstraps, algorithm = "hc",
                     algorithm.args = list(score="bde", iss=10), cluster = cl)

print("Bootstrapping Finished")
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
print("Resample Number Defined")

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
print("Combinations defined")

max_combinations <- nrow(combinations)

probs_list <- vector("list", max_combinations * length(boosts_list))
risk_list <- vector("list", max_combinations * length(boosts_list))

###Function to query the network
compute_gene_stats <- function(gene1, gene2, grain_net, epsilon) {
  
  ##Compute conditional probability
  P_ij <- querygrain(setEvidence(grain_net, nodes = c(gene1), states = c("1"), propagate = TRUE), nodes = gene2)[[1]][2]
  
  ##Compute relative risk
  exposed <- P_ij
  
  unexposed <- querygrain(setEvidence(grain_net, nodes = c(gene1), states = c("0"), propagate = TRUE), nodes = gene2)[[1]][2]
  
  relodds <- (exposed * epsilon) / (unexposed * epsilon)
  
  ##Return a list of dataframes
  list(probs_data = data.frame(Gene_1 = gene1, Gene_2 = gene2, Conditional_Probability = P_ij),
       risk_data = data.frame(Gene_1 = gene1, Gene_2 = gene2, Relative_Risk = relodds))
}

##Execute queries
counter <- 1
for (i in 1:max_combinations) {
  gene1 <- combinations$Gene_1[i]
  gene2 <- combinations$Gene_2[i]
  
  for (bn in boosts_list) {
    
    grain_net <- compile(bnlearn::as.grain(bn), propagate = T)
    
    results <- compute_gene_stats(gene1, gene2, grain_net, epsilon)
    probs_list[[counter]] <- results$probs_data
    risk_list[[counter]] <- results$risk_data
    counter <- counter + 1
  }
}

print("Queries Finished")

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

# Write data frame to a CSV file
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
    return (prob1 / prob2) ###BDPS
  } else {
    return(NA)
  }
}

result <- stats %>%
  distinct(Gene_1, Gene_2) %>%
  rowwise() %>%
  mutate(Score = calculate_ratio(stats, Gene_1, Gene_2)) %>%
  select(Gene_A = Gene_1, Gene_B = Gene_2, Score)

result <- result[!is.na(result$Score),]

write.csv(result, "BDPS.csv", row.names = FALSE)

##Stop cluster
parallel::stopCluster(cl)
print("Stats Created")

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

nodes$group <- lookup[valid_genes]

# Colors
grp <- as.numeric(as.factor(nodes$group))
n <- length(unique(grp))
colors <- brewer.pal.info[brewer.pal.info$colorblind == T, ]
col_vec <- unlist(mapply(brewer.pal, colors$maxcolors, rownames(colors)))
colSide <- sample(col_vec, n)[grp]
nodes$color <- colSide

# Prepare legend
lnodes <- data.frame(color = unique(nodes$color),
                     label = unique(nodes$group),
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