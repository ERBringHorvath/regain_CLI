#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
})

options(repos = c(CRAN = "https://cloud.r-project.org"),
        dplyr.summarise.inform = FALSE)

## ---- CLI ----
opt_list <- list(
  make_option(c("-i","--input"), type="character", help="Input matrix CSV (row.names in col 1)"),
  make_option(c("-M","--metadata"), type="character", help="Metadata CSV (two columns: gene, class/label)"),
  make_option(c("-o","--output_boot"), type="character", help="Output bootstrapped network filename (.rds)"),
  make_option(c("-T","--threads"), type="integer", default=2, help="Threads [default: %default]"),
  make_option(c("-n","--number_of_bootstraps"), type="integer", help="Number of bootstraps (e.g., 300–500)"),
  make_option(c("-r","--resamples"), type="integer", help="Number of resamples for querying"),
  make_option(c("-b","--blacklist"), type="character", default=NULL,
              help="Optional blacklist CSV (no header): from,to"),
  make_option(c("--iss"), type="integer", default=10,
              help="Imaginary sample size for BDe score [default: %default]")
)

parser <- OptionParser(option_list = opt_list,
                       description = "ReGAIN — Bayesian network structure learning + gRain queries")
opt <- parse_args(parser)

required <- c("input","metadata","output_boot","number_of_bootstraps","resamples")
missing <- required[!nzchar(trimws(sapply(required, function(k) as.character(opt[[k]])))) ]
if (length(missing)) {
  print_help(parser); stop(paste("Missing required:", paste(missing, collapse=", ")), call.=FALSE)
}

## ---- dependency check/install (same as before) ----
pkgs <- c('dplyr','parallel','pbapply','BiocManager','RColorBrewer','visNetwork','igraph','reshape2','doParallel','scales')
miss <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(miss) > 0) install.packages(miss)

bioc_pkgs <- c('bnlearn','gRain','progressr','graph')
miss_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[,"Package"])]
if (length(miss_bioc) > 0) BiocManager::install(miss_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(parallel); library(doParallel); library(pbapply)
  library(RColorBrewer); library(visNetwork); library(igraph); library(reshape2)
  library(bnlearn); library(gRain); library(graph); library(scales); library(foreach)
})

## ---- inputs & echo ----
input_file  <- opt$input
metadata_file <- opt$metadata
output_boot <- if (grepl("\\.rds$", opt$output_boot, ignore.case=TRUE)) opt$output_boot else paste0(opt$output_boot, ".rds")
threads     <- opt$threads
nboots      <- opt$number_of_bootstraps
resamples   <- opt$resamples
iss         <- opt$iss
blacklist_path <- opt$blacklist

cat("User Inputs:\n")
cat("  Matrix:               ", input_file, "\n")
cat("  Metadata:             ", metadata_file, "\n")
cat("  Bootstrapped Network: ", output_boot, "\n")
cat("  Threads:              ", threads, "\n")
cat("  # Bootstraps:         ", nboots, "\n")
cat("  # Resamples:          ", resamples, "\n")
cat("  ISS (BDe):            ", iss, "\n")
if (!is.null(blacklist_path)) cat("  Blacklist:            ", blacklist_path, "\n") else cat("  Blacklist:            (none)\n")

stopifnot(file.exists(input_file), file.exists(metadata_file))
data   <- read.csv(input_file, row.names = 1, check.names = TRUE)
d_fact <- data %>% mutate_if(is.numeric, as.factor)
vars   <- colnames(d_fact)

## ---- optional blacklist load & filter ----
blacklist_filtered <- NULL
if (!is.null(blacklist_path)) {
  stopifnot(file.exists(blacklist_path))
  bl <- read.csv(blacklist_path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(bl) < 2) stop("Blacklist must have 2 columns: from,to (no header).")
  colnames(bl)[1:2] <- c("from","to")
  # Keep only arcs whose endpoints are in the dataset
  blacklist_filtered <- subset(bl, from %in% vars & to %in% vars)
  cat("  Blacklist arcs (after filtering to present variables): ", nrow(blacklist_filtered), "\n")
}

## ---- parallel pool ----
n_cores <- max(1, as.integer(threads))
cl <- parallel::makeCluster(n_cores)
on.exit(stopCluster(cl), add = TRUE)
registerDoParallel(cl)

## ---- bootstrap ----
cat("\n \033[32mBootstrapping started.\033[39m\n\n")

algo_args <- list(score = "bde", iss = iss)
if (!is.null(blacklist_filtered) && nrow(blacklist_filtered) > 0) {
  # bnlearn expects a two-column data.frame with names 'from','to'
  algo_args$blacklist <- blacklist_filtered[, c("from","to")]
}

boot <- boot.strength(
  data = d_fact,
  R = nboots,
  algorithm = "hc",
  algorithm.args = algo_args,
  cluster = cl
)

saveRDS(boot, file = output_boot)

avg_boot <- averaged.network(boot, threshold = 0.5)
arcs(avg_boot) <- directed.arcs(avg_boot)

## If blacklist provided, drop those arcs from the averaged network
if (!is.null(blacklist_filtered) && nrow(blacklist_filtered) > 0) {
  avg_boot_clean <- avg_boot
  for (i in seq_len(nrow(blacklist_filtered))) {
    avg_boot_clean <- drop.arc(avg_boot_clean,
                               from = blacklist_filtered$from[i],
                               to   = blacklist_filtered$to[i]
    )
  }
  avg_boot <- avg_boot_clean
}

## ---- prepare to query ----
metadata <- read.csv(metadata_file, check.names = FALSE)
gene_names <- metadata[,1]
lookup <- setNames(as.character(metadata[,2]), metadata[,1])
valid_genes <- intersect(gene_names, colnames(data))

Nlists <- resamples
boosts <- function(d_fact, Nlists, avg_boot) {
  lapply(seq_len(Nlists), function(i) {
    sample_data <- dplyr::slice_sample(d_fact, n = nrow(d_fact), replace = TRUE)
    bnlearn::bn.fit(avg_boot, sample_data, method = 'bayes')
  })
}
boosts_list <- boosts(d_fact, Nlists, avg_boot)

N <- length(valid_genes)
epsilon <- ((N + 0.5) / (N + 1))  # note: currently cancels if multiplied both numerator/denominator

# All ordered pairs, excluding self-pairs
combinations <- expand.grid(Gene_1 = valid_genes, Gene_2 = valid_genes)
combinations <- subset(combinations, Gene_1 != Gene_2)

# If a blacklist exists, remove those directed pairs from combinations
if (!is.null(blacklist_filtered) && nrow(blacklist_filtered) > 0) {
  bl_keys <- paste(blacklist_filtered$from, blacklist_filtered$to, sep = "___")
  combo_keys <- paste(combinations$Gene_1, combinations$Gene_2, sep = "___")
  keep <- !(combo_keys %in% bl_keys)
  combinations <- combinations[keep, , drop = FALSE]
}

cat(paste("\n \033[32mCores registered:\033[39m", n_cores, "\n"))
cat(paste("\n \033[32mNumber of queries:\033[39m", nrow(combinations) * length(boosts_list), "\n"))
cat("\n \033[35mQuerying network. Please be patient.\033[39m\n\n")

compute_gene_stats <- function(gene1, gene2, grain_net, epsilon) {
  P_ij <- querygrain(setEvidence(grain_net, nodes = gene1, states = "1", propagate = TRUE), nodes = gene2)[[1]][2]
  unexposed <- querygrain(setEvidence(grain_net, nodes = gene1, states = "0", propagate = TRUE), nodes = gene2)[[1]][2]
  # NOTE on epsilon: if you intended smoothing, consider (P_ij + eps) / (unexposed + eps)
  relodds <- (P_ij * epsilon) / (unexposed * epsilon)
  list(
    probs_data = data.frame(Gene_1 = gene1, Gene_2 = gene2, Conditional_Probability = P_ij),
    risk_data  = data.frame(Gene_1 = gene1, Gene_2 = gene2, Relative_Risk = relodds)
  )
}

results <- foreach(i = seq_len(nrow(combinations)), .packages = c("bnlearn","dplyr","gRain")) %dopar% {
  gene1 <- combinations$Gene_1[i]
  gene2 <- combinations$Gene_2[i]
  temp_probs <- vector("list", length(boosts_list))
  temp_risk  <- vector("list", length(boosts_list))
  idx <- 1L
  for (bn in boosts_list) {
    grain_net <- compile(as.grain(bn), propagate = TRUE)
    res <- compute_gene_stats(gene1, gene2, grain_net, epsilon)
    temp_probs[[idx]] <- res$probs_data
    temp_risk[[idx]]  <- res$risk_data
    idx <- idx + 1L
  }
  list(probs_data = do.call(rbind, temp_probs), risk_data = do.call(rbind, temp_risk))
}

probs_data <- do.call(rbind, lapply(results, `[[`, "probs_data"))
risk_data  <- do.call(rbind, lapply(results, `[[`, "risk_data"))

probs_stats <- probs_data %>%
  group_by(Gene_1, Gene_2) %>%
  summarise(
    Conditional_Probability_Mean = mean(Conditional_Probability),
    Conditional_Probability_SD   = sd(Conditional_Probability),
    Conditional_Probability_CI_low  = Conditional_Probability_Mean - qt(0.975, n() - 1) * Conditional_Probability_SD / sqrt(n()),
    Conditional_Probability_CI_high = Conditional_Probability_Mean + qt(0.975, n() - 1) * Conditional_Probability_SD / sqrt(n()),
    .groups = "drop"
  )

risk_stats <- risk_data %>%
  group_by(Gene_1, Gene_2) %>%
  summarise(
    Relative_Risk_Mean = mean(Relative_Risk),
    Relative_Risk_SD   = sd(Relative_Risk),
    Relative_Risk_CI_low  = Relative_Risk_Mean - qt(0.975, n() - 1) * Relative_Risk_SD / sqrt(n()),
    Relative_Risk_CI_high = Relative_Risk_Mean + qt(0.975, n() - 1) * Relative_Risk_SD / sqrt(n()),
    .groups = "drop"
  )

stats <- full_join(probs_stats, risk_stats, by = c("Gene_1","Gene_2")) %>% na.omit()
write.csv(stats, "Results.csv", row.names = FALSE)

## ---- BDPS / Fold Change (unchanged) ----
calculate_ratio <- function(stats, gene1, gene2) {
  prob1 <- stats %>% filter(Gene_1 == gene1, Gene_2 == gene2) %>% pull(Conditional_Probability_Mean)
  prob2 <- stats %>% filter(Gene_1 == gene2, Gene_2 == gene1) %>% pull(Conditional_Probability_Mean)
  if (length(prob1) > 0 && length(prob2) > 0) return(prob1 / prob2)
  NA
}

result <- stats %>%
  distinct(Gene_1, Gene_2) %>%
  rowwise() %>%
  mutate(BDPS = calculate_ratio(stats, Gene_1, Gene_2)) %>%
  select(Gene_A = Gene_1, Gene_B = Gene_2, BDPS) %>%
  filter(!is.na(BDPS))

calculate_fold_change <- function(stats, gene1, gene2) {
  fc1 <- stats %>% filter(Gene_1 == gene1, Gene_2 == gene2) %>% pull(Relative_Risk_Mean)
  fc2 <- stats %>% filter(Gene_1 == gene2, Gene_2 == gene1) %>% pull(Relative_Risk_Mean)
  if (length(fc1) > 0 && length(fc2) > 0) return((fc1 / fc2) / 2)
  NA
}

fold_change_results <- stats %>%
  distinct(Gene_1, Gene_2) %>%
  rowwise() %>%
  mutate(Fold_Change = calculate_fold_change(stats, Gene_1, Gene_2)) %>%
  select(Gene_A = Gene_1, Gene_B = Gene_2, Fold_Change) %>%
  filter(!is.na(Fold_Change))

post_hoc <- full_join(result, fold_change_results, by = c("Gene_A","Gene_B")) %>% na.omit()
write.csv(post_hoc, "post_hoc_analysis.csv", row.names = FALSE)

stopImplicitCluster()
cat(" \033[032mStatistics calculated.\033[39m \n")

## ---- network viz (unchanged logic, now using any cleaned avg_boot) ----
net <- graph_from_graphnel(as.graphNEL(avg_boot))
if (vcount(net) < 2 || ecount(net) == 0) stop("Network is too small or has no edges for visualization.")
visZach <- toVisNetworkData(net)
nodes <- visZach$nodes; edges <- visZach$edges

get_width <- function(node_from, node_to, stats) {
  row <- stats %>% filter((Gene_1 == node_from & Gene_2 == node_to) | (Gene_1 == node_to & Gene_2 == node_from))
  if (nrow(row) > 0) {
    width_value <- row$Conditional_Probability_CI_high[1]
    if (width_value <= 0) stop("Encountered non-positive width value.")
    return(width_value)
  }
  0
}

edges$width <- sapply(seq_len(nrow(edges)), function(i) get_width(edges$from[i], edges$to[i], stats))
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
if (length(unique_groups) > length(color_palette)) color_palette <- rep(color_palette, length.out = length(unique_groups))
group_color_lookup <- setNames(color_palette, unique_groups)
nodes$color <- group_color_lookup[nodes$group]

lnodes <- data.frame(color = color_palette, label = unique_groups, font.color = 'white')

network <- visNetwork(nodes = nodes, edges = edges, width = '100%', height = 900) %>%
  visNodes(size = 20, color = list(highlight = 'yellow'), font = list(size = 25)) %>%
  visEdges(smooth = list(enabled = TRUE, type = 'diagonalCross', roundness = 0.1), physics = FALSE, color = "black") %>%
  visIgraphLayout(layout = "layout_with_fr", type = 'full') %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE, labelOnly = TRUE), nodesIdSelection = TRUE) %>%
  visLegend(addNodes = lnodes, width = 0.1, position = 'left', main = 'Gene Class', ncol = 2, useGroups = FALSE)

visSave(network, "Bayesian_Network.html", background = "#F5F4F4")
cat("\n \033[32mAnalysis complete.\033[39m\n\n")
