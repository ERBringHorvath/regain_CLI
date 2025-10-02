#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(optparse)
})

options(repos = c(CRAN = "https://cloud.r-project.org"),
        dplyr.summarise.inform = FALSE)

options(stringsAsFactors = FALSE)

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
              help="Imaginary sample size for BDe score [default: %default]"),
  make_option(c("--no-viz"), action="store_true", default=FALSE, dest = "no_viz",
              help="Skip HTML/PDF visualization"),
  make_option(c("--cp-samples"), type="integer", default=10000,
              help="Monte Carlo samples for cpquery [default: %default]")
)

parser <- OptionParser(option_list = opt_list,
                       description = "ReGAIN — Bayesian network (CPQuery) structure learning + queries")
opt <- parse_args(parser)

no_viz <- isTRUE(opt$no_viz)
cp_samples <- opt$cp_samples

required <- c("input","metadata","output_boot","number_of_bootstraps","resamples")
missing <- required[!nzchar(trimws(sapply(required, function(k) as.character(opt[[k]])))) ]
if (length(missing)) {
  print_help(parser); stop(paste("Missing required:", paste(missing, collapse=", ")), call.=FALSE)
}

pkgs <- c(
  'dplyr','parallel','pbapply','BiocManager','RColorBrewer','visNetwork',
  'igraph','reshape2','doParallel','scales',
  'tidygraph','ggraph','ggplot2','tidyr'   # for PDF plotting
)
miss <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(miss) > 0) install.packages(miss)

bioc_pkgs <- c('bnlearn','graph','progressr')  # gRain not required for cpquery
miss_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[,"Package"])]
if (length(miss_bioc) > 0) BiocManager::install(miss_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(parallel); library(doParallel); library(pbapply)
  library(RColorBrewer); library(visNetwork); library(igraph); library(reshape2)
  library(bnlearn); library(graph); library(scales); library(foreach)
  library(tidygraph); library(ggraph); library(ggplot2); library(tidyr)
})

input_file    <- opt$input
metadata_file <- opt$metadata
output_boot   <- if (grepl("\\.rds$", opt$output_boot, ignore.case=TRUE)) opt$output_boot else paste0(opt$output_boot, ".rds")
threads       <- opt$threads
nboots        <- opt$number_of_bootstraps
resamples     <- opt$resamples
iss           <- opt$iss
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
cat("  Visualization:        ", if (no_viz) "disabled (--no-viz)" else "enabled", "\n")

stopifnot(file.exists(input_file), file.exists(metadata_file))
data   <- read.csv(input_file, row.names = 1, check.names = TRUE)
d_fact <- data %>% mutate_if(is.numeric, as.factor)
vars   <- colnames(d_fact)

## ---- optional blacklist ----
blacklist_filtered <- NULL
if (!is.null(blacklist_path)) {
  stopifnot(file.exists(blacklist_path))
  bl <- read.csv(blacklist_path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(bl) < 2) stop("Blacklist must have 2 columns: from,to (no header).")
  colnames(bl)[1:2] <- c("from","to")
  blacklist_filtered <- subset(bl, from %in% vars & to %in% vars)
  cat("  Blacklist arcs (after filtering to present variables): ", nrow(blacklist_filtered), "\n")
}

## ---- parallel pool ----
n_cores <- max(1, as.integer(threads))
cl <- parallel::makeCluster(n_cores)
on.exit(stopCluster(cl), add = TRUE)
registerDoParallel(cl)
try({ clusterSetRNGStream(cl, 13245) }, silent = TRUE)

cat("\n \033[32mBootstrapping started.\033[39m\n\n")
algo_args <- list(score = "bde", iss = iss)
if (!is.null(blacklist_filtered) && nrow(blacklist_filtered) > 0) {
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

## ---- prepare to query (cpquery) ----
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
combinations <- expand.grid(Gene_1 = valid_genes, Gene_2 = valid_genes)
combinations <- subset(combinations, Gene_1 != Gene_2)

# remove blacklisted directed pairs from query combinations
if (!is.null(blacklist_filtered) && nrow(blacklist_filtered) > 0) {
  bl_keys <- paste(blacklist_filtered$from, blacklist_filtered$to, sep = "___")
  combo_keys <- paste(combinations$Gene_1, combinations$Gene_2, sep = "___")
  keep <- !(combo_keys %in% bl_keys)
  combinations <- combinations[keep, , drop = FALSE]
}

cat(paste("\n \033[32mCores registered:\033[39m", n_cores, "\n"))
cat(paste("\n \033[32mNumber of queries:\033[39m", nrow(combinations) * length(boosts_list), "\n"))
cat("\n \033[35mQuerying network via cpquery. Please be patient.\033[39m\n\n")

# Safely build an expression like `geneName` == '1' / '0'
.make_expr <- function(node, val) {
  parse(text = sprintf("`%s` == '%s'", node, val))[[1]]
}

epsilon = 1e-9

compute_gene_stats <- function(gene1, gene2, bn_fit, epsilon = epsilon, cp_samples = 10000L) {
  # Build expressions (do NOT eval() them here)
  e_event <- .make_expr(gene2, "1")
  e_on    <- .make_expr(gene1, "1")
  e_off   <- .make_expr(gene1, "0")
  
  # Inline the expressions into cpquery() so substitute(event) sees the real exprs
  call_on  <- substitute(bnlearn::cpquery(bn_fit, event = EV, evidence = ED, n = cp_samples),
                         list(EV = e_event, ED = e_on))
  call_off <- substitute(bnlearn::cpquery(bn_fit, event = EV, evidence = ED, n = cp_samples),
                         list(EV = e_event, ED = e_off))
  
  exposed   <- eval(call_on)
  unexposed <- eval(call_off)
  
  # guard against NA from impossible evidence
  if (is.na(exposed))   exposed   <- 0
  if (is.na(unexposed)) unexposed <- 0
  
  rr      <- (exposed + epsilon) / (unexposed + epsilon)
  absdiff <- exposed - unexposed
  
  list(
    probs_data = data.frame(Gene_1 = gene1, Gene_2 = gene2,
                            Conditional_Probability = exposed),
    risk_data  = data.frame(Gene_1 = gene1, Gene_2 = gene2,
                            Relative_Risk = rr),
    abs_data   = data.frame(Gene_1 = gene1, Gene_2 = gene2,
                            Absolute_Risk = absdiff)
  )
}

## ---- parallel query loop ----
results <- foreach(i = seq_len(nrow(combinations)),
                   .packages = c("bnlearn", "dplyr")) %dopar% {
 gene1 <- combinations$Gene_1[i]
 gene2 <- combinations$Gene_2[i]
 
 temp_probs <- vector("list", length(boosts_list))
 temp_risk  <- vector("list", length(boosts_list))
 temp_abs   <- vector("list", length(boosts_list))
 
 for (j in seq_along(boosts_list)) {
   bn_fit <- boosts_list[[j]]
   res <- compute_gene_stats(
     gene1, gene2, bn_fit,
     epsilon    = epsilon,      # pass explicitly
     cp_samples = cp_samples    # for bnL; omit if not defined in bnS
   )
   temp_probs[[j]] <- res$probs_data
   temp_risk[[j]]  <- res$risk_data
   temp_abs[[j]]   <- res$abs_data
 }
 
 list(
   probs_data = do.call(rbind, temp_probs),
   risk_data  = do.call(rbind, temp_risk),
   abs_data   = do.call(rbind, temp_abs)
 )
}

probs_data <- do.call(rbind, lapply(results, `[[`, "probs_data"))
risk_data  <- do.call(rbind, lapply(results, `[[`, "risk_data"))
abs_data   <- do.call(rbind, lapply(results, `[[`, "abs_data"))

## ---- summaries ----
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

abs_stats <- abs_data %>%
  group_by(Gene_1, Gene_2) %>%
  summarise(
    Absolute_Risk_Mean    = mean(Absolute_Risk),
    Absolute_Risk_SD      = sd(Absolute_Risk),
    Absolute_Risk_CI_low  = Absolute_Risk_Mean - qt(0.975, n() - 1) * Absolute_Risk_SD / sqrt(n()),
    Absolute_Risk_CI_high = Absolute_Risk_Mean + qt(0.975, n() - 1) * Absolute_Risk_SD / sqrt(n()),
    .groups = "drop"
  )

stats <- full_join(probs_stats, risk_stats, by = c("Gene_1","Gene_2")) %>%
  left_join(abs_stats, by = c("Gene_1","Gene_2")) %>%
  na.omit()

write.csv(stats, "Query_Results.csv", row.names = FALSE)

## ---- BDPS / Fold Change (same as bnS; no baseline here) ----
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

cat(" \033[032mStatistics calculated.\033[39m \n")

if (no_viz) {
  cat("[INFO] Visualzation disabled (--no-viz).\n")
} else {
  
  ## ---- network viz: HTML + PDF (edge width by Absolute_Risk_CI_high; color by RR) ----
  net <- graph_from_graphnel(as.graphNEL(avg_boot))
  if (vcount(net) < 2 || ecount(net) == 0) stop("Network is too small or has no edges for visualization.")
  
  # shared layout
  set.seed(42)
  layout_fr <- igraph::layout_with_fr(net)
  coords <- data.frame(name = igraph::V(net)$name, x = layout_fr[,1], y = layout_fr[,2], stringsAsFactors = FALSE)
  
  visDat <- toVisNetworkData(net)
  nodes <- visDat$nodes
  edges <- visDat$edges
  
  # keep valid nodes; drop edges referencing removed nodes
  nodes <- nodes[nodes$id %in% valid_genes, ]
  edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, , drop = FALSE]
  
  # node groups/colors
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
  
  # attach coords
  nodes <- nodes %>% left_join(coords, by = c("id" = "name"))
  
  # width/color from stats (unordered pair aggregation)
  stats_pairs <- stats %>%
    transmute(
      Gene_lo = pmin(as.character(Gene_1), as.character(Gene_2)),
      Gene_hi = pmax(as.character(Gene_1), as.character(Gene_2)),
      Absolute_Risk_CI_high,
      Relative_Risk_Mean
    ) %>%
    mutate(pair_key = paste(Gene_lo, Gene_hi, sep = "||"))
  
  pair_summ <- stats_pairs %>%
    group_by(pair_key) %>%
    summarize(
      width_val = suppressWarnings(max(Absolute_Risk_CI_high, na.rm = TRUE)),
      rr_min    = suppressWarnings(min(Relative_Risk_Mean,   na.rm = TRUE)),
      .groups   = "drop"
    ) %>%
    mutate(
      width_val  = ifelse(is.finite(width_val) & width_val > 0, width_val, NA_real_),
      edge_color = ifelse(rr_min < 1, "red", "black")
    )
  
  edges <- edges %>%
    mutate(pair_key = paste(pmin(as.character(from), as.character(to)),
                            pmax(as.character(from), as.character(to)),
                            sep = "||")) %>%
    left_join(pair_summ, by = "pair_key")
  
  edges$width <- scales::rescale(replace_na(edges$width_val, 0), to = c(1, 5))
  edges$color <- ifelse(is.na(edges$edge_color), "black", edges$edge_color)
  
  lnodes <- data.frame(color = color_palette, label = unique_groups, font.color = 'white')
  
  # HTML
  network <- visNetwork(nodes = nodes, edges = edges, width = '100%', height = 900) %>%
    visNodes(size = 20, color = list(highlight = 'yellow'), font = list(size = 25)) %>%
    visEdges(smooth = list(enabled = TRUE, type = 'diagonalCross', roundness = 0.1),
             physics = FALSE, color = "black") %>%
    visIgraphLayout(layout = "layout_with_fr", type = 'full') %>%
    visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE, labelOnly = TRUE),
               nodesIdSelection = TRUE) %>%
    visLegend(addNodes = lnodes, width = 0.1, position = 'left',
              main = 'Variable Type', ncol = 2, useGroups = FALSE)
  visSave(network, "Bayesian_Network.html", background = "#F5F4F4")
  
  # PDF
  tbl <- as_tbl_graph(net) %>%
    activate(nodes) %>%
    left_join(nodes %>% transmute(name = id, x, y, node_color = unname(color)), by = "name")
  
  edge_attrs <- edges %>%
    transmute(pair_key = paste(pmin(as.character(from), as.character(to)),
                               pmax(as.character(from), as.character(to)),
                               sep = "||"),
              width, color)
  
  node_names <- tbl %>% activate(nodes) %>% pull(name)
  tbl <- tbl %>%
    activate(edges) %>%
    mutate(pair_key = paste(pmin(node_names[from], node_names[to]),
                            pmax(node_names[from], node_names[to]),
                            sep = "||")) %>%
    left_join(edge_attrs, by = "pair_key") %>%
    mutate(width = replace_na(width, 1),
           color = replace_na(color, "black"))
  
  p <- ggraph(tbl, layout = "manual", x = x, y = y) +
    geom_edge_link(aes(edge_width = width, edge_colour = color),
                   lineend = "round", show.legend = TRUE) +
    scale_edge_width(range = c(1, 3), guide = "none") +
    scale_edge_colour_identity(name = "Edge (RR<1 = red)",
                               labels = c("black" = "RR >= 1", "red" = "RR < 1"),
                               breaks = c("black", "red"),
                               guide = "legend") +
    geom_node_point(aes(fill = node_color),
                    shape = 21, size = 10, stroke = 0.3, colour = "grey30") +
    scale_fill_identity(name = "Variable Type",
                        breaks = unname(group_color_lookup),
                        labels = names(group_color_lookup),
                        guide = "legend") +
    geom_node_text(aes(label = name), size = 4.8, vjust = -1) +
    theme_void() +
    theme(legend.position = "left")
  
  ggsave("Bayesian_Network.pdf", p, width = 14, height = 10, units = "in")
  cat("\n \033[32mSaved Bayesian_Network.html and Bayesian_Network.pdf.\033[39m\n\n")
  
}

