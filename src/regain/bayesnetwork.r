#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 University of Utah

suppressPackageStartupMessages({ library(optparse) })

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"),
        dplyr.summarise.inform = FALSE)

opt_list <- list(
  make_option(c("-N","--boot"),      type="character", help="Bootstrapped network .rds (from bnS/bnL)"),
  make_option(c("-i","--input"),      type="character", help="Input data matrix CSV (row names in col 1)"),
  make_option(c("-M","--metadata"),  type="character", help="Metadata CSV (two columns: gene, class/label)"),
  make_option(c("-s","--stats"),     type="character", help="Stats CSV from bnS/bnL (has Gene_1,Gene_2, RR/CI cols)"),
  make_option(c("-t","--threshold"), type="double", default=0.5, help="Averaged network threshold [default: %default]"),
  make_option(c("--seed"),           type="integer", default=42,  help="Layout seed [default: %default]"),
  make_option(c("--html-out"),       type="character", default="Bayesian_Network.html", help="HTML output [default: %default]"),
  make_option(c("--pdf-out"),        type="character", default="Bayesian_Network.pdf",  help="PDF output [default: %default]"),
  make_option(c("-b","--blacklist"), type="character", default=NULL, help="Optional blacklist CSV (no header): from,to"),
  make_option(c("--width-metric"),   type="character", default="auto",
              help="Edge-width metric: auto | abs_ci | cp_ci | cp_mean [default: %default]"),
  make_option(c("--rr-threshold"),   type="double", default=1.0, help="RR color threshold [default: %default]")
)

parser <- OptionParser(option_list = opt_list,
                       description = "ReGAIN — Standalone Bayesian network visualization")
opt <- parse_args(parser)

required <- c("boot","input","metadata","stats")
miss <- required[ vapply(required, function(k) { v <- opt[[k]]; is.null(v) || !nzchar(trimws(as.character(v))) }, logical(1)) ]
if (length(miss)) { print_help(parser); stop(paste("Missing required:", paste(miss, collapse=", ")), call.=FALSE) }

cran_pkgs <- c("dplyr","tidyr","igraph","visNetwork","RColorBrewer","tidygraph","ggraph","ggplot2","scales")
bioc_pkgs <- c("bnlearn","graph")
inst <- rownames(installed.packages())
need_cran <- setdiff(cran_pkgs, inst)
if (length(need_cran)) install.packages(need_cran)
need_bioc <- setdiff(bioc_pkgs, inst)
if (length(need_bioc)) { suppressPackageStartupMessages({ if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install(need_bioc, ask=FALSE, update=FALSE) }) }

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(igraph); library(visNetwork); library(RColorBrewer)
  library(tidygraph); library(ggraph); library(ggplot2); library(scales); library(bnlearn); library(graph)
})

# ---- inputs ----
boot_path <- opt$boot
data_path <- opt$input
meta_path <- opt$metadata
stats_path <- opt$stats
thr       <- opt$threshold
seed      <- opt$seed
html_out  <- as.character(opt[["html-out"]]) #safer than opt$`html-out`
pdf_out   <- as.character(opt[["pdf-out"]])
bl_path   <- opt$blacklist
width_mode<- tolower(as.character(opt[["width-metric"]]))
rr_thr    <- as.numeric(opt[["rr-threshold"]])

stopifnot(file.exists(boot_path), file.exists(data_path), file.exists(meta_path), file.exists(stats_path))

cat("Inputs:\n")
cat("  Boot RDS:  ", boot_path, "\n")
cat("  Data CSV:  ", data_path, "\n")
cat("  Metadata:  ", meta_path, "\n")
cat("  Stats CSV: ", stats_path, "\n")
cat("  Threshold: ", thr, "\n")
if (!is.null(bl_path)) cat("  Blacklist: ", bl_path, "\n")

# ---- load data ----cran_pkgs <- c("dplyr","tidyr","igraph","visNetwork","RColorBrewer","tidygraph","ggraph","ggplot2","scales")
bioc_pkgs <- c("bnlearn","graph")
inst <- rownames(installed.packages())
need_cran <- setdiff(cran_pkgs, inst)
if (length(need_cran)) install.packages(need_cran)
need_bioc <- setdiff(bioc_pkgs, inst)
if (length(need_bioc)) { suppressPackageStartupMessages({ if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install(need_bioc, ask=FALSE, update=FALSE) }) }

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(igraph); library(visNetwork); library(RColorBrewer)
  library(tidygraph); library(ggraph); library(ggplot2); library(scales); library(bnlearn); library(graph)
})

# ---- inputs ----
boot_path <- opt$boot
data_path <- opt$input
meta_path <- opt$metadata
stats_path <- opt$stats
thr       <- opt$threshold
seed      <- opt$seed
html_out  <- opt$`html-out`
pdf_out   <- opt$`pdf-out`
bl_path   <- opt$blacklist
width_mode<- tolower(opt$`width-metric`)
rr_thr    <- opt$`rr-threshold`

stopifnot(file.exists(boot_path), file.exists(data_path), file.exists(meta_path), file.exists(stats_path))

cat("Inputs:\n")
cat("  Boot RDS:  ", boot_path, "\n")
cat("  Data CSV:  ", data_path, "\n")
cat("  Metadata:  ", meta_path, "\n")
cat("  Stats CSV: ", stats_path, "\n")
cat("  Threshold: ", thr, "\n")
if (!is.null(bl_path)) cat("  Blacklist: ", bl_path, "\n")

# ---- load data ----
boot <- readRDS(boot_path)
avg_boot <- averaged.network(boot, threshold = thr)
arcs(avg_boot) <- directed.arcs(avg_boot)

if (!is.null(bl_path)) {
  stopifnot(file.exists(bl_path))
  bl <- read.csv(bl_path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(bl) < 2) stop("Blacklist must have 2 columns: from,to (no header).")
  colnames(bl)[1:2] <- c("from","to")
  vars <- nodes(avg_boot)
  bl <- subset(bl, from %in% vars & to %in% vars)
  if (nrow(bl)) {
    for (i in seq_len(nrow(bl))) {
      avg_boot <- drop.arc(avg_boot, from = bl$from[i], to = bl$to[i])
    }
  }
}

data <- read.csv(data_path, row.names = 1, check.names = TRUE)
meta <- read.csv(meta_path, check.names = FALSE)
gene_names <- meta[,1]
lookup <- setNames(as.character(meta[,2]), meta[,1])
valid_genes <- intersect(gene_names, colnames(data))

stats <- read.csv(stats_path, check.names = TRUE)

pick_width_col <- function(stats, mode="auto") {
  mode <- match.arg(mode, c("auto","abs_mean","abs_ci","cp_ci","cp_mean"))
  has <- function(col) col %in% colnames(stats)
  if (mode == "abs_ci" && has("Absolute_Risk_CI_high")) return("Absolute_Risk_CI_high")
  if (mode == "abs_mean" && has("Absolute_Risk_Mean")) return("Absolute_Risk_Mean")
  if (mode == "cp_ci"  && has("Conditional_Probability_CI_high")) return("Conditional_Probability_CI_high")
  if (mode == "cp_mean"&& has("Conditional_Probability_Mean")) return("Conditional_Probability_Mean")
  # auto fallback preference:
  if (has("Absolute_Risk_CI_high")) return("Absolute_Risk_CI_high")
  if (has("Absolute_Risk_Mean")) return("Absolute_Risk_Mean")
  if (has("Conditional_Probability_CI_high")) return("Conditional_Probability_CI_high")
  if (has("Conditional_Probability_Mean")) return("Conditional_Probability_Mean")
  stop("Could not find a suitable width metric column in stats. Expected one of: Absolute_Risk_CI_high, Conditional_Probability_CI_high, Conditional_Probability_Mean")
}
wcol <- pick_width_col(stats, width_mode)

rr_col <- if ("Relative_Risk_Mean" %in% colnames(stats)) "Relative_Risk_Mean" else
          if ("Relative_Risk" %in% colnames(stats)) "Relative_Risk" else
          stop("Stats is missing Relative_Risk_Mean (or Relative_Risk).")

## ---- build igraph ----
net <- graph_from_graphnel(as.graphNEL(avg_boot))
if (vcount(net) < 2 || ecount(net) == 0) stop("Network is too small or has no edges for visualization")

set.seed(seed)
layout_fr <- igraph::layout_with_fr(net)
coords <- data.frame(
  name = igraph::V(net)$name,
  x = layout_fr[,1],
  y = layout_fr[,2],
  stringsAsFactors = FALSE
)

visDat <- toVisNetworkData(net)
nodes <- visDat$nodes
edges <- visDat$edges

# filter to metadata-covered nodes
nodes <- nodes[nodes$id %in% valid_genes, , drop = FALSE]
edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, , drop = FALSE]

# color nodes by metadata group
nodes$group <- lookup[nodes$id]

palette1 <- RColorBrewer::brewer.pal(8, "Set3")
palette2 <- RColorBrewer::brewer.pal(8, "Set2")
palette3 <- RColorBrewer::brewer.pal(12, "Paired")
palette4 <- RColorBrewer::brewer.pal(8, "Dark2")
color_palette <- c(palette1, palette2, palette3, palette4)
unique_groups <- unique(nodes$group)
color_palette <- color_palette[seq_len(length(unique_groups))]
if (length(unique_groups) > length(color_palette)) {
  color_palette <- rep(color_palette, length.out = length(unique_groups))
}
group_color_lookup <- setNames(color_palette, unique_groups)
nodes$color <- group_color_lookup[nodes$group]

# attach coords to nodes
nodes <- nodes %>% left_join(coords, by = c("id" = "name"))

# edge attributes from stats
stats_pairs <- stats %>%
  transmute(
    Gene_lo = pmin(as.character(Gene_1), as.character(Gene_2)),
    Gene_hi = pmax(as.character(Gene_1), as.character(Gene_2)),
    width_metric = .data[[wcol]],
    rr_val = .data[[rr_col]]
  ) %>%
  mutate(pair_key = paste(Gene_lo, Gene_hi, sep = "||"))

pair_summ <- stats_pairs %>%
  group_by(pair_key) %>%
  summarize(
    width_val = suppressWarnings(max(width_metric, na.rm = TRUE)),
    rr_min    = suppressWarnings(min(rr_val,    na.rm = TRUE)),
    .groups   = "drop"
  ) %>%
  mutate(
    width_val = ifelse(is.finite(width_val) & width_val > 0, width_val, NA_real_),
    edge_color = ifelse(rr_min < rr_thr, "red", "black")
  )

edges <- edges %>%
  mutate(pair_key = paste(pmin(as.character(from), as.character(to)),
                          pmax(as.character(from), as.character(to)),
                          sep = "||")) %>%
  left_join(pair_summ, by = "pair_key")

edges$width <- scales::rescale(replace_na(edges$width_val, 0), to = c(1, 5))
edges$color <- ifelse(is.na(edges$edge_color), "black", edges$edge_color)

lnodes <- data.frame(color = color_palette,
                     label = unique_groups,
                     font.color = 'white')

# ---- HTML (visNetwork) ----
network <- visNetwork(nodes = nodes, edges = edges, width = '100%', height = 900) %>%
  visNodes(size = 20, color = list(highlight = 'yellow'), font = list(size = 25)) %>%
  visEdges(smooth = list(enabled = TRUE, type = 'diagonalCross', roundness = 0.1),
           physics = FALSE, color = "black") %>%
  visIgraphLayout(layout = "layout_with_fr", type = 'full') %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE, labelOnly = TRUE),
             nodesIdSelection = TRUE) %>%
  visLegend(addNodes = lnodes, width = 0.1, position = 'left',
            main = 'Variable Type', ncol = 2, useGroups = FALSE)

visSave(network, html_out, background = "#F5F4F4")

# ---- PDF (ggraph) ----
tbl <- as_tbl_graph(net)

node_data <- nodes %>% transmute(name = id, x, y, node_color = unname(color))
tbl <- tbl %>% activate(nodes) %>% left_join(node_data, by = "name")

edge_attrs <- edges %>% transmute(pair_key = paste(pmin(from, to), pmax(from, to), sep = "||"),
                                  width, color)

node_names <- tbl %>% activate(nodes) %>% pull(name)

tbl <- tbl %>%
  activate(edges) %>%
  mutate(pair_key = paste(pmin(node_names[from], node_names[to]),
                          pmax(node_names[from], node_names[to]), sep = "||")) %>%
  left_join(edge_attrs, by = "pair_key") %>%
  mutate(width = replace_na(width, 1),
         color = replace_na(color, "black"))

p <- ggraph(tbl, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(edge_width = width, edge_colour = color),
                 lineend = "round", show.legend = TRUE) +
  scale_edge_width(range = c(1, 3), guide = "none") +
  scale_edge_colour_identity(name = "Edge",
                             labels = c("black" = "RR >= 1", "red" = paste0("RR < ", rr_thr)),
                             breaks = c("black", "red"),
                             guide = "legend") +
  geom_node_point(aes(fill = node_color), shape = 21, size = 10, stroke = 0.3, colour = "grey30") +
  scale_fill_identity(name = "Variable Type",
                      breaks = unname(group_color_lookup),
                      labels = names(group_color_lookup),
                      guide = "legend") +
  geom_node_text(aes(label = name), size = 4.8, vjust = -1) +
  theme_void() +
  theme(legend.position = "left")

ggsave(pdf_out, p, width = 14, height = 10, units = "in")

cat(sprintf("\n\033[32mSaved %s and %s.\033[39m\n\n", html_out, pdf_out))