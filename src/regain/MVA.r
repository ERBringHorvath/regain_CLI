#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 University of Utah

suppressPackageStartupMessages({ library(optparse) })

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"),
        dplyr.summarise.inform = FALSE)

# ---------- CLI ----------
opt_list <- list(
  make_option(c("-i","--input"), type="character", help="Input matrix CSV (first column = feature/gene names)"),
  make_option(c("-m","--method"), type="character", default="euclidean",
              help="Distance: manhattan, euclidean, canberra, clark, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao, chord, hellinger, aitchison, mahalanobis [default: %default]"),
  make_option(c("--k"), type="integer", default=1, help="Number of clusters for k-means (set 0 to auto-select 2..10 by silhouette) [default: %default]"),
  make_option(c("-C","--confidence"), type="double", default=0.95, help="Ellipse confidence level [default: %default]"),
  make_option(c("--seed"), type="integer", default=42, help="Random seed [default: %default]"),
  make_option(c("--label"), type="character", default="auto", help="Sample labels: none | auto | all [default: %default]"),
  make_option(c("--point-size"), type="double", default=3.5, help="Point size [default: %default]"),
  make_option(c("--alpha"), type="double", default=0.75, help="Point alpha [default: %default]"),
  make_option(c("--no-ellipses"), action="store_true", default=FALSE, help="Disable cluster ellipses"),
  make_option(c("--pseudocount"), type="double", default=1e-6, help="Pseudocount for CLR/aitchison [default: %default]"),
  make_option(c("--save-dist"), action="store_true", default=FALSE, help="Export distance matrix as CSV"),
  make_option(c("--png-out"), type="character", default="MVA.png", help="PNG output file [default: %default]"),
  make_option(c("--pdf-out"), type="character", default="MVA.pdf", help="PDF output file [default: %default]"),
  make_option(c("--coords-out"), type="character", default="MVA_coordinates.csv", help="PCoA coordinates CSV [default: %default]"),
  make_option(c("--dist-out"), type="character", default="MVA_distance.csv", help="Distance matrix CSV (with --save-dist) [default: %default]"),
  make_option(c("--pcoa-correction"), type="character", default="auto", help="PCoA negative-eigenvalue correction: auto | none | lingoes | cailliez [default: %default]")
)

parser <- OptionParser(option_list = opt_list,
                       description = "ReGAIN — Multivariate analysis (distance + PCoA + k-means + ellipses)")
opt <- parse_args(parser)

save_dist   <- isTRUE(opt$`save-dist`)   || isTRUE(opt$save_dist)
no_ellipses <- isTRUE(opt$`no-ellipses`) || isTRUE(opt$no_ellipses)

req <- c("input")
missing <- req[ vapply(req, function(k){ v <- opt[[k]]; is.null(v) || !nzchar(trimws(as.character(v))) }, logical(1)) ]
if (length(missing)) { print_help(parser); stop(paste("Missing required:", paste(missing, collapse=", ")), call.=FALSE) }

# ---------- packages ----------
cran_pkgs <- c("dplyr","tibble","vegan","ggplot2","ggrepel","tidyr","ellipse","ape","scales","cluster")
need <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)
# For CLR transform
if (!requireNamespace("compositions", quietly = TRUE)) install.packages("compositions")

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(vegan); library(ggplot2); library(ggrepel)
  library(tidyr); library(ellipse); library(ape); library(scales); library(cluster)
})

set.seed(opt$seed)

## ---- helpers ----
vegdist_methods <- c('manhattan','euclidean','canberra','clark','bray','kulczynski',
                     'jaccard','gower','altGower','morisita','horn','mountford',
                     'raup','binomial','chao','cao')

transform_methods <- c("chord","hellinger","aitchison","mahalanobis")

if (! (tolower(opt$method) %in% c(vegdist_methods, transform_methods)) ) {
  stop(sprintf("Invalid method '%s'. Choose one of: %s",
               opt$method, paste(c(vegdist_methods, transform_methods), collapse=", ")), call.=FALSE)
}

read_matrix <- function(path) {
  dat <- read.csv(path, check.names = FALSE)
  if (ncol(dat) < 2) stop("Matrix must have >= 2 columns (first column = row names).")
  rn <- as.character(dat[[1]])
  mat <- as.data.frame(dat[,-1, drop=FALSE], check.names = FALSE)
  rownames(mat) <- make.unique(rn)
  #We want samples as rows. In ReGAIN matrices, columns are samples -> transpose.
  t(mat) %>% as.data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column(var="sample")
}

clr_transform <- function(X, pseudocount=1e-6) {
  X <- as.matrix(X)
  X <- X + pseudocount
  if (any(X <= 0)) stop("CLR requires strictly positive values after pseudocount.")
  comps <- compositions::clr(X)
  as.data.frame(comps, check.names = FALSE)
}

mahal_dist <- function(X) {
  X <- as.matrix(X)
  S <- stats::cov(X)
  #Guard for singular covariance
  S_inv <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  n <- nrow(X)
  D <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- X[i,]
    for (j in i:n) {
      xj <- X[j,]
      d <- sqrt( t(xi - xj) %*% S_inv %*% (xi - xj) )
      D[i,j] <- d; D[j,i] <- d
    }
  }
  stats::as.dist(D)
}

#Auto-select k via silhouette on 2D PCoA coords
auto_k <- function(coords2d, k_min=2, k_max=10) {
  n <- nrow(coords2d)
  k_max <- max(k_min, min(k_max, n-1))
  best_k <- NA_integer_; best_s <- -Inf
  for (k in k_min:k_max) {
    km <- kmeans(coords2d, centers = k, nstart = 20)
    # skip k if any cluster too small
    if (min(table(km$cluster)) < 3) next
    sil <- cluster::silhouette(km$cluster, dist(coords2d))
    s <- mean(sil[, "sil_width"])
    if (is.finite(s) && s > best_s) { best_s <- s; best_k <- k }
  }
  if (is.na(best_k)) list(k = 1L, sil = NA_real_) else list(k = best_k, sil = best_s)
}

cat("Reading matrix: ", opt$input, "\n")
df <- read_matrix(opt$input)        # rows = samples, cols (2..p) = features
stopifnot(nrow(df) >= 2, ncol(df) >= 3)

X <- df[,-1, drop=FALSE]            # numeric block
# Coerce to numeric safely
X[] <- lapply(X, function(col) {
  if (is.logical(col)) as.numeric(col)
  else if (is.factor(col)) as.numeric(as.character(col))
  else as.numeric(col)
})
if (any(is.na(as.matrix(X)))) {
  warning("NA values detected; replacing with 0.")
  X[is.na(X)] <- 0
}

method <- tolower(opt$method)

dist_obj <- NULL
if (method %in% vegdist_methods) {
  dist_obj <- vegdist(X, method = method)
} else if (method == "chord") {
  Xn <- vegan::decostand(X, method = "normalize")   # chord transform
  dist_obj <- dist(Xn)                              # Euclidean
} else if (method == "hellinger") {
  Xh <- vegan::decostand(X, method = "hellinger")
  dist_obj <- dist(Xh)
} else if (method == "aitchison") {
  Xc <- clr_transform(X, pseudocount = opt$pseudocount)
  dist_obj <- dist(Xc)
} else if (method == "mahalanobis") {
  dist_obj <- mahal_dist(X)
} else {
  stop(sprintf("Unsupported method path for '%s'", method))
}

if (save_dist) {
  cat("Writing distance matrix CSV: ", opt$`dist-out`, "\n")
  D <- as.matrix(dist_obj)
  rownames(D) <- df$sample; colnames(D) <- df$sample
  write.csv(data.frame(sample = rownames(D), D, check.names=FALSE),
            opt$`dist-out`, row.names = FALSE)
}

# ---------- ordination (PCoA with correction) ----------
corr <- tolower(opt$`pcoa-correction`)
if (!corr %in% c("auto","none","lingoes","cailliez")) stop("pcoa-correction must be auto|none|lingoes|cailliez")

needs_correction <- function(p) {
  any(p$values$Eigenvalues < 0)
}

if (corr == "none") {
  pco <- ape::pcoa(dist_obj)
} else if (corr == "lingoes") {
  pco <- ape::pcoa(dist_obj, correction = "lingoes")
} else if (corr == "cailliez") {
  pco <- ape::pcoa(dist_obj, correction = "cailliez")
} else { # auto
  p0 <- ape::pcoa(dist_obj)
  if (needs_correction(p0)) {
    pco <- ape::pcoa(dist_obj, correction = "cailliez")
  } else {
    pco <- p0
  }
}

#Extract first 2 axes
coords <- as.data.frame(pco$vectors[,1:2, drop=FALSE])
colnames(coords) <- c("Axis1","Axis2")
coords <- tibble::add_column(coords, sample = df$sample, .before = 1)

# % variance explained
pct <- round(100 * pco$values$Relative_eig[1:2], 2)
if (any(is.na(pct))) pct <- c(NA, NA)

# ---------- clustering ----------
is_pd_cov <- function(df2) {
  if (nrow(df2) < 3) return(FALSE)
  if (nrow(unique(df2)) < 3) return(FALSE)
  S <- stats::cov(df2)
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  all(is.finite(ev)) && min(ev) > 1e-12
}

km_k <- opt$k
if (km_k == 0) {
  ak <- auto_k(coords[,c("Axis1","Axis2")], k_min = 2, k_max = 10)
  km_k <- ak$k
  cat(sprintf("Auto-selected k = %d (mean silhouette = %.3f)\n", ak$k, ak$sil))
} else if (km_k < 1) {
  stop("--k must be >= 0; use 0 for auto.")
}

clusters <- if (km_k >= 1) {
  kmeans(coords[,c("Axis1","Axis2")], centers = km_k, nstart = 20)$cluster
} else {
  rep(1L, nrow(coords))
}
coords$cluster <- factor(clusters)

# Determine which clusters can draw ellipses (≥3 distinct points & PD covariance)
split_by_cluster <- split(coords[, c("Axis1","Axis2")], coords$cluster)
ok_names <- names(Filter(isTRUE, lapply(split_by_cluster, is_pd_cov)))

# ---------- plotting ----------
lab_mode <- tolower(opt$label)
do_labels <- match.arg(lab_mode, c("none","auto","all"))

p <- ggplot(coords, aes(x = Axis1, y = Axis2)) +
  theme_bw() +
  geom_point(aes(color = cluster), alpha = opt$alpha, size = opt$`point-size`, show.legend = FALSE) +
  labs(x = ifelse(!is.na(pct[1]), sprintf("PCo1 (%.2f%%)", pct[1]), "PCo1"),
       y = ifelse(!is.na(pct[2]), sprintf("PCo2 (%.2f%%)", pct[2]), "PCo2")) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

# labels
if (do_labels == "all") {
  p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 3.5, max.overlaps = Inf, box.padding = 0.4)
} else if (do_labels == "auto") {
  p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 3.2, max.overlaps = 50, box.padding = 0.3)
}

# ellipses (only for drawable clusters; skip if none or --no-ellipses)
if (!no_ellipses && length(ok_names)) {
  coords_ok <- subset(coords, cluster %in% ok_names)
  p <- p +
    stat_ellipse(aes(group = cluster, fill = cluster),
                 data = coords_ok, type = "norm", level = opt$confidence,
                 geom = "polygon", alpha = 0.08, show.legend = FALSE) +
    stat_ellipse(aes(group = cluster, color = cluster),
                 data = coords_ok, type = "norm", level = opt$confidence,
                 linetype = 4, linewidth = 0.6,
                 geom = "path", show.legend = FALSE)
} else if (!no_ellipses) {
  message("Skipping ellipses: no clusters have PD covariance with ≥3 distinct points.")
}

# ---------- outputs ----------
ggsave(filename = opt$`png-out`, plot = p, width = 10, height = 10, dpi = 300)
ggsave(filename = opt$`pdf-out`, plot = p, width = 10, height = 10)
write.csv(coords, opt$`coords-out`, row.names = FALSE)

cat(sprintf("\n\033[32mSaved %s, %s, and %s.\033[39m\n",
            opt$`png-out`, opt$`pdf-out`, opt$`coords-out`))
if (save_dist) cat(sprintf("\033[32mSaved %s.\033[39m\n", opt$`dist-out`))