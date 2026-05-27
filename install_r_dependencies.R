#!/usr/bin/env Rscript

# install_regain_r_deps.R
#
# Install R dependencies required for standalone ReGAIN installations.
#
# Usage:
#   Rscript install_regain_r_deps.R
#
# Notes:
#   - This script installs R packages only.
#   - External command-line tools such as AMRFinderPlus, BLAST, fastANI,
#     Mash, Python, pandas, biopython, and tqdm must be installed separately.
#   - ReGAIN currently expects R >= 4.4.

message("\n=== ReGAIN R dependency installer ===\n")

args <- commandArgs(trailingOnly = TRUE)

install_lib <- NULL

if (length(args) > 0) {
  lib_arg <- grep("^--lib=", args, value = TRUE)

  if (length(lib_arg) == 1) {
    install_lib <- sub("^--lib=", "", lib_arg)
  } else if ("--lib" %in% args) {
    idx <- match("--lib", args)
    if (idx == length(args)) {
      stop("Argument --lib requires a path.", call. = FALSE)
    }
    install_lib <- args[idx + 1]
  }
}

if (!is.null(install_lib)) {
  install_lib <- normalizePath(install_lib, mustWork = FALSE)
  dir.create(install_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(install_lib, .libPaths()))
  message("Installing packages into: ", install_lib)
}

message("Active R library paths:")
message(paste("  -", .libPaths(), collapse = "\n"))

required_r_version <- "4.4.0"

if (getRversion() < required_r_version) {
  stop(
    "ReGAIN requires R >= ", required_r_version,
    ". You are running R ", as.character(getRversion()), ".",
    call. = FALSE
  )
}

# Use a CRAN mirror if one has not already been set.
repos <- getOption("repos")
if (is.null(repos) || is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

cran_packages <- c(
  "optparse",
  "dplyr",
  "tidyr",
  "tibble",
  "igraph",
  "visNetwork",
  "reshape2",
  "scales",
  "ggplot2",
  "ggrepel",
  "tidygraph",
  "ggraph",
  "ellipse",
  "vegan",
  "ape",
  "compositions",
  "cluster",
  "pbapply",
  "foreach",
  "doParallel",
  "progressr",
  "RColorBrewer",
  "bnlearn",
  "gRain"
)

bioc_packages <- c(
  "graph"
)

install_if_missing <- function(pkgs, installer = c("cran", "bioc")) {
  installer <- match.arg(installer)
  
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_pkgs) == 0) {
    message("All ", installer, " packages are already installed.")
    return(invisible(TRUE))
  }
  
  message("\nInstalling missing ", installer, " packages:")
  message(paste("  -", missing_pkgs, collapse = "\n"))
  
  if (installer == "cran") {
    install.packages(missing_pkgs, dependencies = TRUE, lib = .libPaths()[1])
  } else {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      message("\nInstalling BiocManager...")
      install.packages("BiocManager")
    }
    
    BiocManager::install(
	  missing_pkgs,
	  ask = FALSE,
	  update = FALSE,
	  lib = .libPaths()[1]
	)
  }
  
  invisible(TRUE)
}

# Install CRAN deps
install_if_missing(cran_packages, installer = "cran")

# Install Bioconductor dependencies
install_if_missing(bioc_packages, installer = "bioc")

# Final check
all_packages <- c(cran_packages, bioc_packages)

still_missing <- all_packages[
  !vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)
]

message("\n=== Installation summary ===")

if (length(still_missing) == 0) {
  message("Success: all ReGAIN R dependencies are installed.\n")
} else {
  message("Warning: the following packages could not be loaded after installation:")
  message(paste("  -", still_missing, collapse = "\n"))
  message(
    "\nSome packages may require system-level build libraries, especially on Linux.\n",
    "Please inspect the installation output above for compiler or missing-library errors.\n"
  )
  quit(status = 1)
}

