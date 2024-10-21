#!/bin/bash
# bnCPQuery.sh

# Get the directory where the script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Run the R script using its full path
Rscript "$DIR/bayesCPquery.r" "$@"