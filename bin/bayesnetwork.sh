#!/bin/bash
##bayesnetwork.sh

# Get directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

##Run Rscript
Rscript "$DIR/bayesnetwork.r" "$@"