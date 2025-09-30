#!/bin/bash

#MVA

# Get the directory where the script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

Rscript "$DIR/MVA.r" "$@"