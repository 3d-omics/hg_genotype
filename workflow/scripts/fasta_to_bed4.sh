#!/usr/bin/env bash
set -euo pipefail

seqtk seq "$1" \
| paste - - \
| awk '{print $1 "\t" 1 "\t" length($2) "\t"  $1}' \
| tr -d ">"
