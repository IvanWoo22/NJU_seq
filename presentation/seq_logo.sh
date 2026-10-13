#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
[[ $# == 3 && "$2" =~ ^[1-9][0-9]*$ ]] || { echo 'Usage: seq_logo.sh table.tsv sequence_column output.pdf' >&2; exit 2; }
TEMP_DIR="$(mktemp -d)"
trap 'rm -rf -- "$TEMP_DIR"' EXIT
awk -F '\t' -v c="$2" '$0 !~ /^#/ && $1!="Chromosome" && $1!="Position" {print $c}' "$1" >"$TEMP_DIR/sequences"
Rscript "$SCRIPT_DIR/seq_logo.R" "$TEMP_DIR/sequences" "$3"
