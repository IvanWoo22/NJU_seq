#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
[[ $# == 11 && "${11}" =~ ^[1-9][0-9]*$ ]] || { echo 'Require three label/file/score-column triples, output path and positive top N' >&2; exit 2; }
TEMP_DIR="$(mktemp -d)"
trap 'rm -rf -- "$TEMP_DIR"' EXIT
for I in 1 2 3; do
    case "$I" in
      1) INPUT=$2; COL=$3;;
      2) INPUT=$5; COL=$6;;
      3) INPUT=$8; COL=$9;;
    esac
    [[ "$COL" =~ ^[1-9][0-9]*$ ]] || exit 2
    # awk reads all sorted input: no SIGPIPE from head under pipefail.
    awk -F '\t' '$0 !~ /^#/ && $2 ~ /^[0-9]+$/ {print}' "$INPUT" \
        | LC_ALL=C sort -t $'\t' -nrk "$COL,$COL" \
        | awk -F '\t' -v n="${11}" 'NR<=n {print $1 ":" $3 ":" $2}' >"$TEMP_DIR/$I"
done
Rscript "$SCRIPT_DIR/point_venn.R" "$1" "$TEMP_DIR/1" "$4" "$TEMP_DIR/2" "$7" "$TEMP_DIR/3" "${10}"
