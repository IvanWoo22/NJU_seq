#!/usr/bin/env bash
# Mean records per listed site; not per-base sequencing coverage.
set -euo pipefail
[[ $# == 2 ]] || { echo 'Usage: seq_depth.sh read_records.tsv site_records.tsv' >&2; exit 2; }
RECORDS=$(awk 'END {print NR+0}' "$1")
SITES=$(awk 'END {print NR+0}' "$2")
awk -v r="$RECORDS" -v s="$SITES" 'BEGIN {
    printf "Stop times:\t%d\nPositions:\t%d\n",r,s;
    if(s==0) print "Mean records per listed site:\tNA";
    else printf "Mean records per listed site:\t%.2f\n",r/s;
}'
