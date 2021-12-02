#!/usr/bin/env bash
set -e

FP=$(
  cd "$(dirname "${0}")" || exit
  pwd -P
)

FH=$1
COL=$2

cut -f "${COL}" "${FH}" >"${FH}".seqlogo.tmp
Rscript "${FP}"/seq_logo.R "${FH}".seqlogo.tmp "${3}"
