#!/usr/bin/env bash

# Legacy policy: singly aligned reads bypass R1 sequence verification.
set -euo pipefail
[[ $# == 4 ]] || { echo 'Usage: almostunique.sh alignments.tsv R1.fq.gz scratch_dir output.tsv' >&2; exit 2; }
TEMP_DIR="$(mktemp -d "$3/almostunique.XXXXXX")"
trap 'rm -rf -- "$TEMP_DIR"' EXIT

filepath=$(
	cd "$(dirname "${0}")" || exit
	pwd -P
)
"${NJU_PERL:-perl}" "${filepath}"/uniquematch.pl "$1" \
	"$TEMP_DIR"/unique.tmp \
	"$TEMP_DIR"/multi.tmp
"${NJU_PERL:-perl}" "${filepath}"/almostuniquematch.pl "$2" \
	"$TEMP_DIR"/multi.tmp \
	"$TEMP_DIR"/almost.tmp
cat "$TEMP_DIR"/almost.tmp \
	"$TEMP_DIR"/unique.tmp \
	>"$4"
rm "$TEMP_DIR"/multi.tmp
rm "$TEMP_DIR"/almost.tmp
rm "$TEMP_DIR"/unique.tmp
