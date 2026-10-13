#!/usr/bin/env bash
# Historical name: this wrapper DELETES listed IDs; it does not extract them.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ $# != 3 && $# != 5 ]]; then
    echo 'Usage: extract_fastq.sh table.tsv input1.fq output1.fq [input2.fq output2.fq] (delete listed reads)' >&2
    exit 2
fi
TEMP_DIR="$(mktemp -d)"
trap 'rm -rf -- "$TEMP_DIR"' EXIT
awk 'NF {print $1}' "$1" >"$TEMP_DIR/names.txt"
if [[ $# == 3 ]]; then
    "${NJU_PERL:-perl}" "$SCRIPT_DIR/delete_fastq.pl" -n "$TEMP_DIR/names.txt" -i "$2" -o "$3"
else
    "${NJU_PERL:-perl}" "$SCRIPT_DIR/delete_fastq.pl" -n "$TEMP_DIR/names.txt" -i "$2" -o "$3" &
    PID1=$!
    "${NJU_PERL:-perl}" "$SCRIPT_DIR/delete_fastq.pl" -n "$TEMP_DIR/names.txt" -i "$4" -o "$5" &
    PID2=$!
    FAILED=0
    wait "$PID1" || FAILED=1
    wait "$PID2" || FAILED=1
    exit "$FAILED"
fi
