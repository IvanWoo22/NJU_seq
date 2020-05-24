#!/usr/bin/env bash
filepath=$(
  cd "$(dirname "${0}")" || exit
  pwd -P
)
perl "${filepath}"/uniquematch.pl "$1" \
  "$3"/mrna.unique.tmp \
  "$3"/mrna.multi.tmp
perl "${filepath}"/almostuniquematch.pl "$2" \
  "$3"/mrna.multi.tmp \
  "$3"/mrna.almost.tmp
cat "$3"/mrna.almost.tmp \
  "$3"/mrna.unique.tmp \
  >"$4"
rm "$3"/mrna.multi.tmp
rm "$3"/mrna.almost.tmp
rm "$3"/mrna.unique.tmp
