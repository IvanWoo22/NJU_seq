#!/usr/bin/env bash
filepath=$(
  cd "$(dirname "${0}")" || exit
  pwd -P
)
perl "${filepath}"/uniquematch.pl "$1" \
  "$3"/unique.tmp \
  "$3"/multi.tmp
perl "${filepath}"/almostuniquematch.pl "$2" \
  "$3"/multi.tmp \
  "$3"/almost.tmp
cat "$3"/almost.tmp \
  "$3"/unique.tmp \
  >"$4"
rm "$3"/multi.tmp
rm "$3"/almost.tmp
rm "$3"/unique.tmp
