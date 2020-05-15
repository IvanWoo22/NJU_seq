#!/usr/bin/env bash
filepath=$(
  cd "$(dirname "${0}")" || exit
  pwd -P
)

if [[ $# == 5 ]]; then
  awk '{print $1}' "$1" \
  >name.txt
  perl "${filepath}"/delete_fastq.pl \
  -n name.txt \
  -i "$2" \
  -o "$3" &
  perl "${filepath}"/delete_fastq.pl \
  -n name.txt \
  -i "$4" \
  -o "$5" &
  wait
  rm name.txt
elif [[ $# == 3 ]]; then
  awk '{print $1}' "$1" \
  >name.txt
  perl "${filepath}"/delete_fastq.pl \
  -n name.txt \
  -i "$2" \
  -o "$3"
  rm name.txt
else
  echo "Improper Number of Files."
fi
