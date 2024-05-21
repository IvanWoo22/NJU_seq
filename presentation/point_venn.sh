#!/usr/bin/env bash

set -e

filepath=$(
	cd "$(dirname "${0}")" || exit
	pwd -P
)

sort -t $'\t' -nrk "$3","$3" "$2" | head -n "${11}" | awk -F $'\t' '{print $1 $3 $2}' >templist1.tmp
sort -t $'\t' -nrk "$6","$6" "$5" | head -n "${11}" | awk -F $'\t' '{print $1 $3 $2}' >templist2.tmp
sort -t $'\t' -nrk "$9","$9" "$8" | head -n "${11}" | awk -F $'\t' '{print $1 $3 $2}' >templist3.tmp

Rscript "${filepath}"/point_venn.R "$1" templist1.tmp "$4" templist2.tmp "$7" templist3.tmp "${10}"
rm templist1.tmp templist2.tmp templist3.tmp "${10}"*.log
