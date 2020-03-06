#!/usr/bin/env bash
ALL=$(wc -l <"$1")
echo -e "All\t${ALL}"
A=$(awk '$5=="A"' "$1" | wc -l)
G=$(awk '$5=="G"' "$1" | wc -l)
C=$(awk '$5=="C"' "$1" | wc -l)
T=$(awk '$5=="T"' "$1" | wc -l)
echo -e "A\t${A}\nG\t${G}\nC\t${C}\nT\t${T}"

for sig in A G C T; do
  for lef in A G C T; do
    for rgt in A G C T; do
      echo -e "${lef}${sig}${rgt}\t$(awk -va=${sig} -vl=${lef} -vr=${rgt} '$4==l&&$5==a&&$6==r' "$1" | wc -l)"
    done
  done
done
