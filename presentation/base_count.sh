#!/usr/bin/env bash
ALL=$(wc -l <"$1")
A=$(awk '$4=="A"' "$1" | wc -l)
G=$(awk '$4=="G"' "$1" | wc -l)
C=$(awk '$4=="C"' "$1" | wc -l)
U=$(awk '$4=="T"' "$1" | wc -l)
echo -e "All base:\t${ALL}\nA:\t${A}\nG:\t${G}\nC:\t${C}\nU:\t${U}\n"
