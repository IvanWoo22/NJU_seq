#!/usr/bin/env bash
ALL_STOP=$(wc -l <"$1")
ALL_POS=$(wc -l <"$2")
COV=$(echo "scale=2; ${ALL_STOP}/${ALL_POS}" | bc)
echo -e "All stop times:\t${ALL_STOP}\nAll positions:\t${ALL_POS}\nCoverage:\t${COV}\n"
