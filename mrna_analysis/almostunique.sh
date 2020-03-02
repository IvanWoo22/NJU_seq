filepath=$(
  cd "$(dirname "${0}")" || exit
  pwd -P
)
perl "${filepath}"/uniquematch.pl "$1" \
  "$3"/mrna.uniqe.tmp \
  "$3"/mrna.multi.tmp
perl ~/2OMG/mrna_analysis/almostuniquematch.pl "$2" \
  "$3"/mrna.multi.tmp \
  "$3"/mrna.almost.tmp
cat "$3"/mrna.almost.tmp \
  "$3"/mrna.uniqe.tmp \
  >"$4"
