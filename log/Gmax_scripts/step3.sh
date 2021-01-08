PREFIX=$1
THREAD=$2

time pigz -dc output/"${PREFIX}"/mrna_basic.raw.sam.gz |
  parallel --tmpdir /scratch/wangq/wyf/. \
    --pipe --block 2G-1 --no-run-if-empty \
    --linebuffer --keep-order -j "${THREAD}" \
    'awk '\''$6!="*"&&$7=="=" {print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' | perl NJU_seq/mrna_analysis/multimatch_judge.pl' |
  parallel --tmpdir /scratch/wangq/wyf/. \
    --pipe --block 2G-1 --no-run-if-empty \
    --linebuffer --keep-order -j "${THREAD}" \
    'perl NJU_seq/mrna_analysis/multimatch_judge.pl' |
  parallel --tmpdir /scratch/wangq/wyf/. \
    --pipe --block 2G-1 --no-run-if-empty \
    --linebuffer --keep-order -j "${THREAD}" \
    'perl NJU_seq/mrna_analysis/dedup.pl --refstr "Parent=" --transid "ENST" --info data/hsa_exon.info' |
  parallel --tmpdir /scratch/wangq/wyf/. \
    --pipe --block 2G-1 --no-run-if-empty \
    --linebuffer --keep-order -j "${THREAD}" \
    'perl NJU_seq/mrna_analysis/dedup.pl --refstr "Parent=" --transid "ENST" --info data/hsa_exon.info' |
  pigz >temp/"${PREFIX}"/mrna_basic.dedup.tmp.gz
