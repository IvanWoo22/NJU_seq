PREFIX=$1

time bowtie2 -p 72 -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 --score-min C,0,0 \
  --xeq -x index/gmax_protein_coding \
  -1 data/"${PREFIX}"/R1.mrna.fq.gz -2 data/"${PREFIX}"/R2.mrna.fq.gz |
  pigz >output/"${PREFIX}"/mrna_basic.raw.sam.gz
