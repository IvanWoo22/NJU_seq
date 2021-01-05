for PREFIX in NBG00{01..12}; do
bsub -n 24 -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step1.sh ${PREFIX}"
done

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 10 '\''
mkdir temp/{}
time pigz -dcf output/{}/rrna.raw.sam.gz | awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
'\''\'\'''\'' | perl NJU_seq/rrna_analysis/matchquality_judge.pl | perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz \
data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' NBG00{01..12}
'

for PREFIX in NBG00{01..12}; do
bsub -q fat_768 -n 80 -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step2.sh ${PREFIX}"
done