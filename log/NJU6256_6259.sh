for PREFIX in NJU62{56..59} NJU62{68..71}; do
	mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}"
	perl NJU_seq/quality_control/pe_consistency.pl \
		data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R2.fq.gz \
		temp/"${PREFIX}".fq.gz
done

time perl NJU_seq/quality_control/fastq_qc.pl \
	temp/HeLa_RF_NC.fq.gz \
	temp/HeLa_RF_1.fq.gz \
	temp/HeLa_RF_2.fq.gz \
	temp/HeLa_RF_3.fq.gz \
	output \
	HeLa_RF

for PREFIX in NJU62{56..59} NJU62{68..71}; do
	bsub -n 24 -o ${PREFIX}_rrna_alignment.log -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step1.sh ${PREFIX}"
done

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz | awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' | perl NJU_seq/rrna_analysis/matchquality_judge.pl | perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU62{56..59} NJU62{68..71}
'

for PREFIX in NJU62{56..59} NJU62{68..71}; do
	bsub -q fat_768 -n 80 -o ${PREFIX}_mrna_alignment.log -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step2.sh ${PREFIX}"
done

for PREFIX in NJU62{56..59} NJU62{68..71}; do
	bsub -n 24 -o ${PREFIX}_filter.log -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step3.sh ${PREFIX} 16"
done

bsub -n 24 -J "almostunique" '
parallel --keep-order --xapply -j 6 '\''
time pigz -dc temp/{}/mrna_basic.dedup.tmp.gz > temp/{}/mrna_basic.dedup.tmp
time bash NJU_seq/mrna_analysis/almostunique.sh \
  temp/{}/mrna_basic.dedup.tmp \
  data/{}/R1.mrna.fq.gz \
  temp/{} \
  temp/{}/mrna_basic.almostunique.tmp
rm temp/{}/mrna_basic.dedup.tmp
'\'' ::: NBG00{01..12}
'

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 12 '\''
time perl NJU_seq/mrna_analysis/count.pl \
temp/{}/mrna_basic.almostunique.tmp \
>temp/{}/mrna_basic.count.tmp
'\'' ::: NBG00{01..12}
'

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 12 '\''
time pigz -dcf data/gmax.gff3.gz |
awk '\''\'\'''\'' $3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9} '\''\'\'''\'' |
perl NJU_seq/log/Gmax_scripts/merge.pl --refstr "Parent=" --geneid "Glyma." --transid "Glyma." -i temp/{}/mrna_basic.count.tmp -o output/{}/mrna_basic.tsv
'\'' ::: NBG00{01..12}
'

parallel --keep-order -j 4 '
echo {} >>output/{}/mrna_basic.cov
bash NJU_seq/presentation/seq_depth.sh \
temp/{}/mrna_basic.almostunique.tmp \
output/{}/mrna_basic.tsv \
>>output/{}/mrna_basic.cov
' ::: NBG00{01..12}

parallel -j 3 "
perl score_neo.pl \\
output/NBG0001/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NBG00{02..04}

parallel -j 3 "
perl score_neo.pl \\
output/NBG0005/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NBG00{06..08}

parallel -j 3 "
perl score_neo.pl \\
output/NBG0009/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NBG00{10..12}

perl NJU_seq/mrna_analysis/extract_point.pl \
	output/NBG0002/mrna_basic_scored.tsv \
	output/NBG0003/mrna_basic_scored.tsv \
	output/NBG0004/mrna_basic_scored.tsv \
	1000 1 >output/NBG_group1_mrna_scored_1000p.tsv

perl NJU_seq/mrna_analysis/extract_point.pl \
	output/NBG0006/mrna_basic_scored.tsv \
	output/NBG0007/mrna_basic_scored.tsv \
	output/NBG0008/mrna_basic_scored.tsv \
	1000 1 >output/NBG_group2_mrna_scored_1000p.tsv

perl NJU_seq/mrna_analysis/extract_point.pl \
	output/NBG0010/mrna_basic_scored.tsv \
	output/NBG0011/mrna_basic_scored.tsv \
	output/NBG0012/mrna_basic_scored.tsv \
	1000 1 >output/NBG_group3_mrna_scored_1000p.tsv

for PREFIX in group1 group2 group3; do
	perl NJU_seq/presentation/signature_count.pl \
		output/NBG_${PREFIX}_mrna_scored_1000p.tsv \
		output/NBG_${PREFIX}_mrna_signature.pdf
done
