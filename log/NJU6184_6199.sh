for PREFIX in NJU61{84..99} NJU62{20..35}; do
	mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}"
	perl NJU_seq/quality_control/pe_consistency.pl \
		data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R2.fq.gz \
		temp/"${PREFIX}".fq.gz
done

cat NJU_seq/data/ath_rrna/* >data/ath_rrna.fa

pigz -dc data/ath_ncrna.fa.gz \
	| perl NJU_seq/tool/fetch_fasta.pl \
		--stdin -s 'transcript_biotype:rRNA' \
		>>data/ath_rrna.fa
bowtie2-build data/ath_rrna.fa index/ath_rrna
rm data/ath_rrna.fa

for PREFIX in NJU61{84..99} NJU62{20..35}; do
	bsub -n 24 -J "${PREFIX}" -o ${PREFIX}.out "bash NJU_seq/log/Ath_scripts/step1.sh ${PREFIX}"
done

perl test.pl output.36412* | Rscript test.R test.pdf

bsub -n 24 -q largemem -e NJU6184_count.out -J "count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz | awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' | perl NJU_seq/rrna_analysis/matchquality_judge.pl | perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU61{84..99} NJU62{20..35}
'

for RNA in 25s 18s 5-8s; do
	parallel --keep-order --xapply -j 10 "
  perl NJU_seq/rrna_analysis/readend_count.pl NJU_seq/data/ath_rrna/${RNA}.fa temp/{}/rrna.out.tmp ${RNA} >output/{}/rrna_${RNA}.tsv
" ::: NJU61{84..99} NJU62{20..35}
done

for RNA in 25s 18s 5-8s; do
	parallel --keep-order --xapply -j 8 "
  perl NJU_seq/rrna_analysis/score.pl \\
    output/NJU{1}/rrna_${RNA}.tsv \\
    output/NJU{2}/rrna_${RNA}.tsv \\
    output/NJU{3}/rrna_${RNA}.tsv \\
    output/NJU{4}/rrna_${RNA}.tsv \\
      >output/{5}_rrna_${RNA}_scored.tsv
" ::: $(seq 6184 4 6199) $(seq 6220 4 6235) ::: $(seq 6185 4 6199) $(seq 6221 4 6235) ::: $(seq 6186 4 6199) $(seq 6222 4 6235) ::: $(seq 6187 4 6199) $(seq 6223 4 6235) ::: Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF
done

for PREFIX in Ath_sixl Ath_bolt Ath_bloom Ath_fruit Ath_root Ath_stem Ath_flower Ath_leaf; do
	bash NJU_seq/presentation/point_venn.sh \
		Sample1 output/${PREFIX}_rrna_18s_scored.tsv 14 \
		Sample2 output/${PREFIX}_rrna_18s_scored.tsv 15 \
		Sample3 output/${PREFIX}_rrna_18s_scored.tsv 16 \
		output/${PREFIX}_rrna_18s_venn.png 40
done

for PREFIX in NJU61{84..99} NJU62{20..35}; do
	bsub -q largemem -n 24 -e ${PREFIX}_mrna_alignment.log -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step2.sh ${PREFIX}"
done

for PREFIX in NJU61{84..99} NJU62{20..35}; do
	bsub -q largemem -n 24 -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step3.sh ${PREFIX} 16"
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
'\'' ::: NJU61{84..99} NJU62{20..35}
'

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 12 '\''
time perl NJU_seq/mrna_analysis/count.pl \
temp/{}/mrna_basic.almostunique.tmp \
>temp/{}/mrna_basic.count.tmp
'\'' ::: NJU61{84..99} NJU62{20..35}
'

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 8 '\''
time pigz -dcf data/ath.gff3.gz |
awk '\''\'\'''\'' $3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9} '\''\'\'''\'' |
perl NJU_seq/mrna_analysis/merge.pl --refstr "Parent=transcript:" --geneid "AT" --transid "AT" -i temp/{}/mrna_basic.count.tmp -o output/{}/mrna_basic.tsv
'\'' ::: NJU61{84..99} NJU62{20..35}
'

parallel --keep-order -j 4 '
echo {} >>output/{}/mrna_basic.cov
bash NJU_seq/presentation/seq_depth.sh \
temp/{}/mrna_basic.almostunique.tmp \
output/{}/mrna_basic.tsv \
>>output/{}/mrna_basic.cov
' ::: NJU61{84..99} NJU62{20..35}

parallel -j 3 "
perl score_neo.pl \\
output/NJU6184/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NJU61{85..87}

parallel -j 3 "
perl score_neo.pl \\
output/NJU6188/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NJU61{89..91}

parallel -j 3 "
perl score_neo.pl \\
output/NJU6192/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NJU61{93..95}

parallel -j 3 "
perl score_neo.pl \\
output/NJU6196/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NJU61{97..99}

parallel -j 3 "
> perl score_neo.pl \\
> output/NJU6220/mrna_basic.tsv \\
> output/{}/mrna_basic.tsv |
> sort -t $'\t' -nrk 12,12 \\
> >output/{}/mrna_basic_scored.tsv
> " ::: NJU62{21..23}

parallel -j 3 "
perl score_neo.pl \\
output/NJU6224/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NJU62{25..27}

parallel -j 3 "
perl score_neo.pl \\
output/NJU6228/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NJU62{29..31}

parallel -j 3 "
perl score_neo.pl \\
output/NJU6232/mrna_basic.tsv \\
output/{}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_basic_scored.tsv
" ::: NJU62{33..35}

parallel --keep-order --xapply -j 8 "
perl NJU_seq/mrna_analysis/extract_point.pl \
  output/NJU{1}/mrna_basic_scored.tsv \
  output/NJU{2}/mrna_basic_scored.tsv \
  output/NJU{3}/mrna_basic_scored.tsv \
  1000 1 >output/{4}_mrna_scored_1000p.tsv
  " ::: $(seq 6185 4 6199) $(seq 6221 4 6235) ::: $(seq 6186 4 6199) $(seq 6222 4 6235) ::: $(seq 6187 4 6199) $(seq 6223 4 6235) ::: Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF
