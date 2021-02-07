for PREFIX in NJU62{00..15} NJU62{36..51}; do
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

cat NJU_seq/data/ath_rrna/* >data/ath_rrna.fa

pigz -dc data/ath_ncrna.fa.gz |
  perl NJU_seq/tool/fetch_fasta.pl \
    --stdin -s 'transcript_biotype:rRNA' \
    >>data/ath_rrna.fa
bowtie2-build data/ath_rrna.fa index/ath_rrna
rm data/ath_rrna.fa

for PREFIX in NJU62{00..15} NJU62{36..51}; do
  bsub -n 24 -J "${PREFIX}" -e ${PREFIX}.out "bash NJU_seq/log/Ath_scripts/step1.sh ${PREFIX}"
done

perl test.pl output.36412* | Rscript test.R test.pdf

bsub -n 24 -q largemem -J "count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz | awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' | perl NJU_seq/rrna_analysis/matchquality_judge.pl | perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU62{00..15} NJU62{36..51}
'

for RNA in 25s 18s 5-8s; do
  parallel --keep-order --xapply -j 10 "
  perl NJU_seq/rrna_analysis/readend_count.pl NJU_seq/data/ath_rrna/${RNA}.fa temp/{}/rrna.out.tmp ${RNA} >output/{}/rrna_${RNA}.tsv
" ::: NJU62{00..15} NJU62{36..51}
done

for RNA in 25s 18s 5-8s; do
  parallel --keep-order --xapply -j 8 "
  perl NJU_seq/rrna_analysis/score.pl \\
    output/NJU{1}/rrna_${RNA}.tsv \\
    output/NJU{2}/rrna_${RNA}.tsv \\
    output/NJU{3}/rrna_${RNA}.tsv \\
    output/NJU{4}/rrna_${RNA}.tsv \\
      >output/{5}_rrna_${RNA}_scored.tsv
" ::: $(seq 6200 4 6215) $(seq 6236 4 6251) ::: $(seq 6201 4 6215) $(seq 6237 4 6251) ::: $(seq 6202 4 6215) $(seq 6238 4 6251) ::: $(seq 6203 4 6215) $(seq 6239 4 6251) ::: Ath_sixl Ath_bolt Ath_bloom Ath_fruit Ath_root Ath_stem Ath_flower Ath_leaf
done

for PREFIX in Ath_sixl Ath_bolt Ath_bloom Ath_fruit Ath_root Ath_stem Ath_flower Ath_leaf; do
  bash NJU_seq/presentation/point_venn.sh \
    Sample1 output/${PREFIX}_rrna_18s_scored.tsv 14 \
    Sample2 output/${PREFIX}_rrna_18s_scored.tsv 15 \
    Sample3 output/${PREFIX}_rrna_18s_scored.tsv 16 \
    output/${PREFIX}_rrna_18s_venn.png 40
done

bsub -n 24 -q largemem -o NJU6200_count.out -J "count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz | awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' | perl NJU_seq/rrna_analysis/matchquality_judge.pl | perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU62{00..15} NJU62{36..51}
'

for PREFIX in NJU62{00..15} NJU62{36..51}; do
  bsub -q largemem -n 24 -o ${PREFIX}_mrna_alignment.log -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step2.sh ${PREFIX}"
done

for PREFIX in NJU61{84..99} NJU62{20..35}; do
  bsub -q largemem -n 24 -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step3.sh ${PREFIX} 16"
done