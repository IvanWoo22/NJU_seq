for PREFIX in NJU62{52..55} NJU62{64..67}; do
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "bash NJU_seq/log/Sly_scripts/step1.sh ${PREFIX}"
done

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz | awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' | perl NJU_seq/rrna_analysis/matchquality_judge.pl | perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU62{52..55} NJU62{64..67}
'

perl NJU_seq/tool/stat_alignment.pl log/NJU62{52..54,64..67}_rrna_alignment.log | Rscript NJU_seq/tool/draw_table.R output/"${PREFIX}"/rrna.bowtie2.pdf

for PREFIX in NJU62{52..55} NJU62{64..67}; do
  bsub -q fat_768 -n 80 -o log/${PREFIX}_bac_alignment.log -J "${PREFIX}" "bash bac_target.sh ${PREFIX}"
done