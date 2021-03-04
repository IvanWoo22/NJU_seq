for PREFIX in NJU6{276..355}; do
  bsub -n 24 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
    -O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
    arabidopsis/${PREFIX}_*_R1_*.gz arabidopsis/${PREFIX}_*_R2_*.gz -j 20
    "
done

for PREFIX in NJU6{252..255} NJU6{260..267} NJU6{272..275}; do
  mkdir ${PREFIX}
  bsub -n 24 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
    -O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
    plant/${PREFIX}_*_R1_*.gz plant/${PREFIX}_*_R2_*.gz -j 20
    "
done

for PREFIX in NJU6{252..255} NJU6{260..267} NJU6{272..355}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 4 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
  perl NJU_seq/quality_control/pe_consistency.pl \
    data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
    temp/${PREFIX}.fq.gz
    "
done

for PREFIX in NJU6{280..299}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 4 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
  perl NJU_seq/quality_control/pe_consistency.pl \
    data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
    temp/${PREFIX}.fq.gz
    "
done

for PREFIX in NJU6{300..319}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 4 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
  perl NJU_seq/quality_control/pe_consistency.pl \
    data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
    temp/${PREFIX}.fq.gz
    "
done

for PREFIX in NJU6{320..339}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 4 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
  perl NJU_seq/quality_control/pe_consistency.pl \
    data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
    temp/${PREFIX}.fq.gz
    "
done

for PREFIX in NJU6{340..355}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 4 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
  perl NJU_seq/quality_control/pe_consistency.pl \
    data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
    temp/${PREFIX}.fq.gz
    "
done

for PREFIX in NJU6{276..355}; do
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step1.sh ${PREFIX}"
done

bsub -n 24 -o log/NJU6276_6300_rrna_count.log -J "NJU6276_6300_count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz | awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' | perl NJU_seq/rrna_analysis/matchquality_judge.pl | perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU6{276..300}
'

bsub -n 24 -o log/NJU6301_6325_rrna_count.log -J "NJU6301_6325_count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz | awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' | perl NJU_seq/rrna_analysis/matchquality_judge.pl | perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU6{301..325}
'

bsub -n 24 -o log/NJU6326_6355_rrna_count.log -J "NJU6326_6355_count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz | awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' | perl NJU_seq/rrna_analysis/matchquality_judge.pl | perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU6{326..355}
'

for PREFIX in NJU6{276..355}; do
  bsub -q largemem -n 24 -o log/${PREFIX}_mrna_alignment.log -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step2.sh ${PREFIX}"
done

for PREFIX in NJU6{276..355}; do
  bsub -q largemem -n 24 -o log/${PREFIX}_dedup.log -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step3.sh ${PREFIX} 16"
done

for PREFIX in NJU6{184..215}; do
  bsub -q fat_768 -n 80 -o log/${PREFIX}_bac_alignment.log -J "${PREFIX}" "bash bac_target.sh ${PREFIX}"
done

bsub -n 24 -o log/mash.log -J "mash" '
parallel --xapply --keep-order -j 6 "
pigz -dc data/{}/R1.fq.gz >data/{}/R1.fq
perl perbool/filter_fastq.pl --min 21 \
  -i data/{}/R1.fq \
  -o data/{}/R1_filter_min21.fq
rm data/{}/R1.fq
./mash-Linux64-v2.2/mash screen -w -p 4 \
  refseq.genomes.k21s1000.msh \
  data/{}/R1_filter_min21.fq \
  >output/{}/mash_screen.tsv
rm data/{}/R1_filter_min21.fq
" ::: NJU6{184..215}
'

