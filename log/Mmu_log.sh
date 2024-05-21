#!/usr/bin/env bash

THREAD=16

mkdir "20200417"
cd "20200417" || exit
mkdir "data" "index" "temp" "output"

wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gff3.gz -O data/mmu.gff3.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.transcripts.fa.gz -O data/mmu_transcript.fa.gz

cat ~/NJU_seq/data/mmu_rrna/* >data/mmu_rrna.fa
bowtie2-build data/mmu_rrna.fa index/mmu_rrna
rm data/mmu_rrna.fa

pigz -dc data/mmu_transcript.fa.gz \
	| perl ~/NJU_seq/tool/fetch_fasta.pl \
		--stdin -s 'protein_coding' \
		>data/mmu_protein_coding.fa

bowtie2-build --threads "${THREAD}" \
	data/mmu_protein_coding.fa index/mmu_protein_coding
rm data/mmu_protein_coding.fa

############################################################################

ID='NJU61'
PREFIX='Mmu_lung_NC'

mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}"
ln -sf /home/wyf/MgR_data/"${ID}"/R1.fq.gz data/"${PREFIX}"/R1.fq.gz
ln -sf /home/wyf/MgR_data/"${ID}"/R2.fq.gz data/"${PREFIX}"/R2.fq.gz

time bowtie2 -p "${THREAD}" -a -t \
	--end-to-end -D 20 -R 3 \
	-N 0 -L 10 -i S,1,0.50 --np 0 \
	--xeq -x index/mmu_rrna \
	-1 data/"${PREFIX}"/R1.fq.gz -2 data/"${PREFIX}"/R2.fq.gz \
	-S output/"${PREFIX}"/rrna.raw.sam \
	2>&1 \
	| tee output/"${PREFIX}"/rrna.bowtie2.log

time pigz -p "${THREAD}" output/"${PREFIX}"/rrna.raw.sam

time pigz -dcf output/"${PREFIX}"/rrna.raw.sam.gz \
	| parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
		>temp/"${PREFIX}"/rrna.out.tmp

time parallel -j 3 "
  perl ~/NJU_seq/rrna_analysis/readend_count.pl \\
    ~/NJU_seq/data/mmu_rrna/{}.fa temp/${PREFIX}/rrna.out.tmp {} \\
    >output/${PREFIX}/rrna_{}.tsv
  " ::: 28s 18s 5-8s

############################################################################

ID='NJU62'
PREFIX='Mmu_lung_1'

mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}"
ln -sf /home/wyf/MgR_data/"${ID}"/R1.fq.gz data/"${PREFIX}"/R1.fq.gz
ln -sf /home/wyf/MgR_data/"${ID}"/R2.fq.gz data/"${PREFIX}"/R2.fq.gz

time bowtie2 -p "${THREAD}" -a -t \
	--end-to-end -D 20 -R 3 \
	-N 0 -L 10 -i S,1,0.50 --np 0 \
	--xeq -x index/mmu_rrna \
	-1 data/"${PREFIX}"/R1.fq.gz -2 data/"${PREFIX}"/R2.fq.gz \
	-S output/"${PREFIX}"/rrna.raw.sam \
	2>&1 \
	| tee output/"${PREFIX}"/rrna.bowtie2.log

time pigz -p "${THREAD}" output/"${PREFIX}"/rrna.raw.sam

time pigz -dcf output/"${PREFIX}"/rrna.raw.sam.gz \
	| parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
		>temp/"${PREFIX}"/rrna.out.tmp

time parallel -j 3 "
  perl ~/NJU_seq/rrna_analysis/readend_count.pl \\
    ~/NJU_seq/data/mmu_rrna/{}.fa temp/${PREFIX}/rrna.out.tmp {} \\
    >output/${PREFIX}/rrna_{}.tsv
  " ::: 28s 18s 5-8s

############################################################################

ID='NJU63'
PREFIX='Mmu_lung_2'

mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}"
ln -sf /home/wyf/MgR_data/"${ID}"/R1.fq.gz data/"${PREFIX}"/R1.fq.gz
ln -sf /home/wyf/MgR_data/"${ID}"/R2.fq.gz data/"${PREFIX}"/R2.fq.gz

time bowtie2 -p "${THREAD}" -a -t \
	--end-to-end -D 20 -R 3 \
	-N 0 -L 10 -i S,1,0.50 --np 0 \
	--xeq -x index/mmu_rrna \
	-1 data/"${PREFIX}"/R1.fq.gz -2 data/"${PREFIX}"/R2.fq.gz \
	-S output/"${PREFIX}"/rrna.raw.sam \
	2>&1 \
	| tee output/"${PREFIX}"/rrna.bowtie2.log

time pigz -p "${THREAD}" output/"${PREFIX}"/rrna.raw.sam

time pigz -dcf output/"${PREFIX}"/rrna.raw.sam.gz \
	| parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
		>temp/"${PREFIX}"/rrna.out.tmp

time parallel -j 3 "
  perl ~/NJU_seq/rrna_analysis/readend_count.pl \\
    ~/NJU_seq/data/mmu_rrna/{}.fa temp/${PREFIX}/rrna.out.tmp {} \\
    >output/${PREFIX}/rrna_{}.tsv
  " ::: 28s 18s 5-8s

############################################################################

ID='NJU64'
PREFIX='Mmu_lung_3'

mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}"
ln -sf /home/wyf/MgR_data/"${ID}"/R1.fq.gz data/"${PREFIX}"/R1.fq.gz
ln -sf /home/wyf/MgR_data/"${ID}"/R2.fq.gz data/"${PREFIX}"/R2.fq.gz

time bowtie2 -p "${THREAD}" -a -t \
	--end-to-end -D 20 -R 3 \
	-N 0 -L 10 -i S,1,0.50 --np 0 \
	--xeq -x index/mmu_rrna \
	-1 data/"${PREFIX}"/R1.fq.gz -2 data/"${PREFIX}"/R2.fq.gz \
	-S output/"${PREFIX}"/rrna.raw.sam \
	2>&1 \
	| tee output/"${PREFIX}"/rrna.bowtie2.log

time pigz -p "${THREAD}" output/"${PREFIX}"/rrna.raw.sam

time pigz -dcf output/"${PREFIX}"/rrna.raw.sam.gz \
	| parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
		>temp/"${PREFIX}"/rrna.out.tmp

############################################################################

time parallel -j 3 "
  perl ~/NJU_seq/rrna_analysis/readend_count.pl \\
    ~/NJU_seq/data/mmu_rrna/{}.fa temp/${PREFIX}/rrna.out.tmp {} \\
    >output/${PREFIX}/rrna_{}.tsv
  " ::: 28s 18s 5-8s

time parallel -j 3 "
  perl ~/NJU_seq/rrna_analysis/score.pl \\
    output/Mmu_lung_NC/rrna_{}.tsv \\
    output/Mmu_lung_1/rrna_{}.tsv \\
    output/Mmu_lung_2/rrna_{}.tsv \\
    output/Mmu_lung_3/rrna_{}.tsv \\
      >output/Mmu_lung_rrna_{}_scored.tsv
  " ::: 28s 18s 5-8s

############################################################################

PREFIX='Mmu_lung_NC'

time bash ~/NJU_seq/tool/extract_fastq.sh \
	temp/"${PREFIX}"/rrna.out.tmp \
	data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R1.mrna.fq.gz \
	data/"${PREFIX}"/R2.fq.gz data/"${PREFIX}"/R2.mrna.fq.gz

############################################################################

PREFIX='Mmu_lung_1'

time bash ~/NJU_seq/tool/extract_fastq.sh \
	temp/"${PREFIX}"/rrna.out.tmp \
	data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R1.mrna.fq.gz \
	data/"${PREFIX}"/R2.fq.gz data/"${PREFIX}"/R2.mrna.fq.gz

############################################################################

PREFIX='Mmu_lung_2'

time bash ~/NJU_seq/tool/extract_fastq.sh \
	temp/"${PREFIX}"/rrna.out.tmp \
	data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R1.mrna.fq.gz \
	data/"${PREFIX}"/R2.fq.gz data/"${PREFIX}"/R2.mrna.fq.gz

############################################################################

PREFIX='Mmu_lung_3'

time bash ~/NJU_seq/tool/extract_fastq.sh \
	temp/"${PREFIX}"/rrna.out.tmp \
	data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R1.mrna.fq.gz \
	data/"${PREFIX}"/R2.fq.gz data/"${PREFIX}"/R2.mrna.fq.gz

############################################################################

PREFIX='Mmu_lung_NC'

time bowtie2 -p "${THREAD}" -a -t \
	--end-to-end -D 20 -R 3 \
	-N 0 -L 10 --score-min C,0,0 \
	--xeq -x index/mmu_protein_coding \
	-1 data/"${PREFIX}"/R1.mrna.fq.gz -2 data/"${PREFIX}"/R2.mrna.fq.gz \
	-S output/"${PREFIX}"/mrna.raw.sam \
	2>&1 \
	| tee output/"${PREFIX}"/mrna.bowtie2.log

time pigz -p "${THREAD}" output/"${PREFIX}"/mrna.raw.sam

time gzip -dcf output/"${PREFIX}"/mrna.raw.sam.gz \
	| parallel --pipe --block 10M -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl
  ' \
		>temp/"${PREFIX}"/mrna.out.tmp

time gzip -dcf data/mmu.gff3.gz \
	| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
	| perl ~/NJU_seq/mrna_analysis/dedup.pl \
		--refstr "Parent=" \
		--geneid "ENSMUST" \
		-i temp/"${PREFIX}"/mrna.out.tmp \
		-o temp/"${PREFIX}"/mrna.dedup.tmp

time bash ~/NJU_seq/mrna_analysis/almostunique.sh \
	temp/"${PREFIX}"/mrna.dedup.tmp \
	data/"${PREFIX}"/R1.mrna.fq.gz \
	temp/"${PREFIX}" \
	temp/"${PREFIX}"/mrna.almostunique.tmp

time perl ~/NJU_seq/mrna_analysis/count.pl \
	temp/"${PREFIX}"/mrna.almostunique.tmp \
	>temp/"${PREFIX}"/mrna.count.tmp

time gzip -dcf data/mmu.gff3.gz \
	| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
	| perl ~/NJU_seq/mrna_analysis/merge.pl \
		--refstr "Parent=" \
		--geneid "ENSMUSG" \
		--transid "ENSMUST" \
		-i temp/"${PREFIX}"/mrna.count.tmp \
		-o temp/"${PREFIX}"/mrna.position.tmp

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y M; do
	awk -va=${chr} '$1==a&&$3=="+"' \
		temp/"${PREFIX}"/mrna.position.tmp \
		| sort -t $'\t' -nk 2,2 \
			>>output/"${PREFIX}"/mrna.tsv
	awk -va=${chr} '$1==a&&$3=="-"' \
		temp/"${PREFIX}"/mrna.position.tmp \
		| sort -t $'\t' -nrk 2,2 \
			>>output/"${PREFIX}"/mrna.tsv
done

############################################################################

PREFIX='Mmu_lung_1'

time bowtie2 -p "${THREAD}" -a -t \
	--end-to-end -D 20 -R 3 \
	-N 0 -L 10 --score-min C,0,0 \
	--xeq -x index/mmu_protein_coding \
	-1 data/"${PREFIX}"/R1.mrna.fq.gz -2 data/"${PREFIX}"/R2.mrna.fq.gz \
	-S output/"${PREFIX}"/mrna.raw.sam \
	2>&1 \
	| tee output/"${PREFIX}"/mrna.bowtie2.log

time pigz -p "${THREAD}" output/"${PREFIX}"/mrna.raw.sam

time gzip -dcf output/"${PREFIX}"/mrna.raw.sam.gz \
	| parallel --pipe --block 10M -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl
  ' \
		>temp/"${PREFIX}"/mrna.out.tmp

time gzip -dcf data/mmu.gff3.gz \
	| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
	| perl ~/NJU_seq/mrna_analysis/dedup.pl \
		--refstr "Parent=" \
		--geneid "ENSMUST" \
		-i temp/"${PREFIX}"/mrna.out.tmp \
		-o temp/"${PREFIX}"/mrna.dedup.tmp

time bash ~/NJU_seq/mrna_analysis/almostunique.sh \
	temp/"${PREFIX}"/mrna.dedup.tmp \
	data/"${PREFIX}"/R1.mrna.fq.gz \
	temp/"${PREFIX}" \
	temp/"${PREFIX}"/mrna.almostunique.tmp

time perl ~/NJU_seq/mrna_analysis/count.pl \
	temp/"${PREFIX}"/mrna.almostunique.tmp \
	>temp/"${PREFIX}"/mrna.count.tmp

time gzip -dcf data/mmu.gff3.gz \
	| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
	| perl ~/NJU_seq/mrna_analysis/merge.pl \
		--refstr "Parent=" \
		--geneid "ENSMUST" \
		-i temp/"${PREFIX}"/mrna.count.tmp \
		-o temp/"${PREFIX}"/mrna.position.tmp

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y M; do
	awk -va=${chr} '$1==a&&$3=="+"' \
		temp/"${PREFIX}"/mrna.position.tmp \
		| sort -t $'\t' -nk 2,2 \
			>>output/"${PREFIX}"/mrna.tsv
	awk -va=${chr} '$1==a&&$3=="-"' \
		temp/"${PREFIX}"/mrna.position.tmp \
		| sort -t $'\t' -nrk 2,2 \
			>>output/"${PREFIX}"/mrna.tsv
done

############################################################################

PREFIX='Mmu_lung_2'

time bowtie2 -p "${THREAD}" -a -t \
	--end-to-end -D 20 -R 3 \
	-N 0 -L 10 --score-min C,0,0 \
	--xeq -x index/mmu_protein_coding \
	-1 data/"${PREFIX}"/R1.mrna.fq.gz -2 data/"${PREFIX}"/R2.mrna.fq.gz \
	-S output/"${PREFIX}"/mrna.raw.sam \
	2>&1 \
	| tee output/"${PREFIX}"/mrna.bowtie2.log

time pigz -p "${THREAD}" output/"${PREFIX}"/mrna.raw.sam

time gzip -dcf output/"${PREFIX}"/mrna.raw.sam.gz \
	| parallel --pipe --block 10M -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl
  ' \
		>temp/"${PREFIX}"/mrna.out.tmp

time gzip -dcf data/mmu.gff3.gz \
	| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
	| perl ~/NJU_seq/mrna_analysis/dedup.pl \
		--refstr "Parent=" \
		--geneid "ENSMUST" \
		-i temp/"${PREFIX}"/mrna.out.tmp \
		-o temp/"${PREFIX}"/mrna.dedup.tmp

time bash ~/NJU_seq/mrna_analysis/almostunique.sh \
	temp/"${PREFIX}"/mrna.dedup.tmp \
	data/"${PREFIX}"/R1.mrna.fq.gz \
	temp/"${PREFIX}" \
	temp/"${PREFIX}"/mrna.almostunique.tmp

time perl ~/NJU_seq/mrna_analysis/count.pl \
	temp/"${PREFIX}"/mrna.almostunique.tmp \
	>temp/"${PREFIX}"/mrna.count.tmp

time gzip -dcf data/mmu.gff3.gz \
	| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
	| perl ~/NJU_seq/mrna_analysis/merge.pl \
		--refstr "Parent=" \
		--geneid "ENSMUST" \
		-i temp/"${PREFIX}"/mrna.count.tmp \
		-o temp/"${PREFIX}"/mrna.position.tmp

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y M; do
	awk -va=${chr} '$1==a&&$3=="+"' \
		temp/"${PREFIX}"/mrna.position.tmp \
		| sort -t $'\t' -nk 2,2 \
			>>output/"${PREFIX}"/mrna.tsv
	awk -va=${chr} '$1==a&&$3=="-"' \
		temp/"${PREFIX}"/mrna.position.tmp \
		| sort -t $'\t' -nrk 2,2 \
			>>output/"${PREFIX}"/mrna.tsv
done

############################################################################

PREFIX='Mmu_lung_3'

time bowtie2 -p "${THREAD}" -a -t \
	--end-to-end -D 20 -R 3 \
	-N 0 -L 10 --score-min C,0,0 \
	--xeq -x index/mmu_protein_coding \
	-1 data/"${PREFIX}"/R1.mrna.fq.gz -2 data/"${PREFIX}"/R2.mrna.fq.gz \
	-S output/"${PREFIX}"/mrna.raw.sam \
	2>&1 \
	| tee output/"${PREFIX}"/mrna.bowtie2.log

time pigz -p "${THREAD}" output/"${PREFIX}"/mrna.raw.sam

time gzip -dcf output/"${PREFIX}"/mrna.raw.sam.gz \
	| parallel --pipe --block 10M -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl
  ' \
		>temp/"${PREFIX}"/mrna.out.tmp

time gzip -dcf data/mmu.gff3.gz \
	| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
	| perl ~/NJU_seq/mrna_analysis/dedup.pl \
		--refstr "Parent=" \
		--geneid "ENSMUST" \
		-i temp/"${PREFIX}"/mrna.out.tmp \
		-o temp/"${PREFIX}"/mrna.dedup.tmp

time bash ~/NJU_seq/mrna_analysis/almostunique.sh \
	temp/"${PREFIX}"/mrna.dedup.tmp \
	data/"${PREFIX}"/R1.mrna.fq.gz \
	temp/"${PREFIX}" \
	temp/"${PREFIX}"/mrna.almostunique.tmp

time perl ~/NJU_seq/mrna_analysis/count.pl \
	temp/"${PREFIX}"/mrna.almostunique.tmp \
	>temp/"${PREFIX}"/mrna.count.tmp

time gzip -dcf data/mmu.gff3.gz \
	| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
	| perl ~/NJU_seq/mrna_analysis/merge.pl \
		--refstr "Parent=" \
		--geneid "ENSMUST" \
		-i temp/"${PREFIX}"/mrna.count.tmp \
		-o temp/"${PREFIX}"/mrna.position.tmp

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y M; do
	awk -va=${chr} '$1==a&&$3=="+"' \
		temp/"${PREFIX}"/mrna.position.tmp \
		| sort -t $'\t' -nk 2,2 \
			>>output/"${PREFIX}"/mrna.tsv
	awk -va=${chr} '$1==a&&$3=="-"' \
		temp/"${PREFIX}"/mrna.position.tmp \
		| sort -t $'\t' -nrk 2,2 \
			>>output/"${PREFIX}"/mrna.tsv
done

############################################################################

parallel -j 3 "
  bash ~/NJU_seq/presentation/seq_depth.sh \\
    temp/{}/mrna.almostunique.tmp \\
    output/{}/mrna.tsv
  " ::: Mmu_lung_NC Mmu_lung_1 Mmu_lung_2 Mmu_lung_3

#All stop times: 3218328
#All positions:  200922
#Coverage:       16.01
#
#All stop times: 2636715
#All positions:  160061
#Coverage:       16.47
#
#All stop times: 3432243
#All positions:  163911
#Coverage:       20.93
#
#All stop times: 3612607
#All positions:  160264
#Coverage:       22.54

parallel -j 3 "
  perl ~/NJU_seq/mrna_analysis/score.pl 20 \\
    output/Mmu_lung_NC/mrna.tsv \\
    output/{}/mrna.tsv |
      sort -t $'\t' -nrk 10,10 \\
  >output/{}_mrna_scored.tsv
" ::: Mmu_lung_1 Mmu_lung_2 Mmu_lung_3

perl ~/NJU_seq/tool/common.pl \
	output/Mmu_lung_1_mrna_scored.tsv \
	output/Mmu_lung_2_mrna_scored.tsv \
	output/Mmu_lung_3_mrna_scored.tsv \
	>output/Mmu_lung_mrna_scored.tsv

bash ~/NJU_seq/presentation/signature_count.sh \
	output/Mmu_lung_mrna_scored.tsv

pigz -dcf data/mmu.gff3.gz \
	| awk '(($3=="gene")&&($9~/gene_type=protein_coding/)) || (($3=="transcript")&&($9~/transcript_type=protein_coding/)) || ($3=="exon")' \
	| perl ~/NJU_seq/mrna_analysis/filter_nonsenseexon.pl \
		--transwording "transcript" \
		--transid "ID=" \
		--exonid "Parent=" \
	| perl ~/NJU_seq/mrna_analysis/filter_nonsensegene.pl \
	| perl ~/NJU_seq/mrna_analysis/filter_overlapgene.pl \
	| perl ~/NJU_seq/mrna_analysis/judge_altersplice.pl \
		--gene_id "ID=" \
		--trans_wording "transcript" \
		--alter "data/mmu_alter_gene.yml" \
		--unique "data/mmu_unique_gene.yml"

perl ~/NJU_seq/mrna_analysis/stat_altersplice_1.pl \
	data/mmu_alter_gene.yml \
	<output/Mmu_lung_mrna_scored.tsv \
	>temp/Mmu_lung_altergene.tsv
