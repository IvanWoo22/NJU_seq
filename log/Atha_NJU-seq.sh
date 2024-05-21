for SAMPLE in NJUzp{4..10}; do
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${SAMPLE}/R1_clean.fq.gz -p ${SAMPLE}/R2_clean.fq.gz \
		-j 16 ${SAMPLE}/R1.fq.gz ${SAMPLE}/R2.fq.gz \
		2>${SAMPLE}/cutadapt.log

	perl ~/NJU_seq/quality_control/pe_consistency.pl ${SAMPLE}/R1_clean.fq.gz ${SAMPLE}/R2_clean.fq.gz ${SAMPLE}/merged_temp.fq.gz
	time bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x ../index/ath_rrna \
		-1 ${SAMPLE}/R1_clean.fq.gz -2 ${SAMPLE}/R2_clean.fq.gz \
		-S ${SAMPLE}/rrna.raw.sam \
		2>&1 | tee ${SAMPLE}/rrna.bowtie2.log

	time pigz -p 16 ${SAMPLE}/rrna.raw.sam
	time pigz -dcf ${SAMPLE}/rrna.raw.sam.gz \
		| parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j 12 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl' \
			>${SAMPLE}/rrna.out.tmp

	time bash ~/NJU_seq/tool/extract_fastq.sh \
		${SAMPLE}/rrna.out.tmp \
		${SAMPLE}/R1_clean.fq.gz ${SAMPLE}/R1.mrna.fq.gz \
		${SAMPLE}/R2_clean.fq.gz ${SAMPLE}/R2.mrna.fq.gz

	time bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 --score-min C,0,0 \
		--xeq -x ../index/ath_protein_coding \
		-1 ${SAMPLE}/R1.mrna.fq.gz -2 ${SAMPLE}/R2.mrna.fq.gz \
		-S ${SAMPLE}/mrna.raw.sam \
		2>&1 \
		| tee ${SAMPLE}/mrna.bowtie2.log

	parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 12 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' |
    perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl
  ' <${SAMPLE}/mrna.raw.sam | perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl \
		>${SAMPLE}/mrna.out.tmp

	parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 12 '
    /home/linuxbrew/.linuxbrew/bin/perl ~/NJU_seq/mrna_analysis/dedup.pl --refstr "Parent=transcript:" --transid "AT" --info ../data/ath_exon.info
  ' <${SAMPLE}/mrna.out.tmp \
		| /home/linuxbrew/.linuxbrew/bin/perl ~/NJU_seq/mrna_analysis/dedup.pl \
			--refstr "Parent=transcript:" \
			--transid "AT" \
			--info ../data/ath_exon.info \
			>${SAMPLE}/mrna.dedup.tmp

	time perl almostuniquematchneo.pl \
		${SAMPLE}/R1.mrna.fq.gz \
		${SAMPLE}/mrna.dedup.tmp \
		${SAMPLE}/mrna.almostuniquematchneo.tmp

	time perl countneo.pl \
		${SAMPLE}/mrna.almostuniquematchneo.tmp \
		>${SAMPLE}/mrna.count.tmp

	time gzip -dcf ../data/ath.gff3.gz \
		| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
		| perl ~/NJU_seq/mrna_analysis/merge.pl \
			--refstr "Parent=transcript:" \
			--geneid "AT" \
			--transid "AT" \
			-i ${SAMPLE}/mrna.count.tmp \
			-o ${SAMPLE}/mrna.tsv
	pigz ${SAMPLE}/mrna.raw.sam
done

perl ~/NJU_seq/tool/stat_alignment.pl NJUzp{3..10}/rrna.bowtie2.log NJUzp{3..10}/mrna.bowtie2.log \
	| Rscript ~/NJU_seq/tool/draw_table.R bowtie2.stat.pdf
perl /home/ivan/cater_data/ivan/fat/mscore.modified.pl \
	NJUzp6/mrna.tsv NJUzp3/mrna.tsv NJUzp4/mrna.tsv NJUzp5/mrna.tsv \
	>mrna.scored.Nm.MgCl2.tsv
perl /home/ivan/cater_data/ivan/fat/mscore.modified.pl \
	NJUzp10/mrna.tsv NJUzp7/mrna.tsv NJUzp8/mrna.tsv NJUzp9/mrna.tsv \
	>mrna.scored.Nm.Pst.tsv
