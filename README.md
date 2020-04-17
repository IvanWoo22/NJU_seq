# NJU-seq
**N**m site **J**udge **U**niversally **seq**uencing.  
A brand-new version of MgR-seq.  
From clean data to positions with 2-O-methylation.

## 1. Preparation
#### Install packages and software.

Software managed by [brew](https://brew.sh/).
```bash
# generic tools
brew install parallel pigz

# bioinformatics tools
brew install bowtie2
brew install picard-tools samtools
```

Perl packages:
```bash
cpanm YAML::Syck AlignDB::IntSpan PerlIO::gzip
``` 
*To install PerlIO::gzip on [WSL2](https://devblogs.microsoft.com/commandline/announcing-wsl-2/), you might need to install [`zlib.h`](http://www.zlib.net/) manually.*

R packages needed:
`ggplot2` `ggpubr` `gridExtra` `forcats` `dplyr`

Make new folders for analysis.
```bash
# NJU_seq_analysis is the main directory of the following analysis. It can be renamed as you like.
mkdir "NJU_seq_analysis"
cd NJU_seq_analysis
mkdir "data" "index" "temp" "output"
```

## 2. Reference and index
#### Download reference.
Get reference sequence of species from [GENCODE](https://www.gencodegenes.org/) and [Ensembl](http://plants.ensembl.org/Arabidopsis_thaliana/Info/Index?db=core).
```bash
# GENCODE release 33 for human
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz -O data/hsa.gff3.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz -O data/hsa_transcript.fa.gz

# GENCODE release M24 for mouse
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gff3.gz -O data/mmu.gff3.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.transcripts.fa.gz -O data/mmu_transcript.fa.gz

# Ensembl release 46 for Arabidopsis thaliana
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -O data/ath_dna.fa.gz
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz -O data/ath.gff3.gz
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz -O data/ath_transcript.fa.gz
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz -O data/ath_ncrna.fa.gz

# miRBase latest release for miRNA
wget -N ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz -O data/mirna.fa.gz
```

#### Build index
Create index by `bowtie2-build` for mapping.
```bash
THREAD=16

# rRNA index
# 5.8s,18s and 25s rRNA.
cat ~/NJU_seq/data/ath_rrna/* >data/ath_rrna.fa
# rRNA from the reference ncRNA.
pigz -dc data/ath_ncrna.fa.gz |
  perl ~/NJU_seq/tool/fetch_fasta.pl \
  --stdin -s 'transcript_biotype:rRNA' \
  >>data/ath_rrna.fa
bowtie2-build data/ath_rrna.fa index/ath_rrna
rm data/ath_rrna.fa

# mRNA index
# Only protein_coding transcripts will be extract to build index.
pigz -dc data/ath_transcript.fa.gz |
  perl ~/NJU_seq/tool/fetch_fasta.pl \
  --stdin -s 'transcript_biotype:protein_coding' \
  >data/ath_protein_coding.fa

bowtie2-build --threads "${THREAD}" \
  data/ath_protein_coding.fa index/ath_protein_coding
rm data/ath_protein_coding.fa
```

## 3. Data Selection and quality overview
#### Select data for analysing.
Get the sequencing clean data from `MgR_data`.  
*The representation of `ID` can be found in [`sample_list.csv`](/sample_list.csv).*
```bash
ID='NJU45'
PREFIX='Ath_stem_NC'

mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}"
ln -sf /home/wyf/MgR_data/"${ID}"/R1.fq.gz data/"${PREFIX}"/R1.fq.gz
ln -sf /home/wyf/MgR_data/"${ID}"/R2.fq.gz data/"${PREFIX}"/R2.fq.gz
```
*Better to process the other samples in the same group according to the above code box. Here NJU45-48 are in one group.*

#### Quality control for clean data.
Input a `FastQ` file or a `GZ` file of `FastQ`, and then get some quality information.
```bash
PREFIX='Ath_stem_NC'

# For pair-end sequence data, we firstly turn them to single-end file.
perl ~/NJU_seq/quality_control/pe_consistency.pl \
  data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R2.fq.gz \
  temp/"${PREFIX}".fq.gz
# Total:  12627217
# Consistency:  12450207
# Proportion: 98.60%

# PREFIX do as Ath_stem_NC Ath_stem_1 Ath_stem_2 Ath_stem_3.
# $ARGV[-2] should be the directory of output files.
# $ARGV[-1] should be the prefix of the output files.
time perl ~/NJU_seq/quality_control/fastq_qc.pl \
  temp/Ath_stem_NC.fq.gz \
  temp/Ath_stem_1.fq.gz \
  temp/Ath_stem_2.fq.gz \
  temp/Ath_stem_3.fq.gz \
  output \
  Ath_stem
# real  3m28.250s
# user  3m26.499s
# sys   0m0.495s
```
*The quality report created in `/output/Ath_stem.pdf`.*

## 4. Alignment, Count and Score
#### rRNA workflow
Use `bowtie2` to align the data file.
```bash
THREAD=16
PREFIX='Ath_stem_NC'

time bowtie2 -p "${THREAD}" -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 -i S,1,0.50 --np 0 \
  --xeq -x index/ath_rrna \
  -1 data/"${PREFIX}"/R1.fq.gz -2 data/"${PREFIX}"/R2.fq.gz \
  -S output/"${PREFIX}"/rrna.raw.sam \
  2>&1 |
  tee output/"${PREFIX}"/rrna.bowtie2.log
# real  1m53.474s
# user  25m41.846s
# sys   4m15.462s

time pigz -p "${THREAD}" output/"${PREFIX}"/rrna.raw.sam
# real  0m57.962s
# user  3m54.868s
# sys   0m7.795s
```
Filter and count alignment result.
```bash
THREAD=16
PREFIX='Ath_stem_NC'

time pigz -dcf output/"${PREFIX}"/rrna.raw.sam.gz |
  parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
  >temp/"${PREFIX}"/rrna.out.tmp
# real  1m59.185s
# user  20m44.643s
# sys   1m24.060s

time parallel -j 3 "
  perl ~/NJU_seq/rrna_analysis/readend_count.pl \\
    ~/NJU_seq/data/ath_rrna/{}.fa temp/${PREFIX}/rrna.out.tmp {} \\
    >output/${PREFIX}/rrna_{}.tsv
  " ::: 25s 18s 5-8s
# real  0m27.994s
# user  0m59.369s
# sys   0m0.809s
```
Score all sites one by one.
```bash
time parallel -j 3 "
  perl ~/NJU_seq/rrna_analysis/score.pl \\
    output/Ath_stem_NC/rrna_{}.tsv \\
    output/Ath_stem_1/rrna_{}.tsv \\
    output/Ath_stem_2/rrna_{}.tsv \\
    output/Ath_stem_3/rrna_{}.tsv \\
      >output/Ath_stem_rrna_{}_scored.tsv
  " ::: 25s 18s 5-8s
# real  0m0.459s
# user  0m0.661s
# sys   0m0.048s
```

#### Prepare for mRNA.
Extract reads can't be mapped to rRNA.
```bash
PREFIX='Ath_stem_NC'

time bash ~/NJU_seq/tool/extract_fastq.sh \
  temp/"${PREFIX}"/rrna.out.tmp \
  data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R1.mrna.fq.gz \
  data/"${PREFIX}"/R2.fq.gz data/"${PREFIX}"/R2.mrna.fq.gz
# real  1m10.382s
# user  1m45.610s
# sys   0m3.268s
```

#### mRNA workflow
Alignment with protein_coding transcript.
```bash
THREAD=16
PREFIX='Ath_stem_NC'

time bowtie2 -p "${THREAD}" -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 --score-min C,0,0 \
  --xeq -x index/ath_protein_coding \
  -1 data/"${PREFIX}"/R1.mrna.fq.gz -2 data/"${PREFIX}"/R2.mrna.fq.gz \
  -S output/"${PREFIX}"/mrna.raw.sam \
  2>&1 |
  tee output/"${PREFIX}"/mrna.bowtie2.log
# real  6m43.243s
# user  119m34.749s
# sys   8m1.470s

time pigz -p "${THREAD}" output/"${PREFIX}"/mrna.raw.sam
# real  0m22.663s
# user  5m38.943s
# sys   0m9.316s
```
Filter，re-locate and count alignment result.
```bash
THREAD=16
PREFIX='Ath_stem_NC'

time gzip -dcf output/"${PREFIX}"/mrna.raw.sam.gz |
  parallel --pipe --block 10M -j "${THREAD}" '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl
  ' \
  >temp/"${PREFIX}"/mrna.out.tmp
# real  2m31.558s
# user  18m38.617s
# sys   2m8.716s

time gzip -dcf data/ath.gff3.gz |
  awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
    perl ~/NJU_seq/mrna_analysis/dedup.pl \
      --refstr "Parent=transcript:" \
      --geneid "AT" \
      -i temp/"${PREFIX}"/mrna.out.tmp \
      -o temp/"${PREFIX}"/mrna.dedup.tmp
# real  14m32.692s
# user  14m22.635s
# sys   0m10.268s

time bash ~/NJU_seq/mrna_analysis/almostunique.sh \
  temp/"${PREFIX}"/mrna.dedup.tmp \
  data/"${PREFIX}"/R1.mrna.fq.gz \
  temp/"${PREFIX}" \
  temp/"${PREFIX}"/mrna.almostunique.tmp
# real  2m52.775s
# user  2m41.319s
# sys   0m4.881s

time perl ~/NJU_seq/mrna_analysis/count.pl \
  temp/"${PREFIX}"/mrna.almostunique.tmp \
  >temp/"${PREFIX}"/mrna.count.tmp
# real  0m3.072s
# user  0m3.036s
# sys   0m0.036s

time gzip -dcf data/ath.gff3.gz |
  awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
  perl ~/NJU_seq/mrna_analysis/merge.pl \
    --refstr "Parent=transcript:" \
    --geneid "AT" \
    -i temp/"${PREFIX}"/mrna.count.tmp \
    -o temp/"${PREFIX}"/mrna.position.tmp
# real  0m6.802s
# user  0m9.653s
# sys   0m0.227s
 
# The chromosome number changes according to the actual situation.
for chr in 1 2 3 4 5 Mt Pt; do
  awk -va=${chr} '$1==a&&$3=="+"' \
  temp/"${PREFIX}"/mrna.position.tmp |
    sort -t $'\t' -nk 2,2 \
  >>output/"${PREFIX}"/mrna.tsv
  awk -va=${chr} '$1==a&&$3=="-"' \
  temp/"${PREFIX}"/mrna.position.tmp |
    sort -t $'\t' -nrk 2,2 \
  >>output/"${PREFIX}"/mrna.tsv
done
```
Score each covered site.
```bash
parallel -j 3 "
  perl ~/NJU_seq/mrna_analysis/score.pl \\
    output/Ath_stem_NC/mrna.tsv \\
    output/{}/mrna.tsv |
      sort -t $'\t' -nrk 11,11 \\
        >output/{}_mrna_scored.tsv
  " ::: Ath_stem_1 Ath_stem_2 Ath_stem_3

perl ~/NJU_seq/tool/common.pl \
  output/Ath_stem_1_mrna_scored.tsv \
  output/Ath_stem_2_mrna_scored.tsv \
  output/Ath_stem_3_mrna_scored.tsv \
  >output/Ath_stem_mrna_scored.tsv
```
## 5. Statistics and Presentation
Calculate valid sequencing depth (average coverage).
```bash
PREFIX='Ath_stem_NC'

bash ~/NJU_seq/presentation/seq_depth.sh \
  temp/"${PREFIX}"/mrna.almostunique.tmp \
  output/"${PREFIX}"/mrna.tsv
# All stop times: 19227
# All positions:  5632
# Coverage:       3.41shell script
```
See the signature of Nm site.
```bash
bash ~/NJU_seq/presentation/signature_count.sh \
  output/Ath_stem_mrna_scored.tsv
```
Divide annotations into two categories.
```bash
pigz -dcf data/ath.gff3.gz |
  awk '(($3=="gene")&&($9~/biotype=protein_coding/)) \
    || (($3=="mRNA")&&($9~/biotype=protein_coding/)) \
    || ($3=="exon")' |
  perl ~/NJU_seq/mrna_analysis/filter_nonsenseexon.pl |
  perl ~/NJU_seq/mrna_analysis/filter_nonsensegene.pl |
  perl ~/NJU_seq/mrna_analysis/filter_overlapgene.pl |
  perl ~/NJU_seq/mrna_analysis/judge_altersplice.pl \
    data/ath_alter_gene.yml \
    data/ath_unique_gene.yml

perl ~/NJU_seq/mrna_analysis/stat_altersplice_1.pl \
  data/ath_alter_gene.yml \
  <output/Ath_stem_mrna_scored.tsv \
  >temp/Ath_stem_altergene.tsv

perl ~/NJU_seq/mrna_analysis/stat_altersplice_2.pl \
  <temp/Ath_stem_altergene.tsv \
  >output/Ath_stem_altergene_cov.tsv
# Output file is the final result.

pigz -dcf data/ath.gff3.gz |
  awk '($3=="gene") || ($3=="mRNA") || ($3=="CDS") \
    || ($3=="five_prime_UTR") || ($3=="three_prime_UTR")' |
  perl ~/NJU_seq/mrna_analysis/judge_differentregion.pl \
    data/ath_unique_gene.yml \
    >temp/ath_uniquegene_differentregion.yml

perl ~/NJU_seq/mrna_analysis/judge_region.pl \
  data/ath_alter_gene.yml temp/ath_uniquegene_differentregion.yml \
  <output/Ath_stem_mrna_scored.tsv \
  >output/Ath_stem_mrna_scored_region.tsv

perl ~/NJU_seq/mrna_analysis/stat_differentregion_1.pl \
  temp/ath_uniquegene_differentregion.yml \
  <output/Ath_stem_mrna_scored.tsv \
  >temp/Ath_stem_uniquegene.tsv

perl  ~/NJU_seq/mrna_analysis/stat_differentregion_2.pl \
  <temp/Ath_stem_uniquegene.tsv \
  >output/Ath_stem_uniquegene_cov.tsv

perl ~/NJU_seq/mrna_analysis/stat_differentregion_3.pl \
  temp/ath_uniquegene_differentregion.yml \
  <output/Ath_stem_mrna_scored.tsv \
  >output/Ath_stem_uniquegene_distribution.tsv
```

## 6. Motif Found in miRNA
Get miRNA sequence information from ncRNA reference.
```bash
pigz -dc data/mirna.fa.gz |
  perl ~/NJU_seq/tool/fetch_fasta.pl \
    --stdin -s 'thaliana' |
  perl ~/NJU_seq/mrna_analysis/motif_mirna.pl \
    >data/ath_mirna_motif.tsv

perl ~/NJU_seq/mrna_analysis/motif_nm.pl \
  data/ath_dna.fa.gz \
  output/Ath_stem_mrna_scored_sorted.tsv \
  >output/Ath_stem_mrna_motif.tsv

perl ~/NJU_seq/mrna_analysis/motif_compare.pl \
  data/ath_mirna_motif.tsv \
  output/Ath_stem_mrna_motif.tsv
```
