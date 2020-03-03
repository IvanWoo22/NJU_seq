# 2OMG
A brand new version of MgR_seq.

## 1. Preparation
#### Install packages and software.

Software managed by [brew](https://brew.sh/).
```shell script
# generic tools
brew install parallel pigz

# bioinformatics tools
brew install bowtie2
brew install picard-tools samtools
```

Perl packages:
```shell script
cpanm YAML::Syck AlignDB::IntSpan PerlIO::gzip
``` 
*To install PerlIO::gzip on [WSL2](https://devblogs.microsoft.com/commandline/announcing-wsl-2/), you might need to install [`zlib.h`](http://www.zlib.net/) manually.*

R packages needed:
`ggplot2` `ggpubr` `gridExtra` `forcats` `dplyr`

Make new folders for analysis.
```shell script
# 2OMG_analysis is the main directory of the following analysis. It can be renamed as you like.
mkdir "2OMG_analysis"
cd 2OMG_analysis
mkdir "data" "index" "temp" "output"
THREAD=16
```

## 2. Reference and Index
#### Download reference.
Get reference sequence of species from [GENCODE](https://www.gencodegenes.org/) and [Ensembl](http://plants.ensembl.org/Arabidopsis_thaliana/Info/Index?db=core).
```shell script
# GENCODE release 33 for human
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz -O data/hsa.gff3.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz -O data/hsa_transcript.fa.gz

# GENCODE release M24 for mouse
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gff3.gz -O data/mmu.gff3.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.transcripts.fa.gz -O data/mmu_transcript.fa.gz

# Ensembl release 46 for Arabidopsis thaliana
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz -O data/ath.gff3.gz
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz -O data/ath_transcript.fa.gz
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz -O data/ath_ncrna.fa.gz
```

#### Build index
Create index by `bowtie2-build` for mapping.
```shell script
# rRNA index
# 5.8s,18s and 25s rRNA.
cat ~/2OMG/data/ath_rrna/* >data/ath_rrna.fa
# rRNA from the reference ncRNA.
pigz -dc data/ath_ncrna.fa.gz |
  perl ~/2OMG/tool/fetch_fasta.pl \
  --stdin -s "transcript_biotype:rRNA" \
  >>data/ath_rrna.fa
bowtie2-build data/ath_rrna.fa index/ath_rrna
rm data/ath_rrna.fa

# mRNA index
# Only protein_coding transcripts will be extract to build index.
pigz -dc data/ath_transcript.fa.gz |
  perl ~/2OMG/tool/fetch_fasta.pl \
  --stdin -s "transcript_biotype:protein_coding" \
  >data/ath_protein_coding.fa
bowtie2-build --threads ${THREAD} \
  data/ath_protein_coding.fa index/ath_protein_coding
rm data/ath_protein_coding.fa
```

## 3. Data Selection and quality overview
#### Select data for analysing.
Get the sequencing clean data from `MgR_data`.  
*The representation of `ID` can be found in [`sample_list.csv`](/sample_list.csv).*
```shell script
ID=NJU45
PREFIX=Ath_stem_NC

mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}"
ln -sf /home/wyf/MgR_data/${ID}/R1.fq.gz data/${PREFIX}/R1.fq.gz
ln -sf /home/wyf/MgR_data/${ID}/R2.fq.gz data/${PREFIX}/R2.fq.gz
```
*Better to process the other samples in the same group according to the above code box. Here NJU45-48 are in one group.*

#### Quality control for clean data.
Input a `FastQ` file or a `GZ` file of `FastQ`, and then get some quality information.
```shell script
PREFIX=Ath_stem_NC

# For pair-end sequence data, we firstly turn them to single-end file.
perl ~/2OMG/quality_control/pe_consistency.pl \
  data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
  temp/${PREFIX}.fq.gz
# Total:  12627217
# Consistency:  12450207
# Proportion: 98.60%

# PREFIX as Ath_stem_NC Ath_stem_1 Ath_stem_2 Ath_stem_3.
# $ARGV[-2] should be the directory of output files.
# $ARGV[-1] should be the prefix of the output files.
time perl ~/2OMG/quality_control/fastq_qc.pl \
  temp/Ath_stem_NC.fq.gz \
  temp/Ath_stem_1.fq.gz \
  temp/Ath_stem_2.fq.gz \
  temp/Ath_stem_3.fq.gz \
  output \
  Ath_stem
# real    3m28.250s
# user    3m26.499s
# sys     0m0.495s
```
*The quality report created in `/output/Ath_stem.pdf`.*

## 4. Alignment, Count and Score
#### rRNA workflow
Use `bowtie2` to align the data file.
```shell script
time bowtie2 -p ${THREAD} -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 -i S,1,0.50 --np 0 \
  --xeq -x index/ath_rrna \
  -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz \
  -S output/${PREFIX}/rrna.raw.sam \
  2>&1 |
  tee output/${PREFIX}/rrna.bowtie2.log

time pigz -p ${THREAD} output/${PREFIX}/rrna.raw.sam
```
Filter and count alignment result.
```shell script
time gzip -dcf output/${PREFIX}/rrna.raw.sam.gz |
  parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j ${THREAD} '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/2OMG/rrna_analysis/matchquality_judge.pl |
    perl ~/2OMG/rrna_analysis/multimatch_judge.pl
  ' \
  >temp/${PREFIX}/rrna.out.tmp

time parallel -j 3 '
  perl ~/2OMG/rrna_analysis/readend_count.pl \
    ~/2OMG/data/ath_rrna/{}.fa temp/${PREFIX}/rrna.out.tmp {} \
    >output/${PREFIX}/rrna_{}.tsv \
  ' ::: 25s 18s 5-8s
```
Score all sites one by one.
```shell script
time parallel -j 3 '
  perl ~/2OMG/rrna_analysis/score.pl \
    output/Ath_stem_NC/rrna_{}.tsv \
    output/Ath_stem_1/rrna_{}.tsv \
    output/Ath_stem_2/rrna_{}.tsv \
    output/Ath_stem_3/rrna_{}.tsv \
    >output/Ath_stem_rrna_{}.tsv \
  ' ::: 25s 18s 5-8s
```

#### Prepare for mRNA.
Extract reads can't mapped to rRNA.
```shell script
bash ~/2OMG/tool/extract_fastq.sh \
  temp/${PREFIX}/rrna.out.sam \
  data/${PREFIX}/R1.fq.gz data/${PREFIX}/R1.mrna.fq.gz \
  data/${PREFIX}/R2.fq.gz data/${PREFIX}/R2.mrna.fq.gz
```

#### mRNA workflow
Alignment with protein_coding transcript.
```shell script
time bowtie2 -p 8 -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 --score-min C,0,0 \
  --xeq -x index/ath_protein_coding \
  -1 data/${PREFIX}/R1.mrna.fq.gz -2 data/${PREFIX}/R2.mrna.fq.gz \
  -S output/${PREFIX}/mrna.raw.sam \
  2>&1 |
  tee output/${PREFIX}/mrna.bowtie2.log

time pigz -p 8 output/${PREFIX}/mrna.raw.sam
```
Filter，re-locate and count alignment result.
```shell script
time gzip -dcf output/${PREFIX}/mrna.raw.sam.gz |
  parallel --pipe --block 10M -j 6 '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl ~/2OMG/mrna_analysis/multimatch_judge.pl
  ' \
  >temp/${PREFIX}/mrna.out.tmp

time gzip -dcf data/ath.gff3.gz |
  awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
    perl ~/2OMG/mrna_analysis/dedup.pl \
      --refstr "Parent=transcript:" \
      --geneid "AT" \
      -i temp/${PREFIX}/mrna.out.tmp \
      -o temp/${PREFIX}/mrna.dedup.tmp

time bash ~/2OMG/mrna_analysis/almostunique.sh \
  temp/${PREFIX}/mrna.dedup.tmp \
  data/${PREFIX}/R1.mrna.fq.gz \
  temp/${PREFIX} \
  temp/${PREFIX}/mrna.almostunique.tmp

perl ~/2OMG/mrna_analysis/count.pl \
  temp/${PREFIX}/mrna.almostunique.tmp \
  >temp/${PREFIX}/mrna.count.tmp

time gzip -dcf data/ath.gff3.gz |
  awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
  perl ~/2OMG/mrna_analysis/merge.pl \
  --refstr "Parent=transcript:" \
  --geneid "AT" \
  -i temp/${PREFIX}/mrna.count.tmp \
  -o temp/${PREFIX}/mrna.position.tmp

for chr in {1..5} Mt Pt; do
  awk -va=${chr} '$1==a&&$3=="+"' \
    temp/${PREFIX}/mrna.position.tmp |
      sort -t $'\t' -nk 2,2 \
        >>output/${PREFIX}/mrna.tsv
  awk -va=${chr} '$1==a&&$3=="-"' \
    temp/${PREFIX}/mrna.position.tmp |
      sort -t $'\t' -nrk 2,2 \
        >>output/${PREFIX}/mrna.tsv
done

```



