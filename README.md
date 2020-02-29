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
mkdir data
mkdir index
mkdir temp
mkdir output
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
cat ~/2OMG/data/ath_rrna/* >data/ath_rrna.fa
pigz -dc data/ath_ncrna.fa.gz |
  perl ~/2OMG/tool/fetch_fasta.pl \
  --stdin -s "transcript_biotype:rRNA" \
  >>data/ath_rrna.fa
bowtie2-build data/ath_rrna.fa index/ath_rrna
rm data/ath_rrna.fa

# Protein coding mRNA index
pigz -dc data/ath_transcript.fa.gz |
  perl ~/2OMG/tool/fetch_fasta.pl \
  --stdin -s "transcript_biotype:protein_coding" \
  >data/ath_protein_coding.fa
bowtie2-build data/ath_protein_coding.fa index/ath_protein_coding
```

## 3. Data Selection and quality overview
#### Select data for analysing.
Get the sequencing clean data from `MgR_data`.  
*The representation of `ID` can be found in [`sample_list.csv`](/sample_list.csv).*
```shell script
ID=NJU45
PREFIX=Ath_stem_NC
mkdir -p data/${PREFIX}
mkdir -p output/${PREFIX}

ln -sf /home/wyf/MgR_data/${ID}/R1.fq.gz data/${PREFIX}/R1.fq.gz
ln -sf /home/wyf/MgR_data/${ID}/R2.fq.gz data/${PREFIX}/R2.fq.gz
```

#### Quality Control
Input a `FastQ` file or a `GZ` file of `FastQ` and then get some quality information.
```shell script
# For pair-end sequence data, we firstly turn them to single-end file.
perl ~/2OMG/quality_control/pe_consistency.pl \
  data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
  temp/${PREFIX}.fq
perl ~/2OMG/quality_control/fastq_qc.pl \
  temp/Ath_stem_NC.fq \
  temp/Ath_stem_1.fq \
  temp/Ath_stem_2.fq \
  temp/Ath_stem_3.fq \
  output/Ath_stem
```

## 4. Alignment and Count
Use `bowtie2` to align the data file.
```shell script
time bowtie2 -p 8 -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 -i S,1,0.50 --np 0 \
  --xeq -x index/ath_rrna \
  -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz \
  -S output/${PREFIX}/rrna.raw.sam \
  2>&1 |
  tee output/${PREFIX}/rrna.bowtie2.log

time pigz -p 8 data/${PREFIX}/rrna.raw.sam
```
Filter and count alignment result.
```shell script
time gzip -dcf data/${PREFIX}/rrna.raw.sam.gz |
  parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j 6 '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl rrna_analyse_scripts/rrna_judge.pl |
    perl rrna_analyse_scripts/rrna_more_judge.pl
  ' \
  >data/${PREFIX}/rrna.out.sam

time parallel -j 3 \
  perl rrna_analyse_scripts/rrna_count.pl \
    data/ath_rrna/{}.fa data/${PREFIX}/rrna.out.sam {} \
    >data/${PREFIX}/rrna_{}.tsv \
  ::: 25s 18s 5-8s
```




