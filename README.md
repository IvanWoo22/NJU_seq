# 2OMG
 A brand new version of MgR_seq.

## 1. Prepare
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

#### Select data for analysing.
Get the sequencing clean data from `MgR_data`.  
```shell script
ID=NJU45
PREFIX=Ath_stem_NC
ln -sf /home/wyf/MgR_data/${ID}/R1.fq.gz data/${PREFIX}/R1.fq.gz
ln -sf /home/wyf/MgR_data/${ID}/R2.fq.gz data/${PREFIX}/R2.fq.gz
```
*The representation of `ID` can be found in [`sample_list.csv`](/sample_list.csv).*

## 2. Quality Control
#### Input a `FastQ` file or a `GZ` file of `FastQ` and then get some quality information.
```shell script

```

## 3. Build Index
#### Download reference.
Get reference sequence of species from [GENCODE](https://www.gencodegenes.org/) and [Ensembl](http://plants.ensembl.org/Arabidopsis_thaliana/Info/Index?db=core).
```shell script
# GENCODE release 33 for human
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz -O data/hsa.gff3.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz -O data/hsa.fa.gz

# GENCODE release M24 for mouse
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gff3.gz -O data/mmu.gff3.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.transcripts.fa.gz -O data/mmu.fa.gz

# Ensembl release 46 for Arabidopsis thaliana
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz -O data/ath.gff3.gz
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz -O data/ath.fa.gz
```

#### 
Make index for mapping.
```shell script
# rRNA index
cat data/ath_rrna/* > data/ath_rrna.fa
bowtie2-build data/ath_rrna.fa index/ath_rrna

# Protein coding mRNA index
pigz -dc data/gencode.vM24.transcripts.fa.gz |
  perl ~/perbool/fetch_fasta.pl \
  -s "transcript_biotype:protein_coding" --stdin \
  >data/mmu_protein_coding.fa
bowtie2-build data/mmu_protein_coding.fa index/mmu_protein_coding
```




