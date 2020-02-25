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

#### An overview of the data.
Get the sequencing clean data from `MgR_data`.  
The representation of `ID` can be found in [`sample_list.md`](/sample_list.md).
```shell script
ln -sf /home/wyf/MgR_data/${ID}/R1.fq.gz data/${PREFIX}/R1.fq.gz
ln -sf /home/wyf/MgR_data/${ID}/R2.fq.gz data/${PREFIX}/R2.fq.gz
```

#### Get reference and index done.
```shell script
# gencode release 33 for human
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz

# gencode release M24 for mouse
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gff3.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.transcripts.fa.gz
```






## 2. Quality Control
#### Input a `FastQ` file or a `GZ` file of `FastQ` and then get some quality information.
```shell script


```