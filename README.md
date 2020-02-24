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
```shell script

```

#### Get reference and index done.
```shell script

```






## 2. Quality Control
#### Input a `FastQ` file or a `GZ` file of `FastQ` and then get some quality information.
```shell script


```