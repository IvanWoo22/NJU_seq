#!/usr/bin/env bash
CUTADAPT=0
PE=0
R1DEF="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
R2DEF="GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT"
SUFFIX=unset
OUTDIR="."
LOGDIR="."

usage() {
  echo "Usage: Raw_Data_Preparation.sh [ -c | --cutadapt ] [ -p | --pe_consistency ]
                        [ -a R1INDEX ] [ -A R2INDEX ]
                        [ --R1 filename ] [--R2 filename]
                        [ -o | --outdir outdir ] [ --suffix | -s text ]
                        [ -l | --logdir logdir ]"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -a -n Raw_Data_Preparation -o a:A:cdp -l cutadapt,pe_consistency,R1:,R2: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :; do
  case "$1" in
  -c | --cutadapt)
    CUTADAPT=1
    shift
    ;;
  -p | --pe_consistency)
    PE=1
    shift
    ;;
  -a)
    R1INDEX="$2"
    shift 2
    ;;
  -A)
    R2INDEX="$2"
    shift 2
    ;;
  --R1)
    R1="$2"
    shift 2
    ;;
  --R2)
    R2="$2"
    shift 2
    ;;
  -o | --outdir)
    OUTDIR="$2"
    shift 2
    ;;
  --suffix | -s)
    SUFFIX="$2"
    shift 2
    ;;
  -l | --logdir)
    LOGDIR="$2"
    shift 2
    ;;
  --)
    shift
    break
    ;;
  *)
    echo "Unexpected option: $1 - this should not happen."
    usage
    ;;
  esac
done

if [[ -z "$R1" || -z "$R2" ]]; then
  echo "You must provide a pair of fastq files!"
  exit 2
fi

if [ ! -d "$LOGDIR" ]; then
  mkdir "$LOGDIR"
fi

if [ "$CUTADAPT" -ne 0 ]; then
  if [ -z "$R1INDEX" ]; then
    R1INDEX="$R1DEF"
    echo "Use default R1 index: $R1DEF"
  fi
  if [ -z "$R2INDEX" ]; then
    R2INDEX="$R2DEF"
    echo "Use default R2 index: $R2DEF"
  fi
  if [ -z "$R2INDEX" ]; then
    R1outfile="$OUTDIR/R1.fq.gz"
    R2outfile="$OUTDIR/R2.fq.gz"
  else
    R1outfile="${OUTDIR}/${SUFFIX}_R1.fq.gz"
    R2outfile="${OUTDIR}/${SUFFIX}_R2.fq.gz"
  fi
  if [[ -f "$R1outfile" || -f "$R2outfile" ]]; then
    echo "There were files with same output file names!"
    exit 2
  fi
  cutadapt -a $R1INDEX -A $R2INDEX \
    -O 6 -m 10 -e 0.1 --discard-untrimmed -o "${R1outfile}" -p "${R2outfile}" \
    "${R1}" "${R2}" -j 20 2>"${LOGDIR}"/"${SUFFIX}"_cutadapt.log
fi

if [ ! -d "$OUTDIR" ]; then
  mkdir "$OUTDIR"
fi
