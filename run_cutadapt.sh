#!/bin/bash
# export PYTHONPATH=$PYTHONPATH:$HOME/lib/python/:$HOME/lib64/python/:$HOME/lib/python2.6/site-packages/
export PYTHONPATH=$PYTHONPATH:$HOME/lib/python/:$HOME/lib64/python/:/home/pry/lib/python/
#/genomics/share/python2.6.5/bin/python $HOME/bin/cutadapt $@
#python2.6 $HOME/bin/cutadapt $@
$HOME/bin/cutadapt-1.8.1/bin/cutadapt $@

STEPNAME=trim3

# read config
source common2.config
source common.sh
source $CONFIG

echo $FAS_3ADAPTER
TRIMMEDREADS=${DATASET}_trimmed.fastq
UNTRIMMEDREADS=${DATASET}_untrimmed.fastq
TOOSHORTREADS=${DATASET}_tooshort.fastq
RAWREADS=${DATASET}_rawreads.fastq
ADAPTER3SEQ=`tail -n1 $FAS_3ADAPTER`
qsub $QSUBOPTS -N trim.$$ $CODE/run_cutadapt.sh \
  -e $MAX_3ADAPTER_ERROR -o $TRIMMEDREADS --untrimmed-output $UNTRIMMEDREADS \
  --too-short-output $TOOSHORTREADS -a $ADAPTER3SEQ -O $MIN_3ADAPTER_MATCH \
  -m $MIN_TRIMMED_READ_LEN -n 2 --length-tag "length=" $RAWREADS

qsub $QSUBOPTS -hold_jid trim.$$ $CODE/trim_stats_fastq.sh \
  $TOOSHORTREADS $TRIMMEDREADS $UNTRIMMEDREADS ${DATASET}_trim3_stats.txt

