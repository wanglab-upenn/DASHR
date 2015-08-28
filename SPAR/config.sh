#SPAR config file

HOMEDIR=${HOME}

#absolute path to the bin directory
BINDIR=${HOMEDIR}/bin

#absolute path to the SPAR home directory
SPARPATH=${BINDIR}/SPAR

#absolute path to the SPAR output directory
SPARDIR=${HOMEDIR}/SPAR_out

#absolute path to pre-installed STAR, samtools, AWK, etc
STAR=${BINDIR}/STAR-STAR_2.4.0k/bin/Linux_x86_64/STAR # STAR
SAMTOOLS=${BINDIR}/samtools-1.2/samtools # SAMTOOLS
GAWK=${BINDIR}/gawk-4.1.0/gawk
BGTOBIGWIG=${BINDIR}/bedGraphToBigWig

#absolute path to the STAR genome index
genomeDir=${HOMEDIR}/datasets/hg19/star_2.4.0k/  # STAR genome index

#hg 19 chromosome information file
chromInfo=${SPARPATH}/annot/chromInfo.txt

#mapping parameters for STAR
maxMismatchCnt=0
maxMapCnt=100
minMappedLength=14


function printT
{
  echo "`date +'%b %d %H:%M:%S'` ..... $1"
}
