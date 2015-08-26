BGTOBIGWIG=${HOME}/bin/bedGraphToBigWig
chromInfo=/mnt/niagads/users/yyee/paulrnaseqpipe2/coral-1.1.0/bin/chromInfo.txt
BGLIST=$1 # list of begraph files

# process positive strand
cat ${BGLIST} | \
while read bgfile; do
  tissue=${bgfile##*/} # remove directory
  tissue=${tissue%%_*} # remove suffix
  outbigwig=${bgfile/bedgraph/bigWig}
  ${BGTOBIGWIG} ${bgfile} ${chromInfo} ${OUTDIR}/${outbigwig}
  echo ${outbigwig}
done
