
bgfile=$1

source config.sh
outbigwig=${bgfile/bedgraph/bigWig}
${BGTOBIGWIG} ${bgfile} ${chromInfo} ${outbigwig}
