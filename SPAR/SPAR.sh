
source config.sh

set -e

INFASTQ=$1

if [ $# -lt 1 ]
then
  echo "USAGE: `basename $0` reads.fastq"
  exit 1
fi


if [ ! -f "${INFASTQ}" ]; then
  echo -e "*****ERROR: FASTQ file\n${INFASTQ}\nnot found!"
  exit 1
fi


TRIMFASTQ=${INFASTQ}


EXPNAME=`basename ${TRIMFASTQ}`
EXPNAME=${EXPNAME%_*}

OUTDIR=${SPARDIR}/STAR/STAR_m${maxMismatchCnt}_map${maxMapCnt}_${EXPNAME}
OUTDIR=${SPARDIR}/${EXPNAME}_m${maxMismatchCnt}_map${maxMapCnt}
mkdir -p ${OUTDIR}

LOGSPAR=${OUTDIR}/SPAR.log
>${LOGSPAR}

OUTBAM=${OUTDIR}/Aligned.out.filtered.hardClipped.sorted.bam

#echo "EXPNAME=${EXPNAME}"
#echo "${OUTBAM}"

function runScript
{
  bash ${SPARPATH}/scripts/$1
}

function printT
{
  echo "`date +'%b %d %H:%M:%S'` ..... $1" | tee -a ${LOGSPAR}
}

function printL
{
  echo -e "$1" | tee -a ${LOGSPAR}
}

printL "Output directory:\n${OUTDIR}\n"

printT "SPAR run started"

runScript "run_star_smrna2.sh ${TRIMFASTQ} ${maxMismatchCnt} ${maxMapCnt} ${OUTDIR}"

printT "Converting BAM to bedGraph"
runScript "bam_to_bedgraph.sh ${OUTBAM}"

printT "Converting bedGraph to bigWig"
runScript "bedgraph_to_bigwig.sh ${OUTBAM}.pos.bedgraph"
runScript "bedgraph_to_bigwig.sh ${OUTBAM}.neg.bedgraph"

printT "Segmenting [positive strand]"
runScript "segment_bedgraph_entropy.sh ${OUTBAM}.pos.bedgraph pos"
printT "Segmenting [negative strand]"
runScript "segment_bedgraph_entropy.sh ${OUTBAM}.neg.bedgraph neg"

printT "Annotating [positive strand]"
runScript "annotate_segm2.sh ${OUTBAM}.pos.bedgraph.segm"
printT "Annotating [negative strand]"
runScript "annotate_segm2.sh ${OUTBAM}.neg.bedgraph.segm"


printT "DONE."

cat ${OUTBAM}.*.bedgraph.segm.annot.final > ${OUTBAM}.annot.final
cat ${OUTBAM}.*.bedgraph.segm.unannotated.bed > ${OUTBAM}.unannot

awk 'BEGIN{OFS="\t";}{if (NR==FNR) {exprVal=$5; rnaClass=$19; exprPerClass[rnaClass]+=exprVal;classCnt[rnaClass]+=1;totalExprAnnot+=exprVal; totalAnnotPeakCnt+=1}else{exprVal=$5; totalExprUnannot+=exprVal;totalUnannotPeakCnt+=1;}}END{totalPeakCnt=totalAnnotPeakCnt+totalUnannotPeakCnt; totalExpr=totalExprAnnot+totalExprUnannot; for (rnaClass in exprPerClass) print rnaClass, classCnt[rnaClass], exprPerClass[rnaClass], exprPerClass[rnaClass]/totalExpr; print "Unannotated",totalUnannotPeakCnt,totalExprUnannot,totalExprUnannot/totalExpr}' ${OUTBAM}.annot.final ${OUTBAM}.unannot | sort -k1,1 | awk 'BEGIN{OFS="\t"; print "#RNA","Peaks","Reads","Fraction of reads"}{print}' > ${OUTDIR}/mapped_reads_annotation_summary.txt


printL "\n===Output==="
printL "Output directory: ${OUTDIR}"
printL "\nMapping output:"
printL "${OUTBAM}"

printL "Annotation output:"
ls ${OUTBAM}.*.bedgraph.segm.annot.final | tee -a ${LOGSPAR}

printL "Un-annotated output:" 
ls ${OUTBAM}.*.bedgraph.segm.unannotated.bed | tee -a ${LOGSPAR}

printL "\n\n===Run summary==="
printL "FASTQ: ${TRIMFASTQ}"
grep -e "Reads \[all\]" ${OUTDIR}/MAPSTAT.txt | awk 'BEGIN{FS="\t"}{printf "Mapped reads: %d [%.4f%%]\n", $2, $3}' | tee -a ${LOGSPAR}

numAnnot=$(cat ${OUTBAM}.*.bedgraph.segm.annot.final | wc -l)
printL "Annotated loci count: ${numAnnot}"
printL "Annotated loci by RNA class:"
cat ${OUTBAM}.*.bedgraph.segm.annot.final | cut -f 19 | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$1}' | sort -k2,2nr | tee -a ${LOGSPAR}

numUnannot=$(cat ${OUTBAM}.*.bedgraph.segm.unannotated.bed | wc -l)
printL "\nUn-annotated loci count: ${numUnannot}"


echo "Annotation summary:"
cat ${OUTDIR}/mapped_reads_annotation_summary.txt
