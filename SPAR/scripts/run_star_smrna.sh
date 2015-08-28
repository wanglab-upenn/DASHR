# map smRNA reads using STAR
# Usage:
# bash run_star_smRNA.sh <full path to sample FASTQ file>
# Example: 
# bash run_star_smrna.sh /mnt/niagads/users/yyee/tissues/hsa_adpilot_smrna/hsa_adpilot_ctl02058_smrna/1.5/hsa_adpilot_ctl02058_smrna_trimmed.fastq



source config.sh

# STAR output directory (absolute path)
starOutHome=${SPARDIR}/STAR

# Mapping parameters
#multiMapMax=100 
#maxMismatchCnt=0
#minMappedLength=14

FASTQ=$1 # input FASTQ file
         # <tissue>_<SRR>_trimmed.fastq
maxMismatchCnt=$2 #${maxMismatchCnt}
multiMapMax=$3 #${maxMapCnt}
starOutDir=$4
# get tissue and sample name from FASTQ file name
#tissueName=${FASTQ##*/}
#tissueName=${tissueName%%_*}
#sampleName=${FASTQ##*/}
#sampleName=${sampleName#*_}
#sampleName=${sampleName%%_*}

# set output directory as
#   starOutHome/star_m<maxMismatchCnt>_map<multiMapMax>_<tissue>_<sample>
#starOutDir=${starOutHome}/STAR_m${maxMismatchCnt}_map${multiMapMax}_${tissueName}_${sampleName} # output directory (absolute path)
mkdir -p ${starOutDir}

#echo -e "Tissue: ${tissueName}"
#echo -e "Sample: ${sampleName}"
#echo -e "STAR output:\n${starOutDir}"

# run STAR
${STAR} --genomeDir ${genomeDir} --genomeLoad LoadAndKeep  --readFilesIn ${FASTQ} --runThreadN 4 --alignIntronMax 1 --outSAMattributes NH HI NM MD --outFilterMultimapNmax ${multiMapMax} --outReadsUnmapped Fastx --outFilterMismatchNmax ${maxMismatchCnt} --outFilterMatchNmin ${minMappedLength} --outFileNamePrefix ${starOutDir}/


# filter out reads trimmed at 5p end by more than one base
#echo "`date +'%b %d %H:%M:%S'` ..... Removing reads trimmed at 5p"
printT "Removing reads trimmed at 5p"
awk 'BEGIN{OFS="\t"; filtPos=0; filtNeg=0; alnCnt=0; totalAln}
     {
       if (/^@/) {print; next}; # SAM header
       samFlag=$2; cigar=$6; strandFlag=and(samFlag,16);
       totalAln+=1
       if (strandFlag>0) # negative strand
       {
         if (cigar~/[0-9]S$/)
         {
           n=split(cigar,L,/[NMSID]/);
           clip=L[n-1]+0;
           if (clip>1) {filtNeg+=1; next}
         }
       }
       else # positive strand
       {
         if (cigar~/^[0-9]+S/)
         {
           split(cigar,L,/[NMSID]/);
           clip=L[1];
           if (clip>1) {filtPos+=1; next}
         }
       }
       alnCnt+=1;
       print
     }
     END {
      filter5pStatsFile=(FILENAME ".filter_5p_stats");
      print "Alignments before removing 5p clipping", totalAln, totalAln/totalAln*100 > filter5pStatsFile
      print "Removed 5p clipped [positive strand]:",filtPos,filtPos/totalAln*100 > filter5pStatsFile;
      print "Removed 5p clipped [negative strand]:", filtNeg,filtNeg/totalAln*100 > filter5pStatsFile;
      print "Total after removing 5p clipped:",alnCnt,alnCnt/totalAln*100 > filter5pStatsFile;
         }' ${starOutDir}/Aligned.out.sam > ${starOutDir}/Aligned.out.filtered.sam

# convert filtered SAM to BAM file format
#echo "`date +'%b %d %H:%M:%S'` ..... Sorting filtered SAM"
printT "Sorting filtered SAM"
${SAMTOOLS} sort -@ 4 -O bam ${starOutDir}/Aligned.out.filtered.sam -T ${starOutDir}/Aligned.out.filtered.sort > ${starOutDir}/Aligned.out.filtered.sorted.bam

# remove soft-clipped reads
echo "`date +'%b %d %H:%M:%S'` ..... Removing soft-clipped reads"
${SAMTOOLS} view -h ${starOutDir}/Aligned.out.filtered.sorted.bam | awk 'BEGIN{OFS="\t"}{ if ($0~/^@/) {print; next} if ($6~/S/) {softClipCnt+=1; next;} noSoftClipCnt+=1; print }END{print "Alignments [no soft clipping]", noSoftClipCnt > "'${starOutDir}'/softClipStats.txt"; print "Alignments [with soft clipping]", softClipCnt > "'${starOutDir}'/softClipStats.txt"}' | ${SAMTOOLS} view -bS - > ${starOutDir}/Aligned.out.filtered.noS.sorted.bam



# hard-clip SAM file
#echo "`date +'%b %d %H:%M:%S'` ..... Hard-clipping filtered SAM"
printT "Hard-clipping filtered SAM"
awk 'BEGIN {OFS="\t"}{split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' ${starOutDir}/Aligned.out.filtered.sam > ${starOutDir}/Aligned.out.filtered.hardClipped.sam

# hard-clip BAM file
#echo "`date +'%b %d %H:%M:%S'` ..... Sorting hard-clipped filtered SAM"
printT "Sorting hard-clipped filtered SAM"
${SAMTOOLS} sort -@ 4 -O bam ${starOutDir}/Aligned.out.filtered.hardClipped.sam -T ${starOutDir}/Aligned.out.filtered.hardClipped.sort > ${starOutDir}/Aligned.out.filtered.hardClipped.sorted.bam


# index BAM
#echo "`date +'%b %d %H:%M:%S'` ..... Indexing BAM"
printT "Indexing BAM"
${SAMTOOLS} index ${starOutDir}/Aligned.out.filtered.hardClipped.sorted.bam

${SAMTOOLS} index ${starOutDir}/Aligned.out.filtered.noS.sorted.bam

${SAMTOOLS} index ${starOutDir}/Aligned.out.filtered.sorted.bam

# compute stats
#echo "`date +'%b %d %H:%M:%S'` ..... Computing stats"
printT "Computing mapping stats"
statFile=${starOutDir}/MAPSTAT.txt
>${statFile}
nFastQReads=$(grep -e "Number of input reads" ${starOutDir}/Log.final.out | cut -f 2)
nUniqMapReads=$(grep -e "Uniquely mapped reads number" ${starOutDir}/Log.final.out | cut -f 2)
nMultiMapReads=$(grep -e "Number of reads mapped to multiple loci" ${starOutDir}/Log.final.out | cut -f 2)


echo "${nFastQReads} ${nUniqMapReads} ${nMultiMapReads}" | \
   awk 'BEGIN{OFS="\t"}
        {
         print "FASTQ reads:", $1, $1/$1*100
         print "Uniquely-mapped reads:", $2, $2/$1*100
         print "Multi-mapped reads:", $3, $3/$1*100
         print "Total mapped reads:", $2+$3, ($2+$3)/$1*100
        }' >> ${statFile}
cat ${starOutDir}/*.filter_5p_stats ${starOutDir}/softClipStats.txt >> ${statFile}

nNoSReads=$(${SAMTOOLS} view ${starOutDir}/Aligned.out.filtered.noS.sorted.bam | cut -f 1 | sort -u | wc -l)
nSandNoSReads=$(${SAMTOOLS} view ${starOutDir}/Aligned.out.filtered.sorted.bam | cut -f 1 | sort -u | wc -l)
echo -e "${nNoSReads} ${nSandNoSReads}" | \
awk 'BEGIN{OFS="\t"}
     {
       print "Reads [no soft clipping]:", $1, $1/'${nFastQReads}'*100
       print "Reads [with soft clipping]:", $2-$1, ($2-$1)/'${nFastQReads}'*100
       print "Reads [all]:", $2, $2/'${nFastQReads}'*100
     }' >> ${statFile}

# remove SAM files; keep BAMs only
rm ${starOutDir}/*.sam

#printT "DONE."

# ENCODE param
#minMappedLength=16
#${STAR} --genomeDir ${genomeDir}  --readFilesIn ${FASTQ} --runThreadN 4 --alignIntronMax 1 --outSAMattributes All --outFilterMultimapNmax ${multiMapMax} --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin ${minMappedLength} --outFileNamePrefix ${starOutDir}
