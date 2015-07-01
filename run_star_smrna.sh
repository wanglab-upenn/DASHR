# map smRNA reads using STAR
# Usage:
# bash run_star_smRNA.sh <full path to sample FASTQ file>

STAR=/home/pkuksa/bin/STAR-STAR_2.4.0k/bin/Linux_x86_64/STAR
genomeDir=/home/pkuksa/datasets/hg19/star_2.4.0k/  # STAR genome index
SAMTOOLS=/home/pkuksa/bin/samtools-1.1/samtools

# STAR output directory (absolute path)
starOutHome=

# Mapping parameters
multiMapMax=100 #5
maxMismatchCnt=0
minMappedLength=14

# Filtering parameters
filterMismatchCntMax=5

FASTQ=$1 # input FASTQ file

# get tissue and sample name from FASTQ file name
tissueName=${FASTQ##*/}
tissueName=${tissueName%%_*}
sampleName=${FASTQ##*/}
sampleName=${sampleName#*_}
sampleName=${sampleName%%_*}

# set output directory as
#   starOutHome/star_m<maxMismatchCnt>_map<multiMapMax>_<tissue>_<sample>
starOutDir=${starOutHome}/star_m${maxMismatchCnt}_map${multiMapMax}_${tissueName}_${sampleName} # output directory (absolute path)
mkdir -p ${starOutDir}

echo -e "Tissue: ${tissueName}"
echo -e "Sample: ${sampleName}"
echo -e "STAR output:\n${starOutDir}"

# run STAR
${STAR} --genomeDir ${genomeDir} --genomeLoad LoadAndKeep  --readFilesIn ${FASTQ} --runThreadN 4 --alignIntronMax 1 --outSAMattributes NH HI NM MD --outFilterMultimapNmax ${multiMapMax} --outReadsUnmapped Fastx --outFilterMismatchNmax ${maxMismatchCnt} --outFilterMatchNmin ${minMappedLength} --outFileNamePrefix ${starOutDir}/


# filter out reads trimmed at 5p end by more than one base
echo "Removing reads trimmed at 5p"
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
         if (cigar~/^[0-9]S/)
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
      print "Alignments before removing 5p clipping", totalAln, totalAln/totalAln*100 > FILENAME".filter_5p_stats";
      print "Removed 5p clipped [positive strand]:",filtPos,filtPos/totalAln*100 > FILENAME".filter_5p_stats";
      print "Removed 5p clipped [negative strand]:", filtNeg,filtNeg/totalAln*100 > FILENAME".filter_5p_stats";
      print "Total after removing 5p clipped:",alnCnt,alnCnt/totalAln*100 > FILENAME".filter_5p_stats"
         }' ${starOutDir}/Aligned.out.sam > ${starOutDir}/Aligned.out.filtered.sam

# filter out reads mapping to too many places
echo "Removing reads mapped to too many places"
awk 'BEGIN{OFS="\t"; maxHits="'${filterMismatchCntMax}'"+0; alnCnt=0;}
     {
       if ($0~/^@/) {print; next}
       totalAln+=1;
       nHits = substr($12, 6, 3)+0;
       if (nHits <= maxHits) {alnCnt+=1;print;}
     }END{
       print "Alignments before filtering by number of hits:",totalAln,totalAln/totalAln*100 > FILENAME".filter_hit_stats"
       print "Alignments [maxHits="maxHits"]:",alnCnt,alnCnt/totalAln*100 > FILENAME".filter_hit_stats"
         }' ${starOutDir}/Aligned.out.filtered.sam > ${starOutDir}/Aligned.out.filtered.m${filterMismatchCntMax}.sam  

# convert filtered SAM to BAM file format
#${SAMTOOLS} view -bS ${starOutDir}/Aligned.out.filtered.m${filterMismatchCntMax}.sam > ${starOutDir}/Aligned.out.filtered.m${filterMismatchCntMax}.bam
echo "Sorting filtered SAM"
${SAMTOOLS} sort -@ 4 -O bam ${starOutDir}/Aligned.out.filtered.m${filterMismatchCntMax}.sam -T ${starOutDir}/Aligned.out.filtered.m${filterMismatchCntMax}.sort > ${starOutDir}/Aligned.out.filtered.m${filterMismatchCntMax}.sorted.bam
#samtools sort -@ 4 -O bam Aligned.out.filtered.sam -T Aligned.out.filtered.sort > Aligned.out.filtered.sorted.bam

