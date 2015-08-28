# annotate segments / peaks

verbose=0
# input segmentation (both strands with segment IDs)
INSEGM=$1 
# col.5 = segment ID
# col.4 = tissue name
ANNOTS1=annot/hsa19.fulltable.no_mRNA_no_lncRNA_no_piRNA_nomiRP.unique_LOC.bed
ANNOTS2=annot/hsa19.fulltable.miRPonly.unique_LOC.bed
ANNOTS3=annot/hsa19.fulltable.no_mRNA_no_lncRNA_piRNAonly.unique_LOC.bed
ANNOTS4=annot/hsa19.fulltable.no_mRNA_no_lncRNA_rRNAonly.unique_LOC.bed
ANNOTS5=annot/hsa19.fulltable.no_mRNA_no_lncRNA.unique_LOC.bed
#INSEGMBASE=`basename ${INSEGM}`
INSEGMBASE=${INSEGM}
if [ ${verbose} -eq 1 ]; then echo "${INSEGMBASE}"; fi
# step1.
# 1.1 overlap with miRNA-5p/3p/5p3pno, snoRNA, snRNA, scRNA, tRF3, tRF5, tRNA

bedtools intersect -a ${INSEGM}  -b ${ANNOTS1} -wao -s -f 0.8 > ${INSEGMBASE}.step1.bed

# 1.2 col. 8 = -1 put for the next step
numFields=$(awk '{print NF; exit}' ${INSEGM})

if [ ${verbose} -eq 1 ]; then echo "Number of fields: ${numFields}"; fi
overlapField=$(awk 'BEGIN{print '${numFields}'+7;}')
if [ ${verbose} -eq 1 ]; then echo "Overlap field: ${overlapField}"; fi

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.forstep2.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          overlapField=numInputFields+7;
          annotChrStart=numInputFields+2;
          annotChrEnd=numInputFields+3; 
          $overlapField=$overlapField/($annotChrEnd-$annotChrStart);
          print;
       }
     }' ${INSEGMBASE}.step1.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr | \
           awk 'BEGIN{OFS="\t"}{
                  segmID=$5;
                  if (prev_segmID != segmID)
                  {
                     print;
                  }
                  prev_segmID = $5;
                }' > ${INSEGMBASE}.step1.annot.out

# step 2.
# overlap with miRNA precursor

bedtools intersect -a ${INSEGMBASE}.forstep2.bed -b ${ANNOTS2} -wao -s -f 0.8 > ${INSEGMBASE}.step2.bed 

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.forstep3.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          overlapField=numInputFields+7;
          annotChrStart=numInputFields+2;
          annotChrEnd=numInputFields+3; 
          $overlapField=$overlapField/($annotChrEnd-$annotChrStart);
          print;
       }
     }' ${INSEGMBASE}.step2.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr > ${INSEGMBASE}.step2.annot.out 
 
# step 3.
# overlap with piRNA

bedtools intersect -a ${INSEGMBASE}.forstep3.bed -b ${ANNOTS3} -wao -s -f 0.8 > ${INSEGMBASE}.step3.bed 

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.forstep4.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          #overlapField=numInputFields+7;
          #annotChrStart=numInputFields+2;
          #annotChrEnd=numInputFields+3; 
          #$overlapField=$overlapField/($annotChrEnd-$annotChrStart);
          print;
       }
     }' ${INSEGMBASE}.step3.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr | \
     awk 'BEGIN{OFS="\t"}{
                  segmID=$5;
                  if (prev_segmID != segmID)
                  {
                     print;
                  }
                  prev_segmID = $5;
                }' > ${INSEGMBASE}.step3.annot.out

# step 4.
# annotate rRNA

bedtools intersect -a ${INSEGMBASE}.forstep4.bed -b ${ANNOTS4} -wao -s -f 0.8 > ${INSEGMBASE}.step4.bed 

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.forstep5.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          #overlapField=numInputFields+7;
          #annotChrStart=numInputFields+2;
          #annotChrEnd=numInputFields+3; 
          #$overlapField=$overlapField/($annotChrEnd-$annotChrStart);
          print;
       }
     }' ${INSEGMBASE}.step4.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr | \
           awk 'BEGIN{OFS="\t"}{
                  segmID=$5;
                  if (prev_segmID != segmID)
                  {
                     print;
                  }
                  prev_segmID = $5;
                }' > ${INSEGMBASE}.step4.annot.out


# step 5. 
# re-annotate remaining segments/peaks without restriction on the overlapping criteria
bedtools intersect -a ${INSEGMBASE}.forstep5.bed -b ${ANNOTS5} -wao -s > ${INSEGMBASE}.step5.bed 

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.unannotated.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          overlapField=numInputFields+7;
          inputChrStart=2;
          inputChrEnd=3; 
          $overlapField=$overlapField/($inputChrEnd-$inputChrStart);
          print;
       }
     }' ${INSEGMBASE}.step5.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr | \
           awk 'BEGIN{OFS="\t"}{
                  segmID=$5;
                  if (prev_segmID != segmID)
                  {
                     print;
                  }
                  prev_segmID = $5;
                }' > ${INSEGMBASE}.step5.annot.out

cat ${INSEGMBASE}.*.annot.out > ${INSEGMBASE}.annot.final
rm ${INSEGMBASE}*step*.*
