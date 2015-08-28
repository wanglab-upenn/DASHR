INBG=$1 # input bedgraph

strand=$2 # input strand: forward: "pos"; reverse: "neg"

#SEGMBG=${INBG}.segm

processBedgraph="cat ${INBG}"
postprocessBedgraph="cat"

if [ "${strand}" == "neg"  ]; then
  #revBg="'BEGIN{OFS=\"\t\"}{chr = \$1; chrStart =\$2+0; chrEnd=\$3+0; val=\$4; print chr, -chrEnd, -chrStart, val;}' | sort -k1,1 -k2,2n -k3,3n"
  #revBg="'BEGIN{OFS=\"\t\"}{chr = \$1; chrStart=\$2 + 1; chrEnd=\$3 - 1.0; val=\$4; print chr, -chrEnd, -chrStart, val;}' | sort -k1,1 -k2,2n -k3,3n"
  revBg="'BEGIN{OFS=\"\t\"}{chr = \$1; chrStart=\$2; chrEnd=\$3; val=\$4; print chr, -chrEnd, -chrStart, val;}' | sort -k1,1 -k2,2n -k3,3n"
  processBedgraph="cat ${INBG} | awk ${revBg}"
  #postprocessBedgraph="awk 'BEGIN{OFS=\"\t\"}{print \$1, -\$3+0, -\$2+0}' | sort -k1,1 -k2,2n -k3,3n";
  #postprocessBedgraph="awk 'BEGIN{OFS=\"\t\"}{chrStart=-\$3; chrEnd=-\$2+1; print \$1, chrStart, chrEnd, \$4, \$5, \$6}' | sort -k1,1 -k2,2n -k3,3n";
  postprocessBedgraph="awk 'BEGIN{OFS=\"\t\"}{chrStart=-\$3; chrEnd=-\$2+1; \$3=chrEnd; \$2=chrStart; print \$0}' | sort -k1,1 -k2,2n -k3,3n";
fi
 
#echo "${processBedgraph}"
#echo "${postprocessBedgraph}"

eval ${processBedgraph} | \
awk 'function computeEntropy(ent, cntArray)
     {
        numPos=0; 
        total=0;
        for (p in cntArray)
        {
           total+=cntArray[p]; # total number of reads starting
           numPos+=1; # number of distinct starting positions
        }
        if (total == 0 || numPos == 0)
        {
           printf "how did this happen? total=%f, numPos=%d\n", total, numPos;
           exit
        }

        entropy=0;
        maxProb=0;
        for (p in cntArray)
        {
          prob = cntArray[p]/total;
          entropy+=-prob*log(prob);
          if (prob>maxProb)
            maxProb = prob
        }
        entropy = entropy / log(2);
        maxEntropy = log(numPos) / log(2);
        normalizedEntropy = entropy;
        if (numPos>1)
          normalizedEntropy = entropy / maxEntropy;
        end
        split("",ent,":");
        ent[1] = entropy;
        ent[2] = maxEntropy;
        ent[3] = normalizedEntropy;
        ent[4] = total
        ent[5] = numPos
        ent[6] = maxProb
     }
     BEGIN{OFS="\t"; chr="";minHeight=10; minRun=14; minFoldChange=2; peakCnt=0;strand="+"; if ("'${strand}'"=="neg") strand="-"; verbose=0;}
     {
        val = $4;
        #print "BG",$0;
        if (chr!=$1 || ($2>prev_end)) # next chromosome or gap (start > previous end)
        {
            #print "gap,inpeak="inPeak
            if (inPeak==1)
            {
              b = prev_end;
              if ((b-a+1) >= minRun)
              {
                peakCnt+=1;
                peakID=("P" peakCnt);

               # compute entropy
               computeEntropy(ent, cnt5p);
               entropy5p = ent[1];
               normEntropy5p = ent[3];
               maxEntropy5p = ent[2];
               maxProb5p = ent[6]               
               cnt3p[b-a+1]+=val_prev-0;
               computeEntropy(ent, cnt3p);
               entropy3p = ent[1];
               normEntropy3p = ent[3];
               maxEntropy3p = ent[2];
               maxProb3p = ent[6]               
               #print chr,a,b,peakID,exprVal,strand

               print chr,a,b,peakID,exprVal,strand,entropy5p,normEntropy5p,maxEntropy5p, entropy3p, normEntropy3p, maxEntropy3p, maxProb5p, maxProb3p;
               if (verbose==1)
               {
               len = b-a+1+0.0
               for (p in cnt5p)
               { if (p+0.0>len)
                    printf "\nhow did this happen? p=%d, len=%f\n", p,len;
                 printf "%d:%f\t", p, cnt5p[p];
               }
               printf "\n";
               printf "numPos=%d,total=%f\n", ent[5], ent[4]

               for (p in cnt3p)
               {
                 if (p+0.0>len)
                    printf "\nhow did this happen? p=%d,len=%f\n", p,len;
                 printf "%d:%f\t", p, cnt3p[p];
               }
               printf "\n";
               }


               #print chr, a, b, peakID, exprVal, strand;
               #print chr, a, b;
              }
            }
            split("",cnt5p,":");
            split("",cnt3p,":");
            inPeak=-1;
            chr = $1;
        }
        
        fc = 1;
        if (val_prev>0)
          fc = val/val_prev;
        
        if ( (val < minHeight) || (fc <= 1/minFoldChange))
        #if ( val < minHeight )
        {
           #print "here, inPeak="inPeak
           if (inPeak==1)
           {
             b=prev_end; #$2;
             if ((b-a+1) >= minRun)
             {
               #print chr,a,b
               peakCnt+=1;
               peakID=("P" peakCnt);
               # compute entropy
               computeEntropy(ent, cnt5p);
               entropy5p = ent[1];
               normEntropy5p = ent[3];
               maxEntropy5p = ent[2];
               maxProb5p = ent[6];
               cnt3p[$2-a+1]+=val_prev-val;
               computeEntropy(ent, cnt3p);
               entropy3p = ent[1];
               normEntropy3p = ent[3];
               maxEntropy3p = ent[2];
               maxProb3p = ent[6] 
               #print chr,a,b,peakID,exprVal,strand
               print chr,a,b,peakID,exprVal,strand,entropy5p,normEntropy5p,maxEntropy5p, entropy3p, normEntropy3p, maxEntropy3p,maxProb5p,maxProb3p;
               if (verbose==1) {
               len = b-a+1;
               for (p in cnt5p)
               {
                  if (p+0.0>len)
                    printf "\nhow did this happen? p=%d, len=%f\n", p,len;

                 printf "%d:%f\t", p, cnt5p[p];
               }
               printf "\n";
               for (p in cnt3p)
               {
                  if (p+0.0>len)
                    printf "\nhow did this happen? p=%d, len=%f\n", p,len;

                 printf "%d:%f\t", p, cnt3p[p];
               }
               printf "\n";
               }
               #print "peak end:", b, val, val_prev;  
             }
               
           }
           inPeak = -1;
           prev_end = $3;   
           val_prev = $4;
           split("",cnt5p,":");
           split("",cnt3p,":");
           next;
        }
        if (val >= minHeight)
        {  
           if (inPeak<0)
           {
             # fold-change
             fc = minFoldChange;
             if (val_prev > 0)
               fc = val/val_prev;

             if (fc >= minFoldChange)
             {
                a=$2;
                inPeak = 1;
                #print "peak start:", a, val, val_prev
                exprVal = val;
             } 
             #a=$2;
           }
            #inPeak = 1;
        }
        if (val > exprVal) exprVal=val; 
        if (inPeak==1)
        {
          currReadStartCnt=val-val_prev;
          if (currReadStartCnt>0)
            cnt5p[$2-a+1]+=currReadStartCnt; # number of read starts per position
          if (currReadStartCnt<0)
            cnt3p[$2-a+1]+=-currReadStartCnt; # number of read ends per position
        }
        prev_end = $3;   
        val_prev = $4;
     }' | eval ${postprocessBedgraph} > ${INBG}.segm 

