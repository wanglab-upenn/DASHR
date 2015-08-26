INBG=$1 # input bedgraph

cat ${INBG} | \
awk 'BEGIN{OFS="\t"; chr="";minHeight=10; minRun=14; minFoldChange=2}
     {
        val = $4; 
        if (chr!=$1 || ($2>prev_end))
        {
            if (inPeak==1)
            {
              b = prev_end;
              if ((b-a+1) >= minRun)
                print chr, a, b;
            }
            inPeak=-1;
            chr = $1;
        }
        
        fc = 1;
        if (val_prev>0)
          fc = val/val_prev;
        
        if ( (val < minHeight) || (fc <= 1/minFoldChange))
        #if ( val < minHeight )
        {
           if (inPeak==1)
           {
             b= prev_end; #$2;
             if ((b-a+1) >= minRun)
               print chr,a,b
               
           }
           inPeak = -1;
           prev_end = $3;   
           val_prev = $4;
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
             } 
             #a=$2;
           }
            #inPeak = 1;
        }
        prev_end = $3;   
        val_prev = $4;
     }' 
