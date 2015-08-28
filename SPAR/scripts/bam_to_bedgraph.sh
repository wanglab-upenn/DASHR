source config.sh
BAM=$1
${SAMTOOLS} view ${BAM} | \
${GAWK} 'function cmp_num_idx(i1, v1, i2, v2)
      {
       # numerical index comparison, ascending order
       return (i1 - i2)
      }
      function process_chr_coverage( c, n, chr, outfile)
      {
         bs = 1; # block start
         cnt = 0;
         f = 0; # read coverage
         f_prev = 0;
         PROCINFO["sorted_in"]="@ind_num_asc"
         for (i in c)
         {
            ++cnt;
            if (verbose==1 && cnt % 100000 == 0) printf "Processed=%d\n",cnt;
            #if (cnt==1) {bs=i; f_prev = c[i];}
            f+=c[i]
            # print i
            # print i, c[i], f
            if (f != f_prev)
            {
             #if (f_prev > 0.0001)
             if (f_prev >= 1)
               print chr, bs-1, i-1, f_prev > outfile
             bs = i;
            }
            f_prev = f;
         }
      }
    BEGIN{OFS="\t";chrLen=0;verbose=0;}
     {
       if (NR==1) chr_prev = $3;
       chr = $3
       if (chr == chr_prev)
       {
         flag = $2 # SAM flag
         strand = "+";
         if (and(flag,16)>0) strand="-";
         rstart = $4+0 # read start
         rend = $4+length($10)-1+0 # read end
         nHits = substr($12,6);
         w = 1 / nHits;
         if (chrLen<rend) chrLen=rend;
         # print $0
         # print strand, rstart, rend, nHits, w
         if (strand=="+")
         {
             rcov_plus[rstart]+=w;
             rcov_plus[rend]-=w;
         }
         else
         {
             rcov_minus[rstart]+=w;
             rcov_minus[rend]-=w;
         } 
         if ( verbose==1 && (NR % 1000000) == 0 ) printf "Processed %d read alignments\n", NR;
       }
       else # process chromosome
       {
         n = chrLen; #chrLen[chr_prev]
         if (verbose==1) 
           printf "Processing %s [len=%d, Watson strand]\n", chr_prev, n
         outfile = "'${BAM}'.pos.bedgraph"
 
         process_chr_coverage( rcov_plus, n, chr_prev, outfile )     
         split("",rcov_plus,":") # clear coverage array
         if (verbose==1)
           printf "Processing %s [len=%d, Crick strand]\n", chr_prev, n
         outfile = "'${BAM}'.neg.bedgraph"
         process_chr_coverage( rcov_minus, n, chr_prev, outfile )     
         split("",rcov_minus,":") # clear coverage array
         chrLen = 0;
       }
       chr_prev = chr;        
     }
     END{
         # process last chromosome
         n = chrLen; #chrLen[chr_prev]
         if (verbose==1)
           printf "Processing %s [len=%d, Watson strand]\n", chr_prev, n
         outfile = "'${BAM}'.pos.bedgraph"
         process_chr_coverage( rcov_plus, n, chr_prev, outfile )     
         split("",rcov_plus,":") # clear coverage array
         if (verbose==1)
           printf "Processing %s [len=%d, Crick strand]\n", chr_prev, n
         outfile = "'${BAM}'.neg.bedgraph"
         process_chr_coverage( rcov_minus, n, chr_prev, outfile )     
         split("",rcov_minus,":") # clear coverage array
         chrLen = 0;
     }'

printT "BEDGRAPH finished successfuly."
