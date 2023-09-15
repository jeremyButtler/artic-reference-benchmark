refStr="01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta";
ampsStr="tmp.fasta";

while [[ $# -gt 0 ]]; do
# Loop: Get user input
   case $1 in
      -ref) refStr="$2"; shift;;
      -amplicons) ampsStr="$2"; shift;;
   esac

   shift;
done # Loop: Get user input

if [[ ! -f "$refStr" ]]; then
   printf " -ref %s is not a file\n" "$refStr" >&2;
   exit;
fi # check if have a valid reference file

if [[ -f "$2" ]]; then
   printf " -amplicons %s is not a file\n" "$ampsStr" >&2;
   exit;
fi # check if have a valid amplicon file

minimap2 \
    --eqx \
    -a \
    "$refStr" \
    "$ampsStr" |
  00-programs/trimSamFile -stdin |
  awk \
    '
       { # MAIN
          if($0 ~ /^@/){next;}; # sam header

          cigLineStr = $6;
          gsub(/[D]/, ",D,", cigLineStr);
          gsub(/[I]/, ",I,", cigLineStr);
          gsub(/[=]/, ",=,", cigLineStr);
          gsub(/[X]/, ",X,", cigLineStr);
          gsub(/[S]/, ",S,", cigLineStr);

          numCigI = split(cigLineStr, cigAryStr, ",");
          refEndI = $4;   # First base I am working on
          cigStr = "";

          for(iCig = 1; iCig <= numCigI; ++iCig)
          { # loop though each cigar entry
              if(cigAryStr[iCig + 1] == "I")
              { # If this is an indel entry
                 for(iIns=1;iIns <= cigAryStr[iCig];++iIns)
                 { # For all indel positions
                    cigStr = cigStr ",I";
                 } # For all indel positions
              } # If this is an indel entry

              else if(cigAryStr[iCig + 1] == "D")
              { # else if this is a deletion entry
                 for(iDel=1;iDel <= cigAryStr[iCig];++iDel)
                 { # For all indel positions
                    ++refEndI;
                    cigStr = cigStr ",D";
                 } # For all indel positions
              } # else if this is a deletion entry

              else if(cigAryStr[iCig + 1] == "X")
              { # else if this is an snp entry
                 for(iSNP=1;iSNP <= cigAryStr[iCig];++iSNP)
                 { # Loop: unexpand mismtach entries
                    ++refEndI;
                    cigStr = cigStr ",X";
                 } # Loop: unexpand mismtach entries
              } # else if this is an snp entry

              else
              { # else this is a match entry
                 for(iSNP=1;iSNP <= cigAryStr[iCig];++iSNP)
                 { # Loop: unexpand match entries
                    ++refEndI;
                    cigStr = cigStr ",=";
                 } # Loop: unexpand match entries
              } # else this is a match entry

              ++iCig; # Move off the number
          } # loop though each cigar entry

          # Remove the first ","
          cigStr = substr(cigStr, 2);

          printf ">%s %s %s", $1, $4, refEndI;
          printf " %s %s %s\n", length($10), $10, cigStr;
          # Name, first ref position, length, sequence
       } # MAIN
      ' |
  sort -V -k 2,3 |
  awk \
    '
      #****************************************************
      # Sec-0? Sub-02:
      #  - Main bock:
      #  o sec-0? sub-02 cat-01:
      #    - Process the first sequence
      #  o sec-0? sub-02 cat-02:
      #    - 
      #  o sec-0? sub-02 cat-0?:
      #    - Handle gaps between sequences
      #****************************************************

      #++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Sec-0? Sub-02 Cat-01:
      #  - Process the first sequence
      #++++++++++++++++++++++++++++++++++++++++++++++++++++

      { # MAIN
         if(NR == 1)
         { # If this is the frist sequence
            conStartI = $2;
            conRefEndI = $3;
            conLenSeqI = 1;
            conLenCigI = 1;

            # Convert sequence to c-string
            gsub(//, " ", $5);
            tmpSeqI = split($5, tmpSeqCStr, " ");

            # Convert expanded cigar to c-string
            tmpCigI = split($6, tmpCigCStr, ",");

            for(iN = 1; iN < conStartI; ++iN)
            { # Loop: Add ends for missing starting bases
               conSeqCStr[iN] = "N";
               conCigCStr[iN] = "D";
              ++conLenSeqI;
              ++conLenCigI;
            } # Loop: Add ends for missing starting bases

            for(iNt = 1; iNt <= tmpSeqI; ++iNt)
            { # Loop: Add in the first sequence
              conSeqCStr[iNt + conStartI]=tmpSeqCStr[iNt];
              ++conLenSeqI;
            } # Loop: Add in the first sequence

            for(iNt = 1; iNt <= tmpCigI; ++iNt)
            { # Loop: Add in the first cigar entry
              conCigCStr[iNt + conStartI]=tmpCigCStr[iNt];
              ++conLenCigI;
            } # Loop: Add in the first cigar entry

            ++conLenCigI;
            ++conLenSeqI;
            next;
         }; # If this is the frist sequence

         #+++++++++++++++++++++++++++++++++++++++++++++++++
         # Sec-0? Sub-02 Cat-02:
         #  - Process the next amplicon
         #+++++++++++++++++++++++++++++++++++++++++++++++++

         if(conRefEndI > $2)
         { # If I have an overlap
            # Convert the sequence to a c-string
            gsub(//, " ", $5);
            ampLenSeqI = split($5, ampSeqCStr, " ");

            # Convert expanded cigar to c-string
            ampLenCigI = split($6, ampCigCStr, ",");

            ampRefEndI = $3;
            ampStartI = $2;

            overlapI = conRefEndI - ampStartI;
            conSeqOnI = conLenSeqI;
            conCigOnI = conLenCigI;

            #++++++++++++++++++++++++++++++++++++++++++++++
            # Sec-0? Sub-02 Cat-03:
            #  - Find the 1st overlap base in the consensus
            #  - I need to adjust this for indels
            #++++++++++++++++++++++++++++++++++++++++++++++

            refNtI = 1; # Used to detect when at overlap

            while(refNtI <= overlapI)
            { # Loop: Find the ajdusted start (indels)
               if(conCigCStr[iNt] == "I") --conSeqOnI;

               else if(conCigCStr[iNt] == "D") ++refNtI;

               else
               { # Else: Everything matches up
                  ++refNtI;
                  --conSeqOnI;
               } # Else: Everything matches up

               --conCigOnI; # this always has a base
            } # Loop: Find the ajdusted start (indels)

            #++++++++++++++++++++++++++++++++++++++++++++++
            # Sec-0? Sub-02 Cat-04:
            #  - Set up for mergeing overlaps
            #++++++++++++++++++++++++++++++++++++++++++++++

            if($1 ~ /87/)
            { # error round
               for(iCig=conCigOnI;iCig <= conLenCigI;++iCig)
                  printf "%s", conCigCStr[iCig];
               printf "\n\n";

               for(iCig=1; iCig <= overlapI; ++iCig)
               { # For the overlap
                  printf "%s", ampCigCStr[iCig];
                  if(ampCigCStr[iCig] == "I") --iCig;
               } # For the overlap
               printf "\n\n";

                  exit;
            } # error round
            # Find the first reference base in the overlap
            ampSeqOnI = 1;
            ampCigOnI = 1;
            firstTmpSeqI = conSeqOnI + 1;
            firstTmpCigI = conCigOnI + 1;
              # The +1 is to account for undercounting
            lenTmpSeqI = 0;
            lenTmpCigI = 0;

            while(conSeqOnI <= conLenSeqI)
            { # Loop: Trim the new sequence bases

               #+++++++++++++++++++++++++++++++++++++++++++
               # Sec-0? Sub-02 Cat-05:
               #  - Handle insertions in the consensus
               #+++++++++++++++++++++++++++++++++++++++++++

               while(conCigCStr[conCigOnI] == "I")
               { # Loop: Go though all consensus insertions

                  if(ampCigCStr[ampCigOnI] == "I")
                  { # If: this insertion is supported
                     ++lenTmpSeqI;
                     ++lenTmpCigI;
                     ++ampCigOnI;
                     ++ampSeqOnI;

		               # Check if I have a match
                     if(conSeqCStr[conSeqOnI] == ampSeqCStr[ampSeqOnI])
                     { # If: the inerstions are the same
                        tmpSeqCStr[lenTmpSeqCStr] = conSeqCStr[conSeqOnI];

                        tmpCigCStr[lenTmpCigI] = "I";
                     } # If: the inerstions are the same

                     else
                     { # Else: No idea which is correct
                        tmpSeqCStr[lenTmpSeqI] = "N";
                        tmpCigCStr[lenTmpCigI] = "I";
                     } # Else: No idea which is correct
                  } # If: this insertion is supported

                  ++conSeqOnI; # move to next base
                  ++conCigOnI;
               } # Loop: Go though all consensus insertions

               #+++++++++++++++++++++++++++++++++++++++++++
               # Sec-0? Sub-02 Cat-05:
               #  - Handle insertions in the amplicon
               #+++++++++++++++++++++++++++++++++++++++++++

               while(ampCigCStr[ampCigOnI] == "I")
               { # Move past unsupported amplicon inss
                  ++ampCigOnI; # No support from consensus
                  ++ampSeqOnI; # move to next base
               } # Move past unsupported amplicon inss

               #+++++++++++++++++++++++++++++++++++++++++++
               # Sec-0? Sub-02 Cat-06:
               #  - Handle deletions in the consensus
               #+++++++++++++++++++++++++++++++++++++++++++

		         # Check if I have an deletion
		         if(conCigCStr[conCigOnI] == "D")
               { # If: I have a deletion
                  if(ampCigCStr[ampCigOnI] == "D")
                     tmpCigCStr[lenTmpCigI] = "D";

                  else
                  { # Else: This deletion is not supported
                     tmpSeqCStr[lenTmpSeqI] = ampSeqCStr[ampSeqOnI];

                     tmpCigCStr[lenTmpCigI] = "X";
                     ++lenTmpSeqI;
                     ++ampSeqOnI;                
                  } # Else: This deletion is not supported

                  ++lenTmpCigI;
                  ++conCigOnI;
                  ++ampCigOnI;
               } # If: I have a deletion

               #+++++++++++++++++++++++++++++++++++++++++++
               # Sec-0? Sub-02 Cat-07:
               #  - Handle deletions in the amplicon
               #+++++++++++++++++++++++++++++++++++++++++++
		         
               # There is not consensus deletion
               else if(ampSeqCStr[conCigOnI] == "D")
               { # If: I have a deletion
                  tmpSeqCStr[lenTmpSeqI] = conSeqCStr[conSeqOnI];

                  tmpCigCStr[lenTmpCigI] = "X";

                  ++lenTmpSeqI;
                  ++lenTmpCigI;
                  ++conSeqOnI;
                  ++conCigOnI;
                  ++ampCigOnI;
               } # If: I have a deletion

               #+++++++++++++++++++++++++++++++++++++++++++
               # Sec-0? Sub-02 Cat-08:
               #  - Handle matches and mismatches
               #+++++++++++++++++++++++++++++++++++++++++++

               else
               { # Else: Match or mismatch
                 
                  if(conSeqCStr[conSeqOnI] == ampSeqCStr[ampSeqOnI])
                  { # If: my amplicons are a match
                     tmpSeqCStr[lenTmpSeqI] = conSeqCStr[conSeqOnI];

                     tmpCigCStr[lenTmpCigI] = conCigCStr[conCigOnI];
                  } # If: my amplicons are a match

                  else
                  { # Else: there is not agreement
                     tmpSeqCStr[lenTmpSeqI] = "N";
                     tmpCigCStr[lenTmpCigI] = "X";
                  } # Else: there is not agreement

                  ++lenTmpSeqI;
                  ++lenTmpCigI;
                  ++conCigOnI;
                  ++conSeqOnI;
                  ++ampSeqOnI;                
                  ++ampCigOnI;
               } # Else: Match or mismatch
            } # Loop: Trim the new sequence bases

            #++++++++++++++++++++++++++++++++++++++++++++++
            # Sec-0? Sub-02 Cat-09:
            #  - Copy add the merged to the consensus
            #++++++++++++++++++++++++++++++++++++++++++++++

            # Account for overcounting in loop
            --lenTmpCigI;
            --lenTmpSeqI;

            for(iCig = 1; iCig <= lenTmpCigI; ++iCig)
            { # Loop: Copy over the expanded cigar
               conCigCStr[firstTmpCigI] = tmpCigCStr[iCig];
               ++firstTmpCigI;
            } # Loop: Copy over the expanded cigar

            for(iNt = 1; iNt <= lenTmpSeqI; ++iNt)
            { # Loop: Copy over the sequence
               conSeqCStr[firstTmpSeqI] = tmpSeqCStr[iNt];
               ++firstTmpSeqI;
            } # Loop: Copy over the sequence

            while(ampCigOnI <= ampLenCigI)
            { # Loop: Finish copying cigar entries
               conCigCStr[firstTmpCigI] = ampCigCStr[ampCigOnI];
               ++ampCigOnI;
               ++firstTmpCigI;
            } # Loop: Finish copying cigar entries

            while(ampSeqOnI <= ampLenSeqI)
            { # Loop: Finish copying cigar entries
               conSeqCStr[firstTmpSeqI] = ampSeqCStr[ampSeqOnI];
               ++ampSeqOnI;
               ++firstTmpSeqI;
            } # Loop: Finish copying cigar entries
 
            conRefEndI = ampRefEndI;
            conLenSeqI = firstTmpSeqI;
            conLenCigI = firstTmpCigI;
            ++iRnd;
            #if(iRnd > 10) exit;
            #if(iMaskRnd > 3) exit;
         } # If I have an overlap

         #+++++++++++++++++++++++++++++++++++++++++++++++++
         # Sec-0? Sub-02 Cat-11:
         #  - Handle gaps between sequences
         #+++++++++++++++++++++++++++++++++++++++++++++++++

         else
         { # Else I do not have an overlap
            ++iMaskRnd;
            maskI = $2 - conRefEndI;

            for(iN = 1; iN <= maskI; ++iN)
            { # Loop: Add in masking to the start
               conSeqCStr[conLenSeqI] = "N";
               conCigCStr[conLenCigI] = "X";
               ++conLenSeqI;
               ++conLenCigI;
               ++conRefEndI;
            } # Loop: Add in masking to the start

            maskLenI = maskLenI + maskI;

            gsub(//, " ", $5);
            ampLenSeqI = split($5, ampSeqCStr, " ");

            #if($1 ~ /92/) print ">test";
            for(iNt = 1; iNt <= ampLenSeqI; ++iNt)
            { # Loop: add in new bases
               conSeqCStr[conLenSeqI] = ampSeqCStr[iNt];
               #if($1 ~ /92/)printf "%s", ampSeqCStr[iNt];
               ++conLenSeqI;
            } # Loop: add in new bases

            #if($1 ~ /92/){printf "\n\n"; exit;};

            # Convert expanded cigar to c-string
            ampLenCigI = split($6, ampCigCStr, ",");

            for(iCig = 1; iCig <= ampLenCigI; ++iCig)
            { # Loop: add in new cigar entries
               conCigCStr[conLenCigI] = ampCigCStr[iCig];
               ++conLenCigI;

               if(conCigCStr[conLenCigI] != "I")
                  ++conRefEndI;
            } # Loop: add in new cigar entries
            
            #if(iMaskRnd > 3) exit;
         }; # Else I do not have an overlap
      }; # MAIN

      #****************************************************
      # Sec-0? Sub-01:
      #  - END block: print out the consensus sequence
      #****************************************************

      END{
          printf ">con start=%s", conStartI;
          printf " len=%s ref=%s", conLenSeqI, conRefEndI;
          printf " noNs=%s\n", maskLenI;
          
          for(iNt = 1; iNt <= conLenSeqI; ++iNt)
             printf "%s", conSeqCStr[iNt];
          printf "\n";

          #for(iCig = 1; iCig <= conLenCigI; ++iCig)
          #   printf "%s", conCigCStr[iCig];
          #printf "\n";

      }; # END block
    ';
