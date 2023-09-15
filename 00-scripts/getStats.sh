# This script gets some quick stats using alnSeq. It
# expects alnSeqMid (make mid) to be in
# scriptDir/../00-programs/alnSeqMid
# Output:
#  - stdout: Matchs\tSNPs\tInsertions\tDeletions\tNs

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

refStr="$1";
qryStr="$2";
methodStr="$3";
scriptDirStr="$(dirname "$0")";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ "$methodStr" != "Local" ]]; then
    methodStr="-use-needle";
else
    methodStr="-use-water";
fi

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Mask the missing bases at the ends
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Not the best way, but it kinda works
"$scriptDirStr/../00-programs/alnSeqMid" \
    -use-water \
    -format-fasta \
    -line-wrap 0 \
    -ref "$refStr" \
    -query "$qryStr" |
  awk \
      '
        { # Main
           # Check if is reference entry or blank line
           if(NR == 1 || NR == 2 || $0 ~ /^$/) next;
           print $0; # Header for query
           getline;  # sequence
           startGapStr = $0;
           endGapStr = $0;

           # Find the number of gaps at the ends
           sub(/[ATGCN].*/, "", startGapStr);
           gsub(/[-ATGCN]*[ATGCN]/, "", endGapStr);
           startNsI = length(startGapStr);
           endNsI = length(endGapStr);

           # Remove gaps
           gsub(/-/, "", $0);
           startNsStr = "";
           endNsStr = "";

           for(iN = startNsI; iN > 0; --iN)
              startNsStr = "N" startNsStr;
           for(iN = endNsI; iN > 0; --iN)
              endNsStr = endNsStr "N";
           print startNsStr $0 endNsStr;
        } # Main
      ' \
  > "20230914-del-tmp.fasta";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-04:
#  - Get the stats
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

"$scriptDirStr/../00-programs/alnSeqMid" \
      -ref "$refStr" \
      -query "20230914-del-tmp.fasta" \
      $methodStr |
   awk '
      BEGIN{numNsI = 0};
      { # MAIN
         if($2 == "Matches:") {numMatchsI = $3; next;};
         if($2 == "Mismatches:") {numSnpI = $3; next;};
         if($2 == "Insertions:") {numInsI = $3; next;};
         if($2 == "Deletions:") {numDelI = $3; next;};
         if($1 ~ /Qry:/)
         { # If I have a query entry
            sub(/Qry:  */, "", $0);
            gsub(/[ATGC-]/, "", $0);
            numNsI += length($0);
         } # If I have a query entry
      } # MAIN
      END{
         print numMatchsI,numSnpI,numInsI,numDelI,numNsI;
      };
';

rm "20230914-del-tmp.fasta";
