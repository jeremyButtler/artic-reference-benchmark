# This script gets some quick stats using quast. It
# expects quast to be installed in
# script-dir/../00-programs/quast/
# Output:
#  - stdout:
#    o $3=1: query-length\tref-length\tNs\tsnps\tindels
#    o $3=0: query-length\tNs\tsnps\tindels

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

refStr="$1";
qryStr="$2";
#methodStr="$3";
pRefLen=$3; # Print reference length
scriptDirStr="$(dirname "$0")";

if [[ "$pRefLen" == "" ]]; then pRefLen=0;
elif [[ "$pRefLen" -gt 0 ]]; then pRefLen=1;
else pRefLen=0;
fi

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Run quast
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

"$scriptDirStr/../00-programs/quast/quast.py" \
    -r "$refStr" \
    -o "delete-quast" \
    -m 300 \
    -t 3 \
    --no-icarus \
    --no-html \
    --no-plots \
    --silent \
    --space-efficient \
    "$qryStr" \
    > /dev/null &&
  awk \
      -v pRefLen="$pRefLen" \
      '
        BEGIN{FS=OFS="\t"};
        { # MAIN
          if(NR == 1) next;
          if(pRefLen > 0)
             print $16, $17, $41/100, $42/100, $43/100;

          else print $16, $41/100, $42/100, $43/100;
            # $16 is query length
            # $17 is reference length
            # $41/100 is the number of Ns per 1000 bases
            # $42/100 is the number of snps per 1000 bases
            # $43/100 is number insertions per 1000 bases
        } # MAIN
      ' \
      < "delete-quast/transposed_report.tsv";

rm -r "delete-quast" || printf "";

# Old method (needs some debugging
#if [[ "$methodStr" != "Local" ]]; then
#    methodStr="-use-hirschberg";
#else
#    methodStr="-use-water";
#fi
#
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Sec-03:
##  - Mask the missing bases at the ends
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
## Not the best way, but it kinda works
#"$scriptDirStr/../00-programs/alnSeqMid" \
#    -use-water \
#    -format-fasta \
#    -line-wrap 0 \
#    -ref "$refStr" \
#    -query "$qryStr" |
#  awk \
#      '
#        { # Main
#           # Check if is reference entry or blank line
#           if(NR == 1 || NR == 2 || $0 ~ /^$/) next;
#           print $0; # Header for query
#           getline;  # sequence
#           startGapStr = $0;
#           endGapStr = $0;
#
#           # Find the number of gaps at the ends
#           sub(/[ATGCN].*/, "", startGapStr);
#           gsub(/[-ATGCN]*[ATGCN]/, "", endGapStr);
#           startNsI = length(startGapStr);
#           endNsI = length(endGapStr);
#
#           # Remove gaps
#           gsub(/-/, "", $0);
#           startNsStr = "";
#           endNsStr = "";
#
#           for(iN = startNsI; iN > 0; --iN)
#              startNsStr = "N" startNsStr;
#           for(iN = endNsI; iN > 0; --iN)
#              endNsStr = endNsStr "N";
#           print startNsStr $0 endNsStr;
#        } # Main
#      ' \
#  > "20230914-del-tmp.fasta";
#
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Sec-04:
##  - Get the stats
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
#"$scriptDirStr/../00-programs/alnSeqMid" \
#      -ref "$refStr" \
#      -query "20230914-del-tmp.fasta" \
#      $methodStr |
#   awk '
#      BEGIN{numNsI = 0};
#      { # MAIN
#         if($2 == "Matches:") {numMatchsI = $3; next;};
#         if($2 == "Mismatches:") {numSnpI = $3; next;};
#         if($2 == "Insertions:") {numInsI = $3; next;};
#         if($2 == "Deletions:") {numDelI = $3; next;};
#         if($1 ~ /Qry:/)
#         { # If I have a query entry
#            sub(/Qry:  */, "", $0);
#            gsub(/[ATGC-]/, "", $0);
#            numNsI += length($0);
#         } # If I have a query entry
#      } # MAIN
#      END{
#         print numMatchsI,numSnpI,numInsI,numDelI,numNsI;
#      };
#';
#
#rm "20230914-del-tmp.fasta";
