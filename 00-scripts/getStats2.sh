#!/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# getStats TOC: Table Of Contents
# Sec-01:
#  - Variable declerations
# Sec-02:
#  - Get user input
# Sec-04:
#  - Print out the error entry
# Sec-03:
#  - Run quast
# Sec-04:
#  - Run water
# Sec-05:
#  - Print out the error entry
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################################
# Name: getStats.sh
# Use:
#  - This script gets some quick stats using quast or water
#     from emboss.
# Input:
#  -ref: [Required]
#     o Reference sequence (as fasta) to compare query to.
#  -query: [Required]
#     o Query sequence to check accuracy of (as fasta).
#  -p-ref-len: [Yes]
#     o Print the length of the reference sequence.
#     o This can be disabled with -no-ref-len.
#  -p-error: [No]
#     o Print out an error message instead of stats.
#     o This can be disabled with -no-error.
#  -p-header: [No]
#     o Print out a tab delminate header (only).
#     o There will be no newlines printed.
#     o Disabled with -no-head.
#  -quast: [use quast]
#     o Use quast to get stats.
#     o quast.py: $scriptDirStr/../00-programs/quast/
#  -water: [use quast]
#     o Use water from emboss to get stats.
#     o This water should be in your path (type water -h
#       to check if in path).
# Output:
#  - stdout:
#    o Prints query length, reference length (if -p-ref-len
#      used), number of mis-assemblies, N's, snps, indels,
#      insertions, deletions, matches, and program used to
#      stdout in a tab deliminated format.
# Requirments:
#  - This expects water to be in the users path (if using
#    water (EMBOSS) or if using quast, quast to be
#    installed in script-location/../00-programs/quast/
# Note:
#  -  Needle from emboss takes to much memory (12 Gigs for
#     aliging two 29000 nucleotide long genomes.
###########################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#  o sec-01 sub-01:
#    - Input variable declerations
#  o sec-01 sub-02:
#  - Non-user input variables
#  o sec-01 sub-03:
#    - Help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-01 Sub-01:
#  - Input variable declerations
#**********************************************************

refStr="";
qryStr="";
methodStr="quast";
headPrefStr=""; # Prefix to add to header entries
pRefLen=1; # Print reference length
pErrBl=0;  # Print just an error entry
pHeadBl=0; # print the header entry

#**********************************************************
# Sec-01 Sub-02:
#  - Non-user input variable declerations
#**********************************************************

scriptDirStr="$(dirname "$0")";
refLenI=0; # Holds length of the reference sequence
qryLenI=0; # Holds query length (0 for error)

#**********************************************************
# Sec-01 Sub-03:
#  - Help message
#**********************************************************

helpStr="
Name: getStats.sh
Use:
  - Gets stats for a consensus sequence using a reference
    sequence.
Call: $(basename "$0") -ref ref.fa -query query.fa options
Input:
  -ref: [Required]
     o Reference sequence (as fasta) to compare query to.
     o Reference length is the length of all reference
       sequences for quast and the first reference sequence
       for water.
  -query: [Required]
     o Query sequence to check accuracy of (as fasta).
  -p-ref-len: [Yes]
     o Print the length of the reference sequence.
     o This can be disabled with -no-ref-len.
  -p-error: [No]
     o Print out an error message instead of stats.
     o This can be disabled with -no-error.
  -p-header: [No]
     o Print out a tab delminate header (only).
     o There will be no newlines printed.
     o Disabled with -no-head.
  -header-prefix: [None]
     o Prefix to add to each header entry.
  -quast: [use quast]
     o Use quast to get stats (handles multiple query and
       reference sequences).
     o quast.py: $scriptDirStr/../00-programs/quast/
  -water: [use quast]
     o Use water from emboss to get stats.
     o This will only use the first reference sequence in
       the reference fasta.
     o This will merge all query sequences into one entry.
     o This water should be in your path (type water -h
       to check if in path).
Output:
  - stdout:
    o Prints query length, reference length (if -p-ref-len
      used), number of mis-assemblies, N's, snps, indels,
      insertions, deletions, matches, and program used to
      stdout in a tab deliminated format.
Requirments:
  - This expects water to be in the users path (if using
    water (EMBOSS) or if using quast, quast to be
    installed in script-location/../00-programs/quast/

";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Get user input
#  o sec-02 sub-01:
#    - Get user input
#  o sec-02 sub-02:
#    - Check user input
#  o sec-02 sub-03:
#    - Set variables from user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-02 Sub-01:
#  - Get user input
#**********************************************************

while [[ $# -gt 0 ]]; do
# Loop: Get user input
    case $1 in
        -ref) refStr="$2"; shift;;
        -query) qryStr="$2"; shift;;
        -p-ref-len) pRefLen=1;;
        -no-ref-len) pRefLen=0;;
        -p-error) pErrBl=1;;
        -no-error) pErrBl=0;;
        -p-header) pHeadBl=1;;
        -no-head) pHeadBl=0;;
        -header-prefix) headPrefStr="$2"; shift;;

        # Methods
        -quast) methodStr="quast";;
        -water) methodStr="water";;
        # TODO: Add these methods in (first fix alnSeq)
        #-alnseq) medthodStr="alnseq"; shift;;
        #-alnseq-mid) medthodStr="alnseqmid"; shift;;
        #-alnseq-hirsch) medthodStr="hirsch"; shift;;

        # errors and help
        -h) printf "%s\n" "$helpStr"; exit;;
        --h) printf "%s\n" "$helpStr"; exit;;
        -help) printf "%s\n" "$helpStr"; exit;;
        --help) printf "%s\n" "$helpStr"; exit;;
        help) printf "%s\n" "$helpStr"; exit;;
        *) printf "%s\n%s is invalid\n" "$helpStr" "$1";
           exit;;
    esac

    shift;
done # Loop: Get user input

#**********************************************************
# Sec-02 Sub-02:
#  - Check user input
#**********************************************************

if [[ ! -f "$refStr" ]]; then
# Check if the reference is valid
  if [[ "$pErrBl" -lt 1 && "$pHeadBl" -lt 1 ]]; then
  # If this is an actual run
    methodStr="Error";
    pErrBl=1;
    printf "No reference input; printing error\n" >&2;
  fi # If this is an actual run
fi # Check if the reference is valid

if [[ ! -f "$qryStr" ]]; then
# Check if the query is valid
  if [[ "$pErrBl" -lt 1 && "$pHeadBl" -lt 1 ]]; then
  # If this is an actual run
    methodStr="Error";
    pErrBl=1;
    printf "No query input; printing error\n" >&2;
  fi # If this is an actual run
fi # Check if the query is valid

#**********************************************************
# Sec-02 Sub-03:
#  - Set variables from user input
#  o sec-02 sub-03 cat-01:
#    - Set method type for error message and header
#  o sec-02 sub-03 cat-02:
#    - Get the reference sequence length
#  o sec-02 sub-03 cat-03:
#    - Check if have a query sequence
#**********************************************************

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-03 Cat-01:
#  - Set method type for error message and header
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#  o sec-02 sub-03 cat-01:
if [[ "$pErrBl" -gt 0 ]]; then methodStr="Error"; fi
if [[ "$pHeadBl" -gt 0 ]]; then methodStr="Header"; fi

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-03 Cat-02:
#  - Get the reference sequence length
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$pErrBl" -lt 1 && "$pHeadBl" -lt 1 ]]; then
# If I not doing an error or header print (getting stats)
  refLenI="$(\
    awk \
      -v methodStr="$methodStr" \
        '
           BEGIN{headBl = 0; refLenI = 0;};
  
           { # MAIN
             if($0 ~ /^>/)
             { # If I found a header entry
               # Check if using water from emboss (only
               # uses the first reference, so ignore other
               # references
               if(headBl==1 && methodStr=="water") exit;
  
               headBl = 1;
               next;
             } # If I found a header entry
  
             if(headBl == 0) next; # No head found yet
  
             gsub(/[- \t]/, "", $0); # remove white space
             if($0 ~ /^$/) next;   # Skip blank lines
  
             # get length of all sequences in fasta file
             refLenI += length($0);
           } # MAIN
  
           END{print refLenI;};
      ' \
      < "$refStr" \
  )";

  if [[ "$refLenI" -lt 1 ]]; then 
        printf "Reference has no sequence\n" >&2;
        methodStr="Error";
        pErrBl=1;
  fi # If I have no reference sequence
fi # If I not doing an error or header print (getting stats)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-03 Cat-03:
#  - Check if have a query sequence
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$pErrBl" -lt 1 && "$pHeadBl" -lt 1 ]]; then
# If I not doing an error or header print (getting stats)
  qryLenI="$(\
    awk \
        '
           BEGIN{headBl = 0; seqBl = 0;};
  
           { # MAIN
             if($0 ~ /^>/) next; # If this is a header
  
             # Remove whites space an non ATGC characters
             gsub(/[-\t ]/,"",$0); # Remove whitespace & -
             if($0 ~ /^$/) next;  # Skip blank lines
             qryLenI += length($0);

             gsub(/[- Bb D-F d-f H-S h-s U-Z u-z]/, "", $0);
               
             # Determine if have a sequence entry
             if($0 == "") seqBl = 0; # No sequence entry
             else seqBl = 1;
           } # MAIN
  
           END{if(seqBl==0) print 0; else print qryLenI;};
      ' \
      < "$qryStr" \
  )";

  if [[ "$qryLenI" -lt 1 ]]; then 
      printf \
          "One or more query sequences have no sequence\n"\
        >&2;
      methodStr="Error";
      pErrBl=1;
   fi # If the query sequence was invalid
fi # If I not doing an error or header print (getting stats)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Run quast
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ "$methodStr" == "quast" ]]; then
# If using quast
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
   
             # Find the number of convert from err/100kb to
             # acutal number of errors
             convI = $17 / 100000; # $17 = reference length

             entryStr = $16; # Consensus length
             entryStr = entryStr "\t" $45; # aligned length

             if(pRefLen > 0) entryStr = entryStr "\t" $17;

             # TODO: Add in missassemblies stats
             entryStr = entryStr "\t" "TODO:-misAms";

             entryStr = entryStr "\t" $41 * convI; # num Ns
             entryStr = entryStr "\t" $42 * convI; # snps
             entryStr = entryStr "\t" $43 * convI; # indels
             entryStr = entryStr "\t" "NA"; # insertions
             entryStr = entryStr "\t" "NA"; # deletions
             entryStr = entryStr "\t" "NA"; # matches
             entryStr = entryStr "\t" "quast";

             print entryStr;
           } # MAIN
         ' \
         < "delete-quast/transposed_report.tsv";
   
   rm -r "delete-quast" || printf "";
   exit;
fi # If using quast

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-04:
#  - Run water
#  o sec-04 sub-01:
#    - Run water from emboss and awk
#  o sec-04 sub-02:
#    - Get meta data from water header
#  o sec-04 sub-03:
#    - Find the number of insertions
#  o sec-04 sub-04:
#    - Find the number of deletions and Ns
#  o sec-04 sub-05:
#    - Find the reference and query lengths
#  o sec-04 sub-06:
#    - Find snps and print out the data
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-04 Sub-01:
#  - Run water from emboss and awk
#**********************************************************

if [[ "$methodStr" == "water" ]]; then
# If using water from emboss to get stats
   water \
       -asequence "$refStr" \
       -bsequence "$qryStr" \
       -stdout \
       -auto \
       -gapopen 10 \
       -gapextend 2 |
     awk \
         -v pRefLen="$pRefLen" \
         -v refLenI="$refLenI" \
         -v conLenI="$qryLenI" \
         '
           #***********************************************
           # Sec-04 Sub-02:
           #  - Get meta data from water header
           #***********************************************

           BEGIN{OFS="\t"};
           { # MAIN
             # Get to the aligned length
             while(getline) if($0 ~ /^# Length:/) break;
             if($0 !~ /^# Length:/) exit;

             lenAlnI += $3; # length of alignment

             getline; # Identity line
             if($0 == 0) exit;

             getline; # Similarity line
             if($0 == 0) exit;

             sub(/\/.*/, "", $3);
             matchesI += $3; 

             getline; # gaps
             if($0 == 0) exit;

             sub(/\/.*/, "", $3);
             gapI += $3;

             while(getline) if($0 !~ /^#/) break;
             if($0 == 0) exit;

             getline;
             if($0 == 0) exit;

             while($0 !~ /^#/ && $0 !~ /^$/)
             { # Loop: Go though the full alignment

                #******************************************
                # Sec-04 Sub-03:
                #  - Find the number of insertions
                #******************************************

                # Check if at end of alignment

                gsub(/[A-Z a-z]/, "", $3);
                insI += length($3);

                #******************************************
                # Sec-04 Sub-04:
                #  - Find the number of deletions and Ns
                #******************************************

                getline; # move to err entry
                getline; # move to the query sequence
                tmpQryLenI = $4;

                gsub(/[A-M O-Z a-m o-z]/, "", $3);
                nI += length($3);

                # Find the number of deletions
                gsub(/[NnXx]/, "", $3);
                tmpI = length($3);

                # Adjust Ns and add deletions
                nI -= tmpI;
                delI += tmpI;

                getline; # move onto the blank line
                if($0 == 0) exit;

                getline; # move onto the reference entry
                if($0 == 0) exit;
             } # Loop: Go though the full alignment

             #*********************************************
             # Sec-04 Sub-05:
             #  - Find the reference and query lengths
             #*********************************************

             qryLenI += tmpQryLenI;
           } # MAIN

           END{

             #*********************************************
             # Sec-04 Sub-06:
             #  - Find snps and print out the data
             #*********************************************

             # Find the number of snps in the alignment
             snpI = lenAlnI - matchesI - gapsI - nI;

             # add in consensus (full query) length
             entryStr = conLenI;

             # recored the aligned query length
             entryStr = entryStr "\t" qryLenI;
 
             # Add in the consensus length (if requested)
             if(pRefLen >0) entryStr=entryStr "\t" refLenI;

             # Put NA for quast mis-assembly entry
             entryStr = entryStr "\t" "NA";

             entryStr = entryStr "\t" nI;   # number of Ns
             entryStr = entryStr "\t" snpI; # number snps
             entryStr = entryStr "\t" gapI; # number gaps
             entryStr = entryStr "\t" insI; #num insertions
             entryStr = entryStr "\t" delI; # num deletions
             entryStr = entryStr "\t" matchesI; # macthes
             entryStr = entryStr "\t" "water";

             print entryStr;
           } # END
         ';
   
   exit;
fi # If using water from emboss to get stats

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-04:
#  - Print out the error entry
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


if [[ "$pErrBl" -gt 0 ]]; then
  printf "NA"; # Query length
  printf "\tNA"; # aligned length

  if [[ "$pRefLen" -gt 0 ]]; then
     printf "\tNA"; # Reference length
  fi

  { # Print out the rest of the header
    printf "\tNA"; # mis-assemblies
    printf "\tNA"; # Number of Ns
    printf "\tNA"; # Number of snps
    printf "\tNA"; # Number of gaps
    printf "\tNA"; # Number of insertions
    printf "\tNA"; # Number of deletions
    printf "\tNA"; # Number of matches
    printf "\tERROR"; # Number of matches
  } # Print out the rest of the header

  exit;
fi # If the user wanted an error line only

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-05:
#  - Print out the error entry
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ "$pHeadBl" -gt 0 ]]; then
# If I user wanted to print the header
  printf "%sQueryLen" "$headPrefStr";
  printf "%sAlnLen" "$headPrefStr";

  if [[ "$pRefLen" -gt 0 ]]; then
     printf "\t%sRefLen" "$headPrefStr";
  fi

  { # Print out the rest of the header
    printf "\t%sMisAssemblies" "$headPrefStr";
    printf "\t%sNs" "$headPrefStr";
    printf "\t%sSnps" "$headPrefStr";
    printf "\t%sGaps" "$headPrefStr";
    printf "\t%sInss" "$headPrefStr";
    printf "\t%sDels" "$headPrefStr";
    printf "\t%sMatches" "$headPrefStr";
    printf "\t%sStatsProgram" "$headPrefStr";
  } # Print out the rest of the header
fi # If I user wanted to print the header

# Old method (alnSeq needs some debugging). I would also
# need to updated the output stats.
#if [[ "$methodStr" != "Local" ]]; then
#    methodStr="-use-hirschberg";
#else
#    methodStr="-use-water";
#fi
#
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Sec-0?:
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
## Sec-0?:
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
#         if($2 == "Insertions:") {numInI = $3; next;};
#         if($2 == "Deletions:") {numDelI = $3; next;};
#         if($1 ~ /Qry:/)
#         { # If I have a query entry
#            sub(/Qry:  */, "", $0);
#            gsub(/[ATGC-]/, "", $0);
#            numNsI += length($0);
#         } # If I have a query entry
#      } # MAIN
#      END{
#         print numMatchsI,numSnpI,numInI,numDelI,numNsI;
#      };
#';
#
#rm "20230914-del-tmp.fasta";
