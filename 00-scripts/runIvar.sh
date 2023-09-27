#!/usr/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-01 Sub-01:
#  - General input
#**********************************************************

readsStr="";
refStr="";
primFileStr="";
prefixStr="output";
threadsI=3;
keepConQBl=0; # Keep consensus q-score entry

#**********************************************************
# Sec-01 Sub-02:
#  - ivar consensus variables
#**********************************************************

minDepthI=10;    # Min depth to keep an base
minSupDbl="0.5"; # Min percent support to keep an base
minInsSupDbl="0.8"; # Min % support oto keep an insertion
minQI=10;        # min Q-score to keep base

#**********************************************************
# Sec-01 Sub-03:
#  - ivar trim variables
#**********************************************************

minLenI=100; # Min length to keep a trimmed read
minTrimQI=0; # Min q-score to keep a trimmed region
minTrimWindowI=20; # Min window size for trimming
useIvarsTrimBl=0; # 1: use ivars internal trim function

#**********************************************************
# Sec-01 Sub-04:
#  - Script variables (user can not modify)
#**********************************************************

tmpPrimFileStr="$prefixStr-tmp.tsv";
scriptDirStr="$(dirname "$0")";
errBl=0;
isNumericStr="^[0-9][0-9]*$";
isQScoreStr="^[0-9][0-9]?$"; # Limit to 2 digits
isDecPercStr="^[0]*[10]\.?|^[0]*\.[0-9][0-9]*$";
  # ^ means start at first character in string
  # ^[0]* Any number of zeros (including none) can be at
  #    the start
  # [10] means a one or a zero follows (this is 0 or 1)
  # \.? Means that there can be a decimal afterwards
  #   1., 0., 1, 0
  # | means an alternative check
  # ^[0]*\. Means 0's followed by a decimal ex: 0., 00.
  # [0-9][0-9]*$ Means at least one number comes after the
  #   decimal ex: 0.x, 00.xxx

#**********************************************************
# Sec-01 Sub-05:
#  - Help message
#**********************************************************

helpStr="$(basename "$0") -ref ref.fa -fastq reads.fq [...]
Use:
  - Uses ivar to build a consensus
Input:
  General:
     -ref: [Required]
       o Reference to bulid consensus with
     -fastq: [Required]
       o Reads to build consensus with
     -primers: [Required]
       o Bed file with primers to trim reads with
     -ivar-trim: [No]
       o Calls the trimming function in ivar
       o disabled by -no-ivar-trim, which calls
         script-dir/../00-programs/trimPrimers
     -prefix or -p: [$prefixStr]
       o prefix to name the output
     -threads or -t: [$threadsI]
       o Number of threads to use
     -keep-q: [No]
       o Keep the Q-score file for the consensus output by
         ivar
       o -rm-q will disable -keep-q
  ivar consensus:
     -depth: [$minDepthI]
       o Minimum depth to not mask an position
     -min-support: [$minSupDbl]
       o Minimum support to keep an base (0.01 to 1)
     -min-support-ins: [$minInsSupDbl]
       o Minimum support to keep an insertion (0.01 to 1)
     -min-Q: [$minQI]
       o Minimum Q-score to keep a base
  ivar trim:
     -min-len: [$minLenI]
       o Minimum length after trimming to keep a read
     -min-trim-q: [$minTrimQI]
       o Minimum Q-score in window to keep read (??)
     -trim-window: [$minTrimWindowI]
       o Width of window to use in trimming
Output:
  - prefix-con.fa
    o Consensus sequence if no trimming was done
  - prefix-con.qual.txt
    o Q-score entry for prefix-con.fa
    o Only kept if -keep-q is used
  - prefix-con-trim.fa
    o Consensus sequence if trimming was done
  - prefix-con-trim.qual.txt
    o Q-score entry for prefix-con-trim.fa
    o Only kept if -keep-q is used
";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Get and check user input
#  o sec-02 sub-01:
#    - Get user input
#  o sec-02 sub-02:
#    - Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-02 Sub-01:
#  - Get user input
#**********************************************************

while [[ $# -gt 0 ]]; do
# Loop: Get user input
   case $1 in
     -ref) refStr="$2"; shift;;
     -fastq) readsStr="$2"; shift;;
     -primers) primFileStr="$2"; shift;;
     -prefix) prefixStr="$2"; shift;;
     -p) prefixStr="$2"; shift;;
     -threads) threadsI="$2"; shift;;
     -t) threadsI="$2"; shift;;
     -keep-q) keepConQBl=1; shift;;
     -rm-q) keepConQBl=0; shift;;
     -depth) minDepthI="$2"; shift;;
     -min-support) minSupDbl="$2"; shift;;
     -min-support-ins) minInsSupDbl="$2"; shift;;
     -min-Q) minQI="$2"; shift;;
     -min-len) minLenI="$2"; shift;;
     -min-trim-q) minTrimQI="$2"; shift;;
     -trim-window) minTrimWindowI="$2"; shift;;
     -ivar-trim) useIvarsTrimBl=1;;
     -no-ivar-trim) useIvarsTrimBl=0;;
     -h) printf "%s\n" "$helpStr"; exit;;
     --h) printf "%s\n" "$helpStr"; exit;;
     -help) printf "%s\n" "$helpStr"; exit;;
     --help) printf "%s\n" "$helpStr"; exit;;
     help) printf "%s\n" "$helpStr"; exit;;
     *) printf "%s\n%s is invalid\n" "$helpStr" "$1";exit;;
   esac

   shift;
done # Loop: Get user input

#**********************************************************
# Sec-02 Sub-02:
#  - Check user input
#  o sec-02 sub-02 cat-01:
#    - Check if reference is a file
#  o sec-02 sub-02 cat-02:
#    - Check if input reads fastq is a file
#  o sec-02 sub-02 cat-03:
#    - Check if primer file was input/is valid
#  o sec-02 sub-02 cat-04:
#    - Check if a valid prefix was entered
#  o sec-02 sub-02 cat-06:
#    - Check if a valid number of threads was input
#  o sec-02 sub-02 cat-06:
#    - Check if a valid min depth was intput
#  o sec-02 sub-02 cat-07:
#    - Check if a valid min support (snp/match) was intput
#  o sec-02 sub-02 cat-08:
#    - Check if a valid min insertion support
#  o sec-02 sub-02 cat-09:
#    - Check if a valid min quality score was input
#  o sec-02 sub-02 cat-10:
#    - Check if a valid min length was input
#  o sec-02 sub-02 cat-11:
#    - Check if a valid Q-score for trimming was input
#  o sec-02 sub-02 cat-12:
#    - Check if a valid trimming window length was input
#  o sec-02 sub-02 cat-13:
#    - Check if had an error
#**********************************************************

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-01:
#  - Check if input reference is a file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$refStr" == "" ]]; then
   printf "No file provided with -ref\n" >&2;
   errBl=1;
fi # If: The user did not input a reference file

if [[ ! -f "$refStr" ]]; then
   printf "Could not open -ref \"%s\"\n" "$refStr" >&2;
   errBl=1;
fi # If: I could not open the reference file

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-02:
#  - Check if input reads fastq is a file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$readsStr" == "" ]]; then
   printf "No file provided with -fastq\n" >&2;
   errBl=1;
fi # If: The user did not input a reference file

if [[ ! -f "$readsStr" ]]; then
   printf "Could not open -fastq \"%s\"\n" "$readsStr" >&2;
   errBl=1;
fi # If: I could not open the reference file

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-03:
#  - Check if primer file was input/is valid
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$primFileStr" != "" ]]; then
   if [[ ! -f "$primFileStr" ]]; then
      printf \
          "Could not open -primers \"%s\"\n" \
          "$primFileStr" \
        >&2;
      errBl=1;
   fi
fi # If: The user did not input a reference file

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-04:
#  - Check if a valid prefix was entered
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$(basename "$prefixStr")" == "" ]]; then

   if [[ "$(dirname "$prefixStr")" == "" ]]; then
      printf "Blank prefix input; using output\n" >&2;
   else
      printf \
          "Prefix is a directory; using %s/output\n" \
          "$prefixStr" \
        >&2;
   fi # Check if the prefix can be used

   prefixStr="$prefixStr/output";
fi # check if I need to use a default prefix

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-05:
#  - Check if a valid number of threads were input
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! "$threadsI" =~ $isNumericStr ]]; then
   if [[ "$threadsI" == "" ]]; then
      printf "Blank input for -threads; using 3\n" >&2;
   else
      printf \
         "Non-numeric input for -threads; using 3\n" \
         >&2;
   fi # Check if was blank or if invalid

   threadsI=3;
fi # Check if the threads input is valid

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-06:
#  - Check if a valid min depth was intput
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! "$minDepthI" =~ $isNumericStr ]]; then
   if [[ "$minDepthI" == "" ]]; then
      printf "Blank input for -depth\n" >&2;
   else
      printf \
         "Non-numeric input for -depth (%s)\n" \
         "$minDepthI" \
         >&2;
   fi # Check if was blank or if invalid

   errBl=1;
fi # Check if the depth input is valid

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-07:
#  - Check if a valid min support (snp/match) was intput
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! "$minSupDbl" =~ $isDecPercStr ]]; then
   if [[ "$minSupDbl" == "" ]]; then
      printf "Blank input for -min-support\n" >&2;
   else
      printf \
         "Non-decimal input for -min-support (%s)\n" \
         "$minSupDbl" \
         >&2;
   fi # Check if was blank or if invalid

   errBl=1;
fi # Check if the depth input is valid

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-08:
#  - Check if a valid min insertion support
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! "$minInsSupDbl" =~ $isDecPercStr ]]; then
   if [[ "$minInsSupDbl" == "" ]]; then
      printf "Blank input for -min-support-ins\n" >&2;
   else
      printf \
         "Non-numeric input for -min-support-ins (%s)\n" \
         "$minInsSupDbl" \
         >&2;
   fi # Check if was blank or if invalid

   errBl=1;
fi # Check if the depth input is valid

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-09:
#  - Check if a valid min quality score was input
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! "$minQI" =~ $isQScoreStr ]]; then
   if [[ "$minQI" =~ $isNumericStr ]]; then
      printf \
        "-min-q %s is to higher for a Q-score (0 to 93)\n"\
        "$minQI"\
        >&2;
   elif [[ "$minQI" == "" ]]; then
      printf "Blank entry for -min-q\n" >&2;
   else
      printf "-min-q %s is not numeric\n" "$minQI" >&2;
   fi # Find what made the Q-score invalid

   errBl=1;
fi # If this is not a valid q-score range

if [[ "$minQI" -gt 93 ]]; then
   printf
      "-min-q %s is to higher for a Q-score (0 to 93)\n"\
      "$minQI"\
      >&2;
   errBl=1;
fi # Check if the q-score is to high

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-10:
#  - Check if a valid min length was input
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! "$minLenI" =~ $isNumericStr && "$primFileStr" != "" ]];
then
   if [[ "$minLenI" == "" ]]; then
      printf "Blank input for -min-len\n" >&2;
   else
      printf "-min-len %s is not numeric\n" "$minLenI" >&2;
   fi # Check how the input was invalid

   errBl=1;
fi # Check if the user input a valid length

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-11:
#  - Check if a valid Q-score for trimming was input
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! "$minTrimQI" =~ $isQScoreStr && "$primFileStr" != "" ]];
then
   if [[ "$minTrimQI" =~ $isNumericStr ]]; then
      printf \
        "-min-trim-q %s is > max Q-score (0 to 93)\n"\
        "$minTrimQI"\
        >&2;
   elif [[ "$minTrimQI" == "" ]]; then
      printf "Blank entry for -min-trim-q\n" >&2;
   else
      printf \
         "-min-trim-q %s is not numeric\n" \
         "$minTrimQI" \
         >&2;
   fi # Find what made the Q-score invalid

   errBl=1;
fi # If this is not a valid q-score range

if [[ "$minTrimQI" -gt 93 ]]; then
   printf
      "-min-trim-q %s is to > max Q-score (0 to 93)\n"\
      "$minTrimQI"\
      >&2;
   errBl=1;
fi # Check if the q-score for trimming is to high

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-12:
#  - Check if a valid trimming window length was input
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! "$minTrimWindowI" =~ $isNumericStr && "$primFileStr" != "" ]];
then
   if [[ "$minTrimWindowI" == "" ]]; then
      printf "Blank input for -trim-window\n" >&2;
   else
      printf \
         "-trim-window %s is not numeric\n" \
         "$minTrimWindowI" \
         >&2;
   fi # Check how the input was invalid

   errBl=1;
fi # Check if the user input a valid trim window size

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-13:
#  - Check if there was an error
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$errBl" -gt 0 ]]; then exit; fi

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Run ivar
#  o sec-03 sub-01:
#    - Run ivar without trimming
#  o sec-03 sub-02:
#    - Run ivar with trimming
#  o sec-03 sub-03:
#    - Commeted out full trimming notes. Was not worth
#      doing
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-03 Sub-01:
#  - Run ivar without trimming
#**********************************************************

# By having duplicates commands I can avoid temporary files
if [[ "$useIvarsTrimBl" -eq 0 ]]; then
# If: I am using an external primer trim
   if [[ -f "$primFileStr" ]]; then
   # If: I am trimming primers
      awk \
          '{printf ">%s\n%s\n", $4, $7}' \
          < "$primFileStr" \
        > "$prefixStr-tmp-primers.fasta";

      "$scriptDirStr/../00-programs/trimPrimers" \
          -end-trim \
          -primers "$prefixStr-tmp-primers.fasta" \
          -fastq "$readsStr" \
          -out "$prefixStr-tmp.fastq";

      readsStr="$prefixStr-tmp.fastq";
   fi # If: I am trimming primers

   minimap2 \
       -t "$threadsI" \
       --eqx \
       -a \
       -x map-ont \
       "$refStr" \
       "$readsStr" |
     "$scriptDirStr/../00-programs/trimSamFile" -stdin |
     samtools sort -@ "$threadsI" - |
     samtools view \
        -@ "$threadsI" \
        -F 4 \
        -F 256 \
        -F 2048 \
        -b \
        - |
     samtools mpileup -B -aa -A -d 0 -Q 0 - | 
     "$scriptDirStr/../00-programs/ivar" consensus \
       -p "$prefixStr-con" \
       -i "$prefixStr-con" \
       -c "$minInsSupDbl" \
       -t "$minSupDbl" \
       -m "$minDepthI" \
       -n "N" \
       -q "$minQI";
   
   if [[ "$keepConQBl" -lt 1 ]]; then
      rm "$prefixStr-con.qual.txt";
   fi # If I do not want the quality score file

   # If I created a temporary file (readsStr changed to
   # hold name of temporary file)
   if [[ -f "$primFileStr" ]]; then rm "$readsStr"; fi
# If: I am using an external primer trim

#**********************************************************
# Sec-03 Sub-02:
#  - Run ivar without trimming
#**********************************************************

else
# Else: I am trimming primers with ivar
   awk \
       -f "$scriptDirStr/setUpIvarPrimFile.awk" \
       < "$primFileStr"\
     > "$tmpPrimFileStr";

   minimap2 --eqx -a -x map-ont "$refStr" "$readsStr" |
     "$scriptDirStr/../00-programs/trimSamFile" -stdin |
     samtools sort -@ "$threadsI" - |
     samtools view \
        -@ "$threadsI" \
        -F 4 \
        -F 256 \
        -F 2048 \
        -b \
        - |
     "$scriptDirStr/../00-programs/ivar" trim \
       -b "$primFileStr" \
       -f "$tmpPrimFileStr" \
       -m "$minLenI" \
       -s "$minTrimWindowI" \
       -q "$minTrimQI" |
     samtools sort -@ "$threadsI" - |
     samtools mpileup -B -aa -A -d 0 -Q 0 - | 
     "$scriptDirStr/../00-programs/ivar" consensus \
       -p "$prefixStr-con-trim" \
       -i "$prefixStr-con-trim" \
       -c "$minInsSupDbl" \
       -t "$minSupDbl" \
       -m "$minDepthI" \
       -n "N" \
       -q "$minQI";

   if [[ "$keepConQBl" -lt 1 ]]; then
      rm "$prefixStr-con-trim.qual.txt";
   fi # If I do not want the quality score file

   rm "$tmpPrimFileStr";
# Else: I am trimming primers with ivar
fi # Check if I am doing a primer trim run or not

exit;

# file saved as "$prefixStr-con.fa"
# my min depth is scricter than PMC9231111, but I do not
# want a depth of 1
  # -n is mask
  # -q is min quality score
  # -m is min depth
  # -t is min percentate to keep variant
  # -p is the prefix to name everything
  # -i is the header for the file
  # -c is the min insertion frequence (0.8 is default)

#**********************************************************
# Sec-03 Sub-03:
#  - Commeted out full trimming notes. Was not worth doing
#**********************************************************

#This is the full ivar run. I did not find it very usefull.
#  samtools index -@ 3 "$tmpFileStr" &&
#  ivar auto indexs

# Complete protocol
# minimap2 --eqx -a -x map-ont "$refStr" "$readsStr" |
#   samtools sort -@ 3 - |
#   samtools view -F 4 -F 256 -F 2048 -b - |
#   ./src/ivar trim \
#     -q "$minTrimQI" \
#     -m "$minDepthI"0 \
#     -b "$primFileStr" \
#     -f "$tmpPrimFileStr" \
#     -s "$minTrimWindowI" |
#  samtools sort -o delete.bam &&
#  samtools mpileup -B -aa -A -d 0 -Q 0 delete.bam |
#  ./src/ivar variants \
#    -p delet-var \
#    -q "$minQI" \
#    -t "$minSupDbl" \
#    -m "$minDepthI" \
#    -r "$refStr" &&
#  ./src/ivar getmasked \
#    -i delet-var.tsv \
#    -b "$primFileStr" \
#    -f "$tmpPrimFileStr" \
#    -p delet-mask && 
#  /src/ivar removereads \
#    -i delete.bam \
#    -t delet-mask.txt \
#    -p delete-filt \
#    -b "$primFileStr" &&
#  samtools sort delete-filt.bam > tmp.bam &&
#  mv tmp.bam delete-filt.bam &&
#  samtools index delete-filt.bam &&
#  samtools mpileup -B -aa -A -d 0 -Q 0 delete-filt.bam |
#  ./src/ivar consensus \
#    -p delete-filt-con \
#    -i delete \
#    -c "$minInsSupDbl" \
#    -t "$minSupDbl" \
#    -m "$minDepthI" \
#    -n "N" \
#    -q "$minQI";
