#!/usr/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TOC: Table Of contents
#  - Subsample reads from a fastq file
#  o sec-01:
#    - Variable declertions
#  o sec-02:
#    - Get and check user input
#  o sec-03:
#    - Subsample reads
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################################
# Name: subsampleReads.sh
# Use:
#  - Subsamples an input number of reads from a fastq file
# Input:
#  -fastq: [Required]
#    o Fastq file with reads to subsample
#  -depth: [500]
#    o Number of reads to subsample from the fastq file
#  -seed: [1024]
#    o Seed for the random number generator.
# Output:
#  - Prints the subsampled reads to stdout.
###########################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declertions
# o sec-01 sub-01:
#   - User input and script variables
# o sec-01 sub-02:
#   - Help message
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#**********************************************************
# Sec-01 Sub-01:
#  - Usser input and script variables
#**********************************************************

readsStr="";
sampleDepthI="500";
seedI="1024";

numReadsI=0; # Number of reads in file

#**********************************************************
# Sec-01 Sub-02:
#  - Help message
#**********************************************************

helpStr="\
$(basename "$0") -fastq reads.fastq [options...]
Input:
  -fastq: [Required]
    o Fastq file with reads to subsample
  -depth: [$sampleDepthI]
    o Number of reads to subsample from the fastq file
  -seed: [$seedI]
    o Seed for the random number generator.
Output:
  - Prints the subsampled reads to stdout.
";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Get and check user input
# o sec-02 sub-01:
#  - Get user input
# o sec-02 sub-02:
#  - check user input
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#**********************************************************
# Sec-02 Sub-01:
#  - Get user input
#**********************************************************

while [[ $# -gt 0 ]]; do
# loop: get user input
   case $1 in
     -fastq) readsStr="$2"; shift;;
     -depth) sampleDepthI="$2"; shift;;
     -seed) seedI="$2"; shift;;
     -h) printf "%s\n" "$helpStr"; exit;;
     --h) printf "%s\n" "$helpStr"; exit;;
     -help) printf "%s\n" "$helpStr"; exit;;
     --help) printf "%s\n" "$helpStr"; exit;;
     help) printf "%s\n" "$helpStr"; exit;;
     *) printf "%s\n%s is invalid\n" "$helpStr" "$1";exit;;
   esac

   shift;
done # loop: get user input

#**********************************************************
# Sec-02 Sub-02:
#  - Check user input
#**********************************************************

if [[ ! -f "$readsStr" ]]; then
   printf "-fastq %s is not a valid file\n" "$readsStr";
   exit;
fi # If no valid fastq file

if [[ "$seedI" == "\"\"" ]]; then
   seedI="";
fi # If the user wanted a blank seed

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Subsample reads
# o sec-03 sub-01:
#   - Prep (get number of reads) and check if can subsample
# o sec-03 sub-02:
#   - Check if it is better to select the reads I am not
#     going to keep (subsampling > 50% of reads)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#**********************************************************
# Sec-03 Sub-01:
#  - Prep (get number of reads) and check if can subsample
#**********************************************************

numReadsI="$(\
   sed -n 'p;n;n;n;' < "$readsStr" |
      wc -l |
      awk '{print $1}' \
)"; # get the number of reads in the file

if [[ "$sampleDepthI" -ge "$numReadsI" ]]; then
   cat "$readsStr";
   exit;
fi # If subsampling the entire file

#**********************************************************
# Sec-03 Sub-02:
#  - Check if it is better to select the reads I am not
#    going to keep (subsampling > 50% of reads)
#  o sec-03 sub-02 cat-01:
#    - Declare variables and set up for making the list
#  o sec-03 sub-02 cat-02:
#    - Make the do not extract list
#  o sec-03 sub-02 cat-02:
#    - Extract reads of interest
#**********************************************************

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-03 Sub-02 Cat-01:
#  - Declare variables and set up for making the list
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$sampleDepthI" -ge "$((numReadsI / 2))" ]]; then
   awk \
      -v depthI="$sampleDepthI" \
      -v numReadsI="$numReadsI" \
      -v seedI="$seedI" \
      '
         BEGIN{
            readCntI = 1;
            readI = 0;
            del doNotExtAryI[1];

            if(seedI == "") srand();
            else srand(seedI); # Set the random seed

            #++++++++++++++++++++++++++++++++++++++++++++++
            # Sec-03 Sub-02 Cat-02:
            #  - set up the do not extract list
            #++++++++++++++++++++++++++++++++++++++++++++++

            # Get the number of reads to not extract
            depthI = numReadsI - depthI;

            while(readCntI <= depthI)
            { # While I have to build my ignore list
               readI = 1 + int(rand() * numReadsI);

               # Check if the read is already in the list
               if(doNotExtAryI[readI] == 1) continue;

               doNotExtAryI[readI] = 1;
               ++readCntI;
            } # While I have to build my ignore list

            readCntI = 1;
         } # BEGIN BLOCK

         #+++++++++++++++++++++++++++++++++++++++++++++++++
         # Sec-03 Sub-02 Cat-03:
         #  - Extract reads of interest
         #+++++++++++++++++++++++++++++++++++++++++++++++++

         { # MAIN BLOCK (extract reads)
            if(doNotExtAryI[readCntI] == 1)
            { # If the read is on the ignore list
               getline; # sequence
               getline; # + entry
               getline; # Q-score entry

               ++readCntI; # Move to the next read
               next;    # restart on next entries header
            } # If the read is on the ignore list
             
            print $0; # print header entry
            getline;  # Move to sequence entry

            print $0; # print sequence entry
            getline;  # Move to + entry

            print $0; # print + entry
            getline;  # Move to q-score entry

            print $0; # print + entry

            ++readCntI; # move to the next read
         } # MAIN BLOCK (extract reads)
      ' < "$readsStr";
   exit;
fi # If it would be quick to mark reads not extracting

#**********************************************************
# Sec-03 Sub-03:
#  - It is better to select reads to extract
#    (extracting < 50% of reads)
#  o sec-03 sub-03 cat-01:
#    - Declare variables and set up for making the list
#  o sec-03 sub-03 cat-02:
#    - Make the extract list
#  o sec-03 sub-03 cat-02:
#    - Extract reads of interest
#**********************************************************

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-03 Sub-03 Cat-01:
#  - Declare variables and set up for making the list
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

awk \
   -v depthI="$sampleDepthI" \
   -v numReadsI="$numReadsI" \
   -v seedI="$seedI" \
   '
      BEGIN{
         readCntI = 1;
         readI = 0;
         del extractAryI[1];
         srand(seedI);

         #++++++++++++++++++++++++++++++++++++++++++++++
         # Sec-03 Sub-02 Cat-03:
         #  - set up the extract list
         #++++++++++++++++++++++++++++++++++++++++++++++

         while(readCntI <= depthI)
         { # While I have to build my ignore list
            readI = 1 + int(rand() * numReadsI);

            # Check if the read is already in the list
            if(extractAryI[readI] == 1) continue;

            extractAryI[readI] = 1;
            ++readCntI;
         } # While I have to build my ignore list

         readCntI = 1;
      } # BEGIN BLOCK

      #+++++++++++++++++++++++++++++++++++++++++++++++++
      # Sec-03 Sub-02 Cat-03:
      #  - Extract reads of interest
      #+++++++++++++++++++++++++++++++++++++++++++++++++

      { # MAIN BLOCK (extract reads)
         if(extractAryI[readCntI] != 1)
         { # If the read is on the ignore list
            getline; # sequence
            getline; # + entry
            getline; # Q-score entry

            ++readCntI; # Move to the next read
            next;    # restart on next entries header
         } # If the read is on the ignore list
          
         print $0; # print header entry
         getline;  # Move to sequence entry

         print $0; # print sequence entry
         getline;  # Move to + entry

         print $0; # print + entry
         getline;  # Move to q-score entry

         print $0; # print + entry

         ++readCntI; # move to the next read
      } # MAIN BLOCK (extract reads)
   ' < "$readsStr";
exit;
