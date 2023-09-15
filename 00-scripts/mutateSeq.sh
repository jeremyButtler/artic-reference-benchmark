#!/usr/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mutateSeq.sh TOC:
#  o Sec-01:
#    - Variable Declerations
#  o Sec-02:
#    - Read in user input
#  o Sec-03:
#    - Mutate the sequence
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################################
# Name: mutateSeq.sh
# Use:
#   - Randomly mutates bases in the input sequence
# Input:
#   -seq: [Required]
#     o Sequence to mutate (as fasta)
#   -num-mutate: [1000]
#     o Number of mutations to make in sequence
#   -seed: [1]
#     o Seed for random number generator
#   -primers: [NA]
#     o Bed file with 2nd column having the start of
#       each primer and the 3rd column with the end of
#       each primer.
#     o Regions not to mutate
# Output:
#   - Prints mutated sequence to sdtout
###########################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

seqStr="";
mutateI=1000;
seedI=1;
primerStr="";

helpStr="
   $(basename "$0") -seq seq.fasta [options...]
   Use:
     - Randomly mutates bases in the input sequence
   Input:
     -seq: [Required]
       o Sequence to mutate (as fasta)
     -num-mutate: [1000]
       o Number of mutations to make in sequence
     -seed: [1]
       o Seed for random number generator
     -primers: [NA]
       o Bed file with 2nd column having the start of
         each non-mutate region (primer) and the
         3rd column with the end of each non-mutate region
         (primer).
   Output:
     - Mutated sequence to stdout
"
       
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Read in user input
#  o sec-02 sub-01:
#    - Get user input
#  o sec-02 sub-02:
#    - Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-02 Sub-01:
#  - Get user input
#**********************************************************

while [ $# -gt 0 ]; do
# Loop: through all user input
   case $1 in
     -seq) seqStr="$2"; shift;;
     -num-mutate) mutateI="$2"; shift;;
     -seed) seedI="$2"; shift;;
     -primers) primerStr="$2"; shift;;
     *) printf "%s\n%s not valid\n" "$helpStr" "$1"; exit;;
   esac

   shift; # move to next paramater
done # Loop: through all user input

#**********************************************************
# Sec-02 Sub-01:
#  - Check user input
#**********************************************************

if [[ "$seqStr" == "" ]]; then
   printf "No sequence file input (use -seq file.fa)\n";
   exit;
fi # If the user input an invalide sequence

if [[ ! -f "$seqStr" ]]; then
   printf "%s is not a valid file\n" "$seqStr";
   exit;
fi # If the user input an invalide sequence

if [[ ! -f "$primerStr" ]]; then
   if [[ "$primerStr" != "" ]]; then
      printf "%s is not a valid file\n" "$primerStr";
      exit;
   fi
fi # If the user input an invalid primer file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Mutate the sequence
#  o sec-03 sub-01:
#    - Prepare the sequence to be mutated
#  o Sec-03 Sub-02:
#    - Initialize variables and set up table to hold
#      regions not mutating (start of mutate sequence)
#  o sec-03 sub-03:
#    - Mutate the sequence
#  o sec-03 sub-04:
#    - Remove spaces and print out the sequence
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-03 Sub-01:
#  - Prepare the sequence to be mutated
#**********************************************************

awk ' # This is here to convert each sequence to one line
      # and to add spaces between each base
      { # MAIN
         if(NR == 1){print $0; next;}; # Header, 1st seq

         # Check if file has more than one sequence
         if($0 ~ /^>/)
         { # If on a new sequence
            print seqStr; # Print out the last sequence
            seqStr = "";  # Blank for new sequence
            print $0;     # Print out the header
            next;         # Move to the sequence entry
         }; # If on a new sequence

         gsub(//, " ", $0);  # Add space between each base
         seqStr = seqStr $0; # merge sequences
      } # MAIN
      END{print seqStr}; # Print out the final sequence
   ' \
   < "$seqStr" |
   awk \
      -v numMutI="$mutateI" \
      -v seedI="$seedI" \
      -v primersStr="$primerStr" \
      '
         #*************************************************
         # Sec-03 Sub-02:
         #  - Initialize variables and set up table to
         #    hold regions not mutating (start of mutate)
         #*************************************************

         BEGIN{
            FS=" ";       # Separate by spaces

            if(seedI == "") srand();
            else srand(seedI); # Set the random seed

            # Set up variable names                     
            if(numMutI == "") numMutI = 1000;
            iMutate = 0;  # Number of bases mutated
            lenSeqI = 0;  # Length of sequence
            mutateI = 0;  # Base to mutate
            subI = 0;     # What to mutate base to
            delete mutatedAryI[1];
            delete doNotMutateAryI[1];

            if(primersStr != "")
            { # if I have primer coordinates to extract
               while((getline < primersStr) > 0)
               { # Loop: Find primer coordiantes
                  iStart = $2;
                  iEnd = $3;

                  while(iStart <= iEnd)
                  { # Loop: Add each position to the table
                     doNotMutateAryI[iStart] = 1;
                     iStart = iStart + 1;
                  } # Loop: Add each position to the table
               } # Loop: Find primer coordiantes
            } # if I have primer coordinates to extract
         }; # BEGIN block

         #*************************************************
         # Sec-03 Sub-03:
         #  - Mutate the sequence
         #*************************************************

         { # MAIN
            if(NR == 1){print $0; next}; # if on the header
            lenSeqI = length($0) / 2;
            # Find sequence length. The /2 is to account
            # for the spaces between every base

            while(iMutate < numMutI)
            { # Loop: Make mutations
               mutateI = 1 + int(rand() * lenSeqI);
               subI = 1 + int(rand() * 3);

               # Check if I have already mutated this base
               # or if this is a no target region
               if(mutatedAryI[mutateI] == 1) continue;
               if(doNotMutateAryI[mutateI] == 1) continue;
               mutatedAryI[mutateI] = 1;

               if($mutateI == "A" || $mutateI == "a")
               { # if I am mutating an A
                  if(subI == 1) $mutateI = "T";
                  if(subI == 2) $mutateI = "G";
                  if(subI == 3) $mutateI = "C";
               } # if I am mutating an A

               else if($mutateI == "T" || $mutateI == "t")
               { # if I am mutating an T
                  if(subI == 1) $mutateI = "A";
                  if(subI == 2) $mutateI = "C";
                  if(subI == 3) $mutateI = "G";
               } # if I am mutating an T

               else if($mutateI == "G" || $mutateI == "G")
               { # if I am mutating an G
                  if(subI == 1) $mutateI = "C";
                  if(subI == 2) $mutateI = "A";
                  if(subI == 3) $mutateI = "T";
               } # if I am mutating an G

               else if($mutateI == "C" || $mutateI == "C")
               { # if I am mutating an C
                  if(subI == 1) $mutateI = "G";
                  if(subI == 2) $mutateI = "T";
                  if(subI == 3) $mutateI = "A";
               } # if I am mutating an C

               else continue; # Was an N (do not mutate)
               
               iMutate = iMutate + 1;
            }; # Loop Make mutations

            exit;
         }; # MAIN BLOCK

         #*************************************************
         # Sec-03 Sub-04:
         #  - Remove spaces and print out the sequence
         #*************************************************

         END{gsub(/ /, "", $0); print $0};
      ';
exit;
