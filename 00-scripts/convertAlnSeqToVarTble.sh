###########################################################
# Name: convertAlnSeqToVarTble
# Use:
#  - Converts an aligned sequence to a table of differences
# Input:
#  - $1: [Required]
#    o File with both sequences on a single line
#    o The reference should be on the frist line and start
#      with "Ref:  "
#    o The query should be on the frist line and start
#      with "Query:  "
# Output:
#  - Prints a pipe table to stdout
###########################################################

if [[ ! -f "$1" ]]; then
  printf "Input file (%s) can not be opened\n" "$1";
  exit;
fi # if an invalid file was input

# Using sed to separate each base by spaces. This allows
# awk to easily handle the sequences
sed \
    -n \
    '/Ref:/,/Query:/{s/.*:  *//g; s/\(.\)/\1 /g; p;}' \
    < "$1" |
  awk '
     { # MAIN
       # Preparing the table of reference bases to use
       # in the comparison (first entry)

       if(NR == 1)
       { # if setting up the reference table
         printf "| %-9s | Query | Ref |\n", "Base";
         printf "|:---------:|:-----:|:---:|\n", "Base";
         for(baseI = 1; baseI <= NF; ++baseI) refAry[baseI] = $baseI;
         next; # Next line is the query sequence
       } # if setting up the reference table

       # Compare the query to the reference and print
       # differences

       for(baseI = 1; baseI <= NF; ++baseI)
       { # comparing the query base
         if($baseI == "N") continue;
         if(refAry[baseI] == "N") continue;

         else if($baseI == "-")
           printf "| %-9s |   -   |  %s  |\n", baseI, refAry[baseI]

         else if($baseI != refAry[baseI])
           printf "| %-9s |   %s   |  %s  |\n", baseI, $baseI, refAry[baseI];
       } # comparing the query base

       exit;
     } # MAIN
  '
