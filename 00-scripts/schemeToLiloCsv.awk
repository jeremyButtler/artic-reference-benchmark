###########################################################
# Name: schemeToLiloCsv
# call: awk -f schemeToLiloCsv.awk primers.scheme.bed
# Use:
#  - Converts a scheme file from the artic pipeline into
#    a csv file that LILO can read
# Input:
#  - file:
#    o Bed file to convert to csv
# Output:
#  - stdout
#    o Prints converted csv to stdout
###########################################################

{ # MAIN
   # Check if this is an alternate primer
   if($4 ~ /_[aA][Ll][Tt]/) next;

   # Carriage returns '\r' are invisible and mess up the
   # print statement
   gsub(/\r/, "", $0);

   # Get the primer pair number
   tmpStr = $4;
   sub(/[^_]*_/, "", tmpStr);
   sub(/_.*/, "", tmpStr);
   primNumI = int(tmpStr);

   if($4 ~ /LEFT/)
   { # If this is a left primer
      # Make the primer pair name
      nameStr = $4;
      sub(/_.*/, "", nameStr);

      primLeftAryStr[primNumI] = $7;
      nameAryStr[primNumI] = nameStr;
      ++numPrimFoundI;
   } # If this is a left primer

   else primRightAryStr[primNumI] = $7;
} # MAIN

END{ # END
   for(iPrim = 1; iPrim <= numPrimFoundI; ++iPrim)
   { # Loop: Print out all primer pairs
      printf "%s-%03i,%03iF,%s,%03iR,%s\n",
         nameAryStr[iPrim],
         iPrim,
         iPrim,
         primLeftAryStr[iPrim],
         iPrim,
         primRightAryStr[iPrim];
   } # Loop: Print out all primer pairs
} # END
