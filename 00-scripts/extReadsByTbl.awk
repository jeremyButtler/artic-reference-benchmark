###########################################################
# Name: extReadsByTbl
# Use:
#  - Extracts reads from an extraction table and puts the
#    reads for each entry in a separare fastq file.
# Input:
#   - file.sam: [Required]
#     - Sam file to process. Input after you provied all
#       the options
#   -v extTbl=extractionTable.md [Required]
#     - Table to use for read extractiion. This should be
#       made by 00-scripts/readLenPosTbl.sh and have the
#       rows not intrested in deleted.
#   -v prefix="out"
#     - Prefix to name the output files by
#   -v lenBinSizeI=100
#     - Nearest value to round the lengths to. Default is
#       to the 100th value. This should be the same as used
#       for 00-scripts/readLenPosTbl.sh
#   -v posBinSizeI=100
#     - Nearest value to round the position the read mapped
#       to the reference. Default is to the 100th value.
#       This should be the same as used for
#       00-scripts/readLenPosTbl.sh
#   -v help=""
#     - Prints the help message when not ""
# Output:
#  - Prints a fastq file for each entry in the\n";
#    input table (-v extTbl).\n";
#  - Fastq file name is\n";
#    prefix-mappedRefernce-start-readLength\n";
###########################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SOF: Start Of File
#  o fun-01:
#    - Print the help message
#  o begin:
#    - Check user input and convert the extraction table
#      to an array to extract with
#  o main:
#    - Extract reads using the extraction array
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----------------------------------------------------------
# Output:
#  - Prints the help message to the screen
#----------------------------------------------------------
function helpMesg()
{ #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  # Fun-01 TOC:
  #  - Prints out the help message
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/

  printf "awk -f extReadsByTbl.awk -v extTbl=table.md";
  printf " file.sam\n";

  printf "Use:\n";
  printf "  - Extracts reads from a sam file to fastq\n";
  printf "    files using an table from\n";
  printf "    readLenPosTbl.sh\n";

  printf "Input:\n";
  printf "  file.sam: [Required]\n";
  printf "    - Sam file to extract reads from.\n";
  printf "  -v extTbl=table.md [Required]\n";
  printf "    - Table made by readLenPosTbl.sh to\n";
  printf "      extract reads with.\n";
  printf "    - Delete entries you want to ignore.\n";
  printf "  -v prefix=\"out\"\n";
  printf "    - Prefix to add to output file names.\n";
  printf "  -v lenBinSizeI=100\n";
  printf "    - Nearest value to round read lengths\n";
  printf "      down to. This should be the same value\n";
  printf "      as -len-binsize used to make the table\n";
  printf "      with readLenPosTbl.sh.\n";
  printf "  -v posBinSizeI=100\n";
  printf "    - Nearest value to round the read mapping\n";
  printf "      position on the reference down to. \n";
  printf "      This should be the same value as \n";
  printf "      -pos-binsize used to make the table\n";
  printf "      with readLenPosTbl.sh.\n";
  printf "  -v help=\"\"\n";
  printf "    - Will print this help message if the \n";
  printf "      is changed from \"\"\n";
  printf "Output:\n"
  printf "  - Prints a fastq file for each entry in the\n";
  printf "    input table (-v extTbl).\n";
  printf "  - Fastq file name is\n";
  printf "    prefix-mappedRefernce-start-readLength\n";
} # helpMesg

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# BEGIN:
#  - Check user input and convert the extraction table
#    to an array to extract with
#  o begin sub-01:
#    - Check user input
#  o begin sub-02:
#    - Convert input table to an extraction array
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
#  o BEGIN Sub-01:
#    - Check user input
#**********************************************************

BEGIN{
  if(help != "")
  { # If the user wanted the help message
    helpMesg();
    exit;
  } # If the user wanted the help message

  if(extTbl == "")
  { # if no extraction table supplied
    helpMesg();

    printf "\n\n";
    printf "This program needs a table of reads to\n"
    printf "  extract (-v extTabl=table.md). This table\n";
    printf "  is made with readLenPosTbl.sh\n"

    exit;
  } # if no extraction table supplied

  # Set up default values
  if(lenBinSizeI == "") lenBinSizeI = 100;
  if(posBinSizeI == "") posBinSizeI=100;

  #********************************************************
  # BEGIN Sub-02:
  #  - Convert input table to an extraction array
  #********************************************************

  # Move past the header/separator in the pipe table
  getline < extTbl;
  getline < extTbl;

  FS="|";

  while((getline < extTbl) > 0)
  { # While I have entries to add to the extraction array
    gsub(/ /, "", $0); # remove spaces
    extAry[$3][$4][$5] = 1;
  } # While I have entries to add to the extraction array

  FS="\t";
} # BEGIN block

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MAIN:
#  - Extract reads using the extraction array
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

{ # MAIN
  if($1 ~ /^@/){next;}; # Skip headers

  # Get read length
  seqLenI = length($10);

  # Convert the starting position and read length to
  # nearest lenBinSizeI (100th) value. Will discard
  # values less than lenBinSizeI.

  seqLenI = seqLenI - (seqLenI % lenBinSizeI);
  refPosI = $4 - ($4 % posBinSizeI);
  refIdStr = $3;

  #print extAry[refIdStr][refPosI][seqLenI];
  
  if(extAry[refIdStr][refPosI][seqLenI] == 1)
  { # If this is worth printing
    outFile = prefix "-" $3; # add in the reference
    outFile = outFile "-pos" refPosI "-" "len" seqLenI;
    outFile = outFile ".fastq";

    printf "@%s %s %s\n", $1, startPosI, $3 >> outFile;
    printf "%s\n+\n%s\n", $10, $11 >> outFile;
      # $1 is the read id (output as id-positionOnRef)
      # $3 is the reference id
      # $10 is sequence entry
      # $11 is Q-score entry
  } # If this is worth printing
} # MAIN
