###########################################################
# Name: readLenPosTbl.sh
# Use:
#  - Makes a table that has the mapped reference, read
#    position (in hundreds), and read length (in hundreds)
# Input:
#  -sam: [Required]
#    o sam file to get read lengths and read positions for
#  -len-binsize: [100]
#    o Round (floor) read length to nearest bin size
#  -pos-binsize: [100]
#    o Round (floor) read position on reference to nearest
#      bin size
#  -min-length: [200]
#    o Mininmum read length to keep a table entry
#  -min-reads: [200]
#    o Mininmum number of reads needed to keep a table
#      entry
# Output:
#  - Pipe table (markdown format) with reads to stdout
# Requires:
#  - Minimap2
#  - gawk
#  - trimSamFile from find co-infections
#    - https://github.com/jeremybuttler/find--Co-infections
###########################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TOC: Table Of Contents
#   sec-01 Sub-01:
#    - User input variables
#   sec-02:
#    - Read in and check user input
#   sec-03:
#    - Convert sam file to table (one big pipe statment)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
# Sec-01:
#  - Variable declerations
#  o sec-01 sub-01:
#    - User input variables
#  o sec-01 sub-02:
#    - Help message
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>/

#**********************************************************
# Sec-01 Sub-01:
#  - User input variables
#**********************************************************

samStr="";       # Referennce sequence
lenReadBinI=100; # Round read length to nearest 100
posBinI=100;     # Round read position on ref to nearest100
minLenI=200;     # Min read length to keep table entry
minReadsI=200;   # Min numboer of reads to keep table entry

#**********************************************************
# Sec-01 Sub-02:
#  - Help message
#**********************************************************

helpStr="
  bash $(basename "$0") -sam file.sam [options]
  Use:
   - Makes a table that has the mapped reference, read
     position (in hundreds), and read length (in hundreds)
  Input:
   -sam: [Required]
     o sam file to get read lengths and read positions for
   -len-binsize: [100]
     o Round (floor) read length to nearest bin size
   -pos-binsize: [100]
     o Round (floor) read position on reference to nearest
       bin size
   -min-length: [200]
     o Mininmum read length to keep a table entry
   -min-reads: [200]
     o Mininmum number of reads needed to keep a table
       entry
  Output:
   - Pipe table (markdown format) with reads to stdout
  Requires:
   - Minimap2
   - gawk
   - trimSamFile from find co-infections
     - https://github.com/jeremybuttler/find--Co-infections
";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
# Sec-02:
#  - Read in and check user input
#  o sec-02 sub-01:
#    - Read in user input
#  o sec-02 sub-02:
#    - Check user input
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>/

#**********************************************************
# Sec-02 Sub-01:
#  - Read in user input
#**********************************************************

while [ $# -gt 0 ]; do
  case $1 in
    -sam) samStr="$2"; shift;;
    -len-binsize) lenReadBinI="$2"; shift;;
    -pos-binsize) posBinI="$2"; shift;;
    -min-length) minLenI="$2"; shift;;
    -min-reads) minReadsI="$2"; shift;;

    -h) printf "%s\n" "$helpStr"; exit;;
    --h) printf "%s\n" "$helpStr"; exit;;
    -help) printf "%s\n" "$helpStr"; exit;;
    --help) printf "%s\n" "$helpStr"; exit;;

    *) printf "%s\n" "$helpStr";
       printf "%s is an invalid parameter\n" "$1";
       exit;;
  esac

  shift; # move to the next argument
done;  

#**********************************************************
# Sec-02 Sub-02:
#  - Check user input
#**********************************************************

if [[ ! -f "$samStr" ]]; then
  printf "Input file (%s) does not exist\n" "$samStr";
  exit
fi # Check if the sam file exists

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
# Sec-03:
#  - Convert sam file to table (one big pipe statment)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>/

awk \
    -v lenReadBinI="$lenReadBinI" \
    -v posBinI="$posBinI" \
    '
    { # MAIN
      if($1 ~ /^@/){next;}; # Skip headers

      # Get read length
      seqLenInt = length($10);

      # Convert the starting position and read length to
      # nearest lenReadBinI (100th) value. Will discard
      # values less than lenReadBinI.

      seqLenInt = seqLenInt - (seqLenInt % lenReadBinI);
      $4 = $4 - ($4 % posBinI);

      # Print out reference id, starting position, & length
      print $3, $4, seqLenInt;
    }' \
    < "$samStr" |
  sort \
    -n \
    -k 2,3 |
  uniq -c |
  awk \
    -v minReadsI="$minReadsI" \
    -v minLenI="$minLenI" \
    '
      BEGIN{
        OFS="\t";  # ensure tab output

        # print the header
        printf "| No. Reads |";
        printf "        Reference        | Position |"
        printf " Length |\n";

        # print the separator for the header
        printf "|:---------:|";
        printf ":-----------------------:|:--------:|"
        printf ":------:|\n";
      } # BEGIN block

      { # MAIN
        # discard entries with under x (200) reads
        if($1 < minReadsI){next;};

        # discard reads shorter than 200 nucleotides
        if($4 < minLenI){next;};

        # Print out the table entry                 
        printf "| %-9s | %-23s |", $1, $2;
        printf " %-8s | %-6s |\n", $3, $4;
      } # MAIN
'
