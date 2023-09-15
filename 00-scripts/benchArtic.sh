#!/usr/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TOC: Table Of Contents
# - Benchmark the artic pipeline by mutatint the reference
#   by 0%, 1%, 2%, 5%, or 10%
# o sec-01:
#   - Variable declerations
# o sec-02:
#   - Read, check, and set up variables from user input
# o sec-03:
#   - Prepare to benchmark artic
# o sec-04:
#   - Benchmark the artic pipeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################################
# Name: benchArtric.sh
# Use:
#  - Benchmarks the artic pipeline on input reads and
#    scheme. This will back up and mutate the reference
#    sequence.
# Input:
#  -fastq: [Required]
#    o Fastq file of reads to benchmark artic on.
#  -scheme-dir: [Required]
#    o Directory holding your scheme(s).
#    o Path is: scheme-directory/scheme/V-Version/files
#  -scheme: [Required]
#    o Name of the scheme you are using.
#    o Path is: scheme-directory/scheme/V-Version/files
#  -scheme-version: [Required]
#    o Version number of the scheme to use
#    o Path is: scheme-directory/scheme/VVersion/files
#    o This directory should have at least two files.
#      - file-name.scheme.bed: Has the primer information
#      - file-name.reference.fasta: The reference sequence
#  -model: [Optional]
#    o Model to use with Medaka.
#  -prefix: [Optional]
#    o Prefix to add to ouput files.
#  -do-not-mutate:[scheme-dir/scheme/VVersion/*.scheme.bed]
#    o Regions in the reference genome to not mutate.
#    o This file is tab deliminated and has three columns
#      per entry:
#      - Column 1: Anything but a tab or space
#      - Column 2: First base in the region to not mutate
#      - Column 3: Last base in the region to not mutate
#      - EX: \"MN510897 10 23\"
#  -rep: [Optional]
#    o Number of relpicates to run.
#  -seed: [Optional]
#    o Seed to use for the random number generator
#    o The actual seed input is seed * replicate.
#  -t or -threads: [Optional]
#    o Number of threads to use
# Output:
#  - stdout (terminal):
#    o Tsv file with stats from the alignment 
#  - prefix-log.txt:
#    o All the output send out by artic
#  - prefix-percMutX-rep1-scheme-VX.X
#    o Consensus made by artic, whith the prefix being
#      from the user input
#    o percMutX is the percentage the reference genome was
#      muated by.
#    o repX is the repliciate (X = replicate)
#    o scheme is the name of the scheme used
#    o VX.X is the version of the scheme used
###########################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#  o sec-01 sub-01:
#    - Variables that hold user input
#  o sec-01 sub-02:
#    - Variables that are set in the script
#  o sec-01 sub-03:
#    - help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-01 Sub-01:
#  - Variables that hold user input
#**********************************************************

readsStr="";
schemeStr="";
schemeVerStr="";
schemeDirStr=""
maxRepI=5;
prefixStr="out";
threadsI=3;
startSeedI=1024;
medakaModelStr="r941_prom_high_g344";
percMutateI=0; # Percenate to mutate reference by
simpStatsBl=0; # 1: print out only errors and times

#**********************************************************
# Sec-01 Sub-02:
#  - Variables that are set in the script
#**********************************************************

scriptDirStr="$(dirname "$0")";
logStr="";
suffixStr="";
refStr="";
repI=1;
seedI=1;
lenSeqI=0; # Length of reference sequence (found by script)
numMutateI=0;
primersStr="";
  # For some odd reason bash has a ahrd time with the full
  # path for the pimers file

#**********************************************************
# Sec-01 Sub-03:
#  - help message
#**********************************************************

helpStr="bash $(basename "$0") \
   -fastq reads.fastq \
   -scheme name-of-scheme \
   -scheme-dir /path/to/scheme/directory/ \
   -scheme-version X.X;
Input:
   -fastq: [Required]
     o Fastq file of reads to benchmark artic on.
   -perc-mut: [$percMutateI]
     o Percenate to mutate the reference by (whole number)
     o 100 = 1%, 1 = 0.01%
   -scheme-dir: [Required]
     o Directory holding your scheme(s).
     o Path is: scheme-directory/scheme/V-Version/files
   -scheme: [Required]
     o Name of the scheme you are using.
     o Path is: scheme-directory/scheme/V-Version/files
   -scheme-version: [Required]
     o Version number of the scheme to use
     o Path is: scheme-directory/scheme/VVersion/files
     o This directory should have at least two files.
       - file-name.scheme.bed: Has the primer information
         o Make sure there is only one file with a
           .scheme.bed ending.
       - file-name.reference.fasta: The reference sequence
         o Make sure there is only one file with a
           .reference.fasta ending.
   -model: [$medakaModelStr]
     o Model to use with Medaka.
   -prefix: [$prefixStr]
     o Prefix to add to ouput files.
   -do-not-mutate:[scheme-dir/scheme/VVersion/*.scheme.bed]
     o Regions in the reference genome to not mutate.
     o This file is tab deliminated and has three columns
       per entry:
       - Column 1: Anything but a tab or space
       - Column 2: First base in the region to not mutate
       - Column 3: Last base in the region to not mutate
       - EX: \"MN510897 10 23\"
   -rep: [$maxRepI]
     o Number of relpicates to run.
   -seed: [$startSeedI]
     o Seed to use for the random number generator
     o The actual seed input is seed * replicate.
   -t or -threads: [$threadsI]
     o Number of threads to use
Output:
   - stdout
     o Stats to stdout
   - prefix-log.txt
     o All output from artic
     o This file is appended to, so no worries about
       overwriting
   - prefix-ref-repX-X-perc-dif-con.fa
     o Consensus made by artic, whith the prefix being
       from the user input
     o X-perc is the percentage the reference genome was
       muated by.
     o repX is the replicate number.
";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Read, check, and set up variables from user input
#  o sec-02 sub-01:
#    - Read in user input
#  o sec-02 sub-02:
#    - Check user input
#  o sec-02 sub-03:
#    - Set up variables from user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-01 Sub-01:
#  - Read in user input
#**********************************************************

while [ $# -gt 0 ]; do
# Loop: Process user input
   case $1 in
     -fastq) readsStr="$2"; shift;;
     -perc-mut) percMutateI=$2; shift;;
     -do-not-mutate) primersStr="$2"; shift;;
     -seed) startSeedI="$2"; shift;;
     -rep) maxRepI="$2"; shift;;
     -model) medakaModelStr="$2"; shift;;
     -t) threadsI="$2"; shift;;
     -threads) threadsI="$2"; shift;;
     -scheme-version) schemeVerStr="$2"; shift;;
     -scheme-dir) schemeDirStr="$2"; shift;;
     -scheme) schemeStr="$2"; shift;;
     -prefix) prefixStr="$2"; shift;;
     -h) printf "%s\n" "$helpStr"; exit;;
     --h) printf "%s\n" "$helpStr"; exit;;
     -help) printf "%s\n" "$helpStr"; exit;;
     --help) printf "%s\n" "$helpStr"; exit;;
     help) printf "%s\n" "$helpStr"; exit;;
     *) printf "%s\n%s is invalid\n" "$helpStr" "$2";exit;;
   esac   

   shift; # Move to the next parameter
done # Loop: Process user input

#**********************************************************
# Sec-01 Sub-02:
#  - Check user input
#**********************************************************

if [[ ! -f "$readsStr" ]]; then
   printf "No reads (-fastq) were input\n";
   exit;
fi 

if [[ ! -d "$schemeDirStr/$schemeStr/V$schemeVerStr" ]];
then 
   printf \
      "The scheme directory (%s) does not exist\n" \
      "$schemeDirStr/$schemeStr/V$schemeVerStr";
   exit;
fi

#**********************************************************
# Sec-01 Sub-02:
#  - Set up variables from user input
#**********************************************************


logStr="$prefixStr-log.txt";
refStr="$schemeDirStr/$schemeStr/V$schemeVerStr";

if [[ ! -f "$primersStr" ]]; then
   primersStr="$(find "$refStr" -name "*.scheme.bed")";
fi

refStr="$(find "$refStr" -name "*.reference.fasta")";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - prepare to benchmark artic (Pre loop setup). Get
#    reference length, get read count, and activate conda.
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

backUpRefStr="backup-$(basename "$refStr")";
numReadsI="$(
   sed -n 'p;n;n;n;' "$readsStr" | 
      wc -l |
      awk '{print $1}' \
)"; # Get the number of reads in the fastq file
   
lenSeqI="$(\
   awk '
      { # MAIN
         if(NR == 1){next;}
         if($1 ~ /^>/){exit;}
         seqStr = seqStr $0;
      } # MAIN
      END{gsub(/ /, "", seqStr); print length(seqStr);}
      ' "$refStr"\
)"; # find the reference sequence length

# Activate conda
# Make sure conda is active
source /opt/conda/etc/profile.d/conda.sh;

# build consensus with artic
conda activate artic;

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-04:
#  - Benchmark the artic pipeline
#  o sec-03 sub-01:
#    - Back up and mutate reference sequence (if needed)
#  o sec-03 sub-02:
#    - Run artic
#  o sec-03 sub-03:
#    - Clean up and get stats
#  o sec-03 sub-04:
#    - Print out the stats
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-03 Sub-01:
#  - Back up and mutate reference sequence (if needed)
#**********************************************************

while [[ "$repI" -le "$maxRepI" ]]; do
#Loop: though each replicate
   if [[ "$startSeedI" != "" ]]; then
      seedI="$((startSeedI * repI))";
   else
      seedI="";
   fi # Check if using a seed

   suffixStr="percMut$percMutateI-rep$repI-$schemeStr";
   suffixStr="$suffixStr-V$schemeVerStr";

   numMutateI="$(((percMutateI * lenSeqI) / 10000))";
     # allowing for 0.01% to be used
   mv "$refStr" "$backUpRefStr" || exit;
   
   printf "Mutating genome if needed\n" >&2;
   printf "Mutating genome if needed\n" >> "$logStr";

   if [[ "$percMutateI" -gt 0 ]]; then
   # If I am mutated the reference sequence
      "$scriptDirStr/mutateSeq.sh" \
         -seq "$backUpRefStr" \
         -num-mutate "$numMutateI" \
         -primers "$primersStr" \
         -seed "$seedI" \
         > "$refStr";
   
       mutantRefStatsStr="$(
          bash "$scriptDirStr/getStats.sh" \
             "$backUpRefStr" \
             "$refStr" \
       )"; # Get stats for new reference
   else
      cp "$backUpRefStr" "$refStr";
      mutantRefStatsStr="0 0  0  0";
   fi # If I am mutated the reference sequence

   #****************************************************
   # Sec-03 Sub-02:
   #  - Run artic
   #****************************************************
   
   printf "Running artic\n" >&2;
   printf "Running artic\n" >> "$logStr";

   /usr/bin/time \
         -f "%e\t%U\t%S\t%M\t%P\t" \
         -o "tmp.tsv" \
         -a \
      artic \
         minion \
         --skip-nanopolish \
         --medaka \
         --medaka-model "$medakaModelStr" \
         --threads "$threadsI" \
         --read-file "$readsStr" \
         --scheme-directory "$schemeDirStr" \
         --scheme-version "$schemeVerStr" \
         "$schemeStr" \
         "$prefixStr-$suffixStr" >> "$logStr" 2>&1;
         # saving stderr and stdout output to the log file

   #****************************************************
   # Sec-03 Sub-03:
   #  - Clean up and get stats
   #****************************************************
   
   printf "Getting stats\n" >&2;
   printf "Getting stats\n" >> "$logStr";

   timeStr="$(cat tmp.tsv)";
   rm "tmp.tsv";
   rm "$refStr.fai"; # Remove any indexed files
   mv \
      "$prefixStr-$suffixStr.consensus.fasta" \
      "$prefixStr-$suffixStr-con.fa";
   
   # Remove extra files made by artic
   bash "$scriptDirStr/cleanArticOut.sh" \
      "$prefixStr-$suffixStr";
   
   mv "$backUpRefStr" "$refStr";

   statsStr="$(\
      bash "$scriptDirStr/getStats.sh" \
         "$refStr" \
         "$prefixStr-$suffixStr-con.fa" \
   )";
   
   #****************************************************
   # Sec-03 Sub-04:
   #  - Print out the stats
   #****************************************************

   { # Print out the stats
      printf "artic\t%s\t%s\t%s\t%s" \
         "$repI" \
         "$numReadsI" \
         "$percMutateI" \
         "$numMutateI" \

      printf \
        "\t%s\t%s\t%s\t%s\t%s" \
        "$statsStr" \
        "$mutantRefStatsStr" \
        "$timeStr" \
        "$startSeedI" \
        "$seedI";

      printf \
        "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"\
        "$medakaModelStr" \
        "$lenSeqI" \
        "$refStr" \
        "$readsStr" \
        "$primersStr" \
        "$schemeDirStr" \
        "$schemeStr" \
        "$schemeVerStr" \
        "$prefixStr";
   } >&1; # Printing stats to stdout

   repI="$((repI + 1))";
done #Loop: though each replicate

conda deactivate;
