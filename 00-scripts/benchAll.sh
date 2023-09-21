#!/usr/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TOC: Table Of Contents
#  - Benchmark the artic pipeline
#  o sec-01:
#     - Variable declerations
#  o sec-02:
#     - Get and check user input
#  o sec-03:
#     - Pre-loop/benchmarking setup
#  o sec-04:
#     - Run benchmark
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################################
# Name: benchAll.sh
# Use:
#  - Runs benchArtic.sh for all read depths, lengths,
#    percent ref different levels, and replicates
# General Input:
#   -fastq-dir: [05-amplicon-reads]
#     o Directory of fastq files. Each fastq file has a
#       separate amplicon
#   -depth: [30 50 100 200 500 1000 1200]
#     o Read depths to benchmark at.
#     o This can be input multiple times.
#   -perc-mut-ref: [0 100 200 500 1000]
#     o Percentage to mutate the reference by
#       - 100 = 1%, 1 = 0.01%
#     o This can be input multiple times.
#   -model: [r941_prom_high_g344]
#     o Model to use with Medaka.
#   -prefix: [artic]
#     o Prefix to add to ouput files.
#   -do-not-mutate:
#     - [06-alt-schemes/06-no-mutate-regions.tsv]
#     o Regions in the reference genome to not mutate.
#     o This file is tab deliminated and has three columns
#       per entry:
#       - Column 1: Anything but a tab or space
#       - Column 2: First base in the region to not mutate
#       - Column 3: Last base in the region to not mutate
#       - EX: \"MN510897 10 23\"
#   -rep: [5]
#     o Number of relpicates to run.
#   -seed: [1024]
#     o Seed to use for the random number generator
#     o The actual seed input is seed * replicate.
#   -t or -threads: [3]
#     o Number of threads to use
# Artic Unique Input:
#    -scheme-dir: [06-alt-schemes/06-schemes]
#      o Directory holding your scheme for artic.
#      o Path is: scheme-directory/scheme/V-Version/files
#    -scheme: [SARS-CoV-2]
#      o Name of the scheme you are using with artic.
#      o Path is: scheme-directory/scheme/V-Version/files
#    -scheme-version: [V4.1.1500 V4.1.5000 V4.1.11000 V4.1]
#      o Version number of the scheme to use
#      o This can be input multiple times.
#      o Path is: scheme-directory/scheme/VVersion/files
#      o This directory should have at least two files.
#        - file-name.scheme.bed: Has the primer information
#          o Make sure there is only one file with a
#            .scheme.bed ending.
#        - file-name.reference.fasta: The reference sequence
#          o Make sure there is only one file with a
#            .reference.fasta ending.
# Output:
#  - prefix-stats.tsv
#    o Tsv file with stats from the alignment 
#  - prefix-repX-depthX.log
#    - Log holding artic output for each replicate and
#      depth
#  - prefix-repX-depthX-percMutX-rep1-scheme-VX.X-con.fasta
#    o Consensus made by artic, whith the prefix being
#      from the user input
#    o repX is the replicate number (X = replicate)
#    o depthX is the read depth (X = read depth)
#    o percMutX is the percentage the reference genome was
#      muated by.
#    o rep1 is added by benchArtic.sh
#    o scheme is the name of the scheme used
#    o VX.X is the version of the scheme used
###########################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#  o sec-01 sub-01:
#    - Variables that are globally used for input
#  o sec-01 sub-02:
#    - Variables unique to artic
#  o sec-01 sub-03:
#    - Variables unique to buildCon
#  o sec-01 sub-04:
#    - Booleans to disable tests
#  o sec-01 sub-05:
#    - Help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-01 Sub-01:
#  - Variables that are globally used for input
#**********************************************************

ampDirStr="05-amplicon-reads";
  # These are the reads I will subsample from. Each
  # amplicon is split up and only have overlaps at the
  # ends. This is so that when I say 100x read depth, I
  # mean 100x read depth for everything.
seedI=1024;         # Seed for random number generator
repI=5;             # Number of replicates to run
prefixStr="artic";  # What to name output files
threadsI=3;         # Number of threads to use
medakaModelStr="r941_prom_high_g344"; # Model for medaka
depthAryI=(30 50 100 300 500 1000);
  # Read depths to benchmark at
percMutAryI=(0 100 200 500 1000);
  # Amount to mutate the refernce by
  # 0 = 0%, 100 = 1%, 200 = 2%, 500 = 5%, 1000 = 10%
doNotMutateStr="06-alt-schemes/06-no-mutate-regions.tsv";
  # Regions to not mutate

#**********************************************************
# Sec-01 Sub-02:
#  - Variables unique to artic
#**********************************************************

# Settings for Primal scheme
schemeDirStr="06-alt-schemes/06-schemes/";
schemeStr="SARS-CoV-2";
verAryStr=("4.1.1500" "4.1.5000" "4.1.11000" "4.1");
  # V4.1.1500 builds an 1700 nucleotide consensus
  # V4.1.5000 builds an 5700 nucleotide consensus
  # V4.1.11000 builds an 11800 nucleotide consensus
  # V4.1 Builds the full length consensus

#**********************************************************
# Sec-01 Sub-03:
#  - Variables unique to buildCon
#**********************************************************

minLenI=320;
minDepthI=100; # this will be 50% for depths under 100
runDepthI=0;   # acutal depth running
useMedakaBl=1;

#**********************************************************
# Sec-01 Sub-04:
#  - Booleans to disable tests
#**********************************************************

articlBl=1;     # 1: Run the artic test
buildConMajconMedakaBl=1; # 1: run the buildCon test with medaka
buildConMajconBl=1;  # 1: Run the buildCon test without medaka
buildConMedakaBl=1; # 1: Run buildCon test without majcon
liloBl=1;      # 1: Run LILO tests
ivarBl=1;      # run ivar
ivarTrimBl=1;      # run ivar with trimming

#**********************************************************
# Sec-01 Sub-05:
#  - Variables not changed by user
#**********************************************************

scriptDirStr="$(dirname "$0")";
defVerBl=0;     # 0: User has not changed version numbers
defDepthBl=0;   # 0: User has not changed the read depths
defPercMutBl=0; # 0: User has not changed the mutations
statsFileStr=""; # Holds the ouput stats
statsLineStr=""; # line to print to statsFileStr
tmpStr="";      # For holding temporary data/file names
tmpI=0;         # For holding temporary integers
errBl=0;        # 0: no errors
readsStr="";    # Holds the subsampled reads
inputSeedI=0;   # seed given to scripts (multple of seedI)
conPrefStr=""   # prefix for the consensus name
mutantRefStatsStr=""; # Number of differences in mutant ref
subSampDirStr=""; # Directory subsamples saved in
trueRefStr="";

testAleadyRunBl=0;
    # This holds if the test has already been run (1)

#**********************************************************
# Sec-01 Sub-06:
#  - Help message
#**********************************************************

helpStr="bash $(basename "$0") \
General Input:
   -fastq-dir: [$ampDirStr]
     o Directory of fastq files. Each fastq file has a
       separate amplicon
   -depth: [${depthAryI[*]}]
     o Read depths to benchmark at.
     o This can be input multiple times.
   -perc-mut-ref: [${percMutAryI[*]}]
     o Percentage to mutate the reference by
       - 100 = 1%, 1 = 0.01%
     o This can be input multiple times.
   -model: [$medakaModelStr]
     o Model to use with Medaka.
   -prefix: [$prefixStr]
     o Prefix to add to ouput files.
   -do-not-mutate:[$doNotMutateStr]
     o Regions in the reference genome to not mutate.
     o This file is tab deliminated and has three columns
       per entry:
       - Column 1: Anything but a tab or space
       - Column 2: First base in the region to not mutate
       - Column 3: Last base in the region to not mutate
       - EX: \"MN510897 10 23\"
   -rep: [$repI]
     o Number of relpicates to run.
   -seed: [$seedI]
     o Seed to use for the random number generator
     o The actual seed input is seed * replicate.
   -t or -threads: [$threadsI]
     o Number of threads to use
Disable tests:
   -no-artic: [-artic]
     o Disables the artic pipeline tests
     o Artic can be enabled by -artic
   -no-LILO: [-LILO]
     o Disables the LILO tests
     o LILO can be enabled with -LILO
   -no-buildCon-majcon-medaka: [-buildCon-majcon-medaka]
     o Disabled buildCon using majcon and Medaka test
     o This can be enabled with -buildCon-majcon-medaka
   -no-buildCon-majcon: [-buildCon-majcon]
     o Disabled buildCon with only majcon test
     o This can be enabled with -buildCon-majcon
   -no-buildCon-medaka: [-buildCon-medaka]
     o Disable buildCon with only Medaka test
     o This can be enabled with -buildCon-medaka
   -no-ivar: [-ivar]
     o Disable the ivar without trimming test
     o Enable ivar without trimming with -ivar
   -no-ivar-trim: [-ivar-trim]
     o Disable the ivar trimming test
     o Enable ivar with trimming with -ivar-trim
Artic Unique Input:
   -scheme-dir: [$schemeDirStr]
     o Directory holding your scheme for artic.
     o Path is: scheme-directory/scheme/V-Version/files
   -scheme: [$schemeStr]
     o Name of the scheme you are using with artic.
     o Path is: scheme-directory/scheme/V-Version/files
   -scheme-version: [${verAryStr[*]}]
     o Version number of the scheme to use
     o This can be input multiple times.
     o Path is: scheme-directory/scheme/VVersion/files
     o This directory should have at least two files.
       - file-name.scheme.bed: Has the primer information
         o Make sure there is only one file with a
           .scheme.bed ending.
       - file-name.reference.fasta: The reference sequence
         o Make sure there is only one file with a
           .reference.fasta ending.
BuidCon Input:
   -minDepth: [$minDepthI]
     o Min read depth to required for buld con to build a
       consensus
     o If amplicon read depth is beneath this, then it is
       set to amplicon depth / 3.
   -min-len: [$minLenI]
     o Min length to keep reads/consensus.
Output:
    - prefix-stats.tsv
      o Tsv file with stats from the alignment 
    - prefix-repX-depthX.log
      - Log holding artic output for each replicate and
        depth
    -prefix-repX-depthX-percMutX-rep1-scheme-VX.X-con.fasta
      o Consensus made by artic, whith the prefix being
        from the user input
      o repX is the replicate number (X = replicate)
      o depthX is the read depth (X = read depth)
      o percMutX is the percentage the reference genome was
        muated by.
      o rep1 is added by benchArtic.sh
      o scheme is the name of the scheme used
      o VX.X is the version of the scheme used
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

while [ $# -gt 0 ]; do
# Loop: While I have user input to check
   case $1 in
     -fastq-dir) ampDirStr="$2"; shift;;
     -perc-mut-ref)
        if [[ "$defPercMutBl" -eq 0 ]]; then
           percMutAryI=("$2");
           defPercMutBl=1;
        else
           percMutAryI+=("$2"); # append a new depth
        fi # Check if this is the first read depth input
        shift;; # Move to the argument ($2)
     -depth)
        if [[ "$defDepthBl" -eq 0 ]]; then
           depthAryI=("$2");
           defDepthBl=1;
        else
           depthAryI+=("$2"); # append a new depth
        fi # Check if this is the first read depth input
        shift;; # Move to the argument ($2)
     -do-not-mutate) doNotMutateStr="$2"; shift;;
     -seed) seedI="$2"; shift;;
     -rep) repI="$2"; shift;;
     -model) medakaModelStr="$2"; shift;;
     -t) threadsI="$2"; shift;;
     -threads) threadsI="$2"; shift;;
     -minDepth) minDepthI="$2"; shift;;
     -min-len) minLenI="$2"; shift;;
     -scheme-dir) schemeDirStr="$2"; shift;;
     -scheme) schemeStr="$2"; shift;;
     -scheme-version) 
        if [[ "$defVerBl" -eq 0 ]]; then
           verAryStr=("$2");
           defVerBl=1;
        else
           verAryStr+=("$2"); # append to an array
        fi # Check if this is the first version input
        shift;; # Move to the argument ($2)
     -prefix) prefixStr="$2"; shift;;
     -no-artic) articlBl=0;;
     -artic) articlBl=1;;
     -no-LILO) liloBl=0;;
     -LILO) liloBl=1;;
     -no-buildCon-majcon) buildConMajconBl=0;;
     -buildCon-majcon) buildConMajconBl=1;;
     -no-buildCon-majcon-medaka) buildConMajconMedakaBl=0;;
     -buildCon-majcon-medaka) buildConMajconMedakaBl=1;;
     -no-buildCon-medaka) buildConMedakaBl=0;;
     -buildCon-medaka) buildConMedakaBl=1;;
     -ivar) ivarBl=1;;
     -no-ivar) ivarBl=0;;
     -ivar-trim) ivarTrimBl=1;;
     -no-ivar-trim) ivarTrimBl=0;;
     -h) printf "%s\n" "$helpStr"; exit;;
     --h) printf "%s\n" "$helpStr"; exit;;
     -help) printf "%s\n" "$helpStr"; exit;;
     --help) printf "%s\n" "$helpStr"; exit;;
     help) printf "%s\n" "$helpStr"; exit;;
     *) printf "%s\n%s is invalid\n" "$helpStr" "$1";exit;;
   esac

   shift; # Move to next parameter
done # Loop: While I have user input to check

#**********************************************************
# Sec-02 Sub-02:
#  - Check user input
#  o sec-02 sub-02 cat-01:
#    - Check if have amplicons
#  o Sec-02 sub-02 cat-02:
#    - Check if have valid schemes
#  o sec-02 sub-02 cat-03:
#    - Check if I have a valid Medaka model
#  o sec-02 sub-02 cat-04:
#    - Check if user input do not mutate file & if is valid
#  o sec-02 sub-02 cat-05:
#    - Check if user input valid percentages
#  o sec-02 sub-02 cat-06:
#    - Check if user input valid depths
#  o sec-02 sub-02 cat-07:
#    - Check if user input a valid seed
#  o sec-02 sub-02 cat-08:
#    - Check if user input a valid replicate number
#  o sec-02 sub-02 cat-09:
#    - Check if user input a valid number for threads
#  o sec-02 sub-02 cat-10:
#    - Check if user input a valid prefix
#  o sec-02 sub-02 cat-11:
#    - Check if I have invalid input for schemes/amplicons
#**********************************************************

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-01:
#  - Check if have amplicons
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tmpI=0;
errBl=0;

if [[ -d "$ampDirStr" ]]; then
# If the amplicon directory does exist
   for strFq in "$ampDirStr/"*.f*q; do
   # Loop: Check if have fastq files to subsample
      if [[ ! -f "$strFq" ]] ;then
         continue;
      fi # If the null case
   
      tmpI=$((tmpI + 1));
   done # Loop: Check if have fastq files to subsample
   
   if [[ "$tmpI" -le 0 ]]; then
      printf "No fastq files in %s\n" "$ampDirStr" >&2;
      errBl=1; # mark that I have invalid input
   fi # IF there were no fastq files
else
   printf \
       " -fastq-dir %s does not exist\n" \
       "$ampDirStr" \
     >&2;
   errBl=1; # mark that I have invalid input
fi # If the amplicon directory does exist

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-02:
#  - Check if have valid schemes
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for strVer in ${verAryStr[*]}; do
# Loop: Check if all schemes are valid
   if [[ ! -d "$schemeDirStr/$schemeStr/V$strVer" ]]; then
      printf \
          "The scheme directory (%s) does not exist\n" \
          "$schemeDirStr/$schemeStr/V$strVer" \
        >&2;
      errBl=1; # mark that I have invalid input
      continue;
   fi # If I had an invalid scheme
done # Loop: Check if all schemes are valid

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-03:
#  - Check if I have a valid Medaka model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Activate conda
# Make sure conda is active
#source /opt/conda/etc/profile.d/conda.sh;
#  # shellcheck can not follow this, but it is valid
#
## build consensus with artic
#conda activate artic;
#
#tmpStr="$(\
#   medaka tools list_models |
#     grep "$medakaModelStr"
#)";
#
#if [[ "$tmpStr" == "" ]]; then
#   printf \
#       "%s is not a valid Medaka model\n" \
#       "$medakaModelStr" \
#     >&2;
#   errBl=1; # mark that I have invalid input
#fi # If an invalid model was input
#
#conda deactivate;

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-04:
#  - Check if user input do not mutate file & if is valid
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$doNotMutateStr" != "" ]]; then
   if [[ "$doNotMutateStr" == "\"\"" ]]; then
      {
         printf " -do-not-mutate file not input\n";
         printf "The entire reference will be mutated\n"
      } >&2;
   else
      if [[ ! -f "$doNotMutateStr" ]]; then
         {
            printf "file of regions not to mutate";
            printf \
               " (-do-not-mutate %s) does not exist\n" \
               "$doNotMutateStr"; 
            printf "For no file use -do-not-mutate \"\"\n";
         } >&2;

         errBl=1; # mark that I have invalid input
      fi # If do not mutate should exist, but does not
   fi # If is blank input

else
   {
      printf " -do-not-mutate file not input\n";
      printf "The entire reference will be mutated\n"
   } >&2;
fi # check if the user provied a mutate file

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-05:
#  - Check if user input valid percentages
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tmpI=0;
tmpStr="^[0-9][0-9]*$"; # Regular expresion for numeric

for strPerc in ${percMutAryI[*]}; do
   if [[ "$strPerc" == "" ]]; then
      continue;
   fi # If is a blank value

   if ! [[ "$strPerc" =~ $tmpStr ]]; then
      printf \
         " -perc-mut-ref %s is not a number\n" \
         "$strPerc" \
         >&2;
      errBl=1;
      continue;
   fi # If the values is non numeric

   if [[ "$strPerc" -gt 10000 ]]; then
      {
         printf " -perc-mut-ref %s is >" "$strPerc";
         printf " max value of 10000 (100%%)\n"
      } >&2;
      errBl=1;
      continue;
   fi # If the value is to large

   tmpI=$((tmpI + 1));
done

if [[ "$tmpI" -le 0 ]]; then
   printf "Blank values input for -perc-mut-ref\n" >&2;
   errBl=1;
fi # If no mutation percentages were input

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-06:
#  - Check if user input valid depths
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tmpI=0;
tmpStr="^[0-9][0-9]*$"; # Regular expresion for numeric

for strDepth in ${depthAryI[*]}; do
   if [[ "$strDepth" == "" ]]; then
      continue;
   fi # If is a blank value

   if ! [[ "$strDepth" =~ $tmpStr ]]; then
      printf \
         " -depth %s is not a number\n" \
         "$strDepth" \
       >&2;
      errBl=1;
      continue;
   fi # If the values is non numeric

   tmpI=$((tmpI + 1));
done

if [[ "$tmpI" -le 0 ]]; then
   printf "Blank values input for -depth\n" >&2;
   errBl=1;
fi # If no mutation percentages were input

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-07:
#  - Check if user input a valid seed
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tmpStr="^[0-9][0-9]*$"; # Regular expresion for numeric

if ! [[ "$seedI" =~ $tmpStr ]]; then
   if [[ "$seedI" == "" ]]; then
      printf "No seed input, results will be random\n" >&2;
   elif [[ "$seedI" == "\"\"" ]]; then
      printf "No seed input, results will be random\n" >&2;
      seedI="";
   else
      printf " -seed %s is not a number\n" "$seedI" >&2;
      errBl=1;
   fi # If user wanted a blank (no) seed
fi # If the value is non numeric

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-08:
#  - Check if user input a valid replicate number
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tmpStr="^[0-9][0-9]*$"; # Regular expresion for numeric

if ! [[ "$repI" =~ $tmpStr ]]; then
   printf " -rep %s is not a number\n", "$repI" >&2;
   errBl=1;
fi # If the value is non numeric

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-09:
#  - Check if user input a valid number for threads
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tmpStr="^[0-9][0-9]*$"; # Regular expresion for numeric

if ! [[ "$threadsI" =~ $tmpStr ]]; then
   printf \
       " -t %s (threads) is not a number\n", \
       "$threadsI" \
     >&2;
   errBl=1;
fi # If the value is non numeric

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-10:
#  - Check if user input a valid prefix
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$prefixStr" == "" ]]; then
   printf "Blank prefix (-prefix) input\n" >&2;
   printf "Setting prefix to artic-bench\n" >&2;
   prefixStr="artic-bench";
fi # If the prefix was blank

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-11:
#  - Check if I have invalid input for schemes/amplicons
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# If not all schemes were valid
if [[ "$errBl" -gt 0 ]]; then exit; fi

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Pre-loop/benchmarking setup
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

statsFileStr="$prefixStr-stats.tsv";

# This header is large, but also allows complete
# replicatoni of the settings.
if [[ ! -f "$statsFileStr" ]]; then
   { # Print out the header
      printf "program\trep\tdepth\tnoReads\tpercMutate";
      printf "\ttargMut\tlenRef\trefMatch";
      printf "\trefSnp\trefIns\trefDel\trefNs\tmatches";
      printf "\tsnp\tins\tdel\tNs";
      printf "\telpTime\tsysTime\tUserTime\tmemKb";
      printf "\tpercCPU\tseed\tusedSeed\tmedakaModel";
      printf "\tusedMedaka\tusedMajCon\tminLen\tminDepth";
      printf "\tfastq\tschemeDir\tscheme\tschemeVer";
      printf "\tprefix\tconDef\tconScaffoldBuilder";
      printf "\tconStich\n";
   } > "$statsFileStr"
fi

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-04:
#  - Run benchmark
#  o sec-04 sub-01:
#    - Check if using a seed or user wanted no seed
#  o sec-04 sub-02:
#    - Subsample reads for each amplicon
#  o sec-04 sub-03:
#    - Find the reference and primers for each scheme
#  o sec-04 sub-04:
#    - Mutate the reference
#  o sec-04 sub-05:
#    - Get stats for artic
#  o sec-04 sub-06:
#    - Get stats for buildCon
#  o sec-04 sub-07:
#    - Get stats for buildCon no medaka
#  o sec-04 sub-08:
#    - Stats: buildCon, no medaka, ivar polish
#  o sec-04 sub-09:
#    - stats for buildCon just medaka (no majcon)
#  o sec-04 sub-10:
#    - Run LILO
#  o sec-04 sub-11:
#    - Get stats for fixed LILO (auto) + ivar
#  o sec-04 sub-12:
#    - Get stats for fixed LILO (auto)
#  o sec-04 sub-13:
#    - separate scaffold_builer LILO + ivar
#  o sec-04 sub-14:
#    - stats for separate scaffold_builer LILO
#  o sec-04 sub-15:
#    - Get stats for LiloAndStich + ivar
#  o sec-04 sub-16:
#    - Get stats for LiloAndStich
#  o sec-04 sub-17:
#    - Get stats for ivar no trim
#  o sec-04 sub-18:
#    - Get stats for ivar with trim
#  o sec-04 sub-19:
#    - Clean up and move to next round
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-04 Sub-01:
#  - Check if using a seed or user wanted no seed
#**********************************************************

while [[ "$repI" -gt 0 ]]; do
# Loop: Run each replicate

   if [[ "$seedI" != "" ]]; then
      inputSeedI="$((seedI * repI))";
   else
      inputSeedI="";
   fi # If I need to update the seed

   #*******************************************************
   # Sec-04 Sub-02:
   #  - Subsample reads for each amplicon
   #*******************************************************

   for iDepth in ${depthAryI[*]}; do
   # Loop: though all read depths I am testing
      if [[ "$iDepth" -gt "$((minDepthI * 3))" ]]; then
         runDepthI="$minDepthI";
      else
         runDepthI="$((iDepth / 3))";
      fi # Check if I can meet the min read depth

      readsStr="$prefixStr-subSamp-depth$iDepth";
      readsStr="$readsStr-rep$repI.fastq";

      subSampDirStr="$prefixStr-tmp-subsamples";
      mkdir "$subSampDirStr";

      for strAmpFq in "$ampDirStr/"*.f*q; do
      # Loop: Subsample all amplicon fastq files
         if [[ ! -f "$strAmpFq" ]]; then
            continue;
         fi # If is the null case

         # remove fastq ending
         tmpStr="$(
            printf "%s" "$(basename "$strAmpFq")" |
              sed 's/\.fa*s*t*q//;' \
         )";
         tmpStr="$prefixStr-$tmpStr-rep$repI-Depth$iDepth";

         bash "$scriptDirStr/subsampleReads.sh" \
            -fastq "$strAmpFq" \
            -depth "$iDepth" \
            -seed "$inputSeedI" \
           > "$tmpStr.fastq";

         # Merge subsample into a fastq file
         cat "$tmpStr.fastq" >> "$readsStr";
         mv "$tmpStr.fastq" "$subSampDirStr";
      done # Loop: Subsample all amplicon fastq files

      numReadsI="$(
         sed -n 'p;n;n;n;' "$readsStr" | 
            wc -l |
            awk '{print $1}' \
      )"; # Get the number of reads in the fastq file
   
      if [[ ! -f "$readsStr" ]]; then
         continue; # No amplicons were made
      fi # if is the null case

      #****************************************************
      # Sec-04 Sub-03:
      #  - Find the reference and primers for each scheme
      #****************************************************

      for strVer in ${verAryStr[*]}; do
      # Loop: though all scheme versions
         # I need to find the primers and reference I am
         # working with
         primersStr="$schemeDirStr/$schemeStr/V$strVer/";
         primersStr="$(\
             find $primersStr -name *.scheme.bed \
         )";

         refStr="$schemeDirStr/$schemeStr/V$strVer/";
         refStr="$(find $refStr -name *.reference.fasta)";
         backUpRefStr="backup-$(basename "$refStr")";

         lenSeqI="$(\
            awk '
               { # MAIN
                  if(NR == 1){next;}
                  if($1 ~ /^>/){exit;}
                  seqStr = seqStr $0;
               } # MAIN
               END{
                  gsub(/ /, "", seqStr);
                  print length(seqStr);
               }' "$refStr"\
         )"; # find the reference sequence length

         #*************************************************
         # Sec-04 Sub-04:
         #  - Mutate the reference
         #*************************************************

         for iPercMut in ${percMutAryI[*]}; do
         # Loop: though mutating the references
            conPrefStr="$prefixStr-rep$repI-Depth$iDepth";
            conPrefStr="$conPrefStr-percMut$IpercMut";
            conPrefStr="$conPrefStr-$schemeStr-V$strVer";

            numMutateI="$(((iPercMut * lenSeqI) / 10000))";
              # allowing for 0.01% to be used
            mv "$refStr" "$backUpRefStr" || exit;
            trueRefStr="$backUpRefStr";

            if [[ "$iPercMut" -gt 0 ]]; then
            # If I am mutated the reference sequence
               "$scriptDirStr/mutateSeq.sh" \
                  -seq "$backUpRefStr" \
                  -num-mutate "$numMutateI" \
                  -primers "$doNotMutateStr" \
                  -seed "$inputSeedI" \
                  > "$refStr";
            
                mutantRefStatsStr="$lenSeqI  $(
                   bash "$scriptDirStr/getStats.sh" \
                      "$backUpRefStr" \
                      "$refStr" \
                )"; # Get stats for new reference
            else
               cp "$backUpRefStr" "$refStr";
               mutantRefStatsStr="$lenSeqI 0 0  0  0 0";
            fi # If I am mutated the reference sequence

            #**********************************************
            # Sec-04 Sub-05:
            #  - Get stats for artic
            #**********************************************

            testAlreadyRunBl="$( \
              awk \
                  -f "$scriptDirStr/checkIfRunTest.awk" \
                  -v progamStr="artic" \
                  -v depthI="$iDepth" \
                  -v percMutI="$iPercMut" \
                  -v seedI="$inputSeedI" \
                  -v modelStr="$medakaModelStr" \
                  -v usedMedakaBl="true" \
                  -v usedMajConBl="false" \
                  -v minDepthI="NA" \
                  -v minLenI="NA" \
                  -v fqDirStr="$ampDirStr" \
                  -v schemeDirStr="$schemeDirStr" \
                  -v schemeStr="$schemeStr" \
                  -v verStr="$strVer" \
                  -v prefixStr="$prefixStr" \
                  < "$statsFileStr" \
            )";

            if [[ "$articlBl" -gt 0 ]]; then
             if [[ "$testAlreadyRunBl" -eq 0 ]]; then
               /usr/bin/time \
                   -f "%e\t%S\t%U\t%M\t%P" \
                   -o "tmp.tsv" \
                  bash "$scriptDirStr/benchArticNoMut.sh" \
                    -fastq "$readsStr" \
                    -scheme-dir "$schemeDirStr" \
                    -scheme "$schemeStr" \
                    -scheme-version "$strVer" \
                    -model "$medakaModelStr" \
                    -t "$threadsI" \
                    -prefix "$prefixStr";

               if [[ ! -f "$prefixStr-scaffold.fa" ]]; then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-scaffold.fa" \
                  )";
                  rm "$prefixStr-scaffold.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                     ' # Check if this command failed
                       {if($0 ~ /Com/)
                        {
                           print "NA","NA","NA","NA","NA";
                           exit;
                        }
                        else print $0;
                       }
                     ' \
               )";
               rm "tmp.tsv";

               { # Print out the stats
                  printf \
                     "artic\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf "\ttrue\tfalse\tNA\tNA";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\ttrue\tfalse\tfalse\n";
               } >> "$statsFileStr" # Print out the stats
             fi # If I have not already run this test
		      fi # If running artic

            #**********************************************
            # Sec-04 Sub-06:
            #  - Get stats for buildCon
            #**********************************************

            testAlreadyRunBl="$( \
              awk \
                  -f "$scriptDirStr/checkIfRunTest.awk" \
                  -v progamStr="buildCon" \
                  -v depthI="$iDepth" \
                  -v percMutI="$iPercMut" \
                  -v seedI="$inputSeedI" \
                  -v modelStr="$medakaModelStr" \
                  -v usedMedakaBl="true" \
                  -v usedMajConBl="true" \
                  -v minLenI="$minLenI" \
                  -v minDepthI="$runDepthI" \
                  -v fqDirStr="$ampDirStr" \
                  -v schemeDirStr="$schemeDirStr" \
                  -v schemeStr="$schemeStr" \
                  -v verStr="$strVer" \
                  -v prefixStr="$prefixStr" \
                  < "$statsFileStr" \
            )";

            if [[ "$buildConMajconMedakaBl" -gt 0 ]]; then
             if [[ "$testAlreadyRunBl" -eq 0 ]]; then
               /usr/bin/time \
                  -f "%e\t%S\t%U\t%M\t%P" \
                  -o "tmp.tsv" \
                  bash "$scriptDirStr/buildAmpCons.sh" \
                    -fastq "$readsStr" \
                    -ref "$refStr" \
                    -primer-scheme "$primersStr" \
                    -min-length "$minLenI" \
                    -min-depth "$runDepthI" \
                    -model "$medakaModelStr" \
                    -t "$threadsI" \
                    -prefix "$prefixStr";

               if [[ ! -f "$prefixStr-scaffold.fa" ]]; then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-scaffold.fa" \
                  )";
                  rm "$prefixStr-scaffold.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        {if($0 ~ /Com/)
                         {
                            print "NA","NA","NA","NA","NA";
                            exit;
                         }
                         else print $0;
                        }
                      ' \
               )";

               rm "tmp.tsv";
               rm "$prefixStr-amps-cons.fa";

               { # Print out the stats
                  printf \
                     "buildCon\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf \
                    "\t%s\t%s\t%s\t%s" \
                    "true" \
                    "true" \
                    "$minLenI" \
                    "$runDepthI";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\tfalse\ttrue\n";
               } >> "$statsFileStr" # Print out the stats
             fi # If I have not run this test
            fi # If testing buildCon with medaka and majcon

            #**********************************************
            # Sec-04 Sub-07:
            #  - Get stats for buildCon no medaka
            #**********************************************

            testAlreadyRunBl="$( \
              awk \
                  -f "$scriptDirStr/checkIfRunTest.awk" \
                  -v progamStr="buildCon" \
                  -v depthI="$iDepth" \
                  -v percMutI="$iPercMut" \
                  -v seedI="$inputSeedI" \
                  -v modelStr="$medakaModelStr" \
                  -v usedMedakaBl="false" \
                  -v usedMajConBl="true" \
                  -v minLenI="$minLenI" \
                  -v minDepthI="$runDepthI" \
                  -v fqDirStr="$ampDirStr" \
                  -v schemeDirStr="$schemeDirStr" \
                  -v schemeStr="$schemeStr" \
                  -v verStr="$strVer" \
                  -v prefixStr="$prefixStr" \
                  < "$statsFileStr" \
            )";


            if [[ "$buildConMajconBl" -gt 0 ]]; then
             if [[ "$testAlreadyRunBl" -eq 0 ]]; then
               /usr/bin/time \
                  -f "%e\t%S\t%U\t%M\t%P" \
                  -o "tmp.tsv" \
                  bash "$scriptDirStr/buildAmpCons.sh" \
                    -fastq "$readsStr" \
                    -ref "$refStr" \
                    -primer-scheme "$primersStr" \
                    -min-length "$minLenI" \
                    -min-depth "$runDepthI" \
                    -disable-medaka \
                    -t "$threadsI" \
                    -prefix "$prefixStr";

               if [[ ! -f "$prefixStr-scaffold.fa" ]]; then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-scaffold.fa" \
                  )";
                  rm "$prefixStr-scaffold.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        {  if($0 ~ /Com/)
                          {
                            print "NA","NA","NA","NA","NA";
                            exit;
                          }
                          else print $0;
                        }
                     ' \
               )";

               rm "tmp.tsv";
               rm "$prefixStr-amps-cons.fa";

               { # Print out the stats
                  printf \
                     "buildCon\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf \
                    "\t%s\t%s\t%s\t%s" \
                    "false" \
                    "true" \
                    "$minLenI" \
                    "$runDepthI";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\tfalse\ttrue\n";
               } >> "$statsFileStr" # Print out the stats
             fi # If I have not run this test
            fi # If runing buildCon without medaka

            #**********************************************
            # Sec-04 Sub-08:
            #  - Stats: buildCon, no medaka, ivar polish
            #**********************************************

            testAlreadyRunBl="$( \
              awk \
                  -f "$scriptDirStr/checkIfRunTest.awk" \
                  -v progamStr="buildConIvar" \
                  -v depthI="$iDepth" \
                  -v percMutI="$iPercMut" \
                  -v seedI="$inputSeedI" \
                  -v modelStr="$medakaModelStr" \
                  -v usedMedakaBl="false" \
                  -v usedMajConBl="true" \
                  -v minLenI="$minLenI" \
                  -v minDepthI="$runDepthI" \
                  -v fqDirStr="$ampDirStr" \
                  -v schemeDirStr="$schemeDirStr" \
                  -v schemeStr="$schemeStr" \
                  -v verStr="$strVer" \
                  -v prefixStr="$prefixStr" \
                  < "$statsFileStr" \
            )";


            if [[ "$buildConMajconBl" -gt 0 ]]; then
             if [[ "$testAlreadyRunBl" -eq 0 ]]; then
               /usr/bin/time \
                  -f "%e\t%S\t%U\t%M\t%P" \
                  -o "tmp.tsv" \
                  bash "$scriptDirStr/buildAmpCons.sh" \
                    -fastq "$readsStr" \
                    -ref "$refStr" \
                    -primer-scheme "$primersStr" \
                    -min-length "$minLenI" \
                    -min-depth "$runDepthI" \
                    -disable-medaka \
                    -t "$threadsI" \
                    -ivar \
                    -prefix "$prefixStr";

               if [[ ! -f "$prefixStr-scaffold.fa" ]]; then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-scaffold.fa" \
                  )";
                  rm "$prefixStr-scaffold.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        {  if($0 ~ /Com/)
                          {
                            print "NA","NA","NA","NA","NA";
                            exit;
                          }
                          else print $0;
                        }
                     ' \
               )";

               rm "tmp.tsv";
               rm "$prefixStr-amps-cons.fa";

               { # Print out the stats
                  printf \
                     "buildConIvar\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf \
                    "\t%s\t%s\t%s\t%s" \
                    "false" \
                    "true" \
                    "$minLenI" \
                    "$runDepthI";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\tfalse\ttrue\n";
               } >> "$statsFileStr" # Print out the stats
             fi # If I have not run this test
            fi # If runing buildCon without medaka


            #**********************************************
            # Sec-04 Sub-09:
            #  - stats for buildCon just medaka (no majcon)
            #**********************************************

            testAlreadyRunBl="$( \
              awk \
                  -f "$scriptDirStr/checkIfRunTest.awk" \
                  -v progamStr="buildCon" \
                  -v depthI="$iDepth" \
                  -v percMutI="$iPercMut" \
                  -v seedI="$inputSeedI" \
                  -v modelStr="$medakaModelStr" \
                  -v usedMedakaBl="true" \
                  -v usedMajConBl="false" \
                  -v minLenI="$minLenI" \
                  -v minDepthI="$runDepthI" \
                  -v fqDirStr="$ampDirStr" \
                  -v schemeDirStr="$schemeDirStr" \
                  -v schemeStr="$schemeStr" \
                  -v verStr="$strVer" \
                  -v prefixStr="$prefixStr" \
                  < "$statsFileStr" \
            )";

            if [[ "$buildConMedakaBl" -gt 0 ]]; then
             if [[ "$testAlreadyRunBl" -eq 0 ]]; then
               /usr/bin/time \
                  -f "%e\t%S\t%U\t%M\t%P" \
                  -o "tmp.tsv" \
                  bash "$scriptDirStr/buildAmpCons.sh" \
                    -fastq "$readsStr" \
                    -ref "$refStr" \
                    -primer-scheme "$primersStr" \
                    -min-length "$minLenI" \
                    -min-depth "$runDepthI" \
                    -model "$medakaModelStr" \
                    -disable-majcon \
                    -t "$threadsI" \
                    -prefix "$prefixStr";

               if [[ ! -f "$prefixStr-scaffold.fa" ]]; then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-scaffold.fa" \
                  )";
                  rm "$prefixStr-scaffold.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        {  if($0 ~ /Com/)
                          {
                            print "NA","NA","NA","NA","NA";
                            exit;
                          }
                          else print $0;
                          }
                        ' \
               )";

               rm "tmp.tsv";
               rm "$prefixStr-amps-cons.fa";

               { # Print out the stats
                  printf \
                     "buildCon\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf \
                    "\t%s\t%s\t%s\t%s" \
                    "true" \
                    "false" \
                    "$minLenI" \
                    "$runDepthI";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\tfalse\ttrue\n";
               } >> "$statsFileStr" # Print out the stats
             fi # If I have not run this test
            fi # If runing buildCon with only medaka

            #**********************************************
            # Sec-04 Sub-10:
            #  - Run LILO
            #**********************************************

            testAlreadyRunBl="$( \
              awk \
                  -f "$scriptDirStr/checkIfRunTest.awk" \
                  -v progamStr="Lilo" \
                  -v depthI="$iDepth" \
                  -v percMutI="$iPercMut" \
                  -v seedI="$inputSeedI" \
                  -v modelStr="$medakaModelStr" \
                  -v usedMedakaBl="true" \
                  -v usedMajConBl="false" \
                  -v minDepthI="NA" \
                  -v minLenI="NA" \
                  -v fqDirStr="$ampDirStr" \
                  -v schemeDirStr="$schemeDirStr" \
                  -v schemeStr="$schemeStr" \
                  -v verStr="$strVer" \
                  -v prefixStr="$prefixStr" \
                  < "$statsFileStr" \
            )";

            if [[ "$liloBl" -gt 0 ]]; then
             if [[ "$testAlreadyRunBl" -eq 0 ]]; then
               /usr/bin/time \
                  -f "%e\t%S\t%U\t%M\t%P" \
                  -o "tmp.tsv" \
                  bash "$scriptDirStr/benchLILO.sh" \
                    -model "$medakaModelStr" \
                    -t "$threadsI" \
                    -fastq "$(\
                        printf "%s" "$readsStr" |
                        sed 's/\/\//\//g; s/[ \t]//g;' \
                    )" \
                    -ref "$(\
                        printf "%s" "$refStr" |
                        sed 's/\/\//\//g; s/[ \t]//g;' \
                    )" \
                    -primer-scheme "$(\
                        printf "%s" "$primersStr" |
                        sed 's/\/\//\//g; s/[ \t]//g;' \
                    )" \
                    -prefix "$(\
                        printf "%s" "$prefixStr" |
                        sed 's/\/\//\//g; s/[ \t]//g;' \
                    )";

               #*******************************************
               # Sec-04 Sub-11:
               #  - Get stats for fixed LILO (auto) + ivar
               #*******************************************

               bash "$scriptDirStr/runIvar.sh" \
                   -ref  "$prefixStr-scaffold-auto.fa" \
                   -fastq "$readsStr" \
                   -prefix "$prefixStr-lilo-auto-ivar" \
                   -threads "$threadsI";

               if [[ ! -f "$prefixStr-lilo-auto-ivar-con.fa" ]];
               then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                       "$trueRefStr" \
                       "$prefixStr-lilo-auto-ivar-con.fa" \
                  )";
                  rm "$prefixStr-lilo-auto-ivar-con.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        {if($0 ~ /Com/)
                         {
                            print "NA","NA","NA","NA","NA";
                            exit;
                         }
                         else print $0;
                        }
                      ' \
               )";

               rm "$prefixStr-LILO-amps.fa";

               { # Print out the stats
                  printf \
                     "LiloIvar\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf "\ttrue\tfalse\tNA\tNA";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\ttrue\tfalse\tfalse\n";
               } >> "$statsFileStr" # Print out the stats

               #*******************************************
               # Sec-04 Sub-12:
               #  - Get stats for fixed LILO (auto)
               #*******************************************

               if [[ ! -f "$prefixStr-scaffold-auto.fa" ]];
               then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-scaffold-auto.fa" \
                  )";
                  rm "$prefixStr-scaffold-auto.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        {if($0 ~ /Com/)
                         {
                            print "NA","NA","NA","NA","NA";
                            exit;
                         }
                         else print $0;
                        }
                      ' \
               )";

               #rm "$prefixStr-LILO-amps.fa";

               { # Print out the stats
                  printf \
                     "Lilo\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf "\ttrue\tfalse\tNA\tNA";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\ttrue\tfalse\tfalse\n";
               } >> "$statsFileStr" # Print out the stats

               #*******************************************
               # Sec-04 Sub-13:
               #  - separate scaffold_builer LILO + ivar
               #*******************************************

               bash "$scriptDirStr/runIvar.sh" \
                   -ref  "$prefixStr-scaffold-manual.fa" \
                   -fastq "$readsStr" \
                   -prefix "$prefixStr-lilo-man-ivar" \
                   -threads "$threadsI";

               if [[ ! -f "$prefixStr-lilo-man-ivar-con.fa" ]];
               then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                    bash "$scriptDirStr/getStats.sh" \
                      "$trueRefStr" \
                      "$prefixStr-lilo-man-ivar-con.fa"\
                  )";
                  rm "$prefixStr-lilo-man-ivar-con.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        { if($0 ~ /Com/)
                          {
                            print "NA","NA","NA","NA","NA";
                            exit;
                          }
                          else print $0;
                        }
                      ' \
               )";

               { # Print out the stats
                  printf \
                     "LiloIvar\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf "\ttrue\tfalse\tNA\tNA";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\ttrue\tfalse\n";
               } >> "$statsFileStr" # Print out the stats


               #*******************************************
               # Sec-04 Sub-14:
               #  - stats separate scaffold_builer LILO
               #*******************************************

             if [[ ! -f "$prefixStr-scaffold-manual.fa" ]];
               then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-scaffold-manual.fa" \
                  )";
                  rm "$prefixStr-scaffold-manual.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        { if($0 ~ /Com/)
                          {
                            print "NA","NA","NA","NA","NA";
                            exit;
                          }
                          else print $0;
                        }
                      ' \
               )";

               { # Print out the stats
                  printf \
                     "Lilo\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf "\ttrue\tfalse\tNA\tNA";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\ttrue\tfalse\n";
               } >> "$statsFileStr" # Print out the stats

               #*******************************************
               # Sec-04 Sub-15:
               #  - Get stats for LiloAndStich + ivar
               #*******************************************

               bash "$scriptDirStr/runIvar.sh" \
                   -ref  "$prefixStr-scaffold-stich.fa" \
                   -fastq "$readsStr" \
                   -prefix "$prefixStr-lilo-stich-ivar" \
                   -threads "$threadsI";

               if [[ ! -f "$prefixStr-lilo-stich-ivar-con.fa" ]];
               then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                    bash "$scriptDirStr/getStats.sh" \
                      "$trueRefStr" \
                      "$prefixStr-lilo-stich-ivar-con.fa" \
                  )";
                  rm "$prefixStr-lilo-stich-ivar-con.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                        ' # Check if this command failed
                          {  if($0 ~ /Com/)
                             {
                                print "NA","NA","NA","NA","NA";
                                exit;
                             }
                             else print $0;
                          }
                        ' \
               )";

               { # Print out the stats
                  printf \
                   "LiloIvar\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf "\ttrue\tfalse\tNA\tNA";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\tfalse\ttrue\n";
               } >> "$statsFileStr" # Print out the stats

               #*******************************************
               # Sec-04 Sub-16:
               #  - Get stats for LiloAndStich
               #*******************************************

              if [[ ! -f "$prefixStr-scaffold-stich.fa" ]];
               then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                   statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-scaffold-stich.fa" \
                  )";
                  rm "$prefixStr-scaffold-stich.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                        ' # Check if this command failed
                          {  if($0 ~ /Com/)
                             {
                                print "NA","NA","NA","NA","NA";
                                exit;
                             }
                             else print $0;
                          }
                        ' \
               )";

               { # Print out the stats
                  printf \
                   "Lilo\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "$medakaModelStr";
                 printf "\ttrue\tfalse\tNA\tNA";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\tfalse\ttrue\n";
               } >> "$statsFileStr" # Print out the stats

               rm "tmp.tsv";

             fi # If I have not run this test
            fi # If using LILO

            #**********************************************
            # Sec-04 Sub-17:
            #  - Get stats for ivar no trim
            #**********************************************

            testAlreadyRunBl="$( \
              awk \
                  -f "$scriptDirStr/checkIfRunTest.awk" \
                  -v progamStr="ivar" \
                  -v depthI="$iDepth" \
                  -v percMutI="$iPercMut" \
                  -v seedI="$inputSeedI" \
                  -v modelStr="NA" \
                  -v usedMedakaBl="false" \
                  -v usedMajConBl="false" \
                  -v minLenI="NA" \
                  -v minDepthI="NA" \
                  -v fqDirStr="$ampDirStr" \
                  -v schemeDirStr="$schemeDirStr" \
                  -v schemeStr="$schemeStr" \
                  -v verStr="$strVer" \
                  -v prefixStr="$prefixStr" \
                  < "$statsFileStr" \
            )";

            if [[ "$ivarBl" -gt 0 ]]; then
             if [[ "$testAlreadyRunBl" -eq 0 ]]; then
               /usr/bin/time \
                  -f "%e\t%S\t%U\t%M\t%P" \
                  -o "tmp.tsv" \
                  bash "$scriptDirStr/runIvar.sh" \
                    -fastq "$readsStr" \
                    -ref "$refStr" \
                    -t "$threadsI" \
                    -prefix "$prefixStr";

               if [[ ! -f "$prefixStr-con.fa" ]]; then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-con.fa" \
                  )";
                  rm "$prefixStr-con.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        {  if($0 ~ /Com/)
                          {
                            print "NA","NA","NA","NA","NA";
                            exit;
                          }
                          else print $0;
                        }
                     ' \
               )";

               rm "tmp.tsv";

               { # Print out the stats
                  printf \
                     "ivar\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "NA";
                 printf "\tfalse\tfalse\tNA\tNA";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\tfalse\tfalse\n";
               } >> "$statsFileStr" # Print out the stats
             fi # If I have not run this test
            fi # If runing buildCon without medaka

            #**********************************************
            # Sec-04 Sub-18:
            #  - Get stats for ivar with trim
            #**********************************************

            testAlreadyRunBl="$( \
              awk \
                  -f "$scriptDirStr/checkIfRunTest.awk" \
                  -v progamStr="ivarTrim" \
                  -v depthI="$iDepth" \
                  -v percMutI="$iPercMut" \
                  -v seedI="$inputSeedI" \
                  -v modelStr="NA" \
                  -v usedMedakaBl="false" \
                  -v usedMajConBl="false" \
                  -v minLenI="NA" \
                  -v minDepthI="NA" \
                  -v fqDirStr="$ampDirStr" \
                  -v schemeDirStr="$schemeDirStr" \
                  -v schemeStr="$schemeStr" \
                  -v verStr="$strVer" \
                  -v prefixStr="$prefixStr" \
                  < "$statsFileStr" \
            )";

            if [[ "$ivarTrimBl" -gt 0 ]]; then
             if [[ "$testAlreadyRunBl" -eq 0 ]]; then
               /usr/bin/time \
                  -f "%e\t%S\t%U\t%M\t%P" \
                  -o "tmp.tsv" \
                  bash "$scriptDirStr/runIvar.sh" \
                    -fastq "$readsStr" \
                    -ref "$refStr" \
                    -primers "$primersStr" \
                    -t "$threadsI" \
                    -primers "$primersStr" \
                    -prefix "$prefixStr";

               if [[ ! -f "$prefixStr-con-trim.fa" ]]; then
                  statsLineStr="ERR	ERR	ERR	ERR	ERR";
               else
                  statsLineStr="$( \
                     bash "$scriptDirStr/getStats.sh" \
                         "$trueRefStr" \
                         "$prefixStr-con-trim.fa" \
                  )";
                  rm "$prefixStr-con-trim.fa";
               fi

               statsLineStr="$statsLineStr   $(\
                  cat tmp.tsv |
                    awk \
                      ' # Check if this command failed
                        {  if($0 ~ /Com/)
                          {
                            print "NA","NA","NA","NA","NA";
                            exit;
                          }
                          else print $0;
                        }
                     ' \
               )";

               rm "tmp.tsv";

               { # Print out the stats
                  printf \
                     "ivarTrim\t%s\t%s\t%s\t%s\t%s\t%s" \
                     "$repI" \
                     "$iDepth" \
                     "$numReadsI" \
                     "$iPercMut" \
                     "$numMutateI" \
                     "$mutantRefStatsStr";
                  printf \
                     "\t%s\t%s\t%s\t%s" \
                     "$statsLineStr" \
                     "$seedI" \
                     "$inputSeedI" \
                     "NA";
                 printf "\tfalse\tfalse\tNA\tNA";
                 printf \
                     "\t%s\t%s\t%s\t%s\t%s" \
                     "$ampDirStr" \
                     "$schemeDirStr" \
                     "$schemeStr" \
                     "$strVer" \
                     "$prefixStr";
                 # consensus building method
                 printf "\tfalse\tfalse\tfalse\n";
               } >> "$statsFileStr" # Print out the stats
             fi # If I have not run this test
            fi # If runing buildCon without medaka

            #**********************************************
            # Sec-04 Sub-19:
            #  - Clean up and move to next round
            #**********************************************

            mv "$backUpRefStr" "$refStr";
         done # Loop: though mutating the references
      done # Loop: though all scheme versions

      rm -r "$readsStr" "$subSampDirStr"; # No longer need
   done # Loop: though all read depths I am testing


   repI="$((repI - 1))";
done # Loop: Run each replicate

exit;
