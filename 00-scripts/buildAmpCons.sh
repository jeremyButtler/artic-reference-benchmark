###########################################################
# Name: bulidAmpCons.sh
# Use:
#  - Builds a consensus reads produced from tiling primers.
# Input:
#  - See help message
# Output:
#  - prefix-amps-cons.fa
#    o Consensuses for individual amplicons (fasta file)
#  - prefix-scaffold.fa
#    o Stiched together amplicons (scaffold)
###########################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#  o sec-01 sub-01:
#    - Variables holding General user input
#  o sec-01 sub-01:
#    - Variables holding General user input
#  o sec-01 sub-02:
#    - Variables for consensus building or polishing steps
#  o sec-01 sub-03:
#    - Variables specific to ivar
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-01 Sub-01:
#  - Variables holding General user input
#**********************************************************

refStr="";
fqStr="";
readExtractTblStr="";
primersStr="";
inPrefStr="buildcon"; # Input prefix to name stuff
minLenI=300;  # Min read length/consensus length
minDepthI=100; # Min bin depth to build consensus
minPerDepthI="0.1"; # Min percentage of mapped reads
threadsI=3;

#**********************************************************
# Sec-01 Sub-02:
#  - Variables for consensus building or polishing steps
#**********************************************************

medakaModelStr="r941_prom_high_g344"; # Model for medaka
conMethodAryStr=("-use-majCon" "-use-medaka");
   # Holds consensus methods using
lenConMethodI=2; # length of conMethodAryStr

useRefBl=0;   # Use reference in consensus step
ivarBl=0;      # Polish consensus with ivar
polishMedBl=0; # Polish consensus with medakd
minQI=10; # Both majcon and ivar use

#**********************************************************
# Sec-01 Sub-03:
#  - Variables specific to ivar
#**********************************************************

ivarMinSupDbl=0.5;
ivarMinInsSupDbl=0.8;
ivarMinDepthI=10;
minSubDepthI=300;

#**********************************************************
# Sec-01 Sub-04:
#  - Non-user input variables
#**********************************************************


scriptDirStr="$(dirname "$0" | sed 's/^\.\///;')";
numReadsI=0;
prefStr="";         # Temporary prefix
conStr="";
extraOptions=""; # Extra input options
fqIterI=0;       # For checking if I have fastq files
errBl=0;         # To know if I need to quite
tmpFqStr="";
firstConMethodBl=1; #Tells if user modifed consensus method

# path to use to activate conda
condaPathStr="$(\
  conda info | 
    grep -i 'base environment' | 
    sed 's/base.*: //; s/  *.read only.//; s/ //g' \
)";

helpStr="$(basename "$0") \
   -ref ref.fasta \
   -fastq reads.fastq \
   -table read-extract-table.md \
   -primer-scheme primer-scheme.bed;
Use:
   - Builds a consensus with buildcon
   - This trims primers off reads with trimPrimers and 
     then nosiy ends with the reference.
Input:
  General:
    -ref: [Required]
      o Fasta file with the reference sequence
    -fastq: [Required]
      o Fastq file with reads to build amplicons with
    -table: [Required]
      o Table of reads to extract
      o use readLenPosTbl.sh to build this table.
      o Then edit this table to only have the reads you
        want to use to build amplicon consensus with.
    -primer-scheme: [Required]
      o primer scheme used with artic.
      o This will be a bed file with 7 columns
        - Column 4 has the primer name
        - Column 7 has the primer sequence
    -min-length: [$minLenI]
      o Minimum length to map reads and build consensuses
        with.
    -min-depth: [$minDepthI]
      o Minimum depth to build a consensus at.
    -min-perc-depth: [$minPerDepthI]
      o Minimum percentage of total reads needed to build
        a consensus for a read or consensus (0 to 1).
    -prefix: [$inPrefStr]
      o Prefix to name the output amplicon stats file
      o prefix-amp-stats.tsv
    -model: [$medakaModelStr]
      o Model to use with medaka
    -polish-ivar: [No]
      o Use ivar to polish the consensus
      o -polish-ivar is disabled with -no-polish-ivar
    -polish-medaka: [No]
      o Use medaka to polish the consensus
      o -polish-medaka is disabled with -no-polish-medaka
    -use-ref: [No]
      o Use the reference to build the amplicon consensus
      o -use-ref is disabled with -no-ref
    -t: [$threadsI]
      o Number of threads to use
  Method choice:
    - Setting one or more of these options disables the
      default options [use-majcon -use-medaka].
    - Each option can be used multiple times and is run in
      the order submitted. So, order matters.
      o \"-use-majcon -use-medaka\" will run majcon first
        and then medaka.
      o \"-use-medaka -use-majcon\" will run medaka first
        and then majcon.
    -use-medaka: [-use-majcon -use-medaka]
      o Use medaka (Medaka is very slow)
    -use-majcon: [-use-majcon -use-medaka]
      o Use the majority consensus step
    -use-ivar: [-use-majcon -use-medaka]
      o Use ivar to build amplicon consensuses
    -use-racon: [-use-majcon -use-medaka]
      o Use racon to build an amplicon consensus
  Consensus building:
    -min-base-q: [$minQI]
      o Min q-score needed to keep a base (for all methods)
    -ivar-min-sup: [$ivarMinSupDbl]
      o Min support to keep an snp with ivar (0 to 1)
    -ivar-min-ins-sup: [$ivarMinInsSupDbl]
      o Min support to keep insertion with ivar (0 to 1)
    -ivar-min-depth: [$ivarMinDepthI]
      o Min read depth for ivar to not mask a base
    -min-sub-depth: [$minSubDepthI]
      o Min read depth to use subsampled reads from
        buildcon for polishing. Otherwise all reads are
        used.
Output:
   - Prints the consensus to an file
     deletions to stdout. Adds file name to end if
     \"-disable-simp-stats\" is used.
   - Prints the matchs, snps, insertions, deletions, and
     file name for each amplicon to prefix-amp-stats.tsv
";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Get and check user input
#  o sec-02 sub-01:
#    - Get user input
#  o sec-02 sub-02:
#    - Check user input
#  o sec-02 sub-03:
#    - Update variables form user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-02 Sub-01:
#  - Get user input
#**********************************************************

while [[ $# -gt 0 ]]; do
# Loop: Get user input
   case $1 in
       # General input
      -ref) refStr="$2"; shift;;
      -fastq) fqStr="$2"; shift;;
      -table) readExtractTblStr="$2"; shift;;
      -primer-scheme) primersStr="$2"; shift;;
      -min-length) minLenI="$2"; shift;;
      -min-depth) minDepthI="$2"; shift;;
      -min-perc-depth) minPerDepthI="$2"; shift;;
      -prefix) inPrefStr="$2"; shift;;
      -model) medakaModelStr="$2"; shift;;
      -t) threadsI="$2"; shift;;
      -threads) threadsI="$2"; shift;;

      # Polishing steps 
      -no-polish-medaka) polishMedBl=0;;
      -polish-medaka) polishMedBl=1;;
      -polish-ivar) ivarBl=1;;
      -no-polish-ivar) ivarBl=0;;

      -use-ref) useRefBl=1;;
      -no-ref) useRefBl=0;;

      # Consensus building method
      -use-medaka)
         if [[ "$firstConMethodBl" -gt 0 ]]; then
            conMethodAryStr=("-use-medaka");
            firstConMethodBl=0;
            lenConMethodI=1;
         else
            conMethodAryStr+=("-use-medaka");
            lenConMethodI="$((lenConMethodI + 1))";
         fi;; # Check if this is the users first input
      -use-majcon)
         if [[ "$firstConMethodBl" -gt 0 ]]; then
            conMethodAryStr=("-use-majCon");
            firstConMethodBl=0;
            lenConMethodI=1;
         else
            conMethodAryStr+=("-use-majCon");
            lenConMethodI="$((lenConMethodI + 1))";
         fi;; # Check if this is the users first input
      -use-ivar)
         if [[ "$firstConMethodBl" -gt 0 ]]; then
            conMethodAryStr=("-use-ivar");
            firstConMethodBl=0;
            lenConMethodI=1;
         else
            conMethodAryStr+=("-use-ivar");
            lenConMethodI="$((lenConMethodI + 1))";
         fi;; # Check if this is the users first input
      -racon)
         if [[ "$firstConMethodBl" -gt 0 ]]; then
            conMethodAryStr=("-use-racon");
            firstConMethodBl=0;
            lenConMethodI=1;
         else
            conMethodAryStr+=("-use-racon");
            lenConMethodI="$((lenConMethodI + 1))";
         fi;; # Check if this is the users first input

      # consensus settings
      -min-base-q) minQI="$2"; shift;;
      -ivar-min-sup) ivarMinSupDbl="$2"; shift;;
      -ivar-min-ins-sup) ivarMinInsSupDbl="$2"; shift;;
      -ivar-min-depth) ivarMinDepthI="$2"; shift;;
      -min-sub-depth) minSubDepthI="$2"; shift;;

      # help message and errors
      -h) printf "%s\n" "$helpStr"; exit;;
      --h) printf "%s\n" "$helpStr"; exit;;
      -help) printf "%s\n" "$helpStr"; exit;;
      --help) printf "%s\n" "$helpStr"; exit;;
      help) printf "%s\n" "$helpStr"; exit;;
      *)
          printf "%s\n%s is invalid\n" "$helpStr" "$1";
          exit;;
   esac

   shift; # Move to the next parameter
done # Loop: Get user input

#**********************************************************
# Sec-02 Sub-02:
#  - check user input
#  o sec-02 sub-02 cat-01:
#    - Check if the reference is is file
#  o sec-02 sub-02 cat-02:
#    - Check if fastq diretory is a directory with fastqs
#  o sec-02 sub-02 cat-03:
#    - Check if have a valid table to extract reads with
#  o sec-02 sub-02 cat-04:
#    - Check if the primer scheme is a file
#  o sec-02 sub-02 cat-05:
#    - Exit if there are errors
#**********************************************************

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-01:
#  - Check if the reference is is file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! -f "$refStr" ]]; then
   printf \
      " -ref %s is not a file\n" \
      "$refStr" \
      >&2;
   errBl=1;
fi # check if the reference is valid

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-02:
#  - Check if fastq diretory is a directory with fastqs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! -f "$fqStr" ]]; then
# check if the fastq directory has fastqs
   printf \
      " -fastq %s is not a file\n" \
      "$fqStr" \
      >&2;
   errBl=1;
fi # if the fastq directory is a directory

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-03:
#  - Check if have a valid table to extract reads with
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! -f "$readExtractTblStr" ]]; then
# If there is no true reference
   printf \
      " -table %s is not an file\n" \
      "$readExtractTblStr" \
      >&2;
   errBl=1;
fi # If there is no true reference

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-04:
#  - Check if the primer scheme is a file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! -f "$primersStr" ]]; then
   printf \
      " -primer-scheme %s is not an file\n" \
      "$primersStr" \
      >&2;
   errBl=1;
fi # Check if the primers are a file

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-05:
#  - Exit if there are errors
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$errBl" -gt 0 ]]; then exit; fi

#**********************************************************
# Sec-02 Sub-03:
#  - Update variables form user input
#**********************************************************

extraOptions="-threads $threadsI";

if [[ "$medakaModelStr" != "" ]]; then
   extraOptions="$extraOptions -model $medakaModelStr";
fi # If changing the model of medaka I am using
   # This will not tigure medaka to be called, so is safe
   # to always do

# Add in the consensus building methods
for strMethod in ${conMethodAryStr[*]}; do
# Loop: Add all consensus request to buildCon
   extraOptions="$extraOptions $strMethod";
done # Loop: Add all consensus request to buildCon

if [[ "$useRefBl" -gt 0 ]]; then
   extraOptions="$extraOptions -ref $refStr";
fi # Check if using the referenc to build the consensus

ampStatsFileStr="$inPrefStr-amp-stats.tsv";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Build consensuses and get errors
#  o sec-03 sub-01:
#    - split fastq file into amplicons
#  o sec-03 sub-02:
#    - Trim the reads
#  o sec-03 sub-03:
#    - Build consensus
#  o sec-03 sub-04:
#    - Grab a subsample for ivar
#  o sec-03 sub-05:
#    - Build scaffoled (stich amplicons together)
#  o sec-03 sub-06:
#    - run ivar to polish the consensus (if requested)
#  o sec-03 sub-07:
#    - Clean up and exit
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-03 Sub-01:
#  - split fastq file into amplicons
#**********************************************************

awk \
    '{printf ">%s\n%s\n", $4, $7}' \
    < "$primersStr" \
  > "$inPrefStr-tmp-primers.fasta";

mkdir "$inPrefStr-tmp-reads";
mkdir "$inPrefStr-amps-cons";
tmpFqStr="$inPrefStr-tmp-reads/";
tmpFqStr="$tmpFqStr/$(basename "$inPrefStr")-amps";

minimap2 --eqx -a -t "$threadsI" "$refStr" "$fqStr" |
  awk \
      -f "$scriptDirStr/extReadsByTbl.awk" \
      -v extTbl="$readExtractTblStr" \
      -v prefix="$tmpFqStr";

#**********************************************************
# Sec-03 Sub-02:
#  - Trim the reads
#**********************************************************

for strFq in ./"$tmpFqStr"*.fastq;
do # Loop: Build a consensus for each amplicon
   if [[ ! -f "$strFq" ]]; then continue; fi # null case

   prefStr="$inPrefStr-$(\
        printf "%s" "$(basename "$strFq")" |
        sed 's/\.\///; s/.*pos/pos/; s/\.fastq//;' \
   )";

   "$scriptDirStr/../00-programs/trimPrimers" \
       -end-trim \
       -primers "$inPrefStr-tmp-primers.fasta" \
       -fastq "$strFq" \
       -out "$prefStr-tmp.fastq";

   minimap2 \
       --eqx \
       -a \
       "$refStr" \
       "$prefStr-tmp.fastq" |
    "$scriptDirStr/../00-programs/trimSamFile" -stdin |
    samtools view -O fastq -o "$prefStr-tmp-2.fastq";

   rm "$prefStr-tmp.fastq";
   mv "$prefStr-tmp-2.fastq" "$prefStr-tmp.fastq";

   #*******************************************************
   # Sec-03 Sub-03:
   #  - Build consensus
   #*******************************************************

   "$scriptDirStr/../00-programs/buildcon" \
       -min-length "$minLenI" \
       -min-depth "$minDepthI" \
       -min-perc "$minPerDepthI" \
       -prefix "$prefStr" \
       -fastq "$prefStr-tmp.fastq" \
       -ivar-min-depth $ivarMinDepthI \
       -ivar-min-snp-sup $ivarMinSupDbl \
       -ivar-min-ins-sup $ivarMinInsSupDbl \
       -ivar-min-q $minQI \
       -maj-con-min-base-q $minQI \
       -read-con-min-base-q $minQI \
       -read-read-min-base-q $minQI \
       $extraOptions;

   #*******************************************************
   # Sec-03 Sub-04:
   #  - Grab a subsample for ivar
   #*******************************************************

   if [[ "$ivarBl" -gt 0 || "$polishMedBl" -gt 0 ]]; then
      numReadsI="$(\
         wc -l "$prefStr-tmp.fastq" | awk '{print $1/4}' \
      )"; # Find the number of reads in the file

      if [[ "$numReadsI" -gt "$minSubDepthI" ]]; then
         cat \
             "$prefStr--top-reads.fastq" \
           >> "$inPrefStr-ivar.fq";
      else
         cat "$prefStr-tmp.fastq" >>"$inPrefStr-ivar.fq";
      fi # Check if I am keeping a subsample or all reads
   fi # Check if I am keeping reads for ivar 

   rm "$prefStr-tmp.fastq"
   rm "$prefStr--top-reads.fastq";

   mv \
       "$prefStr--clust-0--con.fasta" \
       "$inPrefStr-amps-cons/$(basename "$prefStr")-con.fasta";
done # Loop: Build a consensus for each amplicon

#**********************************************************
# Sec-03 Sub-05:
#  - Build scaffoled (stich amplicons together)
#**********************************************************

cat \
    "$inPrefStr-amps-cons/"*.fasta \
  > "$inPrefStr-amps-cons.fa";

"$scriptDirStr/../00-programs/stich" \
    -ref "$refStr" \
    -amps "$inPrefStr-amps-cons.fa" \
    -out "$inPrefStr-scaffold.fa" \
    -overwrite \
    -prefix "$inPrefStr" \
    -threads "$threadsI";

#**********************************************************
# Sec-03 Sub-06:
#  - run ivar to polish the consensus (if requested)
#**********************************************************

if [[ "$polishMedBl" -gt 0 ]]; then
   . ~/medaka/venv/bin/activate 
   medaka --version || source "$condaPathStr/etc/profile.d/conda.sh";
   medaka --version || conda activate medaka || (printf "No medaka\n" >&2 && exit);

   medaka_consensus \
       -i "$inPrefStr-ivar.fq" \
       -d "$inPrefStr-scaffold.fa" \
       -m "$medakaModelStr" \
       -o "$inPrefStr-medaka-polish-tmp-dir" \
       -t "$threadsI";

   mv \
      "$inPrefStr-medaka-polish-tmp-dir/consensus.fasta" \
      "$inPrefStr-scaffold.fa";
   rm -r "$inPrefStr-medaka-polish-tmp-dir";

   deactivate || conda deactivate;
fi # Check if polishing with medaka

if [[ "$ivarBl" -gt 0 ]]; then
# If using ivar
  minimap2 \
      -t "$threadsI" \
      --eqx \
      -a \
      -x map-ont \
      "$inPrefStr-scaffold.fa" \
      "$inPrefStr-ivar.fq" |
    "$scriptDirStr/../00-programs/trimSamFile" -stdin |
    samtools sort -@ "$threadsI" - |
    samtools view -@ "$threadsI" -F 4 -F 256 -F 2048 -b - |
    samtools mpileup -B -aa -A -d 0 -Q 0 - | 
    "$scriptDirStr/../00-programs/ivar" consensus \
      -p "$inPrefStr-tmp-scaffold" \
      -i "$inPrefStr-scaffold" \
      -c "$ivarMinInsSupDbl" \
      -t "$ivarMinSupDbl" \
      -m "$minDepthI" \
      -n "N" \
      -q "$minQI";

  mv "$inPrefStr-tmp-scaffold.fa" "$inPrefStr-scaffold.fa";
  rm "$inPrefStr-tmp-scaffold.qual.txt";
  rm "$inPrefStr-ivar.fq"; # No longer need
fi # If using ivar

#**********************************************************
# Sec-03 Sub-07:
#  - Clean up and exit
#**********************************************************

rm "$inPrefStr-tmp-primers.fasta";
rm -r "$inPrefStr-tmp-reads";
rm -r "$inPrefStr-amps-cons/";

exit;
