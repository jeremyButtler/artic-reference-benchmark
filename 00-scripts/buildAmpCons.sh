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

refStr="01-input/01-primer-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.reference.fasta";
fqStr="05-amplicon-reads/";
readExtractTblStr="04-read-table/04-buildConTbl.md";
primersStr="01-input/01-primer-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.scheme.bed";
inPrefStr="buildCon"; # Input prefix to name stuff
minLenI=300;  # Min read length/consensus length
minDepthI=20; # Min bin depth to build consensus
threadsI=3;
medakaModelStr="r941_prom_high_g344"; # Model for medaka
useMedakaBl=1;   # 0 is do not use medaka

scriptDirStr="$(dirname "$0" | sed 's/^\.\///;')";
prefStr="";         # Temporary prefix
conStr="";
extraOptions=""; # Extra input options
fqIterI=0;       # For checking if I have fastq files
errBl=0;         # To know if I need to quite
tmpFqStr="";

helpStr="$(basename "$0") \
   -ref ref.fasta \
   -fastq reads.fastq \
   -table read-extract-table.md \
   -primer-scheme primer-scheme.bed;
Use:
   - Builds a consensus with buildCon
   - This trims primers off reads with trimPrimers and 
     then nosiy ends with the reference.
Input:
   -ref: [Required]
     o Fasta file with the reference sequence
   -fastq: [Required]
     o Fastq file with reads to build amplicons with
   -table: [Required]
     o Table of reads to extract
     o use readLenPosTbl.sh to build this table.
     o Then edit this table to only have the reads you want to use
       to build amplicon consensus with.
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
   -prefix: [$inPrefStr]
     o Prefix to name the output amplicon stats file
     o prefix-amp-stats.tsv
   -model: [$medakaModelStr]
     o Model to use with medaka
   -use-medaka: [Yes]
     o Use medaka with buildCon
     o This is disabled by -disable-medaka
   -t: [$threadsI]
     o Number of threads to use
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
      -ref) refStr="$2"; shift;;
      -fastq) fqStr="$2"; shift;;
      -table) readExtractTblStr="$2"; shift;;
      -primer-scheme) primersStr="$2"; shift;;
      -min-length) minLenI="$2"; shift;;
      -min-depth) minDepthI="$2"; shift;;
      -prefix) inPrefStr="$2"; shift;;
      -model) medakaModelStr="$2"; useMedakaBl=1; shift;;
      -use-medaka) useMedakaBl=1;;
      -disable-medaka) useMedakaBl=0;;
      -t) threadsI="$2"; shift;;
      -threads) threadsI="$2"; shift;;
      -h) printf "%s\n" "$helpStr"; exit;;
      --h) printf "%s\n" "$helpStr"; exit;;
      -help) printf "%s\n" "$helpStr"; exit;;
      --help) printf "%s\n" "$helpStr"; exit;;
      help) printf "%s\n" "$helpStr"; exit;;
      *)
          printf "%s\n%s is invalid\n" "$helpStr" "$2";
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

if [[ "$useMedakaBl" -eq 1 ]]; then
   extraOptions="-threads $threadsI -enable-medaka";
   extraOptions="$extraOptions -model $medakaModelStr";
else
   extraOptions="-threads $threadsI";
fi # Check if using medaka

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
#    - Build scaffoled (stich amplicons together)
#  o sec-03 sub-05:
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

   "$scriptDirStr/../00-programs/buildCon" \
       -min-read-read-map-length "$minLenI" \
       -min-read-con-map-length "$minLenI" \
       -min-con-length "$minLenI" \
       -min-reads-per-bin "$minDepthI" \
       -prefix "$prefStr" \
       -fastq "$prefStr-tmp.fastq" \
       $extraOptions;

   rm "$prefStr-tmp.fastq"
   rm "$prefStr--top-reads.fastq";

   mv \
       "$prefStr--clust-0--con.fasta" \
       "$inPrefStr-amps-cons/$(basename "$prefStr")-con.fasta";
done # Loop: Build a consensus for each amplicon

#**********************************************************
# Sec-03 Sub-04:
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
# Sec-03 Sub-04:
#  - Clean up and exit
#**********************************************************

rm "$inPrefStr-tmp-primers.fasta";
rm -r "$inPrefStr-tmp-reads";
rm -r "$inPrefStr-amps-cons/";

exit;
