#!/usr/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TOC: Table Of Contents
#  - Runs a single benchmark for LILO
#  o sec-01:
#     - Variable declerations
#  o sec-02:
#     - Get, check, and apply user input
#  o sec-03:
#     - Run LILO and get stats
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################################
# Name: benchLILO.sh
# Call:
#    bash benchLILO.sh \
#       -ref ref.fasta \
#       -fastq-dir path/to/fastq/dir \
#       -primer-scheme primer-scheme.bed \
#       -true-ref path/to/answer-real-ref.fasta;
# Use:
#  - Benchmarks how accuractly LILO can bulid a consensus
# Input:
#   -ref: [Required]
#     o Fasta file with the reference sequence
#   -fastq: [Required]
#     o Fastq file with reads to build a consensus for
#   -primer-scheme: [Required]
#     o primer scheme used with artic.
#     o This will be a bed file with 7 columns
#       - Column 4 has the primer name
#       - Column 7 has the primer sequence
#   -true-ref: [Required]
#     o The reference for the actual genome. This is used
#       to detect the number of errors.
#   -prefix: [LILO]
#     o Prefix to name the output fasta files
#   -model: [r941_prom_high_g344]
#     o Model to use with medaka
#   -t: [3]
#     o Number of threads to use
# Output:
#  - prefix-scaffold-auto.fa
#    o Fasta file with scaffold produced by a LILO that
#      did not error out.
#  - prefix-scaffold-manual.fa
#    o Fasta file with consensus with a better call to
#      scaffold_builder
#  - prefix-scaffold-stich.fa
#    o Fasta file with consensus built with stich
#  - prefix-LILO-amps.fa
#    o Amplicons consensuses made with LILO
###########################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#  o sec-01 sub-01:
#  - Variables holding user input or this scripts location
#  o sec-01 sub-02:
#    - Script variables
#  o sec-01 sub-03:
#    - Help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-01 Sub-01:
#  - Variables holding user input or this scripts location
#**********************************************************

scriptDirStr="$(dirname "$0")";
fastqStr="test-reads/tmp.fastq";
refStr="$scriptDirStr/../06-alt-schemes/06-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.reference.fasta";
primersStr="$scriptDirStr/../06-alt-schemes/06-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.scheme.bed";
prefStr="LILO";
threadsI=3;
medakaModelStr="r941_prom_high_g344";

#**********************************************************
# Sec-01 Sub-02:
#  - Script variables
#**********************************************************

liloCsvStr="";    # Points to the csv file for LILO
liloConfigStr=""; # Points to the config file for LILO
tmpPrimStr="";    # Points to the primer bed file for LILO
condaPathStr="$(\
  conda info | 
    grep -i 'base environment' | 
    sed 's/base.*: //; s/  *.read only.//; s/ //g' \
)";

#**********************************************************
# Sec-01 Sub-03:
#  - Help message
#**********************************************************

helpStr="$(basename "$0") \
  -ref ref.fasta \
  -fastq-dir path/to/fastq/dir \
  -primer-scheme primer-scheme.bed \
Use:
  - Benchmarks how accuractly LILO can bulid a consensus
Input:
  -ref: [Required]
    o Fasta file with the reference sequence
  -fastq: [Required]
    o Fastq file with reads to build a consensus for
  -primer-scheme: [Required]
    o primer scheme used with artic.
    o This will be a bed file with 7 columns
      - Column 4 has the primer name
      - Column 7 has the primer sequence
  -prefix: [$prefStr]
    o Prefix to name the output fasta files
  -model: [$medakaModelStr]
    o Model to use with medaka
  -t: [$threadsI]
    o Number of threads to use
Output:
  - prefix-LILO-con.fasta
    o Fasta file with the consensus built by LILO
  - prefix-LILO-amps.fasta
    o Fasta file with the amplicons LILO used to build the
      consensus.
";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Get, check, and apply user input
#  o sec-02 sub-01:
#    - Get user input
#  o sec-02 sub-02:
#    - Check user input
#  o sec-02 sub-03:
#    - Set variables using user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-02 Sub-01:
#  - Get user input
#**********************************************************

while [[ $# -gt 0 ]]; do
# Loop: Get user input
   case $1 in
      -ref) refStr="$2"; shift;;
      -fastq) fastqStr="$2"; shift;;
      -primer-scheme) primersStr="$2"; shift;;
      -prefix) prefStr="$2"; shift;;
      -model) medakaModelStr="$2"; useMedakaBl=1; shift;;
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
#    - Check if fastq is a valid file
#  o sec-02 sub-02 cat-03:
#    - Check if the primer scheme is a file
#  o sec-02 sub-02 cat-04:
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
#  - Check fastq is a file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ ! -f "$fastqStr" ]]; then
   printf \
      " -fastq %s is not a file\n" \
      "$fastqStr" \
      >&2;
   errBl=1;
fi # Check if fastq directory had fastq files

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sec-02 Sub-02 Cat-03:
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
# Sec-02 Sub-02 Cat-04:
#  - Exit if there are errors
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [[ "$errBl" -gt 0 ]]; then exit; fi

#**********************************************************
# Sec-02 Sub-03:
#  - Set variables using user input
#**********************************************************

liloCsvStr="$prefStr-Lilo-tmp-primers.csv";
liloConfigStr="$prefStr-Lilo-tmp-config.txt";
tmpPrimStr="$prefStr-Lilo-tmp-primers-scheme.bed";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Run LILO and get stats
#  o sec-03 sub-01:
#    - Prepare files for input into LILO
#  o sec-03 sub-02:
#    - Run LILO
#  o sec-03 sub-03:
#    - Get and print out stats
#  o sec-03 sub-04:
#    - Clean up
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-03 Sub-01:
#  - Prepare files for input into LILO
#**********************************************************

awk \
   '{if($4 ~ /[Aa][Ll][Tt]/){next;} print $0};' \
    < "$primersStr" \
  > "$tmpPrimStr";

awk \
    -f "$scriptDirStr/schemeToLiloCsv.awk" \
    "$primersStr" \
    > "$liloCsvStr";

{ # Build the config file
   printf "scheme: %s\n" "$tmpPrimStr";
   printf "reference: %s\n" "$refStr";
   printf "primers: %s\n" "$liloCsvStr";
   printf "medaka: %s\n" "$medakaModelStr";
} > "$liloConfigStr" # Build the config file

printf "Compressing fastq file for LILO\n" >&2;

mkdir raw;
cp "$fastqStr" "raw/$(basename "$prefStr")-LILO.fastq";
gzip "raw/$(basename "$prefStr")-LILO.fastq";

#**********************************************************
# Sec-03 Sub-02:
#  - Run LILO
#**********************************************************

source "$condaPathStr/etc/profile.d/conda.sh";

conda activate LILO;

  #env PATH="$binPorechopPathStr:$PATH" \
  # This command messes up poreChop. It has a dependecy
  # issue.
env CONDA_PREFIX=$condaPathStr \
  snakemake \
      -k \
      -s \
      "$scriptDirStr/../00-programs/Lilo/LILO" \
      --configfile "$liloConfigStr" \
      --cores "$threadsI";
conda deactivate;

# For some odd reason LILOS scaffold_build call does a poor
# job. So, I am recalling scaffold_builder

mv \
    "$(basename "$prefStr")-LILO_Scaffold.fasta" \
    "$prefStr-scaffold-auto.fa";
# Save the amplicon consensuses
mv \
   "$(basename "$prefStr")-LILO/polished_trimmed.fa" \
   "$prefStr-LILO-amps.fa";


rm -r "$(basename "$prefStr")-LILO_overlap_alignment";
rm -r "$(basename "$prefStr")-LILO.coords";
rm -r "$(basename "$prefStr")-LILO_output.txt";

conda activate scaffold_builder
scaffold_builder.py \
    -i 75 \
    -t 3693 \
    -g 80000 \
    -r "$refStr" \
    -q "$prefStr-LILO-amps.fa" \
    -p "$(basename "$prefStr")-LILO";

rm -r "$(basename "$prefStr")-LILO_overlap_alignment";
rm -r "$(basename "$prefStr")-LILO.coords";
rm -r "$(basename "$prefStr")-LILO_output.txt";

conda deactivate

#**********************************************************
# Sec-03 Sub-04:
#  - Clean up
#**********************************************************

# Remane the output consensus
mv \
   "$(basename "$prefStr")-LILO_Scaffold.fasta" \
   "$prefStr-scaffold-manual.fa";

"$scriptDirStr/../00-programs/stich" \
    -ref "$refStr" \
    -amps "$prefStr-LILO-amps.fa" \
    -out "$prefStr-scaffold-stich.fa" \
    -overwrite \
    -prefix "$prefStr-stich" \
    -threads "$threadsI";

rm -r raw;
rm "$liloCsvStr";
rm "$liloConfigStr";
rm "$tmpPrimStr";
rm -r "$(basename "$prefStr")-LILO";
rm -r "amplicons.bed";
rm -r "porechop";

exit;
