#!/usr/bin/bash

# This script is just a quick script to run the same
# data with different Medaka models. Everything is
# hardcoded.

refStr="01-input/01-primer-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.reference.fasta";
readsStr="01-input/01-SRR21813677-twist.fq";
schemeVersionStr="4.1";
schemeDirStr="01-input/01-primer-schemes/"
schemeStr="SARS-CoV-2";
prefixStr="03-artic/03-modelTest";
threadsI=3;
statsFileStr="03-artic/03-model-stats.tsv";

scriptDirStr="$(dirname "$0")";
statsStr="";

printf "Matches\tSNP\tINS\tDEL\n" > "$statsFileStr";

# Supper accuracy model

artic \
   minion \
   --skip-nanopolish \
   --medaka \
   --medaka-model "r941_prom_sup_g507" \
   --threads "$threadsI" \
   --read-file "$readsStr" \
   --scheme-version "$schemeVersionStr" \
   --scheme-directory "$schemeDirStr" \
   "$schemeStr" \
   "$prefixStr";

statsStr="$(\
   bash "$scriptDirStr/getStats.sh" \
      "$refStr" \
      "$prefixStr.consensus.fasta" \
)";

printf \
   "r941_prom_sup_g507 %s\n" \
   "$statsStr" \
  >> "$statsFileStr";

# Remove extra files made by artic
bash "$scriptDirStr/cleanArticOut.sh" "$prefixStr";
rm "$prefixStr.consensus.fasta";

# McGill documentation model

artic \
   minion \
   --skip-nanopolish \
   --medaka \
   --medaka-model "r941_prom_high_g344" \
   --threads "$threadsI" \
   --read-file "$readsStr" \
   --scheme-version "$schemeVersionStr" \
   --scheme-directory "$schemeDirStr" \
   "$schemeStr" \
   "$prefixStr";

statsStr="$(\
   bash "$scriptDirStr/getStats.sh" \
      "$refStr" \
      "$prefixStr.consensus.fasta" \
)";

printf \
   "r941_prom_high_g344 %s\n" \
   "$statsStr" \
  >> "$statsFileStr"

bash "$scriptDirStr/cleanArticOut.sh" "$prefixStr";
rm "$prefixStr.consensus.fasta";

exit;
