#!/usr/bin/bash

refStr="01-input/01-primer-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.reference.fasta";
readsStr="01-input/01-SRR21813677-twist.fq";
#primersStr="01-input/01-primer-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.scheme.bed ";
primersStr="01-input/01-primers.bed";
  # For some odd reason bash has a ahrd time with the full
  # path for the pimers file

statsFileStr="03-artic/03-artic-stats.tsv";
numMutateI=0;
lenSeqI=0; # Length of reference sequence (found by script)
startSeedI=1;
seedI=1;
repI=1;
maxRepI=5;
medakaModelStr="r941_prom_high_g344";
threadsI=3;
schemeVersionStr="4.1";
schemeDirStr="01-input/01-primer-schemes/"
schemeStr="SARS-CoV-2";
prefixStr="03-artic/03-articBench";

scriptDirStr="$(dirname "$0")";
suffixStr="";

# If I need to make a new stats file
# This header is large, but also allows complete
# replicatoni of the settings.
if [[ ! -f "$statsFileStr" ]]; then
   { # Print out the header
      printf "Program\trep\tNoReads\tPercMutate\ttargMut";
      printf "\tMatches\tsnp\tIns\tDel\trefMatch\trefSnp";
      printf "\trefIns\trefDel\telpTime\tsysTime\tUserTime";
      printf "\tMemKb\tPercCPU\tseed\tusedSeed";
      printf "\tmedakaModel\tlenRef\tRef\tfastq\tprimers";
      printf "\tschemeDir\tscheme\tschemeVersion\tprefix";
      printf "\n";
   } > "$statsFileStr"
fi

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

while [[ "$repI" -le "$maxRepI" ]]; do
#Loop: though each replicate
   seedI="$((startSeedI * repI))";

   for iPercDiff in 0 1 2 5 10; do
   # Loop: Build a consensus for each reference level
      suffixStr="ref-rep$repI-$iPercDiff-perc-dif";
      numMutateI="$(((iPercDiff * lenSeqI) / 100))";
      mv "$refStr" "$backUpRefStr" || exit;
   
      if [[ "$iPercDiff" -gt 0 ]]; then
      # If I am mutated the reference sequence
         "$scriptDirStr/mutateSeq.sh" \
            -seq "$backUpRefStr" \
            -num-mutate "$numMutateI" \
            -seed "$seedI" \
            -primers "$primersStr" \
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
            --scheme-version "$schemeVersionStr" \
            --scheme-directory "$schemeDirStr" \
            "$schemeStr" \
            "$prefixStr-$suffixStr";
   
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
   
      { # Print out the stats
         printf "artic\t%s\t%s\t%s\t%s" \
            "$repI" \
            "$numReadsI" \
            "$iPercDiff" \
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
           "$schemeVersionStr" \
           "$prefixStr";
      } >> "$statsFileStr";
   done # Loop: Build a consensus for each reference level

   repI="$((repI + 1))";
done #Loop: though each replicate
