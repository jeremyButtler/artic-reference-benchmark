# Goal:

The goal of this repository is to see how a poor choice in
  a reference genome affects the accuracy of the consensus
  made by the artic pipeline.

We are also comparing the accuracy of artic to the accuracy
  of LILO, which does not do reference based polishing.

# Methods:

These are very rough, but should give you an idea of what I
  did.

Our goal was to see how accurately the artic pipeline could
  build a consensus as the difference between the
  reference and sequence increased.
We also wanted to see when a method that does not use
  reference polishing, such as LILO, could build more
  accurate consensuses.

For a dataset we used a synthetic SARS-CoV 2 data set
  (twist) that was sequence on a Promethion by the
  McGill Genome Centre (SRA Bio sample SRR21813677).
From McGill's web site we found that this data was likely
  basecalled with Guppy version 3.4.4 using the
  941_prom_high_g344 model.

We were unable to download the complete dataset, but were
  able to download 1.7 million reads (see
  01-input/01-SRR21813677-ids.txt for a list of ids).
We then split these reads into separate amplicons and
  removed any amplicons that had under 1000 read depth.
Also, to avoid having overlapping amplicons (more than just
  the ends), we removed any amplicons that were over 500
  nucleotides or that mapped to a similar position on the
  reference.

  ```
  # Extraction step
  minmap2 -a ref.fasta reads.fastq | 
     awk \
       -f 00-scripts/extReadsByTbl.awk \
       -v extTbl=04-read-table/04-buildConTbl.md \
       -v prefix="amplicons;
  ```

To measure accuracy we looked at the number of SNPs and
  indels using the Needleman alignment from alnSeq
We did not use minimap2, because  minimap2 does not handle
   a large number of Ns very well, which many of our
  consensus had.
We also measured completeness by the number of N's in each
  consensus.
To prevent alignment issues, we ensured that all
  consensuses missing ends filled in with N's before
  aligning.
This was all automated with the getStats.sh script in
  00-scripts.

The pipelines we tested were artic version 1.2.3, LILO
  (Downloaded Sep 5th, 2023), and our own pipeline.
Artic was run with medaka (--medaka --skip-nanopolish) with
  the r941_prom_high_g344 model
  (--medaka-model r941_prom_high_g344) and using the artic
  SARS-CoV 2 V4.1 scheme (--scheme-directory path/to/scheme
  --scheme-version 4.1 SARS-CoV-2).

We ran LILO with three consensus building methods.
The first was using the default consensus output by LILO
  when the conda path error was fixed.
The second was using the polished_trimmed amplicons with
  scaffold_builder using the settings suggested by
  Amanda Warr for manually building the consensus on issue
  #2 in the LILO repository.
The third was with stich, which we used in our own
  pipeline.

We ran LILO with a config file that had the scheme
  (scheme: SARS-CoV-2), reference
  (reference: /path/to/ref.fasta), primers
  (primers: /path/to/primers.scheme.bed), and medaka model
  (medaka: r941_prom_high_g344).
LILO was then called with env to correct the path to
  conda, which causes LILO to crash.
   ```
    env CONDA_PREFIX="$(conda info |
      grep -i 'base environment' | 
      sed 's/base.*: //; s/  *.read only.//; s/ //g' \
   )" -k -s /path/to/Lilo/LILO --configfile config.txt
   ```

For our own pipeline, we first mapped the reads to the
  primers with minimap2 (minimap2 -k5 -w1 -s 20 -P -m 15)
  and then trimmed the mapped primers at the ends of each
  read.
We then trimmed off any junk regions (regions that did not
  map) in the reads using minimap2 and the reference.
After that, we build a consensus for each amplicon using
  buildCon, which uses minimap2
  (see our find--Co-infections repository).
Finally, we stitched the amplicons together using stich
  with minimap2 support (only option currently)
  (see 00-programs).
All of our own code we used for each of these steps can be
  found in 00-programs.

We automated each pipeline using scripts
  (benchArticNoMut.sh, benchLILO.sh, and buildAmpCons.sh),
  which can be found in 00-scripts.

Our levels for this test are read depth, length, and the
  percent difference between the reference and true
  sequence.
For depth, we looked at 30x, 50x, 100x, 200x, 500x, 1000x,
  and 1200x read depth.
We selected the reads for each read depth by randomly
  subsampling each amplicons file.
The seed for each subsample was the same for all depths,
  but changed between replicates.
This step was automated with the subsampleReads.sh script
  in 00-scripts.

For consensuses length we looked at 1700 nucleotides,
  5700 nucleotides, 11800 nucleotides, and 29000
  nucleotides long. 
Scheme files and reference files had entries/parts of 
  sequences removed that did not match these lengths.
To prevent errors, we made sure to only take only the first
  X nucleotides for each sequence.

To make the reference different than the true consensus we
  mutated the reference.
We tested cases were the reference was mutated to 0%, 1%,
  2%, 5%, and 10% different.
When mutating we made sure to not mutate and primer regions
  or areas were high depth consensus had gaps.
Like, subsampling, each level used the same seed, with the
  seed changing between replicates.
This step was automated with a script (mutateSeq.sh).

We ran each test for five replicates.
Also, we automated all of this with the benchAll.sh script
  in 00-scripts.

# Results:

The current graphs is somewhat simple, but shows how the
  reference can change the snp accuracy of the artic
  pipeline. This graph was made using the tsv file
  (filter.tsv) in the 07 folder. There are also some
  initial 30x read depth results in the 09 directory.

![Graph showing how the SNP accuracy of the artic pipeline
  is affected by the reference. Accuracy decreases as the
  reference becomes more distant from the sequence.
  However, even at 10%, the accuracy rarely goes beneath
  99.8% and never beneath 99.75%. Read depth seems to have
  little affect on accuracy. I suspect this is due to
  masking.
](artic-initial-bench.svg)

The figure above shows how the choice of reference can
  change the accuracy of the artic pipeline. However, the
  accuracy is never beneath 99.85%. There is no clear trend
  on the affect of depth.
