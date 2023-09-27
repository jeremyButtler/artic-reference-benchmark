# Goal:

The goal of this repository is to see how a poor choice in
  a reference genome affects the accuracy of the consensus
  made by the artic pipeline.

We are also comparing the accuracy of artic to the accuracy
  of LILO, which does not do reference based polishing.

# found error:

I found that the Neeldeman aligner I built (alnSeq) had an
  issue with adding or misplacing two extra deletions. I
  plan to fix these errors later.

To avoid future errors I am switching to quast and am
  restarting my benchmarking (11-quast-bench).

# License

This repository is under a dual license consisting of the
  CC-0 or MIT license. Take your pick of the one you prefer
  to work with.

The code in 00-programs and 00-scripts is my own and so,
  are under the dual license.
However, the code downloaded with the makefile in
  00-programs (LILO, Porechop, quast, and ivar) are not my
  own and so, are not covered by these licenses.
You can remove all complied programs and downloaded code
  with `cd 00-programs && make removeall`.

# Extra steps

## Installing programs

LILO can be installed with the make file in 00-programs
  `cd 00-programs && make`.
This build will call sudo and require you to enter a
  password (for porechop, required by LILO).

You will need to install the artic pipeline (by conda),
  medaka (conda or source), minimap2, and samtools on your
  own.
I should add these to the Makefile at some point, but will
  probably never bother to.

## Preparing data for graphing

The stats files produced are a mixture of space and tab
  deliminators. Also, there might be some runs that might
  have made an emtpy sequence and thus, never had any
  stats. To do this I use:

```
tr ' ' ',' < file-stats.tsv |
  tr '\t' ',' |
  sed 's/,,*/,/g; s/,$//;' |
  tr ',' '\t' |
  awk '{
      if(NR == 1){lenHeadI = NF; print $0; next;};
      if(NF < lenHeadI){next;}; # if row is missing values
      print $0;
    };' > tmp.tsv &&
  mv tmp.tsv file-stats.tsv;

  # for some odd reason sed has an issue with s/\t\t*/\t/g
```

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

I have not filtered this data set yet with a program like
  filtlong, so my results are based on the raw data on the
  SRA. 
This is a bit of a weak point in my study.

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
  indels using quast.
We also measured completeness by the number of N's in each
  consensus and the number of missing bases (done in R).
The quast call was automated with the getStats.sh script in
  00-scripts.

The pipelines we tested were artic version 1.2.3, ivar 
  version 1.4.2, LILO (Downloaded Sep 5th, 2023).
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

We also ran our own pipeline (buildcon). This is not a 
  pipeline that we expect to be used, since it requires
  the user to identify the reads to use to build each
  amplicons consensus. See buildAmpCon.sh for details.

We automated each pipeline using scripts
  (benchArticNoMut.sh, benchLILO.sh, and runIvar.sh),
  which can be found in 00-scripts.

Our levels for this test are read depth, length, and the
  percent difference between the reference and true
  sequence.
For depth, we looked at 30x, 50x, 100x, 200x, 500x, and
  1000x read depth.
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

## found error:

I found that the Neeldeman aligner I built (alnSeq) had an
  issue with adding or misplacing two extra deletions. I
  plan to fix these errors later.

To avoid future errors I am switching to quast and am
  restarting my benchmarking (11-quast-bench).
