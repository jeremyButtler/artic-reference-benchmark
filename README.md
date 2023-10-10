# Goal:

The goal of this repository is to see how a poor choice in
  a reference genome affects the accuracy of the consensus
  made by the artic pipeline.

We are also comparing the accuracy of artic to the accuracy
  of LILO, which does not do reference based polishing.

# Note about the writing

This is not going to be the most well written repository.
  That is not my goal. This is just here to show what I
  found.

# found error:

I found that the Neeldeman aligner I built (alnSeq) had an
  issue with adding or misplacing two extra deletions. I
  plan to fix these errors later.

To avoid future errors I am switching to quast and am
  restarting my benchmarking (11-quast-bench).
I changed the Ns, snps, and indels from the 100kb used in
  quast to per 1kb.

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

I needed to clean up the stats from the quast benchmarking.
  Something was happening with the medaka polishing runs
  for buildcon. Also I had some errored out entires at low
  read depths or early on.

```
tr '\t' ',' < 11-quast-bench/11-quast-bench-stats.tsv |
  sed 's/[, ][, ]*/,/g; s/[, ]*$//g;' |
  tr ',' '\t' | 
  awk \
      'BEGIN{OFS=FS="\t"};
       {
          if(NR == 1){headI = NF; print $0; next;};
          if($28 == "TRUE") next; 
          if(NF == headI) print $0;
          #print NF;
      };' \
  > 11-quast-bench/11-quast-bench-cleanup.tsv;
```

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
  00-scripts (getStats2.sh gets more data).

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
The third was with stich (in 00-programs).

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
  amplicons consensus. We used only the majcon step (not
  Racon or Medaka), but we tested with and without ivar
  polishing. See buildAmpCon.sh in 00-scripts for details
  on how this pipeline ran.

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
I am currently not planning to graph this data, but it is
  their if you want to try.

To make the reference different than the true consensus we
  mutated the reference.
We tested cases were the reference was mutated to 0%, 1%,
  2%, 5%, and 10% different.
When mutating we made sure to not mutate and primer regions
  or areas were high depth consensus had gaps.
Like, subsampling, each level used the same seed, with the
  seed changing between replicates.
This step was automated with a script (mutateSeq.sh).

We ran each test for ten replicates.
Also, we automated all of this with the benchAll.sh script
  in 00-scripts.

The computer used for benchmarking had an Intel Core
  i5-2500 CPU @ 3.30GHz and 8 gigs of ram.

There is a benchAll2.sh and getStats2.sh, which I never
  used. I made these scripts to get data I missed in my
  first quast run, but found that my results were good
  enough for what I wanted.  I have tested getStats2.sh
  for water from EMBOSS, but not benchAll2.sh. These should
  be better scripts, if they work.

# Results:

Our first step was to narrow down Ivar, LILO, and buildcon
  to a single combination.
We found that Ivar with its built in trimming step was less
  accurate than Ivar with trimPrimers (ivar graphs in
  11-quast-bench), or without primer any trimming (ivar
  graphs in 10-new-bench).
For LILO, we found that LILO with a manual scaffold_builder
  step with Ivar polishing or stich step with Ivar
  polishing was more accurate than the other tested 
  combinations for LILO (lilo graphs in 11-quast-bench).
We found that buildcon with Ivar polishing was the only
  option with any level of accuracy (buildCon graphs in
  11-quast-bench).

![
  Graph showing the number of false snps and indels of
  artic, Ivar, LILO, and buildcon when the read depth is
  over 30x.
](11-quast-bench/all-over30x-snps-indels.svg)

**Figure 1:**
  Graph showing the number of false snps and indels of
  artic, Ivar, LILO, and buildcon when the read depth is
  over 30x.
  
We looked at the indel and snp accuracy of each
  consensus at read depths over 30x and when the reference
  was 0%, 1%, 2%, 5%, and 10% different.
We found that Ivar had the least number of false positive
  snps and indels, followed by either LILO or buildcon
  (figure 1).
Both LILO and buildcon had outliers that were close or had
  over 0.5% false positive snps (figure 1).
However, artic had more outliers then LILO or buildcon
  (figure 1).

We thought that read depth might be a factor, so we made
  graphs for reach read depth (see
  11-quast-bench/all-depthX-snps-indels.svg).
However, we found that in some cases buildcon had less
  outliers at 50x than 500x
  (See 11-quast-bench/all-50X-snps-indels.svg and
   11-quast-bench/all-50X-snps-indels.svg).
We also found that the worst outliers for LILO were at
  1000x read depth
  (see 11-quast-bench/all-1000X-snps-indels.svg).
This could be due to not running enough replicates and so,
  these may be rare cases.
It also could be due to having more complete genomes at 
  deeper read depths, which may include amplicons that
  were more error prone
  (see 11-quast-bench/all-masked.svg).

![
  Graph showing how the number of false positive snps
  changes as the reference becomes more distant from the
  true sequence.
](11-quast-bench/all-over30x-snps-ref.svg)

**Figure 2:**
  Graph showing how the number of false positive snps
  changes as the reference becomes more distant from the
  true sequence. Read depth is over 30x.

Next we looked at how the reference was influencing the
  snp accuracy of our consensus.
We found that reference choice had little affect on
  buildcon, possibly caused a few outliers with over ten
  false positive snps LILO when 10% different, may have
  added an false positive snp for Ivar, and had a large
  increase in false positive snps for artic (figure 2).
For the artic pipeline, we found that even when the
  reference was 1% to 2% different than the true sequence,
  we still had several consensuses with at least one extra 
  false positive snp (figure 2).
We also found less extreme outliers in LILO than buildcon
  or artic when the reference was under 10% different
  (figure 2).

The increase in false positive snps in artic would explain
  many of the outliers in artic in the snp/indel graph
  (figure 1 and 2), which are not explained by read depth
  (11-quast-bench/all-over30x-snps-ref.svg).
Also, we noticed that most consensuses made using Ivar did
  not have the extra false positive snp, even when the
  reference was 10% different
  (see 11-quast-bench/ivar-snps-depth.svg).

![
  Graph showing how depth (30x to 1000x) changes the
  completeness of built consensuses.
](11-quast-bench/all-masked.svg)

**Figure 3:**
  Graph showing how depth (30x to 1000x) changes the
  completeness of built consensuses.

Next we looked at how the read depth changed the
  completeness of each consensus.
We found that Ivar made the most complete consensuses,
  followed by artic, then buildcon, which was closely
  followed by LILO
  (figure 3).
We also found that increasing read depth led to an increase
  in the consensus completeness for every pipeline except
  artic (figure 3).
Also, we noticed that both buildcon and LILO were
  less predictable about how complete a genome was at any
  read depth than artic are Ivar (figure 3).

![
  Graph showing how close the consensus is to the real
  sequence as the reference changes. Reads depths are all
  over 30x.
](11-quast-bench/all-over30x-length-ref.svg)

**Figure 4:**
  Graph showing how close the consensus is to the real
  sequence as the reference changes. Reads depths are all
  over 30x.

We wanted to look at how the length of the consensus
  changed as the reference became more distant to the true
  sequence.
We found that only LILO with scaffold_builder had more than
  one or two outliers for a single reference (figure 4).
In this case scaffold_builder inserted between 100 to 2000
  bases in almost all sequences when the reference was 10%
  different (figure 4).
However, these large insertions were not present when the
  reference was 5% or less different than the true
  sequence (figure 4).
Also, LILO with stich did not suffer from this error, which
  suggests that this is a problem with scaffold_builder or
  our settings for scaffold_builder and not LILO itself
  (figure 4).
There was a similar error, were LILO produced consensuses
  with 40 to 120 false positive indels when the
  scaffold_builder step in LILO did not error out
  (see 11-quast-bench/lilo-snps-indels.svg; this does not
   account for added bases).

![
  Graph showing the time usage of each pipeline as read
  depth increases (starts at 50x read depth).
](11-quast-bench/all-over30x-time.svg)

![
  Graph showing the memory usage of each pipeline as read
  depth increases (starts at 50x read depth).
](11-quast-bench/all-over30x-memory.svg)

**Figure 5:**
  Graph showing the time usage of each pipeline as read
  depth increases (starts at 50x read depth). (a) is time
  usage in minutes. (b) is memory usage in megabytes.

We also wanted to look at the time usage and memory usage
  of each pipeline as read depth increased.
At all read depths, we found that Ivar was the fastest
  program (under two minutes), followed by buildcon
  (under four minutes), then artic (under five minutes),
  and finally LILO (around 20 minutes at 1000x read depth)
  (figure 5a).
Of all programs, artic seemed to have the least increase in
  time usage as read depth increased (figure 5a).

One note I will add is that I know that buildcon, when used
  with Medaka is slower (I think around 30 minutes for
  1000x read depth) than LILO.
So, these time stats only apply to the fastest consensus
  method in buildcon.

For memory usage we noticed that artic used the most
  memory (over 750 mb), followed by LILO (under 600 mb),
  and then followed by either buildcon or Ivar
  (under 200 mb) (figure 5b).
However, the memory usage of buildcon was lower than Ivar
  only after 500x read depth (figure 5b).
This is likely due to buildcon only using only a subsample
  of all reads (300 reads) in its consensus building step.

Again the memory usage for buildcon does not include a
  Medaka step.
The memory usage will be the same as LILO when Medaka is
  used.

# Summary

We set out to determine how the reference affects the
  accuracy of consensus made with artic, Ivar, and LILO.
We found that the reference had little affect on
  consensuses made using Ivar.
That the reference added in many insertions into
  consensuses built by LILO when the reference was 10%
  different than the true sequence, but had little affect
  when the reference was 5% different than the true
  sequence.
That the reference had an affect on the accuracy
  of the artic pipelines consensus, but that it took over
  10% difference to get more than 0.05% increase in snps.

If you are aiming for the most accurate consensus we would
  recommend that you check how close your reads are to the
  reference.
We would then recommend only using artic if the reference
  is 1% or at most 2% different than the possible
  consensus.
The alternatives to artic we would recommend are Ivar or
  LILO combined with Ivar polishing and either a manual
  scaffold_builder call when the reference is under 5%
  different or stich if over 5% different.

We should note that Ivar has not been tested with
  deletions in the reference, which would appear as
  insertions in the reads.
Since Ivar has a stiff insertion threshold (90%), it is
  possible that Ivar would favor deletions in these cases.
However, the most likely targets would be larger
  homopolymers, due to Nanopore reads having more errors
  and thus, less support in larger homopolymer regions.
This error will not likely effect the polishing the
  consensus made by LILO, since Ivar is polishing the
  consensus instead of a reference.
  
Overall, we have no idea what to recommend, except that you
  should make sure that you check the reference when using
  the artic pipeline.
However, if you want to prevent reference errors, we would
  recommend to always check your reference (no mater the
  pipeline), to use multiple references in your check, and
  to select and use the best reference from the check.
This should ensure that the closest reference is always
  used.

# Thanks

- Eric Bortz for encouraging me and pointing me to the type
  of data to look for.
- The Drown lab (Devin Drown and many others) for there
  encouragement.
- My Dad for listening to me talk about this (a lot) and
  finding a computer I could run this on. One thing that
  amazes me is just how much you can do with a little.
