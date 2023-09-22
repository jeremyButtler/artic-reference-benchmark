# Goal:

See how the reference changes the SNP accuracy of the
  artic pipeline.

# Starting out; needed programs and the 00 directories

The 00 directories are not input our a result of processing
  the input. These contain my notes (00-notes), scripts I
  used (00-scripts), and some of the programs that my
  scripts used (00-programs). For the programs I include a
  Makefile to make these programs and ensure they are in
  the correct directory for the scripts that I use.

Before trying to replicate my notes do
     `make -C 00-programs`
  or "`cd 00-programs` and then `make`".
 
After this you will need to install the artic pipeline and
  minimap2.

# Getting data (01-input)

Go to JUMP-RETRY (at bottom) to see how you might do this.
  I acidently wiped the data and had to redownload it.

## Download the data (01)

I download twist (a synthetic SARS-CoV 2 sequence) from
  the SRA. This sample was amplified using the artic-V4.1
  protocol. It was then sequenced on a Nanopore PromethION.

|  Category   |     Information      |
|:-----------:|:--------------------:|
| Bio-project |     SRX17803051      |
| Bio-sample  |     SRR21813677      |
|    Guppy    |   version unknown    |
|    Guppy    |    model unknown     |
|  Published  |     2022-10-05       |
|     ID      |      24726887        |
| Sequencer   | McGill Genome Centre |
|   Primers   |      artic-V4.1      |

The sequences were downloaded from the SRA with
  fasterq-dump version 2.10.9.

`fasterq-dump --split-files SRR21813677`

 fasterq-dump ran out of memory (I think tmp or var) after
  downloading the first 1785053 reads. I used ctrl-C to
  force fasterq-dump to quite before it deleted the
  downloaded reads. Then I removed the last read in every
  split (these reads were not part 1785053 reads) to
  make sure I had a complete fastq entry for each file. I
  then merged the files with
  `cat *SRR21813677*/* > 01-input/01-SRR21813677-twist.fq`.

  **Do not remove the last entry if fasterq-dump downloaded
    all the files. Only do this if fasterq-dump hand an
    error. If fasterq-dump had an error, restart the
    download and hit ctrl-c when you start to see error
    messages. Otherwise fasterq-dump will delete the
    files. I think this error is from having to small of
    a tmp directory**

This next step is not needed, but I did to keep the
  original files. You can do `rm -r *SRR21813677*` instead.

```
# Rename downloaded split files and move to 01-input

mv \
   *SRR21813677* \
   01-input/01-SRR21813677-Partial-twist-split-files;
```

## Extracting the reads I used (01)

I extracted the fastq ids I was able to download to
  01-input/01-SRR21813677-ids.txt with sed.

  **DO NOT RUN THIS COMMAND. THIS IS TO SHOW YOU HOW I
    GOT THE READ IDS. IF YOU USED THIS, YOU WOULD OVERWRITE
    MY READ IDS AND NO LONGER HAVE THE FULL LIST.**

  ```
  sed \
     -n \
     's/^@//;   # Remove @ symbols at start of read ids
      s/ .*//p; # Remove meta data in read ids
      n;n;n;    # Move to the next read id
     ' \
     < 01-input/01-SRR21813677-twist.fq \
    > 01-input/01-SRR21813677-ids.txt;
  ```

## Extracting reads from the fastq file (01)

You can extract the read ids I used from the full twist
  file with fqGetIds or seqkit. fqGetIds is part of the
  programs compiled in 00-programs.

```
# Rename downloaded file (avoid it being used later)
mv \
   01-input/01-SRR21813677-twist.fq \
   01-input/01-SRR21813677-all-twist.fq;

# Extract the read ids I used
00-programs/fqGetIds \
   -f 01-input/01-SRR21813677-ids.txt \
   -fastq 01-input/01-SRR21813677-all-twist.fq \
   -out 01-input/01-SRR21813677-twist.fq;
```


## Find the list of Medaka models (01)

Get models supported by medaka:
   `medaka tools list_models | tr ',' '\n' | grep "_prom"`.

- Possible Models:
  - r103_prom_high_g360
  - r103_prom_snp_g3210
  - r103_prom_variant_g3210
  - r941_prom_fast_g303
  - r941_prom_fast_g507
  - r941_prom_fast_snp_g507
  - r941_prom_fast_variant_g507
  - r941_prom_hac_g507
  - r941_prom_hac_snp_g507
  - r941_prom_hac_variant_g507
  - r941_prom_high_g303
  - r941_prom_high_g330
  - r941_prom_high_g344
  - r941_prom_high_g360
  - r941_prom_high_g4011
  - r941_prom_snp_g303
  - r941_prom_snp_g322
  - r941_prom_snp_g360
  - r941_prom_sup_g507
  - r941_prom_sup_snp_g507
  - r941_prom_sup_variant_g507
  - r941_prom_variant_g303
  - r941_prom_variant_g322
  - r941_prom_variant_g360

## Find model McGill used (01)

Looked at McGill protocls for Guppy
([https://c3g.github.io/covseq_McGill/SARS_CoV2_Sequencing/ONT_overview.html](
  https://c3g.github.io/covseq_McGill/SARS_CoV2_Sequencing/ONT_overview.html))
  and it looks like they are using 9.4.1 flow cells with
  hac (high accuracy) models. The guppy version is 3.4.4.
  The copywrite date on the web page is 2020, so might be
  out of date, though, this might be when my data was
  sequenced.

The model that would be used for R9.4.1 flow cells, with
  Guppy version 3.4.4, and the high accuracy model is
  941_prom_high_g344.

## Download scheme (01)

McGill used the artic SARS-CoV 2 version 4.1 scheme to
  sequence there data. This can be downloaded from github.

```
git clone https://github.com/artic-network/primer-schemes;
mv \
   primer-schemes/SARS-CoV-2/V4.1/* \
   01-input/01-primer-schemes;
```

I also removed and non-SARS-CoV 2 schemes. However, this
  is not needed.

# Read stats (02)

You can probably get a better Q-score stats with nanostat.
  However, I did not have this installed, so I used 
  scoreReads, which is a bit optimistic for Q-scores. It
  looks like it somewhat similar to the spot on the SRA.

Do not trust Q-scores above 13 to much.

If you are looking at the tsv file output by scoreReads,
  then you should know that totalDeletions/totalInsertions
  is the number of deletions/insertions. Insertions and
  Deletions is the number of deletions/insertions
  that are not in homopolymers larger then the input size.
  The default max homopolymer size is 1 (no homopolymer).

If you want to save the sam file, then run minimap2
  separately and use `scoreReads -file` instead of
  `scoreReads -stdin`.

```
# Get stats
minimap2 \
      --eqx \
      -a \
      01-input/01-articFiles/MN908947.fa \
      01-input/01-SRR21813677-twist.fq |
   scoreReads \
      -stdin \
      -min-median-q 0 \
      -min-mean-q 0 \
      -min-aligned-median-q 0 \
      -min-aligned-mean-q 0 \
      -min-read-length 0 \
      -max-read-length 0 \
      -min-q 0 \
   > 02-stats/02-read-stats.tsv;
```

```
# Veiw Median Q-scores
awk \
      '{sub(/\..*/, "", $11); print $11;}' \
      < 02-stats/02-read-stats.tsv |
   sort -n |
   uniq -c |
   less;
```

```
# Veiw Mean Q-scores
awk \
      '{sub(/\..*/, "", $12); print $12;}' \
      < 02-stats/02-read-stats.tsv |
   sort -n |
   uniq -c |
   less;
```

# Run artic (03-artic)

**This set is not needed.**

## Check model (03)

I did a quick check to see if the McGill model was a good
  choice. I just wanted to make sure they did not use the
  supper accuracy model. I hardcoded these two runs in
  `bash 00-scripts/modelTest.sh`. The output file was
  saved as 03-artic/03-model-stats.tsv.

## Run Artic (03)

I did an initial run of artic using an hardcoded bash
  script for the files I downloaded and worked with.
  `bash benchArticHardcode.sh`. This script saved the stats
  to 03-artic/03-artic-stats.tsv. It also saved the
  consensus for each replicate.

I initially ran this hardcoded script with the
  supper-accuracy model. However, I then reran this with
  the correct model. The consensuses were overwritten, but
  the stats are still in 03-artic/03-artic-stats.tsv.

This benchmarking script calls a mutate bash script, which
  mutates any position that a primer does not target. The
  seed in benchArticHarcode.sh is set to 1, so it will
  always give the same results. However, the seed given to
  the mutate script (00-scripts/mutateSeq.sh) is
  1 \* replicate number. So, the results will differ by
  replicate, when the genome is mutated.

I would recommend checking the stats file with sc, if you
  have sc installed (`apt-get install sc`). Sc can be used
  like less (same hot keys) for tsv, csv, and space
  delimited files. It does not handle really large spread
  sheets well, but for small to median size files it works.
  You need psc to translate the files
  (`psc -d ','` for csv).

`cat 03-artic/03-artic-stats.tsv | psc | sc`

## Results (03)

**This was only based on one replicate, so it does not
  mean a lot. I would suggest skipping this step, it is
  not needed.**

From these initial results, it looks like the reference is
  having a small effect on the consensus, but not a large
  effect. With no mutations, there were 4 SNPs and 10
  deletions. At 10% at least 20 bases of roughly 29889
  bases were SNPs, however, there was a large amount of
  variation (max 45 SNPs). Sometimes there were 11 or 12
  deletions or one or two insertions. At 2% there were less
  then 10 bases.

It is interesting that the supper accuracy and high
  accuracy model are giving similar results. I did not have
  this when I was benchmarking the ASHURE data set with
  find--co-infections (this is from memory).

So, the accuracy is not effected by more than 0.2% when the
  reference was 10% different. This is very small, but
  still marks decreased accuracy with poor reference
  choice. I would need to compare this to a non-reference
  corrected method to see if this loss is worth it. I will
  probably build a consensus by hand to do this, using a
  similar workflow to the TBEV data I am working/worked on.

Also, I should improve the benchmarking by checking how
  read depth, read length, and filtering effect this
  result.

# Get read table (04-read-table)

I want to get a read table of the positions of each read.
  This way I can extract the reads mapping to primer
  positions. I am making a temporary sam file and I will
  use this to sift though my reads.

This will give me finer control over read depth.

```
refDirStr="01-input/01-primer-schemes/SARS-CoV-2/V4.1";
minimap2 \
    --eqx \
    -a \
    "$refDirStr/SARS-CoV-2.reference.fasta" \
    01-input/01-SRR21813677-twist.fq \
  > tmp.sam

bash 00-scripts/readLenPosTbl.sh \
   -sam tmp.sam \
   -min-reads 100 \
  > 04-read-table/04-read-table.md;
```

This creates a file with many entries. In this case it
  looks like one entry per 100 base pairs. At this point
  I need to manually shift though these files to find the
  reads I want to keep. 

It looks like most of the kept entries have at least
  1000 reads. There are a few exceptions (> 100, < 200) or
  no reads.

| No. Reads |        Reference        | Position | Length |
|:---------:|:-----------------------:|:--------:|:------:|
| 2388      | MN908947.3              | 0        | 500    |
| 20320     | MN908947.3              | 300      | 500    |
| 10992     | MN908947.3              | 600      | 500    |
| 21451     | MN908947.3              | 900      | 500    |
| 10067     | MN908947.3              | 1200     | 500    |
| 22901     | MN908947.3              | 1500     | 500    |
| 13094     | MN908947.3              | 1800     | 500    |
| 4535      | MN908947.3              | 2100     | 500    |
| 2668      | MN908947.3              | 2400     | 500    |
| 4212      | MN908947.3              | 2800     | 500    |
| 3210      | MN908947.3              | 3300     | 500    |
| 1527      | MN908947.3              | 3600     | 500    |
| 2564      | MN908947.3              | 3900     | 500    |
| 1487      | MN908947.3              | 4300     | 500    |
| 142       | MN908947.3              | 4600     | 400    |
| 535       | MN908947.3              | 5000     | 400    |
| 6365      | MN908947.3              | 5200     | 500    |
| 10969     | MN908947.3              | 5500     | 500    |
| 12911     | MN908947.3              | 5800     | 500    |
| 1622      | MN908947.3              | 6100     | 500    |
| 15463     | MN908947.3              | 6400     | 500    |
| 6145      | MN908947.3              | 6700     | 500    |
| 4714      | MN908947.3              | 7000     | 500    |
| 11465     | MN908947.3              | 7300     | 500    |
| 6404      | MN908947.3              | 7600     | 500    |
| 8029      | MN908947.3              | 7900     | 500    |
| 6464      | MN908947.3              | 8300     | 500    |
| 2247      | MN908947.3              | 8500     | 500    |
| 7861      | MN908947.3              | 8900     | 500    |
| 843       | MN908947.3              | 9100     | 500    |
| 5281      | MN908947.3              | 9400     | 500    |
| 165       | MN908947.3              | 9700     | 300    |
| 26726     | MN908947.3              | 10000    | 500    |
| 8977      | MN908947.3              | 10300    | 500    |
| 8700      | MN908947.3              | 10700    | 500    |
| 9240      | MN908947.3              | 11000    | 500    |
| 15288     | MN908947.3              | 11300    | 500    |
| 3181      | MN908947.3              | 11600    | 500    |
| 7617      | MN908947.3              | 11900    | 500    |
| 5905      | MN908947.3              | 12200    | 500    |
| 7285      | MN908947.3              | 12500    | 500    |
| 2873      | MN908947.3              | 12800    | 500    |
| 11705     | MN908947.3              | 13100    | 500    |
| 2330      | MN908947.3              | 13300    | 500    |
| 2320      | MN908947.3              | 13400    | 500    |
| 5732      | MN908947.3              | 13700    | 500    |
| 4182      | MN908947.3              | 14000    | 500    |
| 12330     | MN908947.3              | 14300    | 500    |
| 2153      | MN908947.3              | 14500    | 500    |
| 121       | MN908947.3              | 14500    | 700    |
| 3004      | MN908947.3              | 15100    | 500    |
| 31300     | MN908947.3              | 15500    | 500    |
| 6796      | MN908947.3              | 15800    | 500    |
| 15152     | MN908947.3              | 16100    | 500    |
| 5280      | MN908947.3              | 16300    | 500    |
| 3759      | MN908947.3              | 16600    | 500    |
| 1646      | MN908947.3              | 16900    | 500    |
| 7175      | MN908947.3              | 17300    | 500    |
| 5042      | MN908947.3              | 17600    | 500    |
| 4092      | MN908947.3              | 17900    | 500    |
| 20664     | MN908947.3              | 18200    | 500    |
| 16081     | MN908947.3              | 18500    | 500    |
| 14701     | MN908947.3              | 18800    | 500    |
| 190       | MN908947.3              | 19100    | 500    |
| 1284      | MN908947.3              | 19400    | 500    |
| 5798      | MN908947.3              | 20000    | 500    |
| 2199      | MN908947.3              | 20300    | 500    |
| 5416      | MN908947.3              | 20600    | 500    |
| 6011      | MN908947.3              | 20900    | 500    |
| 2368      | MN908947.3              | 21200    | 500    |
| 16492     | MN908947.3              | 21500    | 500    |
| 1627      | MN908947.3              | 21800    | 500    |
| 10175     | MN908947.3              | 22000    | 500    |
| 1248      | MN908947.3              | 22400    | 500    |
| 1838      | MN908947.3              | 22600    | 500    |
| 7756      | MN908947.3              | 22700    | 500    |
| 3481      | MN908947.3              | 22900    | 500    |
| 14812     | MN908947.3              | 23200    | 500    |
| 1975      | MN908947.3              | 23500    | 500    |
| 11073     | MN908947.3              | 23800    | 500    |
| 2143      | MN908947.3              | 24100    | 500    |
| 18740     | MN908947.3              | 24400    | 500    |
| 134       | MN908947.3              | 24400    | 900    |
| 4665      | MN908947.3              | 25000    | 500    |
| 5513      | MN908947.3              | 25300    | 500    |
| 4798      | MN908947.3              | 25600    | 500    |
| 13066     | MN908947.3              | 25900    | 500    |
| 15511     | MN908947.3              | 26200    | 500    |
| 22689     | MN908947.3              | 26500    | 500    |
| 2768      | MN908947.3              | 26800    | 500    |
| 7760      | MN908947.3              | 27100    | 500    |
| 8887      | MN908947.3              | 27400    | 500    |
| 8313      | MN908947.3              | 27700    | 500    |
| 11305     | MN908947.3              | 27900    | 500    |
| 14041     | MN908947.3              | 28100    | 500    |
| 18894     | MN908947.3              | 28500    | 500    |
| 5658      | MN908947.3              | 28800    | 500    |
| 12270     | MN908947.3              | 29100    | 500    |
| 1695      | MN908947.3              | 29400    | 500    |

Table: Reads I am am planning to extract. The table is from
  04-read-table/04-read-table-edits-fqExtract.md.

# Extract amplicons (05-amplicon-reads)

## Initial extract (05)

I am using the edited read table to extract the reads that
  I am interested in using for downsampling. These may be
  less, but are at the primer positions or fill in some
  missing regions.

**This is what I start with. However, this will not get
  what I am benchmarking with. I have removed a few of
  the lower read depth files.**
```
awk \
   -f 00-scripts/extReadsByTbl.awk \
   -v extTbl="04-read-table/04-read-table-edits-fqExtract.md"\
   -v prefix="05-amplicon-reads/05-500bp-reads" \
   < tmp.sam;

rm tmp.sam; # Files are large, so try to remove temp files
```

About 772991 reads (797 Mb) were kept from the 1785053
  (2 Gb) of downloaded reads.

## Extract the reads I used for benchmarking (05)

This command will get you the reads I used for my
  benchmarking step. The difference between
  04-read-table-edits-fqExtract.md and 04-buildConTbl.md 
  is that 04-buildConTbl.md has all amplicons with less
  than 1000 reads removed. This is what I benchmarked with.

```
awk \
   -f 00-scripts/extReadsByTbl.awk \
   -v extTbl="04-read-table/04-buildConTbl.md"\
   -v prefix="05-amplicon-reads/05-500bp-reads" \
   < tmp.sam;

rm tmp.sam; # Files are large, so try to remove temp files
```
## Extra removals (05)

**You can skip this if you used 04-buildConTbl.md to
  extract reads**

I wanted to remove any amlicon that had under 1000 reads.
  This way I could sample to 1000 in read depth.

| 142       | MN908947.3              | 4600     | 400    |
| 535       | MN908947.3              | 5000     | 400    |
| 843       | MN908947.3              | 9100     | 500    |
| 165       | MN908947.3              | 9700     | 300    |
| 121       | MN908947.3              | 14500    | 700    |
| 190       | MN908947.3              | 19100    | 500    |
| 134       | MN908947.3              | 24400    | 900    |

```
mv \
   05-500bp-reads-MN908947.3-pos4600-len400.fastq \
   05-500bp-reads-MN908947.3-pos5000-len400.fastq \
   05-500bp-reads-MN908947.3-pos9100-len500.fastq \
   05-500bp-reads-MN908947.3-pos9700-len300.fastq \
   05-500bp-reads-MN908947.3-pos14500-len700.fastq \
   05-500bp-reads-MN908947.3-pos19100-len500.fastq \
   05-500bp-reads-MN908947.3-pos24400-len900.fastq \
   05-under1kDepth;

# Get read counts
wc -l * |
   awk '{print $1 / 4, $2}' |
   sort -n -k 1 |
   less;
```

The lowest read depth is now 1248.

# Set up alternate schemes (06-alt-schemes)

The goal here is to to make some schemes that only use the
  first part of the genome. That way I can test the affect
  of length. I will do this by hand and just keep the first
  X bases for every scheme. I am currently thinking of the
  first 1500 bases, 5000 bases, 11000 bases, and full
  genome.

## Build the no mutate file (06)

I need to make a modified file for mutations that ignores
  all regions that I do not have amplicon depths for and
  primer targets.

I build a consensus using artic and the combined fastq
  files from 05-amplicon-reads. I then used the output
  consensus (06-alt-schemes/06-consenssu-with-Ns.fa) to
  find all N's. These positions were then marked in 
  06-alt-schemes/06-no-mutate-regions.tsv.

```
cat 05-amplicon-reads/*.fastq > tmp.fq;

bash 00-scripts/benchArtic.sh \
   -fastq tmp.fq \
   -scheme-dir 01-input/01-primer-schemes/ \
   -scheme "SARS-CoV-2" \
   -scheme-version "4.1" \
   -prefix test \
   -rep 1;

mv \
   test-ref-rep1-0-perc-dif-con.fa \
   06-alt-schemes/06-consenssu-with-Ns.fa;
rm test-artic-bench-stats.tsv; # I do not care about this
rm tmp.fq;
```

## Build the shorted schemes (06)

My goal is to test lengths of 1500 bases, 5000 bases,
  11000 bases, and the full genome. To do this I need to
  copy and modify the default scheme.

```
cp \
   -r \
   01-input/01-primer-schemes/ \
   06-alt-schemes/06-schemes;

cp \
   -r \
   06-alt-schemes/06-schemes/SARS-CoV-2/V4.1/ \
   06-alt-schemes/06-schemes/SARS-CoV-2/V4.1.1500;

cp \
   -r \
   06-alt-schemes/06-schemes/SARS-CoV-2/V4.1/ \
   06-alt-schemes/06-schemes/SARS-CoV-2/V4.1.5000;


cp \
   -r \
   06-alt-schemes/06-schemes/SARS-CoV-2/V4.1/ \
   06-alt-schemes/06-schemes/SARS-CoV-2/V4.1.11000;
```

### 1500 scheme (06)

06-alt-schemes/06-schemes/SARS-CoV-2/V4.1.1500

For the 1500 scheme I primers 1 to 5. I also left 50 extra
  bases at the end of the sequence, for a total of 1700
  bases.

### 5000 scheme (06)

06-alt-schemes/06-schemes/SARS-CoV-2/V4.1.5000

For the 5000 scheme I kept primers 1 to 18. I also left 50
  extra bases at the end for a total of 5700 bases. The
  one pair of alternate primers (primer pair 10) was
  retained.

### 11000 scheme (06)

06-alt-schemes/06-schemes/SARS-CoV-2/V4.1.11000

For the 11000 scheme I kept primers 1 to 38. I also left
  70 bases at the end of the reference. Making the
  reference 11800 bases long. I also kept the alternate
  primers (both 10's, 23 right, and 27 right).

# Bench artic (07)

I did the first benchmark on artic using my benchmarking
  script. This script has default values set to my files.
  However, all of these can be changed by cmd if needed.

**I did not run this command completely, but there is a 
  tsv file with most of the replicates in 07**

```
bash 00-scripts/benchAllArtic.sh \
   -prefix 07-artic-1st-bench/07-artic-1st-bench \
   -rep 100;

# Remove the consensuses
rm "07-artic-1st-bench/07-artic-1st-bench"*.fasta;
```

The stats will be output to
   07-artic-1st-bench/07-artic-1st-bench-stats.tsv.

I did some quick filtering to remove failures
`awk '{if($0~/Com/){getline; next;}; print $0}' stats.tsv`.

I also had to remove one really odd early entry. The output
   file was saved as
   07-artic-1st-bench/07-artic-1st-bench-stats-filter.tsv.

# Testing bulidCon (08-buildCon)

**This section is me muddling around to get things working.
  This can be skipped and is not needed (move to 09).
  This folder does not even exist in the repository, 
  because it was a dead end.**

I wanted to get the buildCon benchmarking script up and
  running. These are also some initial observations.

First I need to subsample reads

```
mkdir 08-buildCon/08-reads;
iAmp=1;

for strFq in ./05-amplicon-reads/*.fastq; do
   
   bash 00-scripts/subsampleReads.sh \
       -fastq "$strFq" \
       -depth 500 \
     > "08-buildCon/08-reads/amp-$iAmp.fastq";

   iAmp="$((iAmp + 1))";
done
```

The next step was to build the amplicon consensuses.

```
# can not be run with cond activate LILO
# For some odd reason the samtools there errors out
bash 00-scripts/buildAmpCons.sh \
   -ref "01-input/01-primer-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.reference.fasta" \
   -fastq-dir "08-buildCon/08-reads" \
   -primer-scheme "01-input/01-primer-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.scheme.bed" \
   -true-ref "01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta" \
   -prefix "08-buildCon/08" \
   -disable-medaka \
  > 08-buildCon/08-stats.tsv;

rm -r 08-buildCon/08-reads; # No longer need
mkdir 08-buildCon/08-cons;
mv 08-buildCon/*.fasta 08-buildCon/08-cons;
mkdir 08-buildCon/08-alns;
mv 08-buildCon/*aln 08-buildCon/08-alns;
```

I wanted to see how accurate the stitched amplicons were.

```
# Stich together the amplicons
# This is needed to avoid random insertions being input
# into the consensus
cat 08-buildCon/08-cons/*.fasta > tmp.fasta;
awk \
    -f 00-scripts/schemeToLiloCsv.awk \
    06-alt-schemes/06-schemes/SARS-CoV-2/V4.1/SARS-CoV-2.scheme.bed \
  > tmp.csv                        

conda activate LILO
minimap2 \
    01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta \
    tmp.fasta |
  awk \
    -v OFS="\\t" \
    '{if($4-$3>= 320){print $1, $3+10, $4-10}}' \
    - |
  bedtools sort -i - |
  bedtools getfasta -fi tmp.fasta -bed -  |
  awk '{sub(/[ATGC]>.*$/,"",$0);print $0;}' >tmp-bed.fasta;
  # The awk step to to correct and error from bedtools

porechop \
    --adapter_threshold 72 \
    --end_threshold 70 \
    --end_size 30 \
    --extra_end_trim 5 \
    --min_trim_size 3 \
    -f tmp.csv \
    -i tmp-bed.fasta \
    --threads 3 \
    --no_split \
    -o tmp-bed-porechop.fasta;

conda deactivate

# I could not get this step to work always.
# So I never got past this point
conda activate scaffold_builder                                

scaffold_builder.py \
    -i 75 \
    -t 3693 \
    -g 80000 \
    -r 01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta \
     -q tmp-bed-porechop.fasta \
     -p 08-buildCon/08-con;

conda deactivate

mv \
    08-buildCon/08-con_Scaffold.fasta \
    08-buildCon/08-con.fasta;

rm tmp.csv;
rm 08-buildCon/08-con.coords;
rm 08-buildCon/08-con.txt;
rm -r 08-buildCon/08-scaffold_overlap_alignment;
rm tmp.fasta;
rm tmp.fasta.fai;
rm tmp-bed.fasta;
rm tmp-bed-porechop.fasta;

# Get alignments
 
00-programs/alnSeq \
    -use-water \
    -print-positions \
    -print-aligned \
    -query 08-buildCon/08-con.fasta \
    -ref 01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta \
    -out 08-buildCon/08-con.aln;
```

I was trying an alternative way to make a consensus. This
  kinda worked, but used the reference to remove
  insertions.

```
#samtools faidx 01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta;

minimap2 \
    --eqx \
    -a \
    01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta \
    tmp.fasta |
  00-programs/trimSamFile -stdin |
  awk \
    '
       { # MAIN
          if($0 ~ /^@/){next;}; # sam header
          gsub(/[0-9]*[IDS]/, "", $6)
          gsub(/[=X]/, "=", $6)
          numElmI = split($6, lenI, "=");
          lenSeqI = 0;

          for(iElm = 0; iElm <= numElmI; ++iElm)
             lenSeqI = lenSeqI + lenI[iElm];

          print ">" $1, $4, lenSeqI, $10;
          # Name, first ref position, length, sequence
       } # MAIN
      ' |
  sort -V -k 2 |
  awk 
    '
       { # MAIN
          if(NR == 1)
          { # If this is the frist sequence
             startI = $2;
             lenI = $3;
             endI = $2 + $3;
             seqStr = $4;
          } # If this is the frist sequence

          if(endI > $2)
          { # If I have an overlap
             trimI = int(endI - $2);
             seqStr=substr(seqStr,0,length(seqStr-trimI));
             $4 = substr($4, trimI);

             lenI = lenI + $3;
             endI = $2 + $3;
          } # If I have an overlap

          else
          { # Else I do not have an overlap
             maskI = startI - endI;

             for(iN = 0; iN <= maskI; ++iN)
                seqStr = seqStr "N";

             seqStr = seqStr "$4";

             lenI = lenI + $3;
             endI = $2 + $3;
             maskLenI = maskLenI + maskI;
          } # Else I do not have an overlap
       }; # MAIN

       END{
           printf ">con start=%s len=", startI, lenI;
           printf " end= noNs=\n", endI, maskLenI;
           print seqStr;
      }; # END block
    '

lastAmpI="$(\
  awk '
        { # MAIN
           sub(/.*amp-/, "", $0);
           sub(/-.*/, "", $0); 
           if($0 > lastAmpI) lastAmpI = $0;
           getline;
        }; # MAIN
        END{print lastAmpI};
      ' < tmp.fasta \
)";

awk \
    -v lastAmpI="$lastAmpI" \
    '
    { # MAIN
       # Frist amplicon, no triming
       if($0 ~ /amp-1-/)
       { # If on the first amplicon (only trimming end)
          print $0;
          getline;
          seqStr = substr($0, 50, length($0) - 50)
          print $0
       }; # If on the first amplicon (only trimming end)

       ampOnI = $0;
       sub(/.*amp-/, "", ampOnI);
       sub(/-.*/, "", ampOnI); 
       if(ampOnI == lastAmpI)
       print $0;
       getline;
       seqStr = substr($0, 50, length($0) - 50)
       print seqStr;
    }' < tmp.fasta \
  > tmp-trim.fasta;

minimap2 \
   -a 01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta \
    tmp-trim.fasta |
  samtools sort |
  bcftools mpileup \
    -Q 0 \
    -q 0 \
    -f 01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta \
    - \
  > tmp.vcf &&
  bgzip tmp.vcf;

bcftools index tmp.vcf.gz
bcftools consensus -f 01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta tmp.vcf.gz 
```

# Benchmark (09)

At this point I am running my benchmarking command. I have
  set it up so that the defaults all point to my files.

```
bash 00-scripts/benchAll.sh \
    -prefix 09-benchmark/09-benchAll\
    -rep 5;
```

# Redownload data and rerun (JUMP-RETRY)

On Sep, 15th 2023 I had to redownloaad the data. This time
  I did it correctly and was able to get all the data.

`prefetch SRR21813677 && fasterq-dump SRR21813677`

I removed the prefetch download to save some space.

`rm -r SRR21813677`

`mv SRR21813677.fastq 01-input/01-SRR21813677.fastq`;

```
00-programs/fqGetIds \
    -f 01-input/01-SRR21813677-ids.txt \
    -fastq 01-input/01-SRR21813677.fastq \
  > tmp.fq;

minimap2 \
    --eqx \
    -a \
    01-input/01-primer-schemes/SARS-CoV-2/V5.3.2/SARS-CoV-2.reference.fasta \
    tmp.fq \
  > tmp.sam;
awk \
    -f 00-scripts/extReadsByTbl.awk \
    -v extTbl=04-read-table/04-buildConTbl.md \
    -v prefix=05-amplicon-reads/05-500bp-reads \
  < tmp.sam 

rm tmp.fq;
rm tmp.sam;

# Rerunning the benchmark
bash 00-scripts/benchAll.sh \
   -prefix 09-benchmark/09-bench-all-2 \
   -rep 5;
```

# Rerunning each test individually (10)

Running artic and ivar together; both are quick & working.

```
bash 00-scripts/benchAll.sh \
    -no-buildCon-majcon-medaka \
    -no-buildCon-medaka \
    -artic \
    -buildCon-majcon \
    -LILO \
    -ivar \
    -ivar-trim \
    -rep 10 \
    -prefix 10-new-bench/10-bench-artic-ivar;
```
