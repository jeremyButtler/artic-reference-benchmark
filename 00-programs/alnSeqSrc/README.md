# Use
AlnSeq uses a Smith Waterman and Needleman Wunsch alignment
  that, depending on flags used when compiling can run
  with memory usage of O(n \* m / 4) Bytes to O(n \* m)
  Bytes. However, the O(n \* m / 4) is twice as slow as a
  traditional Smith Waterman and Needleman Wunsch
  alignment. AlnSeq also includes an Hirschberg alignment.

There are faster and less memory hungry Waterman Smith
  implementations than alnSeq. One example is the stripped
  Waterman Smith alignment, which I think reduces both
  scoring and direction matrix to just a few rows. See 
  [https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
  for an very fast striped Waterman Smith aligner. This
  one does have an issue when the alignment exceeds 33,000
  bases. However, I think there are alternatives without
  these issues.

This program is dual licensed for MIT and CC0. Pick the
  license which works best for you.

# Current work

I know I said I was finished, but I just thought I would
  try reducing the Waterman directional matrix down to two
  rows and just record the score, start position, and end
  position of the best alignment. My plan would then be to
  use an Hirschberg to find the actual alignment. My logic
  is that a local aligmnent turns into a global alignment
  once you know the start and end. This would reduce memory
  usage down, but would also require running both a
  Waterman and Hirschberg alignment, which takes more time.

I then could also allow recording of the best score for
  each query and reference base or printing of each score
  that passes the min score value to add in additional
  searching. However, I will only print out the score, 
  starting reference base, starting query base, ending
  reference base, and ending query base for alternative
  alignments.

# Building and running alnSeq

## How to build alnSeq

```
# How to install this program

# Install at /usr/local/bin (need root privilage)
make
sudo make install

# Use a different install path
make
make install PREFIX=/path/to/install

# Manual install
make
mv alnSeq /path/to/install
chmod a+x /path/to/install/alnSeq

# Alternative make commands
make fast
  # This command disables the opening gap penalty,
  # directional selection, and two bit arrys. This makes
  # alnSeq behave somewhat similar to bio-alignment.
make mid
  # This command disables directional selection and two
  # bit arrays, but keeps the gap opening penalty.
```

## Extra build options

The flags alnSeq can be compiled with are:

- -DNOGAPOPEN
  - disable gap opening penalty
  - This is faster, but may produce lower quality
    alignments.
- -DBYTEMATRIX
  - Waterman and Needleman use a byte matrix instead of
    a two bit array matrix
  - Alignment takes half the time, but also takes 4x more
    memory.
- -DHIRSCHTWOBIT
  - Have the Hirschberg use two bit arrays instead of byte
    arrays for directions (only relevant if not using
    -DNOGAPOPEN).
  - This doubles the time to make an alignment for a very
    minor decrease in memory usage. 
- -DNOSEQCNVT
  - This prevents the conversion of each base in the query
    and reference sequences to an index for alignment.
  - This option slows down the alignment slightly. The only
    reason to use this option would be if the input case
    of a sequence matters.
- These options force alnSeq to prefer only one direction
  and disables all other options. This does speed up alnSeq
  slightly.
  - It only applies when the insertion (ins), snp, or 
    deletion (del) scores for a base pair are equal.
  - -DSNPINSDEL
  - -DSNPDELINS
  - -DINSSNPDEL
  - -DINSDELSNP
  - -DDELSNPINS
  - -DDELINSSNP

You can compile with these flags using
  `make CFLAGS="flag"`. You can also compile multiple 
  flags with `make CFLAGS="falg1 flag2"`. One example is
  the `make fast` command, which uses make
  `CFLAGS="-DNOGAPOPEN -DINSSNPDEL -DBYTEMATRIX"`.

## How to run alnSeq

```
# help message
alnSeq -h | less

# Get flags set when compiling alnSeq
alnSeq -flags | less

# Alignment formats 

## For a global alignment (Needleman Wunsch)
alnSeq -query query.fasta -ref ref.fasta > alignment.aln

## For a Hirschberg global alignment (slow, but low memory
# usage)
alnSeq -use-hirschberg -query query.fasta -ref ref.fasta > alignment.aln

## For a single local alignment (Waterman Smith)
alnSeq -use-water -query query.fasta -ref ref.fasta > out.aln

# File formatting

## Output an EMBOSS like file
alnSeq -format-emboss -query query.fasta -ref ref.fasta -out out.aln

## Output a clustal file
alnSeq -format-clustal -query query.fasta -ref ref.fasta -out out.aln

## Output a fasta file
alnSeq -format-fasta -query query.fasta -ref ref.fasta -out out.aln

## Trim sequences
alnSeq -print-aligned -query query.fasta -ref ref.fasta -out out.aln

## Print positions for non EMBOSS files
alnSeq -print-positions -query query.fasta -ref ref.fasta -out out.aln
```

# Explaining alnSeq

## What is alnSeq?

AlnSeq is an sequence alignment program that uses either 
  a Needleman Wunsch, Hirschberg, or Smith Waterman
  alignment algorithm. AlnSeq is unique from a traditional
  Needleman Wunsch or Smith Waterman in how it handles the
  scoring matrix and the direction matrix.

One thing I do want to point out is that there are very
  memory efficient algorithms, such as the striped Waterman
  Smith alignment, for optimal local alignments.

## The direction matrix

AlnSeq stores each direction in two bits for the Needleman
  and Waterman alignments. These bits are packed
  into an 8 bit integer. This reduces the directional
  matrix size by 4, but comes at the cost of slower speed.
  This also means that only one direction is stored,
  instead of all possible alternatives.

## The scoring matrix (Needleman/Waterman)

AlnSeq also reduces the scoring matrix down to one row,
  which holds the last scores or the previously updated
  scores. This removes the scoring matrix, but also removes
  the ability to scan the scoring matrix for other high
  scores, which is sometimes done for an Waterman Smith
  alignment. The best score is found while building the
  matrix.

AlnSeq also supports alternative alignments with
  -query-ref-scan-water by storing the best score for each
  reference base and each query base (Starting positions,
  ending positions, score). The score, starting reference
  position, starting query position, ending reference
  position and ending query position are then printed
  out (currently in the same file before the alignment).
  There is no filter, so this will print out everything
  that is at or above -min-score. This includes duplicate
  scores.

## Some light benchmarking

I picked out three programs to compare alnSeq to. The first
  is emboss, which is a more commonly used toolkit. The
  second is bio-alignment (bio or ssw), which had a
  Hirschberg. The last was the
  Complete-Striped-Smith-Waterman-Library (bio or ssw),
  which supports a vectorized, memory efficient, striped
  Smith Wateman alignment. This is not an exhaustive list,
  nor does it include the best programs. It is just a
  simple and quick list.

Also my figures show two different runs of alnSeq. One
  run is the default compile, with the Hirschberg using
  byte arrays and the Waterman and Needleman using two
  bit arrays.
The other run (alnSeqFast) is `make fast` and has the
  gap opening penalty disabled, uses byte arrays, and 
  is hardcoded to prefer insertions, then snps, then
  deletions when scores are equal.
This is the closest you can get to have a similar
  comparison to bio-alignment.

![Memory usage of alnSeq compared to the Waterman Smith,
  Needleman Wunsch, and Hirschberg aligners from
  bio-alignment (bio-align), emobss, and
  the Complete-Striped-Smith-Waterman-Library.
  alnSeq-scan is printing alternative bases.
](
  analysis/20230827-benchmarking/20230827-alnSeq-memory.svg
)

As expected the memory usage was much lower for the
  Hirschberg aligners and striped Smith Waterman
  (local facet; bio_or_ssw), while the non-Hirschberg
  aligners for alnSeq, emboss, and bio-alginment needed
  large amounts of memory.
For the non-Hirschberg's and striped smith waterman
  alignments, alnSeq needed less memory than emboss or
  bio-alignment.
With alnSeq using two bit arrays taking the least amount
  of memory.
When compared to bio-alignments Hirschberg, alnSeq's
  Hirschberg used less memory than bio-alignments
  Hirschberg, however, the memory usage for both
  Hirschbers is small and so the memory saving has no real
  impact.

Bio-alignment's Hirschberg used more memory when two
  different alignments (small-huge, mid-huge, large-huge,
  small-large) were aligned.
From a glance at the code I suspect this was due to 
  bio-alignment allocating memory for its returned scoring
  arrays at each recursion call.
For highly different alignments this would result in
  the midpoint being closer to 1, which would result in
  an nearly identical returned scoring row size in the next
  recursion call.
Since, bio-alignment is not freeing these arrays right
  away, it is possible that these arrays would continue to
  build up at each call, which results in increased memory
  usage.
However, a local alignment should be used instead of an
  global alignment (Needleman/Hirschberg) in these cases.

We found that ssw had lower memory usage when aligning
  similar genomes and a higher memory usage when the
  genomes were very different.
I am not sure why this is happening, but it may be due to
  how much of the matrix it has to construct.
Despite this, ssw still uses less memory than alnSeq's
  Waterman alignment and so, is a better option.

![Time usage of alnSeq compared to the Waterman Smith,
    Needleman Wunsch, and Hirschberg aligners from
    , emobss, and
  Complete-Striped-Smith-Waterman-Library.
](analysis/20230827-benchmarking/20230827-alnSeq-time.svg)

For time usage we found that emboss, alnSeq scan, or
  alnSeqs two bit Hirschberg (not shown) were the slowest
  algorithms.
The fastest programs were ssw (30x faster then alnSeq and
  10x faster than alnSeqFast Waterman), followed by the
  Needleman and Hirschberg from alnSeq-fast, which uses a
  byte matrix and like bio-alignment, ignores gap opening
  penalties.
Some of the extra time needed for Emboss could be due to
  its calculating both a gap opening and gap ending
  penalty.

## Final notes

AlnSeq does not use decimals, so if you want decimals for
  the gap extension penalty, so you will have to multiply
  all scores by 10 to 1000.

# Thanks

- To my dad for being willing to listen and give advice.
- To the people (never met) who coded baba
  [https://baba.sourceforge.net/](https://baba.sourceforge.net/).
  It gave me a great visual on how the Needleman Wunsch
  algorithm worked. Their were many other sources, but this
  was the one that was the most useful to me.
- Wikipedia's entry about the Hirschberg. It is helping me
  understand how the Hirschberg works
- Bio-alignment coded a Hirsberg, which I used to help
  understand how the Hirsberg alignment worked. 
- To all the sights providing guidance on how aligners
  worked. There are to many to list, but each one I used
  helped me understand these algorithms a little better.
