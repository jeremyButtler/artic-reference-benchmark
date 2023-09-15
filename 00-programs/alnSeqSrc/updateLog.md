# Use

This log records how alnSeq has changed between versions.
  Each version is labeled with the version number, year,
  month, and then day. For example version 1.20230615 would
  be version 1, from June 15th, 2023.

# TODO or to fix list

- Set up Smith Waterman Hirschberg combination
  - At this point I am close, since the 20230827 updates
    made my alternative waterman very close to what this
    waterman would be.
- Add in a match matrix. Currently I am using a switch
  table set up for nucleotides to identify matches/snps.
  - This is used for printing the alignment eqx line, not
    for finding an optimal alignment (scoring matrix).
- Need to add a separate file for printing out alternative
  alignments. It currently prints everything in one file.

# Ideas that would be cool, but not worth working on

These ideas are a future vision that is not worth the
  time unless this algorithm becomes useful to others.

- Multithreading & gpu support with openCL
- Get a matrix scan function that recalculates scores on
  fly, so that user can search a completed direction
  matrix.
  - Also allow a function to be input, so user can decide
    if printing or not.
- Add a filter for alternative alignments.

# Log

## 20230827

- Alternative alignment printing now prints out the
  score, starting bases of alignments, and ending bases of
  alignments (Used to be cigars).
  - Still need to add a filter.
- The query-reference scan has been simplified.
- Matrix scan has been removed. It was not that great in
  the first place.
- Issue with Waterman not printing out full length
  alignments has been fixed.
  - It was do to using a `&&` when I should have been
    using `||` in the while loop in printAln
    (alnStruct.c Fun-03 Sec-05 Sub-01).
    ```
    *refAlnStr != defEndAlnFlag ||
    *qryAlnStr != defEndAlnFlag
    ```
- Fixed an issue in -format-emboss, were it would print
  out the wrong alignment length. It now prints out the
  matches + insertions + deletions instead of the longest
  sequence.
- Adding in several compiler time flags
  - -DNOGAPOPEN: Disables gap opening penalties
    - Default gap penatly is set to -6 for this option
    - Removed -gapopen option.
  - -DBYTEMATRIX: Use a character matrix instead of an
    two bit matrix for the Needleman and Waterman
    alignments.
  - Flags to make alnSeq use only one direction
    - Removes the direction selection options.
    - -DSNPINSDEL
    - -DSNPDELINS
    - -DINSSNPDEL
    - -DINSDELSNP
    - -DDELSNPINS
    - -DDELINSSNP
- Added in a `-flag` and `-flag-only` options to print out
  the flags alnSeq was complied with. `-flag` will also
  print out a description of each flag.
- Shortened two bit array function names
  - TwoBitArry -> TwoBit
    - `sed 's/TwoBitArr*y/TwoBit/g'`
  - TwoBitAry  -> TwoBit
    - `sed 's/TwoBitArr*y/TwoBit/g'`
  - twoBitAryMove -> twoBitAryMv
    - `sed 's/twoBitAryMove/twoBitAryMv'`
  - getTwoBitAryIndes -> twoBitGetIndex
- Added to twoBit to two bit function names missing two bit
  - moveToNextLimb    -> mvToNextTwoBitLimb
  - moveToLastLimb    -> twoBitMvToLastLimb
  - moveXElmFromStart -> twoBitMvXElmFromStart
  - lenTwoBit         -> twoBitGetLen
- Updated using this code guide for the new printing and
  alnStruct changes.

## Version 20230811

- Switched from converting my bases to a look up table on
  checks to the read in sequence.
  - (disable with `make CFLAGS="-DNOSEQCNVT"`)
- Improved speed by inlining several functions and adding
  some branchless operations at key points
  (selecting max score, and zeroing waterman values). Other
  places not worth it.
  - All functions in twoBitArrays.c have been inlined and
    are now in twoBitArrays.h.
  - Scoring functions are macros in generalAlnFun.h
- Hirschberg is now compiled with 1 byte arrays by default
  - This is faster, but use a small amount of extra memory.
  - enable twobit version with
    `make CFLAGS=`-DHIRSCHTWOBIT`.
- Made Hirschberg thread safe. This will increased memory
  usage by a very slight amount
  - ~ bytes = 1/4 \* reference length for two bit.
  - ~ bytes = \* reference length for byte arrays
  - This was added in at this point, but I forgot to 
    mention the updated. That is why it is added in at
    a later date.


## Version 20230804

- Fixed base getting removed on printing. It turned out to
  be error in my fasta reader (assigning an index 0
  variable an index 1 value). So, it was not a printing
  error.
  - line 510 in seqStruct.c `spareBuffUL = resBuffUL;` was
    changed to `spareBuffUL = resBuffUL - 1;`.

## Version 20230803

- The Hirschberg now prints using the normal functions.
  - I corrected an error were I was sending in the query as
    the reference to the Hirschberg.
- alnStruct has been changed. It now stores two alignment
  arrays. One for the reference and one for the query.
  - Each array tells if a specific reference base or query
    base mapped to a gap, snp, match, or was soft masked.
  - Also added in are insertion, deletion, snp, and match
    counts.
  - Also includes the first match/snp and last base in the
    alignment
  - The alignment arrays now include matches.
- The Smith Waterman is now printing correctly when a line
  wrap is applied.
- Changed printAln (alnStruct.c/h, function 03, used to be
  04) to take in more input. This allowed me to print
  out more to the header.
  - You can now print the output file in fasta, clustal,
    or an emboss like format. The old expanded cigar
    format is still the default option.
  - I am hoping the change in printing method got rid of a
    rare error I found a few days ago, were one or two
    bases would be replaced with an deletion. I have not
    confirmed if my hope is correct.
- Added in an option `-print-alinged` to only print out the
  aligned bases in the alignment.
- Added in option `-print-position` to include starting and
  ending base numbers at each line.
  - Always set when output is in the EMBOSS like format.
- Updated the help message. I also made a printHelpMesg
  function, which adds in the default values in 
  alnSeqDefualts.h/c to the help message and prints it out.
  This should avoid the help message getting out of date
  in the future.

## Version 20230726

- Hirschberg is working, but the output methods are
  limited.
- Gap opening penalty has been changed to -10. I found this
  worked better when using alnSeq for my projects.

## Version 20230709

- using "-line-wrap 0" will now disable line wrapping
- Worked out some bugs in the Hirschberg alignment
  - Hirschberg still not working

## Version 1.20230630

- Changed sequence lengths and offsets to be longs instead
  of integers
- Put down logic for Hirschberg, which I am currently
  debugging
- Changed the default gap extension and gap opening 
  penalties to be more suited for non-noisy reads
  - Gap opening penalty is now -4 (used to be -1)
  - Gap extension penalty is now -1 (used to be -4)
- Changed
  
## Version 1.20230620

- Fixed an issue were two very different large alignments
  (Tested was: testing/largeTest-query-a.fasta and
  TickBornEncephalitis-reference.fasta) would segfault.
  - This error was due to my inputting unsigned longs
    instead int32_t for twoBitAryMoveForXElm.
  - twoBitAryMoveForXElm now uses unsigned longs instead of
    int32_t for its index

## Version 1.20230619

Other than bug fixes that pop up, this will be my final
  code update until I see that this code is useful to
  others. The only exception is if I find a use for 
  pressing this code further.

- Fixed the printing to stdout on Debain error.
  - This was due to alnST not being allocated properly
    in the cnvtDirMatrixToAlnAry function with calloc. This
    resulted in only a partial alnStruct structure being
    allocated.
  - Also had problem were I freed the alignment matrix to
    soon and so lost the best score. This has been fixed
    by storing the best score in a temporary variable.
- Multi sequence output that prints a best score for each
  query and reference for Waterman Smith alignment no
  longer segfaults
  - Was problem with using double pointers to make an array
    of structures.
  - Setting has been changed to "-query-ref-scan-water".
- Fixed matrix scan for Waterman Smith having reference
  starting base is going out of bounds of unsigned integer.
  - This was caused by my Loop always going negative for a
    counter, so I switched to finding starting base by
    using the ending index (end of the traced alignment
    path or start of alignment) of the two-bit Array to
    find the starting base.
- The reported base range for the alignment is now in index
  0 (first base is position 0) instead of index 1 (first
  base is position 1).
- The min score for matrix scanning and query/reference
  scanning (best score for each reference/query base) has
  been changed to 1000 (at least 200 bases).
  - Score increased for each non-anonymous match is 5.
- Set U and T to be separate scores and updated the default
  EDNA scoring matrix to have U and T be the same.
  - This was set because I built alnSeq to align DNA and
    RNA sequences, which have T and U being the same. This
    would have had no effect on protein alignments, which
    do not use U, but would have had an effect on word
    alignments.
  - AlnSeq ignores case for all alignments.
  - The scoring matrix alnSeq can be set up to have scores
    for the characters a-z.
    
## Version 1.20230616

- The Smith Waterman multi base alignment logic has been
  modified to only record a best score if the alignment
  does not extend to the next row or the current score
  is better then the next score in the alignments path.
  - Currently this segfaults
- Added in a matrix scan feature.
  - This currently is off on the reference starting base

## Version 1.20230615

- Smith Waterman single answer alignment is now working
- Fixed some issues with printing out alignments.
  - Converting matrix index's to base positions was off
  - Alignment printing function was printing extra bases
    or two few bases

## Pre Version 1.20230615

- Logic set up and Needleman Wunsch alignment working
