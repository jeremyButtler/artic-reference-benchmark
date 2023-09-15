/*#########################################################
# Name: trimSam
# Use:
#  - trims soft mask regions off all alignments with
#    sequences in a sam file. Aligments without sequences
#    are ignored & not printed out.
# Libraries:
#   - "samEntryStruct.h"
#   o "cStrToNumberFun.h"
# C Standard Libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>
#########################################################*/

#ifndef TRIMSAM_H
#define TRIMSAM_H

#include "samEntryStruct.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
| trimSam SOH: Start Of Header
|  - fun-01 trimSamReads:
|    o Trims soft mask regions for all reads with a
|      sequence in a sam file
|  - fun-02 trimSamEntry:
|    o Trim soft mask regions off end of sam entry
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: trimSamReads (Fun-01:)
| Use:
|  - Goes though sam file and calls trimSamEntry for each
|    entry
| Input:
|  - samFILE:
|    o Sam file with entries to trim
|  - outFILE:
|    o File to write the trimmed sam file entries to
|  - keepUnmappedReadsBl:
|    o Also print out entries with unmapped reads
| Output:
|  - Prints:
|    o Trimmed sam entries with sequences to outFILE, but
|      ignores sam entries without sequences
\--------------------------------------------------------*/
void trimSamReads(
    FILE *samFILE,               /*Sam file to trim*/
    FILE *outFILE,               /*File to store output*/
    char keepUnmappedReadsBl     /*1: keep unmapped reads*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: trimSamReads
   '  - Goes though sam file and calls trimSamEntry for
   '    each entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: trimSamEntry (Fun-02:)
| Use:
|  - Trims off the soft masked regions of a sam entry
| Input:
|  - samST:
|    o samEntry struct with sam entry to trim
| Output:
|  - Returns:
|    o 0 if suceeded
|    o 2 if header (invalid and ignored)
|    o 4 if an unmapped read (no reference)
|    o 8 if no sequence line
|  - Modifies:
|    o Trims cigar, sequence, & q-score entries in samST.
\--------------------------------------------------------*/
uint8_t trimSamEntry(
    struct samEntry *samST   /*has sam line to trim*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: trimSamEntry
   '  - Trims soft masked regions at start & end of sam
   '    entry
   '  o fun-02 sec-01:
   '    - Variable declerations
   '  o fun-02 sec-02:
   '    - Find how much to trim & trim cigar entry
   '  o fun-02 sec-03:
   '    - Trim the sequence entry
   '  o fun-02 sec-04:
   '    - Trim the q-score entry
   '  o fun-02 sec-05:
   '    - Shift characters for other parts of the cigar
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
