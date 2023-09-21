/*#########################################################
# Name: stichMinimapFun
# Use:
#  - Holds functions used in stiching amplicons into a
#    consensus. These functions use minimap2
# Libraries:
#  - "sitchAmpStruct.h"              (No .c file)
#  - "samFunSrc/trimSam.h"
#  - "samFunSrc/cStrFun.h"           (No .c file)
#  - "stichSetStruct.h"              (No .c file)
#  o "stichDefaults.h"               (No .c file)
#  o "samFunSrc/samEntryStruct.h"    (No .c file)
#  o "samFunSrc/cStrToNumberFun.h"   (No .c file)
#  o "samFunSrc/dataTypeShortHand.h" (No .c file)
#  o "samFunSrc/seqStruct.h"
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdint.h>
#  o <stdio.h>
#  o <string.h>
# Requires:
#  - Minimap2 be in your file path
#########################################################*/

#ifndef STICHMINIMAPFUN_H
#define STICHMINIMAPFUN_H

#include "samFunSrc/trimSam.h"
#include "samFunSrc/dataTypeShortHand.h"
#include "samFunSrc/seqStruct.h"
#include "samFunSrc/cStrFun.h"
#include "stichDefaults.h"
#include "stichAmpStruct.h"
#include "bitTwidle.h"
#include "stichSetStruct.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' stichMinimapFun SOH: Start Of Header
'  - Holds the functions used to stich amplicons together.
'  - This variation uses minimap2 for alignment
'  o fun-01 getAmpPosMinimap:
'    - Gets the position of each amplicon on the reference
'  o fun-02 stichAmpConMinimap:
'    - Stiches together amplicon consensuses to to make
'      an consensus genome
'  o fun-03 stichAmpsMinimap:
'    - Stiches amplions into a consensus
'  o fun-04 samEntryToAlnSeq:
'    - This uses a samEntry struct to get an aligned
'  o fun-05 freeSamEntryAry:
'    - Frees an array of samEntry structures
'  o fun-06 stichAmpConToCStr:
'    - Convert a stichAmpST consensus to a merged, c-string
'      consensus
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: getAmpPosMinmap (Fun-01:)
| Use:
|   - Find the starting positions and ending postions of
|     each amplicon on the reference
| Input:
|  - refST:
|    o seqStruct with reference sequence
|  - qryST:
|    o seqStruct with query sequence
|  - ampFaFILE:
|    o Fasta file handle with amplicons to find positions
|  - numAmpsUL:
|    o This will hold the number of amplicons in ampFaFILE
|  - settings:
|    o stichSet struct with the number of threads to use
|  - outErrUC:
|    o Holds any errors that happened (see output)
| Output:
|   - Modifies:
|     o numAmpsUL to hold the number of amplicons
|     o indexAryUL to hold the index of every sequence.
|     o outErrUC to hold the error, if an error occured
|       - 0: No error
|       - 1: amplicon fasta file error
|       - 2: reference fasta file error
|       - 3: reference fasta has multiple sequences
|       - 4: Could not allocate memory
|       - 5: minimap2 errored out or memory error
|       - 6: No amplicons were kept (none mapped to ref)
|   - Returns:
|     o Pointer to array of samEntry structures
|     o 0 for memory errors
\--------------------------------------------------------*/
struct samEntry * getAmpPosMinimap(
   char *refFileStr,  /*Name of fasta file with reference*/
   char *ampFileStr,  /*Name of fasta file with amplicons*/
   ulong *numAmpsUL,  /*Will have number amplicons kept*/
   struct stichSet *settings, /*Number of threads to use*/
   uchar *outErrUC    /*Reports the error type*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: getAmpPosMinimap
   '  - Gets the position of each amplicon on the reference
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - File checks and memory allocation
   '  o fun-01 Sec-03:
   '     - Align each sequence to the reference
   '  o fun-01 sec-04:
   '    - Reset file pointer and sort the sequences by
   '      reference starting positions
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: stichAmpConMinimap (Fun-02:)
| Use:
|  - Stiches together amplicon consensuses to to make
|    an stichAmpST scaffold with alternative bases
| Input:
|  - ampsAryST:
|    o array of samEntry structs wich contain the sam file
|      entry for all amplicons to stich together
|  - numAmpsUL:
|    o Number of amplicons in ampsAryST
|  - stichSetST:
|    o Settings for stiching amplicons together.
| Output:
|  - Returns:
|    o List of stichAmpST structs that have the
|      uncollapsed (with alternative bases) scaffold
|    o 0 for error
| Note:
|  - This function works with output from getAmpPosMinimap
\--------------------------------------------------------*/
struct stichAmpST * stichAmpConMinimap(
   struct samEntry *ampsAryST,  /*Array of amplicons*/
   ulong numAmpsUL,             /*Number amplicons*/
   struct stichSet *stichSetST  /*Settings for stich*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: stichAmpConMinimap
   '  - Stiches together amplicon consensuses to to make
   '    an consensus genome
   '  o fun-02 sec-01:
   '     - Variable declerations
   '  o fun-02 sec-02:
   '    - Align each amplicon sequence
   '  o fun-02 sec-03:
   '    - Stich together each amplicon sequence
   '  o fun-02 sec-04:
   '    - Deal with disagreements in the amplicons
   '  o fun-02 sec-05:
   '    - Clean up and exit
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: stichAmpsMinimap (Fun-03:)
| Use:
|  - Stiches amplions into a consensus
| Input:
|  - alnSeqStr:
|    o The aligned amplicon sequence to stich into the
|      consensus
|  - ampScoreST:
|    o scoresStruct with the first reference base alnSeqStr
|      starts on (ampScoreST->refStartUL).
|    o ampScoreST->refStartUL should be in index 0
|  - conST:
|    o stichAmpST on the last added base in the consensus.
|      If 0, this will start a new consensus.
|  - refEndUL:
|    o The position of the last added base in conST
|    o This should be index 0;
|  - maskC:
|    o Character to mask amplicon with
| Output:
|  - Modifies:
|    o refEndUL to have to the position last added base
|      added by stichAmpST.
|  - Returns:
|    o A pointer to the last added base in conST.
|      - conST is a double linked list, so you can get to
|        the first base from the last base or the last base
|        from the first base.
|    o 0 if had a memory error (the list in conST is
|      freeded)
|  - Frees:
|    o the list in conST if had a memory error
| Note:
|  - stichAmpsMinimap assumes that the amplicons have been sorted
|    by starting position on reference. With the amplicon
|    mapping to the first reference base coming first.
|    o This can be done with sortScoresStartLen() or
|      sortScoresStartLenIndex() in alnSeqSrc/scoresST.h.
\--------------------------------------------------------*/
struct stichAmpST * stichAmpsMinimap(
   char *alnSeqStr,          /*Amplicon sequence*/
   ulong ampStartUL,        /*First ref base in amplicon*/
   struct stichAmpST *conST, /*Consensus sequence (list)*/
   ulong *refEndUL,          /*Last ref base in consensus*/
   char maskC                /*What to use for maksing*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: stichAmpsMinimap
   '  - Stiches amplions into a consensus
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Check if this is the first amplicon or have a gap
   '      between the consensus and next amplicon 
   '  o fun-03 sec-03:
   '    - Add bases to the overlap
   '  o fun-03 sec-04:
   '    - Add bases to the overlap
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: samEntryToAlnSeq (Fun-04:)
| Use:
|  - This uses samEntry struct to get an aligned sequence.
|  - This function assumes that the samEntry has a sequence
|    and cigar entry. This means the sequence is mapped and
|    is not a secondary or supplemental alignment
| Input:
|  - ampST:
|    o samEntry struct with sequence to convert
|  - errTypeC:
|    o Reports the type of error
| Output:
|  - Returns
|    o A c-string with the aligned query sequence.
|      Insertions are in lower case, with deletions as '-'.
|    o 0 for an error (memory, secondary alignment,
|  - Modifies:
|    o errTypeC to be 0 for no errors
|    o errTypeC to be 1 for invalid sam entry
|      (unmapped, secondary, supplemental, no cigar)
|    o errTypeC to be 2 for a memeory error
\--------------------------------------------------------*/
char * samEntryToAlnSeq(
  struct samEntry *ampST,
  char *errTypeC          /*Holds the error type*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: samEntryToAlnSeq
   '  - This uses a samEntry struct to get an aligned
   '    sequence.
   '  - This function assumes that the samEntry has a
   '    sequence and cigar entry. This means that the
   '    sequence is mapped and is not a secondary or
   '    supplemental alignment
   '  o fun-04 sec-01:
   '     - Variable declerations
   '  o fun-04 sec-02:
   '     - Allocate memory
   '  o fun-04 sec-03:
   '     - Find the first aligned reference base
   '  o fun-04 sec-04:
   '     - Find the first aligned query base
   '  o fun-04 sec-05:
   '     - Build the alignment
   '  o fun-04 sec-06:
   '     - Clean up and return alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: freeSamEntryAry (Fun-06:)
| Use:
|  - Frees an array of samEntry structures
| Input:
|  - samAryST:
|    o Pointer to samEntry array to free
|  - numAmpsUL:
|    o Number of amplicons in the array
| Output:
|  - Frees:
|    - samAryST and its internal heap variables
|  - Modifies:
|    - samAryST to point to 0
\--------------------------------------------------------*/
void freeSamEntryAry(
   struct samEntry **samAryST,
   ulong numAmpsUL
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: freeSamEntryAry
   '  - Frees an array of samEntry structures
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: stichAmpConToCStr (Fun-06:)
| Use:
|  - Converts a list stichAmpST structs (consensus with
|    alternative bases) to a c-string with no alternative
|    bases.
| Input:
|  - conST:
|    o stichAmpST list to merge into a single consensus
|  - settings:
|    o stichSet struct with settings to use to when
|      merging the consensus
| Output:
|  - Returns:
|    o C-string with the merged consensus
|    o 0 for memory errors
\--------------------------------------------------------*/
char * stichAmpConToCStr(
   struct stichAmpST *conST, /*Has consensus to make*/
   struct stichSet *settings  /*Settings for stich*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: stichAmpConToCStr
   '  - Convert a stichAmpST consensus to a merged, c-string
   '    consensus
   '  o fun-06 sec-01:
   '    - Set up for colapsing the consensus list
   '  o fun-06 sec-02:
   '    - Find the possible largest consensus length
   '  o fun-06 sec-03:
   '    - Allocate memory for the consensus
   '  o fun-06 sec-04:
   '    - Merge alternative bases into a single consensus
   '  o fun-06 sec-04:
   '    - Clean up and exit
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
