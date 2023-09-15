/*########################################################
# Name memWater
# Use:
#  o Holds functions doing a memory efficent Smith Waterman
#    pairwise alignments. These aligners return positions,
#    which can be then run through a global/local aligner
#    used to find the actual alignment.
# Includes:
#   - "waterman.h"
#   o "generalAlnFun.h"
#   o "alnStruct.h"
#   o "alnMatrixStruct.h"
#   o "twoBitArrays.h"
#   o "scoresST.h"
#   o "seqStruct.h"
#   o "alnSetStruct.h"
#   o "alnSeqDefaults.h"
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>
#   0 <string.h>
########################################################*/

#ifndef MEMWATER_H
#define MEMWATER_H

#include "waterman.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' memWater SOH: Start Of Header
' o fun-01 memWaterAln:
'   - Run a memory efficent Waterman Smith alignment on
'     input sequences
' o fun-02 mamWaterAltAln:
'   - Run a memory efficent Waterman Smith alignment that
'     returns alternative alignmetns
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: memWaterAln
| Call: memWaterAln(qryST, refST, settings);
| Use:
|   - Performs a memory efficent Smith Waterman alignment
|     on a pair of sequences
| Input;
|   - qryST:
|     o SeqStruct with the query sequence and index 0
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment.
|   - refST:
|     o SeqStruct with the reference sequence and index 0
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment
|   - settings:
|     o alnSet structure with the setttings, such as
|       gap open, gap extend, scoring matrix, and preffered
|       direction.
| Output:
|  - Returns:
|    o scoresStruct with the best score
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
struct scoresStruct * memWaterAln(
    struct seqStruct *qryST, /*query sequence*/
    struct seqStruct *refST, /*ref sequence*/
    struct alnSet *settings  /*Settings for alignment*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: WatermanAltAln
   '  - Run a Waterman Smith alignment on input sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in initial negatives for ref
   '  o fun-01 sec-04:
   '    - Fill the matrix with scores
   '  o fun-01 sec-05:
   '    - Set up for returing the matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: memWaterAltAln
| Call: memWaterAltAln(qryST, refST, settings);
| Use:
|   - Performs a memory efficent Smith Waterman alignment
|     on a pair of sequences
| Input;
|   - qryST:
|     o SeqStruct with the query sequence and index 0
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment
|   - refST:
|     o SeqStruct with the reference sequence and
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment
|   - settings:
|     o alnSet structure with the setttings, such as
|       gap open, gap extend, min score, scoring matrix,
|       and preffered direction.
| Output:
|  - Returns:
|    o alnMatrixStruct with the best score and alternative
|      alignments (includes best score)
|      - This does not return a direction matrix
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
struct alnMatrixStruct * memWaterAltAln(
    struct seqStruct *qryST,   /*query sequence and data*/
    struct seqStruct *refST,   /*ref sequence and data*/
      /* both qryST and refST have the sequence,
      `  they also have the point to start the alignment
      `  seqST->offsetUL (index 0) and the point to end
      `  the alignment seqST->endAlnUL (index 0).
      */
    struct alnSet *settings /*Settings for the alignment*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: WatermanAltAln
   '  - Run a Waterman Smith alignment on input sequences
   '  o fun-02 sec-01:
   '    - Variable declerations
   '  o fun-02 sec-02:
   '    - Allocate memory for alignment
   '  o fun-02 sec-03:
   '    - Fill in initial negatives for ref
   '  o fun-02 sec-04:
   '    - Fill the matrix with scores
   '  o fun-02 sec-05:
   '    - Set up for returing the matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
