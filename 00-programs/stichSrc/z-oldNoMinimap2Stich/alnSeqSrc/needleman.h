/*#########################################################
# Name needleman
# Use:
#  o Holds functions for doing a pairwise Needleman Wunsch
#    alignment
# Libraries:
#   - "generalAlnFun.h"
#   - "alnStruct.h"
#   - "alnMatrixStruct.h"
#   o "twoBitArrays.h"
#   o "scoresST.h"
#   o "seqStruct.h"
#   o "alnSetStruct.h"
#   o "alnSeqDefaults.h"
# C Standard libraries Used:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  - fun-01 NeedleManWunschAln:
'    o Perform a Needleman-Wunsch alignment on the two
'      input sequences
'  - fun-02 updateDirAndScore:
'    o Picks best score and direction for the current base
'      pairs being compared in a Needleman Wunsch alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef NEEDLEMAN_H
#define NEEDLEMAN_H

#include "generalAlnFun.h"
#include "alnStruct.h"
#include "alnMatrixStruct.h"

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnMatrixStruct with the directoin matrix and best
|      score struct pointing to the connor left cell in
|      the direction matrix
\--------------------------------------------------------*/
struct alnMatrixStruct * NeedlemanAln(
    struct seqStruct *queryST, // query sequence and data
    struct seqStruct *refST,  // ref sequence and data
      // both queryST and refST have the sequence,
      // they also have the point to start the alignment
      // seqST->offsetUI (index 0) and the point to end
      // the alignment seqST->endAlnUI (index 0).
      // CURRENTLY THIS ONLY WORKS FOR FULL ALIGNMENTS
    struct alnSet *settings // Settings for the alignment
    // *startI and *endI paramaters should be index 1
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: NeedlemanAln
   '  - Perform a Needleman-Wunsch alignment on the two
   '    input sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in the initial negatives for the reference
   '  o fun-01 sec-04:
   '    - Fill the matrix with scores
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o scoreOnL and dirOnUC to hold best score & direction
\--------------------------------------------------------*/
void updateDirAndScoreNeedle(
    struct twoBitAry *dirOnST,// Cell to update
    struct alnSet *alnSetST,  // for score selection
    long *scoreTopL,          // Score for an insertion
    long *scoreDiagnolL,      // Score for an match/snp
    long *scoreLeftL,         // The score for an deletion
    long *scoreOnL            // Score to update
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: updateDirAndScoreNeedle
   '  - Picks the best score and direction for the current
   '    base pairs being compared in a Needleman
   '    Wunsch alignment
   '  o fun-02 sec-01:
   '    - Matches->insertions->deletions
   '  o fun-02 sec-02:
   '    - Matches->deletions->insertions
   '  o fun-02 sec-03:
   '    - Insertions->matches->deletions
   '  o fun-02 sec-04:
   '    - Deletions->matches->insertions
   '  o fun-02 sec-05:
   '    - Insertions->deletions->matches
   '  o fun-02 sec-06:
   '    - Deletions->insertions->matches
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void selectMatchInsLeft(
    struct twoBitAry *dirOnST,// Cell to update
    long *scoreTopL,          // Score for an insertion
    long *scoreDiagnolL,      // Score for an match/snp
    long *scoreLeftL,         // The score for an deletion
    long *scoreOnL            // Score to update
);

void selectInsMatchLeft(
    struct twoBitAry *dirOnST,// Cell to update
    long *scoreTopL,          // Score for an insertion
    long *scoreDiagnolL,      // Score for an match/snp
    long *scoreLeftL,         // The score for an deletion
    long *scoreOnL            // Score to update
);

#endif
