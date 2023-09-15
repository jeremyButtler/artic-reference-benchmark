/*########################################################
# Name waterman
# Use:
#  o Holds functions doing a Waterman-Smith pairwise
#    alignments.
# Includes:
#   - "generalAlnFun.h"
#   - "alnStruct.h"
#   - "alnMatrixStruct.h"
#   o "twoBitArrays.h"
#   o "scoresST.h"
#   o "seqStruct.h"
#   o "alnSetStruct.h"
#   o "alnSeqDefaults.h"
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>  // by alnSetStructure.h
#   - <string.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o fun-01 WatermanSmithAln:
'    - Perform a Waterman Smith alignment on input
'      sequences
'  o fun-02 addBestBaseScore:
'    - Adds a score and index to the kept scores list
'  o fun-03 printMatrixCig:
'    - Prints out a cigar for an single path in a
'      direction matrix
'  o fun-04 printAltWaterAlns:
'    - Prints out the best aligment and the saved
'       alterantive alignments  (best alignment for each
'       base) to a file
'  o fun-06 updateDirScoreWaterSingle:
'    - Picks the best score and direction for the current
'      base pairs being compared in a Waterman-Smith
'      alignment
'    - Inlined function is at bottom of this file
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef WATERMAN_H
#define WATERMAN_H

#include <string.h>

#include "generalAlnFun.h"
#include "alnStruct.h"
#include "alnMatrixStruct.h"

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnMatrixStruct with the direction matrix and scores
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
struct alnMatrixStruct * WatermanAln(
    struct seqStruct *qryST, /*query sequence and data*/
    struct seqStruct *refST,   /*ref sequence and data*/
      /* both qryST and refST have the sequence,
      `  they also have the point to start the alignment
      `  seqST->offsetUL (index 0) and the point to end
      `  the alignment seqST->endAlnUL (index 0).
      */
    struct alnSet *settings /*Settings for the alignment*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: WatermanAln
   '  - Run a Waterman Smith alignment on input sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in initial negatives for ref
   '  o fun-01 sec-04:
   '    - Fill the matrix with scores
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnMatrixStruct with the direction matrix, best
|      score, and alternative alignments (includes best
|      score)
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
struct alnMatrixStruct * WatermanAltAln(
    struct seqStruct *qryST, /*query sequence and data*/
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

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints out all saved alternative aliginments as
|      cigars to prefix--alt.aln
|  - Returns
|    o 1 for success
|    o 2 for file error
|    o 64 for memory error
\--------------------------------------------------------*/
void printAltWaterAlns(
  struct alnMatrixStruct *alnMtxST,
     /*Has alternative alignments*/
  long minScoreL, /*Min score to keep an alternative*/
  FILE *outFILE   /*File to write alignments to*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: printAltWaterAlns
   '  - Prints out all saved alternatives alignments
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Open the output file
   '  o fun-04 sec-03:
   '    o Print out the reference alignments
   '  o fun-04 sec-04:
   '    - Print out the query aligments
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o scoreOnL to hold the best score
|    o dirOnUC to hold the best direction
\--------------------------------------------------------*/
static inline void waterTwoBitMaxScore(
    struct twoBitAry *dirOnST, /*Holds best direction*/
    struct alnSet *alnSetST,   /*for score selection*/
    long *insScL,              /*Insertion Score*/
    long *snpScL,              /*SNP/match score*/
    long *delScL,              /*Deletion score*/
    long *retScL               /*Holds best score*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: updateDirScoreWaterSingle
   '  - Picks the best score and direction for the current
   '    base pairs being compared in a Waterman Smith
   '    alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    uint8_t dirUC = 0;
    unsigned long maskUL = 0;

    #if defined SNPINSDEL
       snpInsDel(dirUC,*retScL,*insScL,*snpScL,*delScL);
    #elif defined SNPDELINS
       snpDelIns(dirUC,*retScL,*insScL,*snpScL,*delScL);
    #elif defined INSSNPDEL 
       insSnpDel(dirUC,*retScL,*insScL,*snpScL,*delScL);
    #elif defined INSDELSNP
       insDelSnp(dirUC,*retScL,*insScL,*snpScL,*delScL);
    #elif defined DELSNPINS
       delSnpIns(dirUC,*retScL,*insScL,*snpScL,*delScL);
    #elif defined DELINSSNP
       delInsSnp(dirUC,*retScL,*insScL,*snpScL,*delScL);
    #else
      switch(alnSetST->bestDirC)
      { /*Switch; get an snp/match priority*/
        case defSnpInsDel:
          snpInsDel(dirUC,*retScL,*insScL,*snpScL,*delScL);
          break;
        case defSnpDelIns:
          snpDelIns(dirUC,*retScL,*insScL,*snpScL,*delScL);
          break;
        case defInsSnpDel: 
          insSnpDel(dirUC,*retScL,*insScL,*snpScL,*delScL);
          break;
        case defInsDelSnp:
          insDelSnp(dirUC,*retScL,*insScL,*snpScL,*delScL);
          break;
        case defDelSnpIns:
          delSnpIns(dirUC,*retScL,*insScL,*snpScL,*delScL);
          break;
        case defDelInsSnp:
          delInsSnp(dirUC,*retScL,*insScL,*snpScL,*delScL);
          break;
      } /*Switch; get an snp/match priority*/
    #endif

    /*Faster than an if check (almost 30% faster). This
    ` makes it pretty close to the needleman times
    */
    maskUL = 0 - (*retScL > 0);
    dirUC = dirUC & maskUL;
    *retScL = *retScL & maskUL;
    changeTwoBitElm(dirOnST, dirUC);

    return;
} /*waterTwoBitMaxScore*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o scoreOnL to hold the best score
|    o dirOnUC to hold the best direction
\--------------------------------------------------------*/
static inline void waterByteMaxScore(
    char *dirC,                /*Holds best direction*/
    struct alnSet *alnSetST,   /*for score selection*/
    long *insScL,              /*Insertion Score*/
    long *snpScL,              /*SNP/match score*/
    long *delScL,              /*Deletion score*/
    long *retScL               /*Holds best score*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: updateDirScoreWaterSingle
   '  - Picks the best score and direction for the current
   '    base pairs being compared in a Waterman Smith
   '    alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   unsigned long maskUL = 0;

   #if defined SNPINSDEL
      snpInsDel(*dirC,*retScL,*insScL,*snpScL,*delScL);
   #elif defined SNPDELINS
      snpDelIns(*dirC,*retScL,*insScL,*snpScL,*delScL);
   #elif defined INSSNPDEL 
      insSnpDel(*dirC,*retScL,*insScL,*snpScL,*delScL);
   #elif defined INSDELSNP
      insDelSnp(*dirC,*retScL,*insScL,*snpScL,*delScL);
   #elif defined DELSNPINS
      delSnpIns(*dirC,*retScL,*insScL,*snpScL,*delScL);
   #elif defined DELINSSNP
      delInsSnp(*dirC,*retScL,*insScL,*snpScL,*delScL);
   #else
     switch(alnSetST->bestDirC)
     { /*Switch; get an snp/match priority*/
       case defSnpInsDel:
         snpInsDel(*dirC,*retScL,*insScL,*snpScL,*delScL);
         break;
       case defSnpDelIns:
         snpDelIns(*dirC,*retScL,*insScL,*snpScL,*delScL);
         break;
       case defInsSnpDel: 
         insSnpDel(*dirC,*retScL,*insScL,*snpScL,*delScL);
         break;
       case defInsDelSnp:
         insDelSnp(*dirC,*retScL,*insScL,*snpScL,*delScL);
         break;
       case defDelSnpIns:
         delSnpIns(*dirC,*retScL,*insScL,*snpScL,*delScL);
         break;
       case defDelInsSnp:
         delInsSnp(*dirC,*retScL,*insScL,*snpScL,*delScL);
         break;
     } /*Switch; get an snp/match priority*/
   #endif

   /*Faster than an if check (almost 30% faster). This
   ` makes it pretty close to the needleman times
   */
   maskUL = 0 - (*retScL > 0);
   *dirC = *dirC & maskUL;
   *retScL = *retScL & maskUL;

   return;
} /*waterByteMaxScore*/

/*--------------------------------------------------------\
| Name: updateStartePos
| Fun-06 TOC:
| Call: updateStartPos(
|          dir,
|          lastRefStart,
|          lastQryStart,
|          refStartPtr,
|          qryStartPtr
|       );
| Input:
|   - dir:
|     o Direction the current cell points to
|   - lastSnpDir:
|     o The previous direction a snp would move to
|     o This is here to allow stops to update
|   - lastRefStart:
|     o Starting point of the previous reference base
|   - lastQryStart:
|     o Starting point of the previous query base
|   - refStartPtr:
|     o Pointer to the reference position to update
|   - qryStartPtr:
|     o Pointer to the query position to update
|   - refBaseStr
|     o Base on in the reference sequence
|   - refStr
|     o First base in the reference sequence
|     o refBaseStr - refStr = reference position (index 0)
|   - qryBaseStr
|      o Base on in the query sequence
|   - qryStr 
|     o First base in the query sequence
|     o qryBaseStr - qryStr = query position (index 0)
| Output:
|  - Modifies:
|    o lastRefStart, lastQryStart to hold the previous
|      starting values for this cell (this is for SNPs)
|    o refStartPtr and qryStartPtr to hold the first base
|      in the alignment.
| Note:
|  - This assumes you are only using one array to keep
|    track of each query and reference starting position.
\--------------------------------------------------------*/
#define updateStartPos( \
   dir, \
   lastSnpDir,   /*Direction of previous snp*/\
   lastRefStart, \
   lastQryStart, \
   refStartPtr, \
   qryStartPtr, \
   refBaseStr, /*Base on in the reference sequence*/ \
   refStr,     /*First base in the reference sequence*/ \
   qryBaseStr, /*Base on in the query sequence*/ \
   qryStr      /*First base in the query sequence*/ \
){\
   unsigned long swapUL = 0; \
   switch(dir) \
   { /*Switch: Find the starting index*/ \
      case defMvStop: \
      /*Case: Last alignment ends (reset start)*/ \
         (lastRefStart) = *(refStartPtr); \
         (lastQryStart) = *(qryStartPtr); \
         break; \
      /*Case: Last alignment ends (reset start)*/ \
      case defMvDel: \
      /*Case: Added an deletion*/ \
         (lastRefStart) = *(refStartPtr); \
         (lastQryStart) = *(qryStartPtr); \
         \
         *(refStartPtr) = *((refStartPtr) - 1); \
         *(qryStartPtr) = *((qryStartPtr) - 1); \
         break; \
      /*Case: Added an deletion*/ \
      case defMvIns: \
      /*Case: Added an insertion*/ \
         (lastQryStart) = *(qryStartPtr); \
         (lastRefStart) = *(refStartPtr); \
         break;  /*Already have starting bases*/ \
      /*Case: Added an insertion*/ \
      case defMvSnp: \
      /*Case: Added an snp/match*/ \
         if(lastSnpDir == defMvStop)\
         { /*If: A snp/match would move to a stop*/\
            (lastRefStart) = *(refStartPtr); \
            (lastQryStart) = *(qryStartPtr); \
            *(refStartPtr) = (refBaseStr) - (refStr); \
            *(qryStartPtr) = (qryBaseStr) - (qryStr); \
            break; \
         } /*If: A snp/match would move to a stop*/\
         \
         swapUL = (lastRefStart); \
         (lastRefStart) = *(refStartPtr); \
         *(refStartPtr) = swapUL; \
         \
         swapUL = (lastQryStart); \
         (lastQryStart) = *(qryStartPtr); \
         *(qryStartPtr) = swapUL; \
         break; \
      /*Case: Added an snp/match*/ \
   } /*Switch: Find the starting index*/ \
} /*updateStartPos*/

/*--------------------------------------------------------\
| Name: keepAltScore
| Fun-07 TOC:
| Call: keepAltScore(
|          score,
|          qryAlt,
|          refAlt,
|          refStart,
|          refEnd,
|          qryStart,
|          qryEnd
|       );
| Input:
|   - score:
|     o Score to check if worth keeping as an alternative
|   - qryAlt:
|     o scoresStruct pointer with previous best query score
|   - refAlt:
|     o scoresStruct pointer with previous best ref score
|   - refStart:
|     o First reference base in the alignment
|   - refEnd:
|     o Last (current) reference base in the alignment
|   - qryStart:
|     o First query base in the alignment
|   - qryEnd:
|     o Last (current) query base in the alignment
| Output:
|  - Modifies:
|    o qryAlt to hold starting points, ending points, and
|      the score of the alignment (only if as good or
|      better alignment).
|    o refAlt to hold starting points, ending points, and
|      the score of the alignment (only if not chosen for
|      query and is as good or better alignment).
\--------------------------------------------------------*/
#define keepAltScore( \
   score, \
   minScore, \
   qryAlt, \
   refAlt, \
   refStart, \
   refEnd, \
   qryStart, \
   qryEnd \
){\
   if((score) >= (minScore) && (qryAlt)->scoreL < (score))\
   { /*If I am keeping score as a best query score*/ \
      (qryAlt)->scoreL = (score); \
      (qryAlt)->refStartUL = (refStart); \
      (qryAlt)->refEndUL = (refEnd); \
      (qryAlt)->qryStartUL = (qryStart); \
      (qryAlt)->qryEndUL = (qryEnd); \
   } /*If I am keeping score as a best query score*/ \
   \
   if((score) >= (minScore) && (refAlt)->scoreL < (score))\
   /*if(refAlt->scoreL < score & score >= minScore) \*/\
   { /*If I am keeping score as a best reference score*/ \
      (refAlt)->scoreL = (score); \
      (refAlt)->refStartUL = (refStart); \
      (refAlt)->refEndUL = (refEnd); \
      (refAlt)->qryStartUL = (qryStart); \
      (refAlt)->qryEndUL = (qryEnd); \
   } /*If I am keeping score as a best reference score*/ \
}

#endif
