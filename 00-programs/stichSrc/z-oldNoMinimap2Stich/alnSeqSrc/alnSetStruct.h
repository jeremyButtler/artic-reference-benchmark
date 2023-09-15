/*######################################################################
# Name alignmentsFun
# Use:
#  o Holds functions for doing pairwise alignments (Needleman/waterman)
# Includes:
#   - "alnSeqDefaults.h"
#   - "cStrToNumberFun.h"
#   - "twoBitArrays.h"
#   - "seqStruct.h"
# C Standard libraries:
#   o <stdio.h>
#   o <stdlib.h>
#   o <stdint.h>
######################################################################*/

#ifndef ALNSETSTRUCT_H
#define ALNSETSTRUCT_H

#include "alnSeqDefaults.h"
#include "cStrToNumberFun.h"
#include "twoBitArrays.h"
#include "seqStruct.h"

/*#include <stdio.h>*/

#define defClearNonAlph (1 | 2 | 4 | 8 | 16) // clear 64 bit and case
#define defToUper (1 | 2 | 4 | 8 | 16 | 64)
    // Clear 32nd bit (marks lower case)

#define defMoveStop 0    // Do not move
#define defMoveLeft 1    // Move left (deletion) in alignment matrix
#define defMoveUp 2      // Move up (insertion) in alignment matrix
#define defMoveDiagnol 3 // Move on a diagnol (snp/match) in alignment

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o st-01 alnSet:
'     o Holds settings for my alignment program
'  o fun-01 initAlnSet:
'    - Set all values in altSet (alingment settings)
'      structure to defaults
'  o fun-02 freeAlnSet:
'    o Frees and alnSet (alignment settings) structure
'  o fun-03 setBasePairScore:
'    - Changes SNP/Match penalty for one query/reference
'      combination
'  o fun-04 getBaseScore:
'    - Get the score for a pair of bases from an alignment
'  - fun-05 readInScoreFile
'     o Reads in a file of scores for a scoring matrix
'  o fun-12 seqToLookupIndex:
'    - Converts a sequence to a look up table index
'      (table is in alnSetStruct.c/h)
'  o fun-13 lookupIndexToSeq:
'    - Converts a sequence of lookup indexs back into
'      uppercase characters (a-z)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| ST-01: alnSet
| Use: Holds settings for my alignment program
\--------------------------------------------------------*/
typedef struct alnSet
{ /*alnSet*/
   // Line wrap for printing out an alignment
   unsigned short lineWrapUS;
   unsigned short lenFileNameUS;
   char pBasePosBl; /*1 Print out base numbers*/
   char pFullAlnBl;
     /*1: Print out the full alignmnet
     ` 0: Print out the aligned region
     */
   char formatFlag;
     /*defExpandCig: is default format (S D I = X)
     ` defEMBOSS: is EMBOSS format (| space)
     ` defClustal: is clustal format (* space)
     */
   char justScoresBl;
     /*1: print coordiantes and scores only
     ` 0: Print the alignment
     */

   /*Preference for alignment algorithim used*/
   char useNeedleBl;
   char useWaterBl;
   char useHirschBl;
   char memWaterBl;

   /*Directional priorities (see alnSeqDefualts.h for
   ` options)
   */
   char bestDirC;

   /*General alignment variables*/
   int32_t gapOpenI;    /*Penalty for starting an indel*/
   int32_t gapExtendI;  /* Penalty for extending an indel*/

   /*If the user does not want to convert the sequence*/
   #ifdef NOSEQCNVT
       int16_t snpPenaltyC[26][26]; /*Scoring matrix*/
       /*No null, subtract 1*/
   #else
       int16_t snpPenaltyC[27][27]; /*Scoring matrix*/
       /*Do not subtract 1, that way 0 is null*/
   #endif
     /* Size is due to wanting a look up table that can
     `  handle anonymous bases. Most cells will be set to
     `  0. 
     `  How to get score
     `  score =
     `   snpPenaltyC[(uint8_t) (base1 & defClearNonAlph)-1]
     `              [(uint8_t) (base2 & defClearNonAlph)-1]
     */

   /*Waterman smith specific variables*/
   char refQueryScanBl;   /*Keep best score for ref/query*/
     /* If set to 1: Recored a best score for each base
     `  in the reference and query in a Waterman alignment
     */
   long minScoreL;  /*Min score to keep alignment*/
}alnSet;

/*---------------------------------------------------------------------\
| Output: Modifies: alnSetST to have default alignment settings values
\---------------------------------------------------------------------*/
void initAlnSet(
    struct alnSet *alnSetST // Alinment settings structure to initialize
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: initAlnSet
   '  - Set values in altSet (alingment settings) structure to defaults
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Frees the alnSet structure (does not set pointer to 0)
\---------------------------------------------------------------------*/
void freeAlnSet(
    struct alnSet *alnSetST,  // Alignment settings structure to free
    char stackBl              // 1: alnSetSt on stack; 0: on heap
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-1 Sub-1: freeAlnSet
   '  - Frees and alnSet (alignment settings) structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: one score in an snp/match scoring matrix
\---------------------------------------------------------------------*/
void setBasePairScore(
    const char *queryBaseC, // Query base to change score for
    const char *refBaseC,   // Reference base to change score for
    int16_t newScoreC,      // New value for [query][ref] combination
    struct alnSet *alnSetST // structure with scoring matrix to change
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-1 Sub-1: setBasePairScore
   '  - Changes SNP/Match penalty for one query/reference combination
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o score of a single pair of bases
\--------------------------------------------------------*/
static inline int16_t getBaseScore(
    const char *queryBaseC,
    const char *refBaseC,
       /*Both query and reference bases should be
       ` converted to an look up table index
       ` (seqToLookupIndex). Using -DNOSEQCNVT at compile
       ` time will have this function do the look up.
       */
    struct alnSet *alnSetST /*has scoring matrix*/
){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  ' Fun-04 TOC: getBaseScore:
  '  - Get the score for a pair of bases from an alignment
  \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   #ifdef NOSEQCNVT
      return
          alnSetST->snpPenaltyC
              [(uint8_t) (*queryBaseC & defClearNonAlph)-1]
              [(uint8_t) (*refBaseC & defClearNonAlph) -1];
   #else /*Else the base is pre-converted*/
      return
         alnSetST->snpPenaltyC
            [(uint8_t) *queryBaseC][(uint8_t) *refBaseC];
   #endif

} // getBaseScore

unsigned long readInScoreFile(
    struct alnSet *alnSetST,  // structure with scoring matrix to change
    FILE *scoreFILE           // File of scores for a scoring matrix
);  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-05 TOC: readInScoreFile
    '  - Reads in a file of scores for a scoring matrix
    '  o fun-05 sec-1: Variable declerations and buffer set up
    '  o fun-05 sec-2: Read in line and check if comment
    '  o fun-05 sec-3: Get score, update matrix, & move to next line
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o seqST->seqCStr to have look up table indexs (1-27, 
|      with null as 0) instead of bases
\--------------------------------------------------------*/
void seqToLookupIndex(
   struct seqStruct *seqST /*Sequence to convert*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: seqToLookupIndex
   '  - Converts a sequence to a look up table index
   '    (table is in alnSetStruct.c/h)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o seqST->seqCStr to have bases instead of look up
|      table indexs
\--------------------------------------------------------*/
void lookupIndexToSeq(
   struct seqStruct *seqST /*Sequence to convert*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: lookupIndexToSeq
   '  - Converts a sequence of lookup indexs back into
   '    uppercase characters (a-z)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
