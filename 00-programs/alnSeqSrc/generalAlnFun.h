/*########################################################
# Name generalAlnFun
# Use:
#  o Holds general functions used in my Needleman Wunsch
#    or Waterman Smith alignment.
# Libraries:
#   - "alnSetStruct.h"
#   o "alignmentSettings.h"
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>  // Used by alnSetStructures.h
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
' o fun-01 twoBitMaxScore:
'   - Picks the best score and direction for the current
'     base pairs being compared in an alignment. This
'     function is set up for two bit arrays.
'   - Function is inlied and 2nd to last function in header
' o fun-02 charMaxScore:
'   - Picks the best score and direction for the current
'     base pairs being compared in an alignment. This
'     function is set up for charters
'   - Function is inlined and is last function in header
'  o fun-03 alnMaxScore:
'    - Picks the best score for the current base pairs
'      being compared in an alignment.
'  o fun-04 checkIfBasesMatch:
'    - Are two bases are same? (includes anonymous bases)
'  o macro-01 indelScore:
'    - Calculates the score for an indel
'  o macro-02 max:
'    - Find the maximum value (branchless)
'  o macro-03 ifMax:
'    - Set a value (ret) to a value based on which value
'      is greater.
'  o macro-04 snpInsDel:
'    - Selects the max score and direction for max score.
'    - snps and then insertions when max scores are equal.
'  o macro-05 snpDelIns:
'    - Selects the max score and direction for max score.
'    - perfers snps and then deletions when max scores are
'      equal.
'  o macro-06 insSnpDel:
'    - Selects the max score and direction for max score.
'    - Perfers insertions and then snps when max scores
'      are equal.
'  o macro-07 insDelSnp:
'    - Selects the max score and direction for max score.
'    - Perfers indels and then deletions when max scores
'      are equal.
'  o macro-08: delSnpIns
'    - Selects the max score and direction for max score.
'    - perfers deletions and then snps when max scores are
'      equal.
'  o macro-09: delInsSnp
'    - Selects the max score and direction for max score.
'    - Perfers insertions and then deletions when max
'      scores are equal.
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef GENERALLALNFUN_H
#define GENERALLALNFUN_H

#include "alnSetStruct.h"

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o 1: if bases were a match
|    o 0 if bases do not mach
\--------------------------------------------------------*/
char checkIfBasesMatch(
    char *queryBaseC,// Query base to compare to reference
    char *refBaseC   // Reference base to compare to query
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-1 Sub-1: checkIfBasesMatch
   '  - Are two bases are same? (includes anonymous bases)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Macro-01: indelScore
| Use:
|  - Calculates the score for an indel
| Input:
|   - retScore:
|     o Holds the score for the indel
|   - dirC:
|     o Direction of the last cell
|   - oldScore:
|     o The previous score to add to the gap penalty
|   - settings:
|     o Has the gap opeing and extension penalties.
| Output:
|  - Modifies:
|    o retScore to hold the score for an indel
\--------------------------------------------------------*/
#define indelScore(retScore, dirC, oldScore, settings){\
   (retScore) = 0 - ((dirC) == defMvSnp); \
   (retScore) = \
        (oldScore) \
      + ((settings)->gapOpenI & (retScore)) \
      + ((settings)->gapExtendI & (~(retScore))); \
    /* Branchless, apply a gap opening and extension
    `   penalty. This is faster then my old switch
    `   statment method.
    ` tmpDirUL = 0 - (direction != diagnol)
    `   direction != diagnol
    `     Becomes 1 if the move results in a gap
    `   0 - direction != diagnol
    `     Becomes 0 - 1 if have a gap. This results in
    `       all bits being set to 1.
    `     Becomes 0 - 0 if have no gaps. This is 0.
    ` score + (gapOpen & (!dir)) + (gapExtend & dir)
    `   gapOpen & (~dir)
    `     ~dir gives zero if all bits are set to one
    `       (gap), but all bits set to one if 0.
    `     If ~dir is zero, gapOpen = 0, else applies the
    `       gap opening penalty.
    `   gapExtend & (dir)
    `     Becomes zero when there are no gaps (dir = 0).
    */ \
}

/*--------------------------------------------------------\
| Macro-02: max
| Use:
|  - Find the maximum value (branchless)
| Input:
|  o ret:
|    - Value to hold the maximum value
|  o x:
|    - First value to compare, ret is set to x if x >= y
|  o y:
|    - second value to compare, ret is set to y if y > x
| Output:
|  - Sets:
|    - Sets ret to x if x is greater than or equal to y
|    - Sets ret to y if y is greather than x
| From:
|  - https://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
\--------------------------------------------------------*/
#define max(ret, x, y){\
   (ret) = (x) ^ (((x) ^ (y)) & (-((x) < (y)))); \
   /*Logic:
   `  x < y:
   `    if x < y I get 0 (x > y); else I get 1
   `  -(x < y):
   `    If x < y I get -1 (111...)
   `    If x >= y I get 0
   `  x ^ ((x ^ y) & 111...): (only when x < y)
   `    This gives me x
   `  x ^ (x ^ y) & 0: (only when x >= y)
   `    This gives me y
   */ \
}

/*--------------------------------------------------------\
| Macro-03: ifMax
| Use:
|  - Set a value (ret) to a value based on which value
|    is greater.
| Input:
|  o ret:
|    - This will hold the return value
|  o x:
|    - First value to compare, (if x >= y)
|  o y:
|    - second value to compare, (if y > x)
|  o xRet:
|    - Value to set ret of x is >= y
|  o yRet:
|    - Value to set ret of y is > x
| Output:
|  - Sets:
|    - ret to xRet if x is greater than or equal to y
|    - ret to yRet if y is greater than x
\--------------------------------------------------------*/
#define ifMax(ret, x, y, xRet, yRet){\
   (ret) = (xRet) ^ (((xRet) ^ (yRet)) & (-((x) < (y)))); \
   /*This follows the same logic as max(ret, x, y), except
   ` instead of setting ret to the highest value, I set
   ` ret to xRet if x is >= to y or yRet if y > x.
   */ \
}

/*--------------------------------------------------------\
| Macro-04: snpInsDel
| Use:
|  - Selects the max score and direction for max score.
|  - This macro perfers snps and then insertions when max
|    scores are equal.
| Input:
|  o maxDir
|    - This will hold the direction of the max score
|  o maxSc:
|    - This will hold the max score
|  o insSc:
|    - Score for having an insertion at this position
|  o snpSc
|    - Score for having an SNP/match at this position
|  o delSc:
|    - Score for having an deletion at this position
| Output:
|  - Sets:
|    o Sets maxDir to the direction of the max score
|    - Sets maxSc to the max score
\--------------------------------------------------------*/
#define snpInsDel(maxDir, maxSc, insSc, snpSc, delSc){\
   max(maxSc, (snpSc), (insSc)); \
   ifMax((maxDir), (snpSc), (insSc), defMvSnp, defMvIns);\
   ifMax((maxDir), (maxSc), (delSc), (maxDir), defMvDel);\
   max((maxSc), (maxSc), (delSc)); \
}

/*--------------------------------------------------------\
| Macro-05: snpDelIns
| Use:
|  - Selects the max score and direction for max score.
|  - This macro perfers snps and then deletions when max
|    scores are equal.
| Input:
|  o maxDir
|    - This will hold the direction of the max score
|  o maxSc:
|    - This will hold the max score
|  o insSc:
|    - Score for having an insertion at this position
|  o snpSc
|    - Score for having an SNP/match at this position
|  o delSc:
|    - Score for having an deletion at this position
| Output:
|  - Sets:
|    o Sets maxDir to the direction of the max score
|    - Sets maxSc to the max score
\--------------------------------------------------------*/
#define snpDelIns(maxDir, maxSc, insSc, snpSc, delSc){\
   max(maxSc, (snpSc), (delSc)); \
   ifMax((maxDir), (snpSc), (delSc), defMvSnp, defMvDel);\
   ifMax((maxDir), (maxSc), (insSc), (maxDir), defMvIns);\
   max((maxSc), (maxSc), (insSc)); \
}

/*--------------------------------------------------------\
| Macro-06: insSnpDel
| Use:
|  - Selects the max score and direction for max score.
|  - This macro perfers insertions and then snps when max
|    scores are equal.
| Input:
|  o maxDir
|    - This will hold the direction of the max score
|  o maxSc:
|    - This will hold the max score
|  o insSc:
|    - Score for having an insertion at this position
|  o snpSc
|    - Score for having an SNP/match at this position
|  o delSc:
|    - Score for having an deletion at this position
| Output:
|  - Sets:
|    o Sets maxDir to the direction of the max score
|    - Sets maxSc to the max score
\--------------------------------------------------------*/
#define insSnpDel(maxDir, maxSc, insSc, snpSc, delSc){\
   max(maxSc, (insSc), (snpSc)); \
   ifMax((maxDir), (insSc), (snpSc), defMvIns, defMvSnp);\
   ifMax((maxDir), (maxSc), (delSc), (maxDir), defMvDel);\
   max((maxSc), (maxSc), (delSc)); \
}

/*--------------------------------------------------------\
| Macro-07: insDelSnp
| Use:
|  - Selects the max score and direction for max score.
|  - This macro perfers indels and then deletions when max
|    scores are equal.
| Input:
|  o maxDir
|    - This will hold the direction of the max score
|  o maxSc:
|    - This will hold the max score
|  o insSc:
|    - Score for having an insertion at this position
|  o snpSc
|    - Score for having an SNP/match at this position
|  o delSc:
|    - Score for having an deletion at this position
| Output:
|  - Sets:
|    o Sets maxDir to the direction of the max score
|    - Sets maxSc to the max score
\--------------------------------------------------------*/
#define insDelSnp(maxDir, maxSc, insSc, snpSc, delSc){\
   max(maxSc, (insSc), (delSc)); \
   ifMax((maxDir), (insSc), (delSc), defMvIns, defMvDel);\
   ifMax((maxDir), (maxSc), (snpSc), (maxDir), defMvSnp);\
   max((maxSc), (maxSc), (snpSc)); \
}

/*--------------------------------------------------------\
| Macro-08: delSnpIns
| Use:
|  - Selects the max score and direction for max score.
|  - This macro perfers deletions and then snps when
|    max scores are equal.
| Input:
|  o maxDir
|    - This will hold the direction of the max score
|  o maxSc:
|    - This will hold the max score
|  o insSc:
|    - Score for having an insertion at this position
|  o snpSc
|    - Score for having an SNP/match at this position
|  o delSc:
|    - Score for having an deletion at this position
| Output:
|  - Sets:
|    o Sets maxDir to the direction of the max score
|    - Sets maxSc to the max score
\--------------------------------------------------------*/
#define delSnpIns(maxDir, maxSc, insSc, snpSc, delSc){\
   max(maxSc, (delSc), (snpSc)); \
   ifMax((maxDir), (delSc), (snpSc), defMvDel, defMvSnp);\
   ifMax((maxDir), (maxSc), (insSc), (maxDir), defMvIns);\
   max((maxSc), (maxSc), (insSc)); \
}

/*--------------------------------------------------------\
| Macro-09: delInsSnp
| Use:
|  - Selects the max score and direction for max score.
|  - This macro perfers insertions and then deletions when
|    max scores are equal.
| Input:
|  o maxDir
|    - This will hold the direction of the max score
|  o maxSc:
|    - This will hold the max score
|  o insSc:
|    - Score for having an insertion at this position
|  o snpSc
|    - Score for having an SNP/match at this position
|  o delSc:
|    - Score for having an deletion at this position
| Output:
|  - Sets:
|    o Sets maxDir to the direction of the max score
|    - Sets maxSc to the max score
\--------------------------------------------------------*/
#define delInsSnp(maxDir, maxSc, insSc, snpSc, delSc){\
   max(maxSc, (delSc), (insSc)); \
   ifMax((maxDir), (delSc), (insSc), defMvDel, defMvIns);\
   ifMax((maxDir), (maxSc), (snpSc), (maxDir), defMvSnp);\
   max((maxSc), (maxSc), (snpSc)); \
}

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o scoreOnL and dirOnST to hold best score & direction
\--------------------------------------------------------*/
static inline void twoBitMaxScore(
    struct twoBitAry *dirOnST, /*Holds best direction*/
    struct alnSet *alnSetST,   /*for score selection*/
    long *insScL,              /*Insertion Score*/
    long *snpScL,              /*SNP/match score*/
    long *delScL,              /*Deletion score*/
    long *retScL               /*Holds best score*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: twoBitMaxScore
   '  - Picks the best score and direction for the current
   '    base pairs being compared in an alignment. This
   '    function is set up for two bit arrays.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   uint8_t dirUC = 0;

   #if defined SNPINSDEL
      snpInsDel(dirUC,*retScL,*insScL,*snpScL,*delScL);
      changeTwoBitElm(dirOnST, dirUC);
      return;
   #elif defined SNPDELINS
      snpDelIns(dirUC,*retScL,*insScL,*snpScL,*delScL);
      changeTwoBitElm(dirOnST, dirUC);
      return;
   #elif defined INSSNPDEL 
      insSnpDel(dirUC,*retScL,*insScL,*snpScL,*delScL);
      changeTwoBitElm(dirOnST, dirUC);
      return;
   #elif defined INSDELSNP
      insDelSnp(dirUC,*retScL,*insScL,*snpScL,*delScL);
      changeTwoBitElm(dirOnST, dirUC);
      return;
   #elif defined DELSNPINS
      delSnpIns(dirUC,*retScL,*insScL,*snpScL,*delScL);
      changeTwoBitElm(dirOnST, dirUC);
      return;
   #elif defined DELINSSNP
      delInsSnp(dirUC,*retScL,*insScL,*snpScL,*delScL);
      changeTwoBitElm(dirOnST, dirUC);
      return;

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

    changeTwoBitElm(dirOnST, dirUC);
    return;
  #endif
} /*twoBitMaxScore*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o scoreOnL and dirOnC to hold best score & direction
\--------------------------------------------------------*/
static inline void charMaxScore(
    char *dirOnC,              /*Holds best direction*/
    struct alnSet *alnSetST,   /*for score selection*/
    long *insScL,              /*Insertion Score*/
    long *snpScL,              /*SNP/match score*/
    long *delScL,              /*Deletion score*/
    long *retScL               /*Holds best score*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: charMaxScore
   '  - Picks the best score and direction for the current
   '    base pairs being compared in an alignment. This
   '    function is set up for charters
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   #if defined SNPINSDEL
     snpInsDel(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
     return;
   #elif defined SNPDELINS
     snpDelIns(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
     return;
   #elif defined INSSNPDEL 
     insSnpDel(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
     return;
   #elif defined INSDELSNP
     insDelSnp(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
     return;
   #elif defined DELSNPINS
     delSnpIns(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
     return;
   #elif defined DELINSSNP
     delInsSnp(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
     return;

   #else
   switch(alnSetST->bestDirC)
   { /*Switch; get an snp/match priority*/
      case defSnpInsDel:
        snpInsDel(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
        return;
      case defSnpDelIns:
        snpDelIns(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
        return;
      case defInsSnpDel: 
        insSnpDel(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
        return;
      case defInsDelSnp:
        insDelSnp(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
        return;
      case defDelSnpIns:
        delSnpIns(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
        return;
      case defDelInsSnp:
        delInsSnp(*dirOnC,*retScL,*insScL,*snpScL,*delScL);
        return;
   } /*Switch; get an snp/match priority*/

   return;
   #endif
} /*charMaxScore*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o scoreOnL and dirOnC to hold best score & direction
\--------------------------------------------------------*/
static inline void alnMaxScore(
    struct alnSet *alnSetST,   /*for score selection*/
    long *insScL,              /*Insertion Score*/
    long *snpScL,              /*SNP/match score*/
    long *delScL,              /*Deletion score*/
    long *retScL               /*Holds best score*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: alnMaxScore
   '  - Picks the best score for the current base pairs
   '    being compared in an alignment.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   #if defined SNPINSDEL
      max(*retScL, *snpScL, *insScL);
      max(*retScL, *retScL, *delScL);
      return;
   #elif defined SNPDELINS
      max(*retScL, *snpScL, *delScL);
      max(*retScL, *retScL, *insScL);
      return;
   #elif defined INSSNPDEL
      max(*retScL, *insScL, *snpScL);
      max(*retScL, *retScL, *delScL);
      return;
   #elif defined INSDELSNP
      max(*retScL, *insScL, *delScL);
      max(*retScL, *retScL, *snpScL);
      return;
   #elif defined DELSNPINS
      max(*retScL, *delScL, *snpScL);
      max(*retScL, *retScL, *insScL);
      return;
   #elif defined DELINSSNP
      max(*retScL, *delScL, *insScL);
      max(*retScL, *retScL, *snpScL);
      return;
   #else
      switch(alnSetST->bestDirC)
      { /*Switch; get an snp/match priority*/
         case defSnpInsDel:
           max(*retScL, *snpScL, *insScL);
           max(*retScL, *retScL, *delScL);
           return;
         case defSnpDelIns:
           max(*retScL, *snpScL, *delScL);
           max(*retScL, *retScL, *insScL);
           return;
         case defInsSnpDel: 
           max(*retScL, *insScL, *snpScL);
           max(*retScL, *retScL, *delScL);
           return;
         case defInsDelSnp:
           max(*retScL, *insScL, *delScL);
           max(*retScL, *retScL, *snpScL);
           return;
         case defDelSnpIns:
           max(*retScL, *delScL, *snpScL);
           max(*retScL, *retScL, *insScL);
           return;
         case defDelInsSnp:
           max(*retScL, *delScL, *insScL);
           max(*retScL, *retScL, *snpScL);
           return;
      } /*Switch; get an snp/match priority*/

      return;
   #endif
} /*alnMaxScore*/

#endif

