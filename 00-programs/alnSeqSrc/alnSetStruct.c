/*#########################################################
# Name: alnSetStructure
# Use:
#   - Holds the alignment settings structure and the
#     functions need to initalize, free, and modify
#     the score matrix
# Libraries:
#  - "alnSeqDefaults.h"
#  - "seqStruct.h"
# C Standard Librarys Used:
#  - <stdint.h>
#  - <stdlib.h>
#  - <stdio.h>
\########################################################*/

#include "alnSetStruct.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  o fun-01 initAlnSet:
'    - Set all values in altSet (alingment settings)
'      structure to defaults
'  o fun-02 freeAlnSet:
'    o Frees and alnSet (alignment settings) structure
'  o fun-03 setBasePairScore:
'    - Changes SNP/Match penalty for one query/reference
'      combination
'  o fun-04 getBasePairScore:
'    - Get the score for a pair of bases from an alignment
'    - Inlined function (in header file only)
'  o fun-05 readInScoreFile
'     o Reads in a file of scores for a scoring matrix
'  o fun-06 seqToLookupIndex:
'    - Converts a sequence to a look up table index
'      (table is in alnSetStruct.c/h)
'  o fun-07 lookupIndexToSeq:
'    - Converts a sequence of lookup indexs back into
'      uppercase characters (a-z)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  o Modifies:
|    - alnSetST to have default alignment settings values
\--------------------------------------------------------*/
void initAlnSet(
    struct alnSet *alnSetST // Has settings to initialize
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: initAlnSet
   '  - Set values in altSet (alingment settings)
   '    structure to defaults
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*Set the default output settings*/
   alnSetST->lineWrapUS = defLineWrap;
   alnSetST->lenFileNameUS = 1024;

   alnSetST->pBasePosBl = defPPos;
   alnSetST->pFullAlnBl = defPAln;
   alnSetST->formatFlag = defFormat;
   alnSetST->justScoresBl = defJustScoresBl;

   /*Choose the aligment style*/
   alnSetST->useNeedleBl = defUseNeedle;
   alnSetST->useWaterBl = defUseWater;
   alnSetST->useHirschBl = defUseHirsch;
   alnSetST->memWaterBl = defUseMemWater;

   /*Select direction to keep if everything is equal*/
   alnSetST->bestDirC = defBestDir;

   /*General alignment variables*/
   alnSetST->gapOpenI = defGapOpen;
   alnSetST->gapExtendI = defGapExtend;

   /*Waterman specific variables*/
   alnSetST->refQueryScanBl = defQueryRefScan;
   alnSetST->minScoreL = defMinScore;
  
   /*score matrix Initalize*/
   for(uint8_t colUC = 0; colUC < 26; ++colUC)
   { // loop for all columns in the comparison matrix
       for(uint8_t rowUC = 0; rowUC < 26; ++rowUC)
           alnSetST->snpPenaltyC[colUC][rowUC] = 0;
           // Most of these cells will never be used
           // But are needed to build the table
   } // loop for all columns in the comparison matrix

   // Set up scores for non-anonmyous base pairs
   setBasePairScore("a", "a", defAToA, alnSetST);
   setBasePairScore("a", "t", defAToT, alnSetST);
   setBasePairScore("a", "u", defAToT, alnSetST);
   setBasePairScore("a", "g", defAToG, alnSetST);
   setBasePairScore("a", "c", defAToC, alnSetST);

   setBasePairScore("t", "a", defTToA, alnSetST);
   setBasePairScore("t", "t", defTToT, alnSetST);
   setBasePairScore("t", "g", defTToG, alnSetST);
   setBasePairScore("t", "c", defTToC, alnSetST);

   setBasePairScore("u", "a", defTToA, alnSetST);
   setBasePairScore("u", "u", defTToT, alnSetST);
   setBasePairScore("u", "g", defTToG, alnSetST);
   setBasePairScore("u", "c", defTToC, alnSetST);

   setBasePairScore("g", "a", defGToA, alnSetST);
   setBasePairScore("g", "t", defGToT, alnSetST);
   setBasePairScore("g", "u", defGToT, alnSetST);
   setBasePairScore("g", "g", defGToG, alnSetST);
   setBasePairScore("g", "c", defGToC, alnSetST);

   setBasePairScore("c", "a", defCToA, alnSetST);
   setBasePairScore("c", "t", defCToT, alnSetST);
   setBasePairScore("c", "u", defCToT, alnSetST);
   setBasePairScore("c", "g", defCToG, alnSetST);
   setBasePairScore("c", "c", defCToC, alnSetST);

   // non-anonymous base and anonymous base pairs
   setBasePairScore("a", "w", defAToW, alnSetST);
   setBasePairScore("a", "s", defAToS, alnSetST);
   setBasePairScore("a", "m", defAToM, alnSetST);
   setBasePairScore("a", "k", defAToK, alnSetST);
   setBasePairScore("a", "r", defAToR, alnSetST);
   setBasePairScore("a", "y", defAToY, alnSetST);
   setBasePairScore("a", "b", defAToB, alnSetST);
   setBasePairScore("a", "d", defAToD, alnSetST);
   setBasePairScore("a", "h", defAToH, alnSetST);
   setBasePairScore("a", "v", defAToV, alnSetST);
   setBasePairScore("a", "n", defAToN, alnSetST);
   setBasePairScore("a", "x", defAToX, alnSetST);

   setBasePairScore("w", "a", defWToA, alnSetST);
   setBasePairScore("s", "a", defSToA, alnSetST);
   setBasePairScore("m", "a", defMToA, alnSetST);
   setBasePairScore("k", "a", defKToA, alnSetST);
   setBasePairScore("r", "a", defRToA, alnSetST);
   setBasePairScore("y", "a", defYToA, alnSetST);
   setBasePairScore("b", "a", defBToA, alnSetST);
   setBasePairScore("d", "a", defDToA, alnSetST);
   setBasePairScore("h", "a", defHToA, alnSetST);
   setBasePairScore("v", "a", defVToA, alnSetST);
   setBasePairScore("n", "a", defNToA, alnSetST);
   setBasePairScore("x", "a", defXToA, alnSetST);

   setBasePairScore("c", "w", defCToW, alnSetST);
   setBasePairScore("c", "s", defCToS, alnSetST);
   setBasePairScore("c", "m", defCToM, alnSetST);
   setBasePairScore("c", "k", defCToK, alnSetST);
   setBasePairScore("c", "r", defCToR, alnSetST);
   setBasePairScore("c", "y", defCToY, alnSetST);
   setBasePairScore("c", "b", defCToB, alnSetST);
   setBasePairScore("c", "d", defCToD, alnSetST);
   setBasePairScore("c", "h", defCToH, alnSetST);
   setBasePairScore("c", "v", defCToV, alnSetST);
   setBasePairScore("c", "n", defCToN, alnSetST);
   setBasePairScore("c", "x", defCToX, alnSetST);

   setBasePairScore("w", "C", defWToC, alnSetST);
   setBasePairScore("s", "C", defSToC, alnSetST);
   setBasePairScore("m", "C", defMToC, alnSetST);
   setBasePairScore("k", "C", defKToC, alnSetST);
   setBasePairScore("r", "C", defRToC, alnSetST);
   setBasePairScore("y", "C", defYToC, alnSetST);
   setBasePairScore("b", "C", defBToC, alnSetST);
   setBasePairScore("d", "C", defDToC, alnSetST);
   setBasePairScore("h", "C", defHToC, alnSetST);
   setBasePairScore("v", "C", defVToC, alnSetST);
   setBasePairScore("n", "C", defNToC, alnSetST);
   setBasePairScore("x", "C", defXToC, alnSetST);

   setBasePairScore("g", "w", defGToW, alnSetST);
   setBasePairScore("g", "s", defGToS, alnSetST);
   setBasePairScore("g", "m", defGToM, alnSetST);
   setBasePairScore("g", "k", defGToK, alnSetST);
   setBasePairScore("g", "r", defGToR, alnSetST);
   setBasePairScore("g", "y", defGToY, alnSetST);
   setBasePairScore("g", "b", defGToB, alnSetST);
   setBasePairScore("g", "d", defGToD, alnSetST);
   setBasePairScore("g", "h", defGToH, alnSetST);
   setBasePairScore("g", "v", defGToV, alnSetST);
   setBasePairScore("g", "n", defGToN, alnSetST);
   setBasePairScore("g", "x", defGToX, alnSetST);

   setBasePairScore("w", "g", defWToG, alnSetST);
   setBasePairScore("s", "g", defSToG, alnSetST);
   setBasePairScore("m", "g", defMToG, alnSetST);
   setBasePairScore("k", "g", defKToG, alnSetST);
   setBasePairScore("r", "g", defRToG, alnSetST);
   setBasePairScore("y", "g", defYToG, alnSetST);
   setBasePairScore("b", "g", defBToG, alnSetST);
   setBasePairScore("d", "g", defDToG, alnSetST);
   setBasePairScore("h", "g", defHToG, alnSetST);
   setBasePairScore("v", "g", defVToG, alnSetST);
   setBasePairScore("n", "g", defNToG, alnSetST);
   setBasePairScore("x", "g", defXToG, alnSetST);

   setBasePairScore("t", "w", defTToW, alnSetST);
   setBasePairScore("t", "s", defTToS, alnSetST);
   setBasePairScore("t", "m", defTToM, alnSetST);
   setBasePairScore("t", "k", defTToK, alnSetST);
   setBasePairScore("t", "r", defTToR, alnSetST);
   setBasePairScore("t", "y", defTToY, alnSetST);
   setBasePairScore("t", "b", defTToB, alnSetST);
   setBasePairScore("t", "d", defTToD, alnSetST);
   setBasePairScore("t", "h", defTToH, alnSetST);
   setBasePairScore("t", "v", defTToV, alnSetST);
   setBasePairScore("t", "n", defTToN, alnSetST);
   setBasePairScore("t", "x", defTToX, alnSetST);

   setBasePairScore("w", "t", defWToT, alnSetST);
   setBasePairScore("s", "t", defSToT, alnSetST);
   setBasePairScore("m", "t", defMToT, alnSetST);
   setBasePairScore("k", "t", defKToT, alnSetST);
   setBasePairScore("r", "t", defRToT, alnSetST);
   setBasePairScore("y", "t", defYToT, alnSetST);
   setBasePairScore("b", "t", defBToT, alnSetST);
   setBasePairScore("d", "t", defDToT, alnSetST);
   setBasePairScore("h", "t", defHToT, alnSetST);
   setBasePairScore("v", "t", defVToT, alnSetST);
   setBasePairScore("n", "t", defNToT, alnSetST);
   setBasePairScore("x", "t", defXToT, alnSetST);

   // Set u and t to same scores (U is RNA version of T)
   setBasePairScore("u", "w", defTToW, alnSetST);
   setBasePairScore("u", "s", defTToS, alnSetST);
   setBasePairScore("u", "m", defTToM, alnSetST);
   setBasePairScore("u", "k", defTToK, alnSetST);
   setBasePairScore("u", "r", defTToR, alnSetST);
   setBasePairScore("u", "y", defTToY, alnSetST);
   setBasePairScore("u", "b", defTToB, alnSetST);
   setBasePairScore("u", "d", defTToD, alnSetST);
   setBasePairScore("u", "h", defTToH, alnSetST);
   setBasePairScore("u", "v", defTToV, alnSetST);
   setBasePairScore("u", "n", defTToN, alnSetST);
   setBasePairScore("u", "x", defTToX, alnSetST);

   setBasePairScore("w", "u", defWToT, alnSetST);
   setBasePairScore("s", "u", defSToT, alnSetST);
   setBasePairScore("m", "u", defMToT, alnSetST);
   setBasePairScore("k", "u", defKToT, alnSetST);
   setBasePairScore("r", "u", defRToT, alnSetST);
   setBasePairScore("y", "u", defYToT, alnSetST);
   setBasePairScore("b", "u", defBToT, alnSetST);
   setBasePairScore("d", "u", defDToT, alnSetST);
   setBasePairScore("h", "u", defHToT, alnSetST);
   setBasePairScore("v", "u", defVToT, alnSetST);
   setBasePairScore("n", "u", defNToT, alnSetST);
   setBasePairScore("x", "u", defXToT, alnSetST);

   // anonymous base and anonymous base pairs
   setBasePairScore("w", "w", defWToW, alnSetST);
   setBasePairScore("w", "s", defWToS, alnSetST);
   setBasePairScore("w", "m", defWToM, alnSetST);
   setBasePairScore("w", "k", defWToK, alnSetST);
   setBasePairScore("w", "r", defWToR, alnSetST);
   setBasePairScore("w", "y", defWToY, alnSetST);
   setBasePairScore("w", "b", defWToB, alnSetST);
   setBasePairScore("w", "d", defWToD, alnSetST);
   setBasePairScore("w", "h", defWToH, alnSetST);
   setBasePairScore("w", "v", defWToV, alnSetST);
   setBasePairScore("w", "n", defWToN, alnSetST);
   setBasePairScore("w", "x", defWToX, alnSetST);

   setBasePairScore("s", "w", defSToW, alnSetST);
   setBasePairScore("s", "s", defSToS, alnSetST);
   setBasePairScore("s", "m", defSToM, alnSetST);
   setBasePairScore("s", "k", defSToK, alnSetST);
   setBasePairScore("s", "r", defSToR, alnSetST);
   setBasePairScore("s", "y", defSToY, alnSetST);
   setBasePairScore("s", "b", defSToB, alnSetST);
   setBasePairScore("s", "d", defSToD, alnSetST);
   setBasePairScore("s", "h", defSToH, alnSetST);
   setBasePairScore("s", "v", defSToV, alnSetST);
   setBasePairScore("s", "n", defSToN, alnSetST);
   setBasePairScore("s", "x", defSToX, alnSetST);

   setBasePairScore("m", "w", defMToW, alnSetST);
   setBasePairScore("m", "s", defMToS, alnSetST);
   setBasePairScore("m", "m", defMToM, alnSetST);
   setBasePairScore("m", "k", defMToK, alnSetST);
   setBasePairScore("m", "r", defMToR, alnSetST);
   setBasePairScore("m", "y", defMToY, alnSetST);
   setBasePairScore("m", "b", defMToB, alnSetST);
   setBasePairScore("m", "d", defMToD, alnSetST);
   setBasePairScore("m", "h", defMToH, alnSetST);
   setBasePairScore("m", "v", defMToV, alnSetST);
   setBasePairScore("m", "n", defMToN, alnSetST);
   setBasePairScore("m", "x", defMToX, alnSetST);

   setBasePairScore("k", "w", defKToW, alnSetST);
   setBasePairScore("k", "s", defKToS, alnSetST);
   setBasePairScore("k", "m", defKToM, alnSetST);
   setBasePairScore("k", "k", defKToK, alnSetST);
   setBasePairScore("k", "r", defKToR, alnSetST);
   setBasePairScore("k", "y", defKToY, alnSetST);
   setBasePairScore("k", "b", defKToB, alnSetST);
   setBasePairScore("k", "d", defKToD, alnSetST);
   setBasePairScore("k", "h", defKToH, alnSetST);
   setBasePairScore("k", "v", defKToV, alnSetST);
   setBasePairScore("k", "n", defKToN, alnSetST);
   setBasePairScore("k", "x", defKToX, alnSetST);

   setBasePairScore("r", "w", defRToW, alnSetST);
   setBasePairScore("r", "s", defRToS, alnSetST);
   setBasePairScore("r", "m", defRToM, alnSetST);
   setBasePairScore("r", "k", defRToK, alnSetST);
   setBasePairScore("r", "r", defRToR, alnSetST);
   setBasePairScore("r", "y", defRToY, alnSetST);
   setBasePairScore("r", "b", defRToB, alnSetST);
   setBasePairScore("r", "d", defRToD, alnSetST);
   setBasePairScore("r", "h", defRToH, alnSetST);
   setBasePairScore("r", "v", defRToV, alnSetST);
   setBasePairScore("r", "n", defRToN, alnSetST);
   setBasePairScore("r", "x", defRToX, alnSetST);

   setBasePairScore("y", "w", defYToW, alnSetST);
   setBasePairScore("y", "s", defYToS, alnSetST);
   setBasePairScore("y", "m", defYToM, alnSetST);
   setBasePairScore("y", "k", defYToK, alnSetST);
   setBasePairScore("y", "r", defYToR, alnSetST);
   setBasePairScore("y", "y", defYToY, alnSetST);
   setBasePairScore("y", "b", defYToB, alnSetST);
   setBasePairScore("y", "d", defYToD, alnSetST);
   setBasePairScore("y", "h", defYToH, alnSetST);
   setBasePairScore("y", "v", defYToV, alnSetST);
   setBasePairScore("y", "n", defYToN, alnSetST);
   setBasePairScore("y", "x", defYToX, alnSetST);

   setBasePairScore("b", "w", defBToW, alnSetST);
   setBasePairScore("b", "s", defBToS, alnSetST);
   setBasePairScore("b", "m", defBToM, alnSetST);
   setBasePairScore("b", "k", defBToK, alnSetST);
   setBasePairScore("b", "r", defBToR, alnSetST);
   setBasePairScore("b", "y", defBToY, alnSetST);
   setBasePairScore("b", "b", defBToB, alnSetST);
   setBasePairScore("b", "d", defBToD, alnSetST);
   setBasePairScore("b", "h", defBToH, alnSetST);
   setBasePairScore("b", "v", defBToV, alnSetST);
   setBasePairScore("b", "n", defBToN, alnSetST);
   setBasePairScore("b", "x", defBToX, alnSetST);

   setBasePairScore("d", "w", defDToW, alnSetST);
   setBasePairScore("d", "s", defDToS, alnSetST);
   setBasePairScore("d", "m", defDToM, alnSetST);
   setBasePairScore("d", "k", defDToK, alnSetST);
   setBasePairScore("d", "r", defDToR, alnSetST);
   setBasePairScore("d", "y", defDToY, alnSetST);
   setBasePairScore("d", "b", defDToB, alnSetST);
   setBasePairScore("d", "d", defDToD, alnSetST);
   setBasePairScore("d", "h", defDToH, alnSetST);
   setBasePairScore("d", "v", defDToV, alnSetST);
   setBasePairScore("d", "n", defDToN, alnSetST);
   setBasePairScore("d", "x", defDToX, alnSetST);

   setBasePairScore("h", "w", defHToW, alnSetST);
   setBasePairScore("h", "s", defHToS, alnSetST);
   setBasePairScore("h", "m", defHToM, alnSetST);
   setBasePairScore("h", "k", defHToK, alnSetST);
   setBasePairScore("h", "r", defHToR, alnSetST);
   setBasePairScore("h", "y", defHToY, alnSetST);
   setBasePairScore("h", "b", defHToB, alnSetST);
   setBasePairScore("h", "d", defHToD, alnSetST);
   setBasePairScore("h", "h", defHToH, alnSetST);
   setBasePairScore("h", "v", defHToV, alnSetST);
   setBasePairScore("h", "n", defHToN, alnSetST);
   setBasePairScore("h", "x", defHToX, alnSetST);

   setBasePairScore("v", "w", defVToW, alnSetST);
   setBasePairScore("v", "s", defVToS, alnSetST);
   setBasePairScore("v", "m", defVToM, alnSetST);
   setBasePairScore("v", "k", defVToK, alnSetST);
   setBasePairScore("v", "r", defVToR, alnSetST);
   setBasePairScore("v", "y", defVToY, alnSetST);
   setBasePairScore("v", "b", defVToB, alnSetST);
   setBasePairScore("v", "d", defVToD, alnSetST);
   setBasePairScore("v", "h", defVToH, alnSetST);
   setBasePairScore("v", "v", defVToV, alnSetST);
   setBasePairScore("v", "n", defVToN, alnSetST);
   setBasePairScore("v", "x", defVToX, alnSetST);

   setBasePairScore("n", "w", defNToW, alnSetST);
   setBasePairScore("n", "s", defNToS, alnSetST);
   setBasePairScore("n", "m", defNToM, alnSetST);
   setBasePairScore("n", "k", defNToK, alnSetST);
   setBasePairScore("n", "r", defNToR, alnSetST);
   setBasePairScore("n", "y", defNToY, alnSetST);
   setBasePairScore("n", "b", defNToB, alnSetST);
   setBasePairScore("n", "d", defNToD, alnSetST);
   setBasePairScore("n", "h", defNToH, alnSetST);
   setBasePairScore("n", "v", defNToV, alnSetST);
   setBasePairScore("n", "n", defNToN, alnSetST);
   setBasePairScore("n", "x", defNToX, alnSetST);

   setBasePairScore("x", "w", defXToW, alnSetST);
   setBasePairScore("x", "s", defXToS, alnSetST);
   setBasePairScore("x", "m", defXToM, alnSetST);
   setBasePairScore("x", "k", defXToK, alnSetST);
   setBasePairScore("x", "r", defXToR, alnSetST);
   setBasePairScore("x", "y", defXToY, alnSetST);
   setBasePairScore("x", "b", defXToB, alnSetST);
   setBasePairScore("x", "d", defXToD, alnSetST);
   setBasePairScore("x", "h", defXToH, alnSetST);
   setBasePairScore("x", "v", defXToV, alnSetST);
   setBasePairScore("x", "n", defXToN, alnSetST);
   setBasePairScore("x", "x", defXToX, alnSetST);

   return;
} // initAlnSet

/*--------------------------------------------------------\
| Output: Frees alnSet structure (pointer not set to 0)
\--------------------------------------------------------*/
void freeAlnSet(
    struct alnSet *alnSetST,  // structure to free
    char stackBl       // 1: alnSetSt on stack; 0: on heap
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-1 Sub-1: freeAlnSet
   '  - Frees and alnSet (alignment settings) structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(stackBl & 1) return; // Nothing to do
   if(alnSetST != 0) free(alnSetST); // No heap variables
   return; 
} // freeAlnSet
    
/*--------------------------------------------------------\
| Output:
|  o Modifies:
|    - one score in an snp/match scoring matrix
\--------------------------------------------------------*/
void setBasePairScore(
    const char *queryBaseC, // Query base to change score
    const char *refBaseC, // Reference base to change score
    int16_t newScoreC,  // New score for [query][ref] pair
    struct alnSet *alnSetST // has scoring matrix to change
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: setBasePairScore
   '  - Sets the SNP/Match penalty for a query/reference
   '    base pair
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   #ifdef NOSEQCNVT
      alnSetST->snpPenaltyC
          [(uint8_t) (*queryBaseC & defClearNonAlph) - 1]
          [(uint8_t) (*refBaseC & defClearNonAlph) - 1]
          =
          newScoreC;
   #else
      alnSetST->snpPenaltyC
          [(uint8_t) (*queryBaseC & defClearNonAlph)]
          [(uint8_t) (*refBaseC & defClearNonAlph)]
          =
          newScoreC;
   #endif

   return;
} // setBasePairScore

/*--------------------------------------------------------\
| Output:
|  o Modifies:
|    - Score matrix in alngSetST to hold the scores from
|      the file (scoreFILE)
\--------------------------------------------------------*/
unsigned long readInScoreFile(
    struct alnSet *alnSetST,// has scoring matrix to change
    FILE *scoreFILE    // File of scores for scoring matrix
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: readInScoreFile
   '  - Reads in a file of scores for a scoring matrix
   '  o fun-05 sec-1: Variable declerations & set up
   '  o fun-05 sec-2: Read in line and check if comment
   '  o fun-05 sec-3: Convert score & add to matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-1: Variable declerations and set up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint16_t lenBuffUS = 1024;
   char buffCStr[lenBuffUS];
   char *tmpCStr = 0;
   int16_t scoreS = 0;

   buffCStr[lenBuffUS - 1] = '\0';
   buffCStr[lenBuffUS - 2] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-2: Read in line and check if comment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(fgets(buffCStr, 1024, scoreFILE))
   { // While I have scores to read in
       
       if(buffCStr[0] == '/' && buffCStr[1] == '/')
       { // On a comment, move onto the next line
           while(
               buffCStr[lenBuffUS - 2] != '\0' &&
               buffCStr[lenBuffUS - 2] != '\n'
           ) { // While have more buffer to read in
               buffCStr[lenBuffUS - 2] = '\0';
               fgets(buffCStr, 1024, scoreFILE);
           } // While have more buffer to read in

           // Reset the buffer
           buffCStr[lenBuffUS - 2] = '\0';

           continue;
       } // On a comment, move onto the next line

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-05 Sec-3: Convert score & add to matrix
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       if(buffCStr[0] == '\n')
           continue;                        // Blank line

       if(buffCStr[0] < 64 && buffCStr[2] < 64)
           return ftell(scoreFILE);  // Invalid character
       
       tmpCStr = cStrToInt16(&buffCStr[4], &scoreS);

       setBasePairScore(
         &buffCStr[0],
         &buffCStr[2],
         scoreS,
         alnSetST
       ); // Add the score to the matrix

       if(tmpCStr == &buffCStr[3])
           return ftell(scoreFILE);         // No score

       while(
           buffCStr[lenBuffUS - 2] != '\0' &&
           buffCStr[lenBuffUS - 2] != '\n'
       ) { // While have more buffer to read in
           buffCStr[lenBuffUS - 2] = '\0';
           fgets(buffCStr, 1024, scoreFILE);
       } // While have more buffer to read in

       // Reset the buffer
       buffCStr[lenBuffUS - 2] = '\0';
   } // While I have scores to read in

   return 0;
} // readInScoreFile

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o seqST->seqCStr to have look up table indexs (1-27, 
|      with null as 0) instead of bases
\--------------------------------------------------------*/
void seqToLookupIndex(
   struct seqStruct *seqST /*Sequence to convert*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: seqToLookupIndex
   '  - Converts a sequence to a look up table index
   '    (table is in alnSetStruct.c/h)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   #ifdef NOSEQCNVT
      return;
      /*Do nothing,
      ` look up functions will handle conversion.
      */
   #endif

   char *tmpStr = seqST->seqCStr;

   while(*tmpStr != '\0')
   { /*Loop: convert bases to lookup table values*/
      *tmpStr &= defClearNonAlph;
      ++tmpStr;
   } /*Loop: convert bases to lookup table values*/

   return;
} /*seqToLookupIndex*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o seqST->seqCStr to have bases instead of look up
|      table indexs
\--------------------------------------------------------*/
void lookupIndexToSeq(
   struct seqStruct *seqST /*Sequence to convert*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: lookupIndexToSeq
   '  - Converts a sequence of lookup indexs back into
   '    uppercase characters (a-z)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   #ifdef NOSEQCNVT
      return; /*Do nothing, never did a conversion*/
   #endif

   char *tmpStr = seqST->seqCStr;

   while(*tmpStr != '\0')
   { /*Loop: convert bases to lookup table values*/
      *tmpStr |= 64;
      ++tmpStr;
   } /*Loop: convert bases to lookup table values*/

   return;
} /*lookupIndexToSeq*/
