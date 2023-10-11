/*#########################################################
# Name: scoreST
# Use:
#  - Hols the score structure,which holds the index and
#    score of a single cell in a Needleman-Wunsch or
#    Smith-Waterman direction matrix
# Libraries:
# C Standard Libraries:
#  - <stdlib.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o struct-01: scoresStruct
'    - Holds the score for a single direction matrix
'      positon
'  o fun-01 sortScores:
'    - Sorts an array of scores structs by score using
'      shell short.  Order is greatest to least.
'  o fun-02 swapScoreSTs:
'    - Swaps values in two score structures
'    - Macro (in scoresST.h)
'  o fun-03 initScoresST:
'    - Sets scores in a scores struture to 0
'    - Macro (in scoresST.h)
'  o fun-04 freeScoresST:
'    - Frees score structure if on heap, else does nothing
'    - Macro (in scoresST.h)
'  o fun-05 sortScoresStartAndLen:
'    - Sorts an array of scores structs by score using
'      shell short. Order is least to greatest.
'  o fun-06 cpScoreST:
'    - Copies the values from one scoreStruct to another
'      socre struct
'  o fun-07 freeScoresSTAry:
'    - Frees an array of scoresStruct's
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef SCORESST_H
#define SCORESST_H

#include <stdlib.h>

/*--------------------------------------------------------\
| Struct-01: scoresStruct
|  - Holds the score for a single direction matrix positon
\--------------------------------------------------------*/
typedef struct scoresStruct
{ /*scoresStruct*/
  unsigned long refStartUL;/*1st ref base in alignment*/
  unsigned long refEndUL;  /*last ref base in alignment*/
  unsigned long qryStartUL;/*1st query base in alignment*/
  unsigned long qryEndUL;  /*last query base in alingment*/

  long scoreL;             /*Score of alignment*/
}scoresStruct;

/*--------------------------------------------------------\
| Outptu:
|  - Modifies
|    o sorts the scoresST array between the index of the
|      first element to the index of the last element
\--------------------------------------------------------*/
void sortScores(
   struct scoresStruct **scoresST, /*scores to sort*/
   unsigned long firstElmUL,       /*1st element to sort*/
   unsigned long lastElmUL         /*last element to sort*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: sortScores
   '  - Sorts an array of scores structs by score using
   '    shell short.  Order is greatest to least.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-01 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: swapScoreSTs
| Call: swapScoreSTs(firstScore, secondScore)
| Fun-02 TOC:
| Use:
|  - Swap values in firstScore and secondScore
| Input:
|  - firstScore
|    o A pointer to a scoresStruct to swap values in
|  - secondScore
|    o A pointer to a scoresStruct to swap values in
| Output:
|  - Modifies:
|    o firstScore to have the values in secondScore
|    o secondScore to have the values in firstScore
\--------------------------------------------------------*/
#define swapScoreSTs(firstScore, secondScore){ \
  unsigned long swapUL = 0; \
  \
  swapUL = (firstScore)->refStartUL; \
  (firstScore)->refStartUL = (secondScore)->refStartUL; \
  (secondScore)->refStartUL = swapUL; \
  \
  swapUL = (firstScore)->refEndUL; \
  (firstScore)->refEndUL = (secondScore)->refEndUL; \
  (secondScore)->refEndUL = swapUL; \
  \
  swapUL = (firstScore)->qryStartUL; \
  (firstScore)->qryStartUL = (secondScore)->qryStartUL; \
  (secondScore)->qryStartUL = swapUL; \
  \
  swapUL = (firstScore)->qryEndUL; \
  (firstScore)->qryEndUL = (secondScore)->qryEndUL; \
  (secondScore)->qryEndUL = swapUL; \
  \
  swapUL = (unsigned long) (firstScore)->scoreL; \
  (firstScore)->scoreL = (secondScore)->scoreL; \
  (secondScore)->scoreL = (long) swapUL; \
} /*swapScoreSTs*/

/*--------------------------------------------------------\
| Name: initScoresST
| Call: initScoresST(scoreST)
| Fun-03 TOC:
| Use:
|  - Initializes (sets to 0) all values in a scoresStruct
| Input:
|  - scoreST:
|    o A scoresStruct to initialize
| Output:
|  - Modifies:
|    o scoreST to have values set to 0
\--------------------------------------------------------*/
#define initScoresST(scoreST){\
   (scoreST)->refStartUL = 0; \
   (scoreST)->refEndUL = 0; \
   (scoreST)->qryStartUL = 0; \
   (scoreST)->qryEndUL = 0; \
   (scoreST)->scoreL = 0; \
} /*initScoresST*/

/*--------------------------------------------------------\
| Name: freeScoresST
| Fun-04 TOC:
| Call: freeScoresST(scoresPtrST, stackBl)
| Use:
|  - Frees score structure and any variables insdie the
|    structure
| Input:
|  - scoresPtrST:
|    o Pointer to scoresST structer to free
|  - stackBl:
|    o Tells if freeing a heap or stack allocated variable
|    o 1: Freeing a heap variable
|    o 0: Freeing a stack variable (do nothing)
| Output:
|  - Frees
|    o 1: frees the score structure
|    o 0: Does nothing currently
\--------------------------------------------------------*/
#define freeScoresST(scoresPtrST, heapBl){ \
  if((heapBl)) free((scoresPtrST)); \
} /*freeScoresST*/

/*--------------------------------------------------------\
| Name: sortScoresStartAndLen
| Call: sortScoresStartAndLen(&scoresArray, 1, arrayLength)
| Use:
|  - Sorts the array of scoresStructs by the starting
|    position of the reference and length.
| Input:
|  - scoresST:
|    o Array of positions (in scoresStructs) to sort
|    o Each struct has the reference start and end
|  - firstElmUL:
|    o 1st element to sort in scoresST
|  - lastElmUL:
|    o last element to sort in scoresST
| Output:
|  - Modifies
|    o sorts the scoresST array by reference postion first
|      and alignment length second.
\--------------------------------------------------------*/
void sortScoresStartAndLen(
   struct scoresStruct **scoresST, /*scores to sort*/
   unsigned long firstElmUL,       /*1st element to sort*/
   unsigned long lastElmUL         /*last element to sort*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: sortScoresStartAndLen
   '  - Sorts an array of scores structs by score using
   '    shell short. Order is least to greatest.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-01 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: cpScoreST
| Call: cpScoreST(pasteScore, cpScore)
| Fun-06 TOC:
| Use:
|  - Copies the values from one scoreStruct to another
|    socre struct
| Input:
|  - pasteScore:
|    o Pointer to struct to copy scores to
|  - cpScore
|    o Pointer to struct to copy scores from
| Output:
|  - Modifies:
|    o pasteScore to have the values in cpScore
\--------------------------------------------------------*/
#define cpScoreST(pasteScore, cpScore){ \
  (pasteScore)->refStartUL = (cpScore)->refStartUL; \
  (pasteScore)->refEndUL = (cpScore)->refEndUL; \
  (pasteScore)->qryStartUL = (cpScore)->qryStartUL; \
  (pasteScore)->qryEndUL = (cpScore)->qryEndUL; \
  (pasteScore)->scoreL = (cpScore)->scoreL; \
} /*swapScoreSTs*/

/*--------------------------------------------------------\
| Name: freeScoresSTAry
| Call: freeScoreSTAry(scoreAryST, lenAryUL, 0)
| Fun-07 TOC:
| Use:
|  - Frees an array of score structures
| Input:
|  - scoreAryST:
|    o Start of an array (not list) of scoreStructs to free
|  - LenAryUL
|    o Number of scoreStruts in the array
|  - onStackBl:
|    o 1: Treat the array as if on stack. Free pointers in
|      structer, but not structure
|    o 0: Treat array as if on heap (free array)
| Output:
|  - Frees:
|    o scoreAryST
\--------------------------------------------------------*/
static inline void freeScoresSTAry(
   struct scoresStruct *scoreAryST,
   unsigned long lenAryUL,
   char onStackBl
){if(!onStackBl) free(scoreAryST); return;}

#endif
