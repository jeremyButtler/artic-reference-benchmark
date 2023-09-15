/*#########################################################
# Name: scoreST
# Use:
#  - Hols the score structure,which holds the index and
#    score of a single cell in a Needleman-Wunsch or
#    Smith-Waterman direction matrix
# Libraries:
#  - dataTypeShortHand.h"
# C Standard Libraries:
#  - <stdlib.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o struct-01: scoresStruct
'    - Holds the score for a single direction matrix
'      positon
'  o fun-01 initScoresST:
'    - Sets scores in a scores struture to 0
'    - Macro (in scoresST.h)
'  o fun-02 freeScoresST:
'    - Frees score structure if on heap, else does nothing
'    - Macro (in scoresST.h)
'  o fun-03 freeScoresSTAry:
'    - Frees an array of scoresStruct's
'  o fun-04 swapScoreSTs:
'    - Swaps values in two score structures
'    - Macro (in scoresST.h)
'  o fun-05 cpScoreST:
'    - Copies the values from one scoreStruct to another
'      socre struct
'  o fun-06 sortScores:
'    - Sorts an array of scores structs by score using
'      shell short.  Order is greatest to least.
'  o fun-07 sortScoresStartLen:
'    - Sorts an array of scores structs by score using
'      shell short. Order is least to greatest.
'  o fun-08 sortScoresStartLenIndex:
'    - Gets the number of rounds to run shell sort for
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef SCORESST_H
#define SCORESST_H

#include <stdlib.h>
#include "dataTypeShortHand.h"
  /*I am using 59 char line wraps and the unsigned takes
  ` 9 characters. So, I need something shorter. In this
  ' case I renamed unsigned variables to: ulonglong, ulong,
  ' uint, ushort, and uchar.
  */

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
| Name: initScoresST
| Call: initScoresST(scoreST)
| Fun-01 TOC:
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
| Fun-02 TOC:
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
| Name: freeScoresSTAry
| Call: freeScoreSTAry(scoreAryST, lenAryUL, 0)
| Fun-03 TOC:
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
   ulong lenAryUL,
   char onStackBl
){if(!onStackBl) free(scoreAryST); return;}

/*--------------------------------------------------------\
| Name: swapScoreSTs
| Call: swapScoreSTs(firstScore, secondScore)
| Fun-04 TOC:
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
  ulong swapUL = 0; \
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
  swapUL = (ulong) (firstScore)->scoreL; \
  (firstScore)->scoreL = (secondScore)->scoreL; \
  (secondScore)->scoreL = (long) swapUL; \
} /*swapScoreSTs*/

/*--------------------------------------------------------\
| Name: cpScoreST
| Call: cpScoreST(pasteScore, cpScore)
| Fun-05 TOC:
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
| Output:
|  - Modifies
|    o sorts the scoresST array between the index of the
|      first element to the index of the last element
\--------------------------------------------------------*/
static inline void sortScores(
   struct scoresStruct **scoresST, /*scores to sort*/
   ulong firstElmUL,       /*1st element to sort*/
   ulong endElmUL         /*last element to sort*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: sortScores
   '  - Sorts an array of scores structs by score using
   '    shell short.  Order is greatest to least.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-06 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*Number of elements to sort*/
  ulong numElmUL = (ulong) endElmUL - (ulong) firstElmUL;

  /*Number of sorting rounds*/
  ulong subUL = 0;
  ulong nextElmUL = 0;
  ulong lastElmUL = 0;
    /*Array of sizes for the sub arrays in shell sort*/
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-06 Sec-02:
  ^  - Find the max search value
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] + 1*/
  subUL = 1; /*Initialzie first array*/
  while(subUL < numElmUL) subUL = (3 * subUL) + 1;
  subUL = (subUL - 1) / 3; /*Account for overshooting*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-06 Sec-03:
  ^  - Sort the scores structure array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  while(subUL > 0)
  { /*loop trhough all sub arrays sort the subarrays*/
    for(ulong indexUL = 0; indexUL <= subUL; ++indexUL)
    { /*For each element in the subarray*/
      for(ulong elmUL = indexUL;
          elmUL + subUL <= (ulong) endElmUL;
          elmUL += subUL
      ) { /*Loop; swap each nth element of the subarray*/
        nextElmUL = elmUL + subUL;

        if(
              (*scoresST + elmUL)->scoreL
            > (*scoresST + nextElmUL)->scoreL
        ) { /*If I need to swap an element*/
          swapScoreSTs(
            (*scoresST + elmUL),
            (*scoresST + nextElmUL)
          ); /*Swap scores around*/

          nextElmUL = elmUL;
          lastElmUL = elmUL;

          while(lastElmUL >= subUL)
          { /*loop; move swapped element back*/
             lastElmUL -= subUL;

             if(
                 (*scoresST + nextElmUL)->scoreL
               < (*scoresST + lastElmUL)->scoreL
             ) break; /*If this element is positioned*/

             swapScoreSTs(
               (*scoresST + nextElmUL),
               (*scoresST + lastElmUL)
             ); /*Swap scores around*/

             nextElmUL = lastElmUL;
          } /*loop; move swapped element back*/
        } /*If I need to swap elements*/
      } /*Loop; swap each nth element of the subarray*/
    } /*For each element in the subarray*/

    subUL = (subUL - 1) / 3; /*Move to the next round*/
  } /*loop through all sub arrays to sort the subarrays*/

  return; /*Finshed sorting the array*/
} /*sortScores*/

/*--------------------------------------------------------\
| Name: sortScoresStartLen
| Call: sortScoresStartLen(&scoresArray, 1, arrayLength)
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
static inline void sortScoresStartLen(
   struct scoresStruct **scoresST, /*scores to sort*/
   ulong firstElmUL,       /*1st element to sort*/
   ulong endElmUL          /*last element to sort*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: sortScoresStartLen
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

  /*Number of elements to sort*/
  ulong numElmUL = endElmUL - firstElmUL;

  /*Number of sorting rounds*/
  ulong subUL = 0;
  ulong nextElmUL = 0;
  ulong lastElmUL = 0;
    /*Array of sizes for the sub arrays in shell sort*/
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-07 Sec-02:
  ^  - Find the max search value
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] + 1*/
  subUL = 1; /*Initialzie first array*/
  while(subUL < numElmUL) subUL = (3 * subUL) + 1;
  subUL = (subUL - 1) / 3; /*Account for overshooting*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-07 Sec-03:
  ^  - Sort the scores structure array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  while(subUL > 0)
  { /*loop trhough all sub arrays sort the subarrays*/
    for(ulong indexUL = 0; indexUL <= subUL; ++indexUL)
    { /*For each element in the subarray*/
      for(ulong elmUL = indexUL;
          elmUL + subUL <= endElmUL;
          elmUL += subUL
      ) { /*Loop; swap each nth element of the subarray*/
        nextElmUL = elmUL + subUL;

        if(
          (/*Check if the later position has a lower start*/ 
              (*scoresST + elmUL)->refStartUL
            > (*scoresST + nextElmUL)->refStartUL
          ) || 
          ((  /*If both scores have an equal start*/
                (*scoresST + elmUL)->refStartUL
             == (*scoresST + nextElmUL)->refStartUL
           ) &&
           (  /*Check to find the shorter lengths*/
               (*scoresST + elmUL)->refEndUL
             > (*scoresST + nextElmUL)->refEndUL
          ))
        ) { /*If I need to swap an element*/
          swapScoreSTs(
            (*scoresST + elmUL),
            (*scoresST + nextElmUL)
          ); /*Swap scores around*/

          nextElmUL = elmUL;
          lastElmUL = elmUL;

          while(lastElmUL >= subUL)
          { /*loop; move swapped element back*/
             lastElmUL -= subUL;

             if(
               (/*Check if later position has lower start*/ 
                   (*scoresST + nextElmUL)->refStartUL
                 > (*scoresST + lastElmUL)->refStartUL
               ) || 
               ((  /*If both scores have an equal start*/
                     (*scoresST + nextElmUL)->refStartUL
                  == (*scoresST + lastElmUL)->refStartUL
                ) &&
                (  /*Check to find the shorter lengths*/
                    (*scoresST + nextElmUL)->refEndUL
                  > (*scoresST + lastElmUL)->refEndUL
               ))
             ) break; /*If this element is positioned*/

             swapScoreSTs(
               (*scoresST + nextElmUL),
               (*scoresST + lastElmUL)
             ); /*Swap scores around*/

             nextElmUL = lastElmUL;
          } /*loop; move swapped element back*/
        } /*If I need to swap elements*/
      } /*Loop; swap each nth element of the subarray*/
    } /*For each element in the subarray*/

    subUL = (subUL - 1) / 3; /*Move to the next round*/
  } /*loop through all sub arrays to sort the subarrays*/

  return; /*Finshed sorting the array*/
} /*sortScoresStartLen*/

/*--------------------------------------------------------\
| Name: sortScoresStartLenIndex
| Call: sortScoresStartLenIndex(
|          &scoresArray,
|          1,
|          arrayLength,
|          indexsToMatch)
| Use:
|  - Sorts the array of scoresStructs by the starting
|    position of the reference and length.
|  - This also makes sure indexsToMatch stays in sync with
|    the sorted scores structures.
| Input:
|  - scoresST:
|    o Array of positions (in scoresStructs) to sort
|    o Each struct has the reference start and end
|  - firstElmUL:
|    o 1st element to sort in scoresST
|  - lastElmUL:
|    o last element to sort in scoresST
|  - altAryUL:
|    o Alternate array of unsigned longs (ulong) to keep
|      sorted with scoresST.
| Output:
|  - Modifies
|    o sorts the scoresST array by reference postion first
|      and alignment length second.
|    o altAryUL to stay in sync with scoresST
\--------------------------------------------------------*/
static inline void sortScoresStartLenIndex(
   struct scoresStruct **scoresST, /*scores to sort*/
   ulong firstElmUL,       /*1st element to sort*/
   ulong endElmUL,         /*last element to sort*/
   ulong *altAryUL  /*Needs to stay with scoresST*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: sortScoresStartLenIndex
   '  - Sorts an array of scores structs by score using
   '    shell short. Order is least to greatest.
   '  - This also makes sue that altAryUL stays in sync
   '    with the scoresST array.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-08 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   '  o fun-08 sec-02:
   '    - Find the max search value
   '  o fun-08 sec-03:
   '    - Sort the scores structure array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*Number of elements to sort*/
  ulong numElmUL = (ulong) endElmUL - (ulong) firstElmUL;

  /*Number of sorting rounds*/
  ulong subUL = 0;
  ulong nextElmUL = 0;
  ulong lastElmUL = 0;
  ulong swapUL = 0;
    /*Array of sizes for the sub arrays in shell sort*/
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-08 Sec-02:
  ^  - Find the max search value
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] + 1*/
  subUL = 1; /*Initialzie first array*/
  while(subUL < numElmUL) subUL = (3 * subUL) + 1;
  subUL = (subUL - 1) / 3; /*Account for overshooting*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-08 Sec-03:
  ^  - Sort the scores structure array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  while(subUL > 0)
  { /*loop trhough all sub arrays sort the subarrays*/
    for(ulong indexUL = 0; indexUL <= subUL; ++indexUL)
    { /*For each element in the subarray*/
      for(ulong elmUL = indexUL;
          elmUL + subUL <= (ulong) endElmUL;
          elmUL += subUL
      ) { /*Loop; swap each nth element of the subarray*/
        nextElmUL = elmUL + subUL;

        if(
          (/*Check if the later position has a lower start*/ 
              (*scoresST + elmUL)->refStartUL
            > (*scoresST + nextElmUL)->refStartUL
          ) || 
          ((  /*If both scores have an equal start*/
                (*scoresST + elmUL)->refStartUL
             == (*scoresST + nextElmUL)->refStartUL
           ) &&
           (  /*Check to find the shorter lengths*/
               (*scoresST + elmUL)->refEndUL
             > (*scoresST + nextElmUL)->refEndUL
          ))
        ) { /*If I need to swap an element*/
          swapScoreSTs(
            (*scoresST + elmUL),
            (*scoresST + nextElmUL)
          ); /*Swap scores around*/

          /*Keep the alternate array in sync*/
          swapUL = *(altAryUL + elmUL);
          *(altAryUL + elmUL) = *(altAryUL + nextElmUL);
          *(altAryUL + nextElmUL) = swapUL;

          nextElmUL = elmUL;
          lastElmUL = elmUL;

          while(lastElmUL >= subUL)
          { /*loop; move swapped element back*/
             lastElmUL -= subUL;

             if(
               (/*Check if later position has lower start*/ 
                   (*scoresST + nextElmUL)->refStartUL
                 > (*scoresST + lastElmUL)->refStartUL
               ) || 
               ((  /*If both scores have an equal start*/
                     (*scoresST + nextElmUL)->refStartUL
                  == (*scoresST + lastElmUL)->refStartUL
                ) &&
                (  /*Check to find the shorter lengths*/
                    (*scoresST + nextElmUL)->refEndUL
                  > (*scoresST + lastElmUL)->refEndUL
               ))
             ) break; /*If this element is positioned*/

             swapScoreSTs(
               (*scoresST + nextElmUL),
               (*scoresST + lastElmUL)
             ); /*Swap scores around*/

             /*Keep the alternate array in sync*/
             swapUL = *(altAryUL + nextElmUL);
             *(altAryUL+nextElmUL) = *(altAryUL+lastElmUL);
             *(altAryUL + lastElmUL) = swapUL;

             nextElmUL = lastElmUL;
          } /*loop; move swapped element back*/
        } /*If I need to swap elements*/
      } /*Loop; swap each nth element of the subarray*/
    } /*For each element in the subarray*/

    subUL = (subUL - 1) / 3; /*Move to the next round*/
  } /*loop through all sub arrays to sort the subarrays*/

  return; /*Finshed sorting the array*/
} /*sortScoresStartLenIndex*/


#endif
