/*#########################################################
# Name: scoreST
# Use:
#  - Hols the score structure,which holds the index and
#    score of a single cell in a Needleman-Wunsch or
#    Smith-Waterman direction matrix
# Libraries:
# C Standard Libraries:
#########################################################*/

#include "scoresST.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
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
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: sortScores
   '  - Sorts an array of scores structs by score using
   '    shell short. Order is greatest to least.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-01 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*Number of elements to sort*/
  unsigned long numElmUL = lastElmUL - firstElmUL;

  /*Number of sorting rounds*/
  unsigned long rndsUL = (numElmUL / 3);
  unsigned long subAryUL[rndsUL];
    /*Array of sizes for the sub arrays in shell sort*/

  unsigned long lastScoreUL = 0;
  unsigned long tmpScoreUL = 0;
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-01 Sec-02:
  ^  - Build the sub array sizes for each round of shell
  ^     sort (sub array)
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] + 1*/
  subAryUL[0] = 1; /*Initialzie first array*/

  for(unsigned long cntUL = 1; cntUL < rndsUL; ++cntUL)
    subAryUL[cntUL] = (3 * subAryUL[cntUL - 1]) + 1;

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-01 Sec-03:
  ^  - Sort the scores structure array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  while(rndsUL > 0)
  { /*loop trhough all sub arrays sort the subarrays*/
    --rndsUL; /*Move to the next subarray round*/
      /*Account for index 1*/

    for(
      unsigned long elmUL = 0;
      elmUL < subAryUL[rndsUL];
      ++elmUL
    ) { /*For each element in the subarray*/

      for(
        unsigned long subAryOnUL = firstElmUL + elmUL;
        subAryOnUL + subAryUL[rndsUL] <= lastElmUL;
        subAryOnUL += subAryUL[rndsUL]
      ) { /*Loop; swap each nth element of the subarray*/

        if(
          (*(scoresST + subAryOnUL))->scoreL <
          (*(scoresST+subAryOnUL+subAryUL[rndsUL]))->scoreL
        ) { /*If I need to swap elements*/

          swapScoreSTs(
            *(scoresST + subAryOnUL),
            *(scoresST + subAryOnUL + subAryUL[rndsUL])
          ); /*Swap scores around*/

          lastScoreUL = subAryOnUL;
          tmpScoreUL = lastScoreUL - subAryUL[rndsUL];

          while(
            tmpScoreUL >= firstElmUL &&
            (*(scoresST + lastScoreUL))->scoreL >
              (*(scoresST + tmpScoreUL))->scoreL
          ) { /*loop; move swapped element back*/
            swapScoreSTs(
              *(scoresST + lastScoreUL),
              *(scoresST + lastScoreUL - subAryUL[rndsUL])
            ); /*Swap socres around*/

            tmpScoreUL -= subAryUL[rndsUL];
          } /*loop; move swapped element back*/
        } /*If I need to swap elements*/
      } /*Loop; swap each nth element of the subarray*/
    } /*For each element in the subarray*/
  } /*loop through all sub arrays to sort the subarrays*/

  return; /*Finshed sorting the array*/
} /*sortScores*/

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
|    o Frees the score struct if requested
\--------------------------------------------------------*/
