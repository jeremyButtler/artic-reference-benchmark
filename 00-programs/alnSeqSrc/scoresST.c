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
'  o fun-05 sortScoresStartAndLen:
'    - Sorts an array of scores structs by score using
'      shell short. Order is least to greatest.
'  o fun-06 cpScoreST:
'    - Copies the values from one scoreStruct to another
'      socre struct
'  o fun-07 freeScoresSTAry:
'    - Frees an array of scoresStruct's
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
        unsigned long subOnUL = firstElmUL + elmUL;
        subOnUL + subAryUL[rndsUL] <= lastElmUL;
        subOnUL += subAryUL[rndsUL]
      ) { /*Loop; swap each nth element of the subarray*/

        if(
          (*(scoresST + subOnUL))->scoreL <
          (*(scoresST+subOnUL+subAryUL[rndsUL]))->scoreL
        ) { /*If I need to swap elements*/

          swapScoreSTs(
            *(scoresST + subOnUL),
            *(scoresST + subOnUL + subAryUL[rndsUL])
          ); /*Swap scores around*/

          lastScoreUL = subOnUL;
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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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

  /*Number of elements to sort*/
  unsigned long numElmUL = lastElmUL - firstElmUL;

  /*Number of sorting rounds*/
  unsigned long rndsUL = (numElmUL / 3);
  unsigned long subAryUL[rndsUL];
    /*Array of sizes for the sub arrays in shell sort*/

  unsigned long lastScoreUL = 0;
  unsigned long tmpScoreUL = 0;
  unsigned long cmpOneUL = 0;
  unsigned long cmpTwoUL = 0;
  unsigned long lenOneUL = 0;
  unsigned long lenTwoUL = 0;
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-05 Sec-02:
  ^  - Build the sub array sizes for each round of shell
  ^     sort (sub array)
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] + 1*/
  subAryUL[0] = 1; /*Initialzie first array*/

  for(unsigned long cntUL = 1; cntUL < rndsUL; ++cntUL)
    subAryUL[cntUL] = (3 * subAryUL[cntUL - 1]) + 1;

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-05 Sec-03:
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
        unsigned long subOnUL = firstElmUL + elmUL;
        subOnUL + subAryUL[rndsUL] <= lastElmUL;
        subOnUL += subAryUL[rndsUL]
      ) { /*Loop; swap each nth element of the subarray*/

        cmpOneUL = (*(scoresST + subOnUL))->refStartUL;
        cmpTwoUL =
           (  *(scoresST+subOnUL+subAryUL[rndsUL])
           )->refStartUL;

        lenOneUL =
            (*(scoresST + subOnUL))->refEndUL - cmpOneUL;
        lenTwoUL =
           (*(scoresST+subOnUL+subAryUL[rndsUL]))->refEndUL
         - cmpTwoUL;

        if(
             cmpOneUL > cmpTwoUL
          || (cmpOneUL == cmpTwoUL && lenOneUL == lenTwoUL)
        ) { /*If I need to swap an element*/
          swapScoreSTs(
            *(scoresST + subOnUL),
            *(scoresST + subOnUL + subAryUL[rndsUL])
          ); /*Swap scores around*/

          lastScoreUL = subOnUL;
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
} /*sortScoresStartAndLen*/

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

