/*#########################################################
# Name: hirschberg
# Use:
#  - Holds functions for doing a hirschberg global
#    alignment
# Libraries:
#  - "alnStruct.h"
#  o "alnStruct.h"
#  o "generalAlnFun.h"
#  o "twoBitArrays.h"
#  o "scoresST.h"
#  o "seqStruct.h"
#  o "alnSetStruct.h"
#  o "alnSeqDefaults.h"
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdint.h>
#  o <stdio.h>  // by alnSetStructure.h
#  o <string.h>
#########################################################*/

#include "hirschberg.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  o fun-01 Hirschberg:
'    - Sets up for and calls the recursvie function to
'      run a Hirschberg alignment
'  o fun-02 HirschbergFun:
'    - Does the recursive part of a Hirschberg alignment
'  o fun-03 scoreForwardHirsch:
'    - Does a single round of scoring for a hirschberg
'      alignment (forward direction)
'    - INLING HAS NO EFFECT
'  o fun-04 scoreReverseHirsch:
'    - Does a single round of scoring for a hirschberg
'      alignment (reverse direction)
'    - INLING HAS NO EFFECT
'  o fun-05 positionSingleRefBase:
'    - Align a single reference base to a query sequence
'    - INLING HAS NO EFFECT
'  o fun-06 twoBitAlnToAlnST:
'    - Converts a two bit array with an alignment to an
'      alnStruct structure
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o A alignment structure with the alignment.
|    o 0 For memory errors
\--------------------------------------------------------*/
struct alnStruct * Hirschberg(
  struct seqStruct *refST, /*Reference sequence to align*/
  struct seqStruct *qryST, /*Qeury sequence to align*/
    /* For refST and qryST, use seqStruct->offsetUL to set
    `  the starting point for the alignmnet and
    `  seqStruct->endAlnUL to set the ending point
    */
  struct alnSet *settings /*Settings to use for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Hirschberg
   '  - Sets up for and calls the recursvie function to
   '    run a Hirschberg alignment
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Memory allocation (set up for Hirschberg)
   '  o fun-01 sec-03:
   '    - Run the hirschberg alignment
   '  o fun-01 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   unsigned long lenRefUL =
     refST->endAlnUL - refST->offsetUL + 1;
     /*+1 to convert to index 1 (values are index 0)*/
   unsigned long lenQryUL =
     qryST->endAlnUL - qryST->offsetUL + 1;
     /*+ 1 to convert to index 1 (values are index 0)*/

   long *forwardScoreRowL = 0;
   long *reverseScoreRowL = 0;

   struct alnStruct *alnST = 0;

   #if defined HIRSCHTWOBIT
      struct twoBitAry *refAln = 0;
      struct twoBitAry *qryAln = 0;
   #else
      char *refAln = 0;
      char *qryAln = 0;
    #endif

   #if defined HIRSCHTWOBIT
      struct twoBitAry *dirRow = 0;
   #else
      char *dirRow = 0;
    #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-02:
   ^   - Memory allocation (set up for Hirschberg)
   ^   o fun-01 sec-02 sub-01:
   ^     - Initalize the ouput alignment structure 
   ^   o fun-01 sec-02 sub-02:
   ^     - Initalize the scoring rows
   ^   o fun-01 sec-02 sub-03:
   ^     - Initalize the direction rows
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-02 Sub-01:
   *  - Initalize the ouput alignment structure 
   \******************************************************/

   #ifdef HIRSCHTWOBIT
      refAln = makeTwoBit(lenRefUL + 1, 0);
       /*+ 2 for index 0*/

      if(refAln == 0) return 0;

      /*Mark the end of the alignment*/
      twoBitMvXElmFromStart(refAln, lenRefUL);
      changeTwoBitElm(refAln, defEndAlnFlag);
      twoBitMvXElmFromStart(refAln, 0);

   #else
      refAln = calloc(lenRefUL + 1, sizeof(char));
      if(refAln == 0) return 0;
      *(refAln + lenRefUL) = defEndAlnFlag;
   #endif 

   #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
      dirRow = makeTwoBit(lenRefUL + 1, 0);

      if(dirRow == 0)
      { /*If I could not make another direction row*/
         freeTwoBit(refAln, 0, 0);
         return 0;
      } /*If I could not make another direction row*/

    #elif !defined NOGAPOPEN
      dirRow = calloc(lenRefUL + 1, sizeof(char));

      if(dirRow == 0)
      { /*If I could not make another direction row*/
         free(refAln);
         refAln = 0;
         return 0;
      } /*If I could not make another direction row*/
   #endif

   #ifdef HIRSCHTWOBIT
      qryAln = makeTwoBit(lenQryUL, 0);

      if(qryAln == 0)
      { /*If had a memroy allocation error*/
        freeTwoBit(refAln, 0, 0);
        freeTwoBit(dirRow, 0, 0);
        return 0;
      } /*If had a memroy allocation error*/

      /*Mark the end of the alignment*/
      twoBitMvXElmFromStart(qryAln, lenQryUL);
      changeTwoBitElm(qryAln, defEndAlnFlag);
      twoBitMvXElmFromStart(qryAln, 0);

   #else
      qryAln = calloc(lenQryUL + 1, sizeof(char));

      if(qryAln == 0)
      { /*If I could not make another direction row*/
         free(refAln);
         refAln = 0;
         #ifndef NOGAPOPEN
            free(dirRow);
            dirRow = 0;
         #endif
         return 0;
      } /*If I could not make another direction row*/
    #endif

   /******************************************************\
   * Fun-01 Sec-02 Sub-02:
   *  - Initalize the scoring rows
   \******************************************************/

   /* I am using full length arrays to make the later
   `  steps eaiser. This takes more memory, but makes life
   `  nicer
   */

   forwardScoreRowL = malloc(sizeof(long) * lenRefUL);

   if(forwardScoreRowL == 0)
   { /*If had a memory allocatoin error*/
     #ifdef HIRSCHTWOBIT
        freeTwoBit(refAln, 0, 0);
        freeTwoBit(qryAln, 0, 0);
     #else
        free(refAln);
        refAln = 0;
        free(qryAln);
        qryAln = 0;
     #endif
 
     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        freeTwoBit(dirRow, 0, 0);
     #elif !defined NOGAPOPEN
        free(dirRow);
        dirRow = 0;
     #endif

     return 0;
   } /*If had a memory allocatoin error*/

   reverseScoreRowL = malloc(sizeof(long) * lenRefUL);

   if(reverseScoreRowL == 0)
   { /*If had a memory allocatoin error*/
     #ifdef HIRSCHTWOBIT
        freeTwoBit(refAln, 0, 0);
        freeTwoBit(qryAln, 0, 0);
     #else
        free(refAln);
        refAln = 0;
        free(qryAln);
        qryAln = 0;
     #endif
 
     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        freeTwoBit(dirRow, 0, 0);
     #elif !defined NOGAPOPEN
        free(dirRow);
        dirRow = 0;
     #endif

     free(forwardScoreRowL);
     return 0;
   } /*If had a memory allocatoin error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-03:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Sening in offset values, because alignment array is
   ` sized to the alignmnet region
   */
   HirschbergFun(
     refST->seqCStr + refST->offsetUL,
     0,                /*1st reference base to align*/
     lenRefUL,         /*Length of ref region to align*/
     qryST->seqCStr + qryST->offsetUL,
     0,                /*1st query base to align*/
     lenQryUL,         /*length of query target region*/
     forwardScoreRowL, /*For scoring*/
     reverseScoreRowL, /*For scoring*/
     refAln,      /*Holds the reference alignment*/
     qryAln,      /*Holds the query alignment*/
     dirRow,      /*Direction row for thread safe scoring*/
     settings     /*Settings for the alignment*/
   );
     /*dirRow becomes a dummy variable for -DNOGAPOPEN*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-04:
   ^    - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   free(forwardScoreRowL);
   free(reverseScoreRowL);

   alnST = twoBitAlnToAlnST(refST, qryST, refAln, qryAln);

   #ifdef HIRSCHTWOBIT
      freeTwoBit(refAln, 0, 0);
      freeTwoBit(qryAln, 0, 0);
   #else
      free(refAln);
      refAln = 0;
      free(qryAln);
      qryAln = 0;
   #endif
 
   #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
      freeTwoBit(dirRow, 0, 0);
   #elif !defined NOGAPOPEN
      free(dirRow);
      dirRow = 0;
   #endif

   return alnST; /*Is 0 if twoBitAlnToAlnSt failed*/
} /*Hirschberg*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the output alignment
\--------------------------------------------------------*/
void HirschbergFun(
  char *refSeqCStr,          /*Reference sequence*/
  unsigned long refStartUL,
    /*index 0: 1st bast to algin in reference*/
  unsigned long refLenUL,
    /*index 1 length of region to Align*/

  char *qrySeqCStr,          /*Query sequence*/
  unsigned long qryStartUL,/*index 0 Starting query base*/
  unsigned long qryLenUL,/*index 1 Length of query region*/

  long *forwardScoreRowL,   /*Holds final forward row*/
  long *reverseScoreRowL,   /*For finding reverse scores*/
    /* both the forward and reverse scoring rows must be
    `  the size of the full length reference.
    */
  
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *refAlnST,/*Holds ref alignment*/
     struct twoBitAry *qryAlnST,/*Holds query alignment*/
     struct twoBitAry *dirRow,
       /*spare direction row (makes thread safe)*/
  #else
     char *refAlnST,   /*Holds output reference alignment*/
     char *qryAlnST,   /*Holds the output query alignment*/
     char *dirRow,/*direction row for thread safe scoring*/
  #endif
  struct alnSet *settings /*Settings to use for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: HirschbergFun
   '  - Does the recursive part of a Hirschberg alignment
   '  o fun-02 sec-01:
   '    - Variable declerations
   '  o fun-02 sec-02:
   '    - Check if on a leaf (final part of alignment
   '  o fun-02 sec-03:
   '    - Get scores
   '  o fun-02 sec-04:
   '    - Find the midpoint
   '  o fun-02 sec-05:
   '    - Run the next hirschberg alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint8_t bitC = 0;
   long forwardIndelColL = 0;
   long reverseIndelColL = 0;
   unsigned long midPointUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-02:
   ^  - Check if on a leaf (final part of alignment
   ^  o fun-02 sec-02 sub-01:
   ^    - Handle cases were I have just insertions
   ^  o fun-02 sec-02 sub-02:
   ^    - Handle cases were I have just deletions
   ^  o fun-02 sec-02 sub-03:
   ^    - Handle cases were I have to align last ref base
   ^  o fun-02 sec-02 sub-04:
   ^  - Handle cases were I have to align last query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-01:
   *  - Handle cases were I have just insertions
   \******************************************************/

   if(refLenUL == 0)
   { /*If all remaing bases are query insertions*/
     #ifdef HIRSCHTWOBIT
        twoBitMvXElmFromStart(qryAlnST, qryStartUL);
     #else
        qryAlnST += qryStartUL;
     #endif

     for(
       unsigned long numInsUI = 0;
       numInsUI < qryLenUL;
       ++numInsUI
     ){ /*Loop: fill in the insertions*/
       #ifdef HIRSCHTWOBIT
          changeTwoBitElm(qryAlnST, defGapFlag);
          twoBitMvToNextElm(qryAlnST);
       #else
          *qryAlnST = defGapFlag;
          ++qryAlnST;
       #endif
     } /*Loop: fill in the insertions*/

     return; /*Nothing else to do*/
   } /*If all remaing bases are query insertions*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-02:
   *  - Handle cases were I have just deletions
   \******************************************************/

   if(qryLenUL == 0)
   { /*If all remaing bases are query deletions*/
     #ifdef HIRSCHTWOBIT
        twoBitMvXElmFromStart(refAlnST, refStartUL);
     #else
        refAlnST += refStartUL;
     #endif

     for(
       unsigned long numInsUI = 0;
       numInsUI < refLenUL;
       ++numInsUI
     ){ /*Loop: fill in the insertions*/
       #ifdef HIRSCHTWOBIT
          changeTwoBitElm(refAlnST, defGapFlag);
          twoBitMvToNextElm(refAlnST);
       #else
          *refAlnST = defGapFlag;
          ++refAlnST;
       #endif
     } /*Loop: fill in the insertions*/

     return; /*Nothing else to do*/
   } /*If all remaing bases are query deletions*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-03:
   *  - Handle cases were I have to align last ref base
   \******************************************************/

   if(refLenUL == 1)
   { /*If I have to align the last reference base*/
     #ifdef HIRSCHTWOBIT
        twoBitMvXElmFromStart(qryAlnST, qryStartUL);
        twoBitMvXElmFromStart(refAlnST, refStartUL);
     #endif

     if(qryLenUL == 0)
     { /*If bases are aligned (one reference & one query)*/
        #ifdef HIRSCHTWOBIT
           changeTwoBitElm(refAlnST, defGapFlag);
        #else
           refAlnST += refStartUL;
           *refAlnST = defGapFlag;
        #endif
        return; // Finished
     } /*If bases are aligned (one reference & one query)*/

     if(qryLenUL == 1)
     { /*If bases are aligned (one reference & one query)*/
       qrySeqCStr += qryStartUL;
       refSeqCStr += refStartUL;

       if(checkIfBasesMatch(qrySeqCStr, refSeqCStr) == 1)
         bitC = defMatchFlag;

        else bitC = defSnpFlag;

        #ifdef HIRSCHTWOBIT
           changeTwoBitElm(qryAlnST, bitC);
           changeTwoBitElm(refAlnST, bitC);
        #else
           qryAlnST += qryStartUL;
           *qryAlnST = bitC;

           refAlnST += refStartUL;
           *refAlnST = bitC;
        #endif
        return; /*Finished*/
     } /*If bases are aligned (one reference & one query)*/

     positionSingleBase(
       *(refSeqCStr + refStartUL),/*ref base*/
       refStartUL,              /*Position of ref base*/
       qrySeqCStr,              /*first base of query*/
       qryStartUL,              /*positoin of query*/
       qryLenUL,                /*Length of the query*/
       refAlnST,                /*Array to hold alignment*/
       qryAlnST,                /*Array to hold alignment*/
       settings                 /*Has Scoring variables*/
     );

     return; /*This base is now aligned*/
   } /*If I have to align the last reference base*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-04:
   *  - Handle cases were I have to align last query base
   \******************************************************/

   if(qryLenUL == 1)
   { /*If I have to align the last query base*/
     #ifdef HIRSCHTWOBIT
        twoBitMvXElmFromStart(qryAlnST, qryStartUL);
        twoBitMvXElmFromStart(refAlnST, refStartUL);
     #endif

     if(refLenUL == 0)
     { /*If bases are aligned (one reference & one query)*/
        #ifdef HIRSCHTWOBIT
           changeTwoBitElm(qryAlnST, defGapFlag);
        #else
           qryAlnST += qryStartUL;
           *qryAlnST = defGapFlag;
        #endif
        return; /*Finished*/
     } /*If bases are aligned (one reference & one query)*/

     if(refLenUL == 1)
     { /*If bases are aligned (one reference & one query)*/
       qrySeqCStr += qryStartUL;
       refSeqCStr += refStartUL;

       if(checkIfBasesMatch(qrySeqCStr, refSeqCStr) == 1)
         bitC = defMatchFlag;

        else bitC = defSnpFlag;

        #ifdef HIRSCHTWOBIT
           changeTwoBitElm(qryAlnST, bitC);
           changeTwoBitElm(refAlnST, bitC);
        #else
           qryAlnST += qryStartUL;
           *qryAlnST = bitC;
           refAlnST += refStartUL;
           *refAlnST = bitC;
        #endif
        return; /*Finished*/
     } /*If bases are aligned (one reference & one query)*/

     positionSingleBase(
       *(qrySeqCStr + qryStartUL),/*ref base*/
       qryStartUL,              /*Position of ref base*/
       refSeqCStr,              /*first base of reference*/
       refStartUL,              /*positoin of query*/
       refLenUL,                /*Length of the query*/
       qryAlnST,                /*Array to hold alignment*/
       refAlnST,                /*Array to hold alignment*/
       settings                 /*Has Scoring variables*/
     );

     return; /*Finshed aligning this query base*/
   } /*If I have to align the last query base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-03:
   ^  - Get scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    forwardIndelColL = 
      scoreForwardHirsch(
        refSeqCStr,       /*Entire reference sequence*/
        refStartUL,       /*Starting base of ref target*/
        refLenUL,         /*length of ref target region*/
        qrySeqCStr,       /*Query seq with coordinates*/
        qryStartUL,       /*Starting base of query target*/
        qryLenUL / 2,     /*Length of query target region*/
        forwardScoreRowL, /*Array of scores to fill*/
        refAlnST,         /*direction row for gap extend*/
        settings          /*setttings to use*/
    ); /*Get the scores for the forward direction*/
    /*For -DNOGAPOPEN, refAlnST is ignored*/

    reverseIndelColL = 
      scoreReverseHirsch(
        refSeqCStr,         /*Entire reference sequence*/
        refStartUL,         /*Starting base of ref target*/
        refLenUL,           /*length of ref target region*/
        qrySeqCStr,         /*Query seq with coordinates*/
        qryStartUL + (qryLenUL / 2),/*new query start*/
        qryLenUL - (qryLenUL / 2),  /*New query length*/
        reverseScoreRowL,  /*Array of scores to fill*/
        dirRow,            /*direction row for gap extend*/
        settings           /* setttings to use*/
      ); /*Get the scores for the reverse direction*/
      /*For -DNOGAPOPEN, dirRow is ignored*/
      /* I can get away with queryLen/2 here, because 
      `  queryLen is index 1 and the function takes in
      `  an lenth 1 argument
      `  I made this section thread safe by using refAlnST
      `    for the forward score and dirRow for the
      `    reverse row.
      */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-04:
   ^   - Find the midpoint
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *(forwardScoreRowL + refStartUL + refLenUL - 1) +=
     reverseIndelColL;

   midPointUL = refLenUL;

   for(
     unsigned long baseUL = 0;
     baseUL < refLenUL - 1;
     ++baseUL
   ) { /*Loop; add up all scores*/
     *(forwardScoreRowL + refStartUL + baseUL) += 
       *(reverseScoreRowL + refStartUL + baseUL + 1);
       /*The reverse row is already reversed*/

     if(
       *(forwardScoreRowL + refStartUL + baseUL) >
       *(forwardScoreRowL + refStartUL + midPointUL -1)
     ) midPointUL = baseUL + 1;
   } /*Loop; add up all scores*/

   *(reverseScoreRowL + refStartUL) += forwardIndelColL;

   if(
     *(reverseScoreRowL + refStartUL) >
     *(forwardScoreRowL + refStartUL + midPointUL - 1)
   ) midPointUL = 0;


   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-05:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   HirschbergFun(
     refSeqCStr,      /*Full reference sequence*/
     refStartUL,      /*Full reference sequence*/
     midPointUL,      /*Length of new reference sequence*/
     qrySeqCStr,      /*Full query sequence*/
     qryStartUL,      /*Start of queyr target region*/
     qryLenUL / 2,    /*Length of query target region*/
     forwardScoreRowL,/*For scoring*/
     reverseScoreRowL,/*Has last line of scores*/
     refAlnST,        /*direction row for gap extend*/
     qryAlnST,        /*Holds the alignment codes*/
     dirRow,          /*For threadsafe scoreing*/
     settings         /*Settings for the alignment*/
   );

   HirschbergFun(
     refSeqCStr,              /*Full reference sequence*/
     refStartUL + midPointUL, /*New reference start*/
     refLenUL - midPointUL,   /*New reference end*/
     qrySeqCStr,              /*Full query sequence*/
     qryStartUL + (qryLenUL / 2), /*New query start*/
     qryLenUL - (qryLenUL / 2),/*New query length*/
     forwardScoreRowL,        /*For scoring*/
     reverseScoreRowL,        /*Has last line of scores*/
     refAlnST,                /*Holds reference alingment*/
     qryAlnST,                /*Holds query alingment*/
     dirRow,                  /*For threadsafe scoring*/
     settings                 /*Settings for alignment*/
   );

   return;
} /*HirschbergFun*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o The indel column score
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      Needleman Wunsch / Smith Waterman alignment
|    o dirRowSt to hold the last row of directions in a
|      Needleman Wunsch / Smith Waterman alignment
\--------------------------------------------------------*/
long scoreForwardHirsch(
  char *refSeqCStr,          /*Reference sequence*/
  unsigned long refStartUL,  /*index 0 starting ref base*/
  unsigned long refLenUL,    /*index 1 Length of target*/

  char *qrySeqCStr,          /*Query sequence*/
  unsigned long qryStartUL, /*Index 0 Starting query base*/
  unsigned long qryLenUL,    /*index 1 length of target*/

  long *scoreRowPtrL,        /*Array of scores to fill*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *dirRowST,/*direction row*/
  #else
     char *dirRowST,      /*direction row, for gap extend*/
  #endif
  struct alnSet *settings    /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: scoreForwardHirsch
   '  - Does a single round of scoring for a hirschberg
   '    alignment (forward direction)
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-03 sec-03:
   '    - Score till on the last row
   '  o fun-03 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *qryStr = 0;
   char *refStr = 0;

   long *scoreOnPtrL = scoreRowPtrL + refStartUL;
   long indelColL = 0;
   long insScoreL = 0;
   long delScoreL = 0;
   long matchScoreL = 0;
   long nextMatchScoreL = 0;

   #if !defined HIRSCHTWOBIT && !defined NOGAPOPEN
      char *dirRowStart = 0;
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Set up the first row (indel row) of scores
   ^  o fun-03 sec-02 sub-01:
   ^    - Set up the first two elements (no gat extend)
   ^  o fun-03 sec-02 sub-02:
   ^    - Set up remaing elements (indel) (uses gap extend)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-01:
   *  - Set up the first two elements (no gat extend)
   \******************************************************/

  #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
     twoBitMvXElmFromStart(dirRowST, refStartUL);
     changeTwoBitElm(dirRowST, defMvDel);
     twoBitMvToNextElm(dirRowST);
  #elif !defined NOGAPOPEN
     dirRowST += refStartUL;
     dirRowStart = dirRowST;
     *dirRowST = defMvDel;
     ++dirRowST;
  #endif

   indelColL = 0;

   #if !defined NOGAPOPEN
      *scoreOnPtrL = settings->gapOpenI;
   #else
      *scoreOnPtrL = settings->gapExtendI;
   #endif

   ++scoreOnPtrL;

   /******************************************************\
   * Fun-03 Sec-02 Sub-02:
   *  - Set up remaing elements (indels) (uses gap extend)
   \******************************************************/

   refStr = refSeqCStr + refStartUL + 1;
   while(refStr < refSeqCStr + refStartUL + refLenUL)
   { /*Loop:Set the initial blank scores*/
     *scoreOnPtrL = *(scoreOnPtrL-1) +settings->gapExtendI;

     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        changeTwoBitElm(dirRowST, defMvDel);
        twoBitMvToNextElm(dirRowST);
     #elif !defined NOGAPOPEN
        *dirRowST = defMvDel;
        ++dirRowST;
     #endif

     ++scoreOnPtrL;
     ++refStr;
   } /*Loop:Set the initial blank scores*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - Score till on the last row
   ^  o fun-04 sec-03 sub-01:
   ^    - Set up the first element (indel, has gap open)
   ^  o fun-04 sec-03 sub-02:
   ^    - Find the scores for the next row
   ^  o fun-04 sec-03 sub-03:
   ^    - Select best score for direction
   ^  o fun-04 sec-03 sub-04:
   ^    - Find the scores for the next indels
   ^  o fun-04 sec-03 sub-05:
   ^    - Find the best score for last base pair in row
   ^  o fun-04 sec-03 sub-06:
   ^    - Move to start of row (direction matrix)
   ^  o fun-04 sec-03 sub-07:
   ^    - Find the scores for the 1st baise pair in the row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-03 Sub-01:
   *  - Set up the first element (indel, has gap open)
   \******************************************************/

   /*Move back to the start of row*/
   #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
      twoBitMvXElmFromStart(dirRowST, refStartUL);
   #elif !defined NOGAPOPEN
      dirRowST = dirRowStart;
   #endif

   scoreOnPtrL = scoreRowPtrL + refStartUL;

   /*Fine the first match/snp score (first ref base)*/
   nextMatchScoreL =
     getBaseScore(
       qrySeqCStr + qryStartUL,   /*first query base*/
       refSeqCStr + refStartUL,   /*first ref base*/
       settings                   /*Has score matrix*/
   ); /*Get the score for the enxt base*/

   nextMatchScoreL += indelColL;
     /*This is here so I can overwrite the array with the
     ` new scores
     */

   #if !defined NOGAPOPEN
      indelColL += settings->gapOpenI;
   #else
      indelColL += settings->gapExtendI;
   #endif

   /*Find the first insertion and deletion scores*/
   insScoreL =  *scoreOnPtrL + settings->gapExtendI;
   delScoreL = indelColL + settings->gapExtendI;

   /******************************************************\
   * Fun-03 Sec-03 Sub-02:
   ^  - Find the scores for the next row
   \******************************************************/

   refStr = refSeqCStr + refStartUL + 1;
     /*Already found the score for the first base pair*/
   qryStr = qrySeqCStr + qryStartUL;

   while(qryStr < qrySeqCStr + qryStartUL + qryLenUL)
   { /*Loop: score all query bases (rows)*/

     while(refStr < refSeqCStr + refStartUL + refLenUL)
     { /*Loop:score all query bases (columns)*/

       /*Get the score for the next match. This allows me
       ` to overwite the last diagnol and thus, use only
       ` one row of scoring
       */

       matchScoreL = nextMatchScoreL;

       nextMatchScoreL =
         getBaseScore(
           qryStr,     /*Current query base*/
           refStr, /*next reference base*/
           settings         /*Has score matrix*/
       ); /*Get the score for the next base*/

       nextMatchScoreL += *scoreOnPtrL;

       /**************************************************\
       * Fun-03 Sec-03 Sub-03:
       *  - Select best score for direction
       \**************************************************/

       #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
          twoBitMaxScore(
            dirRowST,        // Direction matrix
            settings,        // has direction preference
            &insScoreL,      // Score for an insertion
            &matchScoreL,    // Score for an deletion
            &delScoreL,      // Score for an match/snp
            scoreOnPtrL        // Score position to update
          ); /*Update the score and direction*/
       #elif !defined NOGAPOPEN
          charMaxScore(
            dirRowST,        // Direction matrix
            settings,        // has direction preference
            &insScoreL,      // Score for an insertion
            &matchScoreL,    // Score for an deletion
            &delScoreL,      // Score for an match/snp
            scoreOnPtrL        // Score position to update
          ); /*Update the score and direction*/
       #else
          alnMaxScore(
             settings,
             &insScoreL,
             &matchScoreL,
             &delScoreL,
             scoreOnPtrL
          );
       #endif

       /**************************************************\
       * Fun-03 Sec-03 Sub-04:
       *  - Find the scores for the next indels
       \**************************************************/

       /* The deletion scores are based on the found base,
       `  So I can find the next score before moving
       `  Get the deletion score
       */
       #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
          indelScore(
             delScoreL,
             getTwoBitElm(dirRowST),
             *scoreOnPtrL,
             settings
          ); /*Macro from generalAlnFun.h*/

          twoBitMvToNextElm(dirRowST);
      #elif !defined NOGAPOPEN
         indelScore(
             delScoreL,
             *dirRowST, /*For no gap open, is a dummy*/
             *scoreOnPtrL,
             settings
          ); /*Macro from generalAlnFun.h*/

          ++dirRowST;
      #else
         delScoreL = *scoreOnPtrL + settings->gapExtendI;
      #endif

       ++scoreOnPtrL;

       /* Finding indel scores at end, so that I can keep
       `  the indel column in a separate variable
       `  Get the insertion score
       */
       #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
          indelScore(
             insScoreL,
             getTwoBitElm(dirRowST),
             *scoreOnPtrL,
             settings
          ); /*Macro from generalAlnFun.h*/
       #elif !defined NOGAPOPEN
          indelScore(
              insScoreL,
              *dirRowST, /*For no gap open, is a dummy*/
              *scoreOnPtrL,
              settings
           ); /*Macro from generalAlnFun.h*/
       #else
          insScoreL = *scoreOnPtrL + settings->gapExtendI;
       #endif

       ++refStr;
     } /*Loop:score all query bases (columns)*/

     /****************************************************\
     * Fun-03 Sec-03 Sub-05:
     *  - Find the best score for last base pair in row
     \****************************************************/

     /*Update the final score in the row*/
     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        twoBitMaxScore(
          dirRowST,
          settings,
          &insScoreL,
          &nextMatchScoreL,
          &delScoreL,
          scoreOnPtrL
        ); /*Update the score and direction*/
     #elif !defined NOGAPOPEN
        charMaxScore(
          dirRowST,
          settings,
          &insScoreL,
          &nextMatchScoreL,
          &delScoreL,
          scoreOnPtrL
        ); /*Update the score and direction*/
     #else
       alnMaxScore(
          settings,
          &insScoreL,
          &nextMatchScoreL,
          &delScoreL,
          scoreOnPtrL
        );
     #endif

     /****************************************************\
     * Fun-03 Sec-03 Sub-06:
     *  - Move to start of row (direction matrix)
     \****************************************************/

     /*Move back to the start of row*/
     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        twoBitMvXElmFromStart(dirRowST, refStartUL);
     #elif !defined NOGAPOPEN
        dirRowST = dirRowStart;
     #endif

     scoreOnPtrL = scoreRowPtrL + refStartUL;

     /****************************************************\
     * Fun-03 Sec-03 Sub-07:
     *  - Find the scores for the 1st baise pair in the row
     \****************************************************/

     /*I need to refind the insertion and deletion scores*/
     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        indelScore(
           insScoreL,
           getTwoBitElm(dirRowST),
           *scoreOnPtrL,
           settings
        ); /*Macro from generalAlnFun.h*/
     #elif !defined NOGAPOPEN
        indelScore(
            insScoreL,
            *dirRowST, /*For no gap open, is a dummy*/
            *scoreOnPtrL,
            settings
         ); /*Macro from generalAlnFun.h*/
     #else
        insScoreL = *scoreOnPtrL + settings->gapExtendI;
     #endif

     /*Move to the next base/restart reference base*/
     refStr = refSeqCStr + refStartUL + 1;
     ++qryStr;

     /*Find the first match/snp score (first ref base)*/
     nextMatchScoreL =
       getBaseScore(
         qryStr,                    /*Next query base*/
         refStr - 1,                /*1st reference base*/
         settings                   /*Has score matrix*/
     ); /*Get the score for the next base*/

     nextMatchScoreL += indelColL;
     indelColL += settings->gapExtendI;
     delScoreL = indelColL+settings->gapExtendI;
   } /*Loop: score all query bases (rows)*/
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-04:
  ^  - Clean up
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Correct for being on the next row*/
   indelColL -= settings->gapExtendI;
   return indelColL;
} /*scoreForwardHirsch*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o The indel column score
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      backwards Needleman Wunsch /Smith Waterman alignment
|    o dirRowST to hold the last row of directions in a
|      backwards Needleman Wunsch /Smith Waterman alignment
\--------------------------------------------------------*/
long scoreReverseHirsch(
  char *refSeqCStr,          /*Reference sequence*/
  unsigned long refStartUL,  /*index 0 starting ref base*/
  unsigned long refLenUL,    /*index 1 Length of target*/

  char *qrySeqCStr,          /*Query sequence*/
  unsigned long qryStartUL, /*Index 0 Starting query base*/
  unsigned long qryLenUL,    /*index 1 length of target*/

  long *scoreRowPtrL,        /*Array of scores to fill*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *dirRowST,/*direction row*/
  #else
     char *dirRowST,      /*direction row, for gap extend*/
  #endif
  struct alnSet *settings    /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: scoreReverseHirsch
   '  - Does a single round of scoring for a hirschberg
   '    alignment (reverse direction)
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-04 sec-03:
   '    - Score till on the last row
   '  o fun-04 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *qryStr = 0;
   char *refStr = 0;

   long *scoreOnPtrL = 0;
     /*refLenUL is index 1, so  need -1 to get index 0*/

   long insScoreL = 0;
   long delScoreL = 0;
   long matchScoreL = 0;
   long nextMatchScoreL = 0; /*Score for the next base*/
   long indelColL = 0;       /*Holds indel column values*/

   #if !defined HIRSCHTWOBIT && !defined NOGAPOPEN
      char *endDirRowStr = 0;
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-02:
   ^  - Set up the first row (indel row) of scores
   ^  o fun-04 sec-02 sub-01:
   ^    - Set up the first two elements (no gat extend)
   ^  o fun-04 sec-02 sub-02:
   ^    - Set up remaing elements (indel) (uses gap extend)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-04 Sec-02 Sub-01:
   *  - Set up the first two elements (no gat extend)
   \******************************************************/

   #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
     twoBitMvXElmFromStart(dirRowST,refStartUL+refLenUL-1);
     changeTwoBitElm(dirRowST, defMvDel);
     twoBitMvBackOneElm(dirRowST);
   #elif !defined NOGAPOPEN
     dirRowST += refStartUL + refLenUL - 1;
     endDirRowStr = dirRowST;
     *dirRowST = defMvDel;
     --dirRowST;
   #endif

   scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL - 1;
     /*- 1 to account for refLenUL being index 1*/
   indelColL = 0;

   #if !defined NOGAPOPEN
      *scoreOnPtrL = settings->gapOpenI;
   #else
      *scoreOnPtrL = settings->gapExtendI;
   #endif

   --scoreOnPtrL;

   /******************************************************\
   * Fun-04 Sec-02 Sub-02:
   *  - Set up remaing elements (indels) (uses gap extend)
   \******************************************************/

   /* Loop from the second to last base till the start,
   `  Starting at second to last, because already did the
   `  first base
   */
   refStr = refSeqCStr + refStartUL + refLenUL - 2;

   while(refStr > refSeqCStr + refStartUL - 1)
   { /*Loop:Set the initial blank scores*/
     *scoreOnPtrL = *(scoreOnPtrL+1) +settings->gapExtendI;

     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        changeTwoBitElm(dirRowST, defMvDel);
        twoBitMvBackOneElm(dirRowST);
     #elif !defined NOGAPOPEN
        *dirRowST = defMvDel;
        --dirRowST;
     #endif

     --scoreOnPtrL;
     --refStr;
   } /*Loop:Set the initial blank scores*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-03:
   ^  - Score till on the last row
   ^  o fun-04 sec-03 sub-01:
   ^    - Set up the first element (indel, has gap open)
   ^  o fun-04 sec-03 sub-02:
   ^    - Find the scores for the next row
   ^  o fun-04 sec-03 sub-03:
   ^    - Select best score for direction
   ^  o fun-04 sec-03 sub-04:
   ^    - Find the scores for the next indels
   ^  o fun-04 sec-03 sub-05:
   ^    - Find the best score for last base pair in row
   ^  o fun-04 sec-03 sub-06:
   ^    - Move to start of row (direction matrix)
   ^  o fun-04 sec-03 sub-07:
   ^    - Find the scores for the 1st baise pair in the row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-04 Sec-03 Sub-01:
   *  - Set up the first element (indel, has gap open)
   \******************************************************/

   /*Move back to the start of row (refLenUL is index 1)*/
   #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
     twoBitMvXElmFromStart(dirRowST,refStartUL+refLenUL-1);
   #elif !defined NOGAPOPEN
     dirRowST = endDirRowStr;
   #endif

   scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL - 1;

   /*Find the first match/snp score (first ref base)*/
   nextMatchScoreL =
     getBaseScore(
       qrySeqCStr + qryStartUL + qryLenUL - 1,
       refSeqCStr + refStartUL + refLenUL - 1,
       settings
   ); /*Get the score for the next base*/

   nextMatchScoreL += indelColL;
     /* This is here so I can overwrite the array with the
     `  new scores
     */

   /*Find the first insertion and deletion scores*/

   #if !defined NOGAPOPEN
      indelColL += settings->gapOpenI;
   #else
      indelColL += settings->gapExtendI;
   #endif

   insScoreL = indelColL + settings->gapExtendI;
   delScoreL = indelColL + settings->gapExtendI;

   /******************************************************\
   * Fun-04 Sec-03 Sub-02:
   ^  - Find the scores for the next row
   \******************************************************/

   refStr = refSeqCStr + refStartUL + refLenUL - 2;
     /*Start one base back to account for the
     ` pre-calcualted scores for the first base
     */
   qryStr = qrySeqCStr + qryStartUL + qryLenUL - 1;

   while(qryStr > qrySeqCStr + qryStartUL - 1)
   { /*Loop: score all query bases (rows)*/

     while(refStr > refSeqCStr + refStartUL - 1)
     { /*Loop:score all query bases (columns)*/

       /* Get the score for the next match. This allows me
       `  to overwite the last diagnol and thus, use only
       `  one row of scoring
       */

       matchScoreL = nextMatchScoreL;

       nextMatchScoreL =
         getBaseScore(
           qryStr,     /*Current query base*/
           refStr,     /*next reference base*/
           settings    /*Has score matrix*/
       ); /*Get the score for the enxt base*/
          /*At worst case refStr will be '\0'*/

       nextMatchScoreL += *scoreOnPtrL;

       /**************************************************\
       * Fun-04 Sec-03 Sub-03:
       *  - Select best score for direction
       \**************************************************/

       #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
          twoBitMaxScore(
            dirRowST,
            settings,
            &insScoreL,
            &matchScoreL,
            &delScoreL,
            scoreOnPtrL
          ); /*Update the score and direction*/
       #elif !defined NOGAPOPEN
          charMaxScore(
            dirRowST,        /*Direction matrix*/
            settings,        /*has direction preference*/
            &insScoreL,      /*Score for an insertion*/
            &matchScoreL,    /*Score for an deletion*/
            &delScoreL,      /*Score for an match/snp*/
            scoreOnPtrL      /*Score position to update*/
          ); /*Update the score and direction*/
       #else
          alnMaxScore(
             settings,
             &insScoreL,
             &matchScoreL,
             &delScoreL,
             scoreOnPtrL
          );
       #endif

       /**************************************************\
       * Fun-04 Sec-03 Sub-04:
       *  - Find the scores for the next indels
       \**************************************************/

       /* The deletion scores are based on the found base,
       `  So I can find the next score before moving
       `  Get the deletion score
       */
       #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
          indelScore(
             delScoreL,
             getTwoBitElm(dirRowST),
             *scoreOnPtrL,
             settings
          ); /*Macro from generalAlnFun.h*/

          twoBitMvBackOneElm(dirRowST);
       #elif !defined NOGAPOPEN
         indelScore(
             delScoreL,
             *dirRowST, /*For no gap open, is a dummy*/
             *scoreOnPtrL,
             settings
          ); /*Macro from generalAlnFun.h*/

          --dirRowST;
       #else
         delScoreL = *scoreOnPtrL + settings->gapExtendI;
      #endif

       --scoreOnPtrL;

       /* Finding indel scores at end, so that I can keep
       `  the indel column in a separate variable
       `  Get the insertion score
       */
       #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
          indelScore(
             insScoreL,
             getTwoBitElm(dirRowST),
             *scoreOnPtrL,
             settings
          ); /*Macro from generalAlnFun.h*/
       #elif !defined NOGAPOPEN
          indelScore(
              insScoreL,
              *dirRowST, /*For no gap open, is a dummy*/
              *scoreOnPtrL,
              settings
           ); /*Macro from generalAlnFun.h*/
       #else
          insScoreL = *scoreOnPtrL + settings->gapExtendI;
       #endif

       --refStr;
     } /*Loop:score all query bases (columns)*/

     /****************************************************\
     * Fun-04 Sec-03 Sub-05:
     *  - Find the best score for last base pair in row
     \****************************************************/

     /*Update the final score in the row*/
     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        twoBitMaxScore(
          dirRowST,
          settings,
          &insScoreL,
          &nextMatchScoreL,
          &delScoreL,
          scoreOnPtrL
        ); /*Update the score and direction*/
     #elif !defined NOGAPOPEN
        charMaxScore(
          dirRowST,
          settings,
          &insScoreL,
          &nextMatchScoreL,
          &delScoreL,
          scoreOnPtrL
        ); /*Update the score and direction*/
     #else
       alnMaxScore(
          settings,
          &insScoreL,
          &nextMatchScoreL,
          &delScoreL,
          scoreOnPtrL
        );
     #endif

     /****************************************************\
     * Fun-04 Sec-03 Sub-06:
     *  - Move to start of row (direction matrix)
     \****************************************************/

     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
       twoBitMvXElmFromStart(
          dirRowST,
          refStartUL + refLenUL - 1
       );
     #elif !defined NOGAPOPEN
        dirRowST = endDirRowStr;
     #endif

     scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL -1;

     /****************************************************\
     * Fun-04 Sec-03 Sub-07:
     *  - Find the scores for the 1st baise pair in the row
     \****************************************************/

     /*I need to refind the insertion and deletion scores*/
     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        indelScore(
           insScoreL,
           getTwoBitElm(dirRowST),
           *scoreOnPtrL,
           settings
        ); /*Macro from generalAlnFun.h*/
     #elif !defined NOGAPOPEN
        indelScore(
            insScoreL,
            *dirRowST, /*For no gap open, is a dummy*/
            *scoreOnPtrL,
            settings
         ); /*Macro from generalAlnFun.h*/
     #else
        insScoreL = *scoreOnPtrL + settings->gapExtendI;
     #endif

     /*Reset the base on*/
     refStr = refSeqCStr + refStartUL + refLenUL - 2;
     --qryStr;

     /*Find the first match/snp score (first ref base)*/
     nextMatchScoreL =
       getBaseScore(
         qryStr,             /*Next query base*/
         refStr + 1,
         settings            /*Has score matrix*/
     ); /*Get the score for the enxt base*/

     nextMatchScoreL += indelColL;
     indelColL += settings->gapExtendI;
     delScoreL = indelColL+settings->gapExtendI;
   } /*Loop: score all query bases (rows)*/
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-04 Sec-04:
  ^  - Clean up
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Correct for being on the next row*/
   indelColL -= settings->gapExtendI;
   return indelColL;
} /*scoreReverseHirsch*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the alignment for the single
|      base aligned to the sequence
\--------------------------------------------------------*/
void positionSingleBase(
  char baseC,             /*Single base to align to a seq*/
  unsigned long baseIndexUL,/*Index base is at*/
  char *seqCStr,            /*Sequence position on*/
  unsigned long startOfSeqUL,
    /*Index 0 of first base to align bascC to in seqCStr*/
  unsigned long lenSeqUL,
    /*Index 1; Length of the aligned region in seqCStr*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *baseCAlnST,
       /* Two bit alingment array for the sequence having
       `  baseC
       */
     struct twoBitAry *seqAlnST, 
       /* Two bit alignment array for the sequence alinging
       `  baseC to
       */
  #else
     char *baseCAlnST,
     char *seqAlnST, 
  #endif
  struct alnSet *settings        /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: positionSingleRefBase
   '  - Align a single base to a sequence
   '  o fun-05 sec-01:
   '    - Variable declerations
   '  o fun-05 sec-02:
   '    - Find the reference bases position on the query
   '  o fun-05 sec-03:
   '    - Fill in insertions and reference base position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   long insScoreL = 0;      /*The score of an insertion*/
   long delScoreL = 0;      /*the score of a deleton*/
   long matchScoreL = 0;    /*The score of a match*/
   long curScoreL = 0;      /*The score of the last row*/

   char *seqBaseCStr = seqCStr + startOfSeqUL;
   char *endSeqCStr = seqCStr +startOfSeqUL +lenSeqUL-1;
   char *lastMatchCStr = 0;
   char dirC = 0; /*Removes some #if defined statments*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-02:
   ^  - Find the reference bases position on the query
   ^  o fun-05 sec-02 sub-01:
   ^    - Find the first scores for the loop
   ^  o fun-05 sec-02 sub-02:
   ^    - Find the remaing scores
   ^  o fun-05 sec-02 sub-03:
   ^     - Figure out which of hte final scores to keep
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-05 Sec-02 Sub-01:
   *  - Find the first scores for the loop
   \******************************************************/

   #if defined HIRSCHTWOBIT
      twoBitMvXElmFromStart(baseCAlnST, baseIndexUL);
      twoBitMvXElmFromStart(seqAlnST, startOfSeqUL);
   #else
      baseCAlnST += baseIndexUL;
      seqAlnST += startOfSeqUL;
   #endif

   matchScoreL = getBaseScore(seqBaseCStr,&baseC,settings);

   #if !defined NOGAPOPEN
      insScoreL = settings->gapOpenI;
   #else
      insScoreL = settings->gapExtendI;
   #endif

   delScoreL = insScoreL;
   ++seqBaseCStr;

   /******************************************************\
   * Fun-05 Sec-02 Sub-02:
   *  - Find the remaing scores
   \******************************************************/
 
   do { /*While I have bases to compare*/
     charMaxScore(
       &dirC,
       settings,        /*has direction preference*/
       &insScoreL,      /*Score for an insertion*/
       &matchScoreL,    /*Score for an deletion*/
       &delScoreL,      /*Score for an match/snp*/
       &curScoreL       /*Score position to update*/
     ); /*Update the score*/
     /*Directionalty determines which direction to select
     ` when one or more directions are equal (ins, snp,del)
     */

     matchScoreL =
         insScoreL
       + getBaseScore(seqBaseCStr,&baseC,settings);

     insScoreL += settings->gapExtendI;

     switch(dirC)
     { /*Switch: Check if keeping the score*/
       case defMvSnp:
         lastMatchCStr = seqBaseCStr - 1;
         #if !defined NOGAPOPEN
            delScoreL = curScoreL + settings->gapOpenI;
         #else
            delScoreL = curScoreL + settings->gapExtendI;
         #endif
         break;

       case defMvDel:
       case defMvIns:
         delScoreL = curScoreL + settings->gapExtendI;
         break;
       case defMvStop: break; /*Never fires*/
     } /*Switch: Check if keeping the score*/

     ++seqBaseCStr;
   } while(seqBaseCStr <= endSeqCStr);
   /*While I have bases to compare*/

   /******************************************************\
   * Fun-05 Sec-02 Sub-03:
   *  - Figure out which of the final scores to keep
   \******************************************************/

   charMaxScore(
     &dirC,
     settings,        /*has direction preference*/
     &insScoreL,      /*Score for an insertion*/
     &matchScoreL,    /*Score for an deletion*/
     &delScoreL,      /*Score for an match/snp*/
     &curScoreL       /*Score position to update*/
   ); /*Update the final score*/
   /*Directionalty determines which direction to select
   ` when one or more directions are equal (ins, snp,del)
   */

   switch(dirC)
   { /*Switch: Check if keeping the score*/
     case defMvSnp:
       lastMatchCStr = seqBaseCStr - 1;
       break;

     case defMvDel: break;
     case defMvIns: break;
     case defMvStop: break; // Never fires
   } /*Switch: Check if keeping the score*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-03:
   ^  - Fill in the insertions and reference base position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   /* A series of deletions and insertions are prefered
   `  over matches and smps. In this case put the base
   `  at the first base. There is no good position
   */
   if(lastMatchCStr == 0)
     lastMatchCStr = seqCStr + startOfSeqUL;

   seqBaseCStr = seqCStr + startOfSeqUL;
   /* No need to change query position since previouis
   `  loop only usedo one direction position
   */

   /*Add in the insertions at the start*/
   while(seqBaseCStr < lastMatchCStr)
   { /*While I have insertions to fill*/
     #ifdef HIRSCHTWOBIT
        changeTwoBitElm(seqAlnST, defGapFlag);
        twoBitMvToNextElm(seqAlnST);
     #else
        *seqAlnST = defGapFlag;
        ++seqAlnST;
     #endif
     ++seqBaseCStr;
   } /*While I have insertions to fill*/
   
   /*Add in the position of the base*/
   if(checkIfBasesMatch(seqBaseCStr, &baseC))
   { /*IF the bases matched*/
     #ifdef HIRSCHTWOBIT
        changeTwoBitElm(baseCAlnST, defMatchFlag);
        changeTwoBitElm(seqAlnST, defMatchFlag);
     #else
        *baseCAlnST = defMatchFlag;
        *seqAlnST = defMatchFlag;
     #endif
   } /*IF the bases matched*/

    else
    { /* Else this was a SNP*/
     #ifdef HIRSCHTWOBIT
        changeTwoBitElm(baseCAlnST, defSnpFlag);
        changeTwoBitElm(seqAlnST, defSnpFlag);
     #else
        *baseCAlnST = defSnpFlag;
        *seqAlnST = defSnpFlag;
     #endif
    } /*Else this was a SNPS*/

   #ifdef HIRSCHTWOBIT
      twoBitMvToNextElm(seqAlnST);
   #else
      ++seqAlnST;
   #endif

   ++seqBaseCStr;

   /*Finish adding in the insertions at the end*/
   while(seqBaseCStr <= endSeqCStr)
   { /*While I have insertions to fill*/
     #ifdef HIRSCHTWOBIT
        changeTwoBitElm(seqAlnST, defGapFlag);
        twoBitMvToNextElm(seqAlnST);
     #else
        *seqAlnST = defGapFlag;
        ++seqAlnST;
     #endif
     ++seqBaseCStr;
   } /*While I have insertions to fill*/

   return;
} /*positionSingleRefBase*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o 0: for error
|    o pointer to alnStruct with alignment
\--------------------------------------------------------*/
struct alnStruct * twoBitAlnToAlnST(
  struct seqStruct *refST,
   /*Has reference alignment start/end & reference length*/
  struct seqStruct *qryST,
   /*Has query alignment start/end and query length*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *refAlignment,
       /*Two bit array with the reference alignment*/
     struct twoBitAry *qryAlignment
       /*Two bit array with the query alignment*/
  #else
     char *refAlignment, /*has reference alignment*/
     char *qryAlignment /*has query alignment*/
  #endif
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: twoBitAlnToAlnST
   '  - Converts a two bit array with an alignment to an
   '    alnStruct structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   char *refAlnStr = 0;
   char *qryAlnStr = 0;
   uint8_t bitUC = 0;

   long refIndexUL = 0;
   long qryIndexUL = 0;

   long refFirstAlnBaseL = -1;
   long refLastAlnBaseL = 0;
   long qryFirstAlnBaseL = -1;
   long qryLastAlnBaseL = 0;

   long numDelsL = 0;
   long numInssL = 0;
   long numSnpsL = 0;
   long numMatchesL = 0;

   struct alnStruct *alnST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-02:
   ^  - Allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   alnST = malloc(sizeof(alnStruct));
   if(alnST == 0) return 0; /*Memory error*/
   initAlnST(alnST);

   refAlnStr = calloc(refST->lenSeqUL + 1, sizeof(char));

   if(refAlnStr == 0)
   { /*If had a memor error*/
      freeAlnST(alnST, 1); /*1 for freeing heap memory*/
      alnST = 0;
      return 0;
   } /*If had a memor error*/

   alnST->refAlnStr = refAlnStr;

   qryAlnStr = calloc(qryST->lenSeqUL + 1, sizeof(char));

   if(qryAlnStr == 0)
   { // IF had a memory allocation error
     freeAlnST(alnST, 1); /*1 for freeing heap memory*/
     alnST = 0;
     refAlnStr = 0;
     return 0;
   } // IF had a memory allocation error

   alnST->qryAlnStr = qryAlnStr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-03:
   ^  - Add in the alignment
   ^  o fun-06 sec-03 sub-01:
   ^    - Add in the reference alignment
   ^  o fun-06 sec-03 sub-02:
   ^    - Add in the query alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Add in the reference alignment
   \******************************************************/

   #ifdef HIRSCHTWOBIT
      twoBitMvXElmFromStart(refAlignment, 0);
   #endif

   /*the refAlignment array starts at refST->offsetUL*/
   refIndexUL = refST->offsetUL;
   refAlnStr += refIndexUL;

   #ifdef HIRSCHTWOBIT
      bitUC = getTwoBitElm(refAlignment);
   #else
      bitUC = (uint8_t) *refAlignment;
   #endif

   while(bitUC != defEndAlnFlag)
   { /*Loop: Add reference aligned bases to alnStruct*/
      switch(bitUC)
      { /*Switch: Check the error type*/
         case 0: break;

         case defSnpFlag:
         /*Case: snp*/
            ++numSnpsL;

            if(refFirstAlnBaseL < 0) /*Start of alignment*/
               refFirstAlnBaseL = refIndexUL;

            refLastAlnBaseL=refIndexUL;/*end of alignment*/
            break;
         /*Case: snp*/

         case defMatchFlag:
         /*Case: matches*/
            ++numMatchesL;

            if(refFirstAlnBaseL < 0) /*Start of alignment*/
               refFirstAlnBaseL = refIndexUL;

            refLastAlnBaseL=refIndexUL;/*end of alignment*/
            break;
         /*Case: matches*/

         case defGapFlag:
         /*Case: Deletions*/
            ++numDelsL;
            break;
         /*Case: Deletions*/
      } /*Switch: Check the error type*/

      *refAlnStr = bitUC;

      #ifdef HIRSCHTWOBIT
         twoBitMvToNextElm(refAlignment);
         bitUC = getTwoBitElm(refAlignment);
      #else
         ++refAlignment;
         bitUC = (uint8_t) *refAlignment;
      #endif

      ++refAlnStr;
      ++refIndexUL;
   } /*Loop: Add reference aligned bases to alnStruct*/

   if(refFirstAlnBaseL >= 0) /*Start of alignment*/
      alnST->refStartAlnUL = refFirstAlnBaseL;
   else alnST->refStartAlnUL = refST->lenSeqUL;
     /*One base of last base*/

   if(refLastAlnBaseL > 0)
      alnST->refEndAlnUL = refLastAlnBaseL;
   else alnST->refEndAlnUL = refST->lenSeqUL;
     /*One base of last base*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Add in the query alignment
   \******************************************************/

   #ifdef HIRSCHTWOBIT
      twoBitMvXElmFromStart(qryAlignment, 0);
   #endif
      
   /*The qryAlignment array starts at qryST->offset*/
   qryIndexUL = qryST->offsetUL;
   qryAlnStr += qryIndexUL;

   #ifdef HIRSCHTWOBIT
      bitUC = getTwoBitElm(qryAlignment);
   #else
      bitUC = (uint8_t) *qryAlignment;
   #endif

   while(bitUC != defEndAlnFlag)
   { /*Loop: Add query aligned bases to alnStruct*/
      switch(bitUC)
      { /*Switch: Check the error type*/
         case 0: break;

         case defSnpFlag:
         /*Case: snp*/
            if(qryFirstAlnBaseL < 0) /*Start of alignment*/
               qryFirstAlnBaseL = qryIndexUL;

            qryLastAlnBaseL=qryIndexUL;/*end of alignment*/
            break;
         /*Case: snp*/

         case defMatchFlag:
         /*Case: matches*/
            if(qryFirstAlnBaseL < 0) /*Start of alignment*/
               qryFirstAlnBaseL = qryIndexUL;

            qryLastAlnBaseL=qryIndexUL;/*end of alignment*/
            break;
         /*Case: matches*/

         case defGapFlag:
         /*Case: Insertions*/
            ++numInssL;
            break;
         /*Case: inerstions*/
      } /*Switch: Check the error type*/

      *qryAlnStr = bitUC;

      #ifdef HIRSCHTWOBIT
         twoBitMvToNextElm(qryAlignment);
         bitUC = getTwoBitElm(qryAlignment);
      #else
         ++qryAlignment;
         bitUC = (uint8_t) *qryAlignment;
      #endif

      ++qryIndexUL;
      ++qryAlnStr;
   } /*Loop: Add query aligned bases to alnStruct*/

   if(qryFirstAlnBaseL >= 0) /*Start of alignment*/
      alnST->qryStartAlnUL = qryFirstAlnBaseL;
   else alnST->qryStartAlnUL = qryST->lenSeqUL;
     /*One base of last base*/

   if(qryLastAlnBaseL > 0)
      alnST->qryEndAlnUL = qryLastAlnBaseL;
   else alnST->qryEndAlnUL = qryST->lenSeqUL;
     /*One base of last base*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Add in the alignment stats
   \******************************************************/

   alnST->numInssUL = (unsigned long) numInssL;
   alnST->numDelsUL = (unsigned long) numDelsL;
   alnST->numSnpsUL = (unsigned long) numSnpsL;
   alnST->numMatchesUL = (unsigned long) numMatchesL;

   alnST->lenAlnUL =
      (unsigned long)
      (numInssL + numDelsL + numSnpsL + numMatchesL);

   alnST->refLenUL = refST->lenSeqUL;
   alnST->qryLenUL = qryST->lenSeqUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-04:
   ^  - Add in softmasking
   ^  o fun-06 sec-03 sub-01:
   ^    - Add soft masking to the un-aligned ending bases
   ^  o fun-06 sec-03 sub-02:
   ^    - Add soft masking to the un-aligned starting bases
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Add soft masking to the un-aligned ending bases
   \******************************************************/

    while(refIndexUL < refST->lenSeqUL)
    { /*Loop: Apply mask to starting reference bases*/
       *refAlnStr = defSoftMaskFlag;
       ++refAlnStr;
       ++refIndexUL;
    } /*Loop: Apply mask to starting reference bases*/

    while(qryIndexUL < qryST->lenSeqUL)
    { /*Loop: Apply mask to starting query bases*/
       *qryAlnStr = defSoftMaskFlag;
       ++qryAlnStr;
       ++qryIndexUL;
    } /*Loop: Apply mask to starting query bases*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-02:
   *  - Add soft masking to the un-aligned starting bases
   \******************************************************/

    refAlnStr = alnST->refAlnStr; 
    while(*refAlnStr == 0)
    { /*Loop: Apply mask to starting reference bases*/
       *refAlnStr = defSoftMaskFlag;
       ++refAlnStr;
    } /*Loop: Apply mask to starting reference bases*/

    qryAlnStr = alnST->qryAlnStr; 
    while(*qryAlnStr == 0)
    { /*Loop: Apply mask to starting query bases*/
       *qryAlnStr = defSoftMaskFlag;
       ++qryAlnStr;
    } /*Loop: Apply mask to starting query bases*/

    return alnST;
} /*twoBitAlnToAlnST*/
