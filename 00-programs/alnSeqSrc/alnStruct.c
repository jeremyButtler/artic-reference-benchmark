/*######################################################### # Name: alnStruct
# Use:
#  - Holds the alingment structure that stores the
#    alignment. This also includes the functions needed
#    to maintain and print out the alignment
# Libraries:
#  - "twoBitArrays.h"
#  - "seqStruct.h"
#  - "scoresST.h"
#  - "alnSetStruct.h"
#  o "alnSeqDefaults.h"
# C Standard Libraries:
#  - <time.h>
#  - <string.h>
#  o <stdlib.h>
#  o <stdio.h>
#  o <stdint.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  o fun-01 alnSTToSeq:
'    - Makes an alignment for an single sequence
'      (reference or query)
'  o fun-02 dirMatrixToAlnST:
'    - Builds an anlignment array for the input direction
'  o fun-03 printAln:
'    - Prints out an alignment
'  o fun-04 initAlnST:
'    - Initalize all values in alnST to 0
'  o fun-05 freeAlnST:
'    - Frees alnST and all variables in alnST
'  o fun-06 pEMBOSSHead:
'    - Prints out the EMBOSS header to a file
'  o fun-07 pExpandCigHead:
'    - Prints out the expanded cigar header entry
'  o fun-08 capIdLen:
'    - Caps id length in seqST first white space or the
'      character at maxIdLenI - 1.
'  o fun-09 addPosToBuff:
'    - Adds a base position to a buffer. This uses the
'      printing format and position as a guide to
'      determine if printing.
'  o fun-10 eqxAddGap
'    - Adds a gap to an eqx buffer
'  o fun-11 eqxAddSMask
'    - Adds an soft mask to an eqx buffer
'  o fun-12 eqxAddSnp
'    - Adds an SNP entry to an eqx buffer
'  o fun-13 eqxAddMatch
'    - Adds an match entry to an eqx buffer
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "alnStruct.h"

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Heap alloacted C-string with alignment for the
|      input sequence
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
char * alnSTToSeq(
    struct seqStruct *seqST,/*Has sequence to work with*/
    char qryBl,             /*1: working on query; 0: ref*/
    struct alnStruct *alnST,/*Has alignment array*/
    char extAlnRegionBl     /*1: Only keep aligned region*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: alnSTToSeq
   '  - Makes an alignment for an single sequence
   '    (reference or query)
   '  o fun-01 Sec-01:
   '    - Variable declerations 
   '  o fun-01 Sec-02:
   '    - Allocate memory & identify if ref or query seq
   '  o fun-01 Sec-03:
   '    - Add softmasking to the start
   '  o fun-01 Sec-04:
   '    - Make the aligned sequence
   '  o fun-01 Sec-05:
   '    - Add soft masking to the end
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-01:
   ^  - Variable declerations 
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *baseStr = seqST->seqCStr;
   char *tmpRetStr = 0;
   char *alignedSeqStr = 0;/*Returned aligned sequence*/
   char *seqAlnStr = 0;  /*Alignment array for sequence*/
   char *otherAlnStr = 0;/*Other sequence in alignment*/

   unsigned long seqStartUL = 0;
   unsigned long otherStartUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory & identify if ref or query sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   alignedSeqStr =
     calloc(
       alnST->qryLenUL + alnST->refLenUL + 2,
       sizeof(char)
   ); /*Allocate and initialize the alignment array*/

   if(seqAlnStr == 0) return 0;  // memory error

   tmpRetStr = alignedSeqStr;

   if(qryBl)
   { /*If I am building a query sequence*/
     seqAlnStr = alnST->qryAlnStr;
     otherAlnStr = alnST->refAlnStr;

     seqStartUL = alnST->qryStartAlnUL;
     otherStartUL = alnST->refStartAlnUL;
   } /*If I am building a query sequence*/

   else seqAlnStr = alnST->qryAlnStr;
   { /*If I am building an reference sequence*/
     seqAlnStr = alnST->refAlnStr;
     otherAlnStr = alnST->qryAlnStr;

     seqStartUL = alnST->refStartAlnUL;
     otherStartUL = alnST->qryStartAlnUL;
   } /*If I am building an reference sequence*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Add softmasking to the start
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(otherStartUL > seqStartUL)
   { /*While I have softmasked bases only for reference*/
     if(!extAlnRegionBl)
     { /*If I am returing the complete sequence*/
       *tmpRetStr = '-';
       ++tmpRetStr;
     } /*If I am returing the complete sequence*/

     ++otherAlnStr;
     --otherStartUL;
   } /*While I have softmasked bases only for reference*/

   while(*seqAlnStr == defSoftMaskFlag)
   { /*While I have softmasking at the start*/
     if(!extAlnRegionBl)
     { /*If I am returing the complete sequence*/
       *tmpRetStr = *baseStr;
       ++tmpRetStr;
     } /*If I am returing the complete sequence*/

     ++baseStr;
     ++seqAlnStr;

     if(*otherAlnStr == defSoftMaskFlag) ++otherAlnStr;
   } /*While I have softmasking at the start*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Make the aligned sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(*seqAlnStr != defEndAlnFlag)
   { /*While I have query bases to add*/
     if(*otherAlnStr == defGapFlag)
     { /*If I have an gap in the other sequence*/
       ++otherAlnStr;
       continue;
     } /*If I have an gap in the other sequence*/

     if(*otherAlnStr == defSoftMaskFlag) break;
       
     switch(*seqAlnStr)
     { /*Switch; check the error type*/
       case defEndAlnFlag:   goto makeAlnSeqCleanUp;
       case defSoftMaskFlag: goto makeAlnSeqCleanUp;

       case defSnpFlag:
       case defMatchFlag:
       /*Cases (part): Match and SNP*/
         ++otherAlnStr;
         /*Case defGapFlag: has the rest of the case*/
       /*Cases (part): Match and SNP*/

       case defGapFlag:
       /*Cases: Match, SNP, Gap*/
         ++seqAlnStr;
         *tmpRetStr = *baseStr;
         ++tmpRetStr;
         ++baseStr;
         break;
       /*Cases: Match, SNP, Gap*/
     } /*Switch; check the error type*/
   } /*While I have bases or an alignment to copy over*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Add soft masking to the end
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  makeAlnSeqCleanUp:

  /*Check if I am only extracting the aligned region*/
  if(extAlnRegionBl) goto makeAlnSeqEnd;

  while(*seqAlnStr == defEndAlnFlag)
  { /*While I have query softmaksing to add to the end*/
    if(*otherAlnStr != defEndAlnFlag) ++otherAlnStr;
    *tmpRetStr = *baseStr;
    ++tmpRetStr;
    ++baseStr;
    ++seqAlnStr;
  } /*While I have query softmaksing to add to the end*/

  while(*otherAlnStr != defEndAlnFlag)
  { /*While I have reference softmaksing to add to end*/
    *tmpRetStr = *baseStr;
    ++tmpRetStr;
    ++baseStr;
    ++otherAlnStr;
  } /*While I have softmaksing to add to end*/

  makeAlnSeqEnd:

  *tmpRetStr = '\0'; // Make into a c-string

  return alignedSeqStr;
} /*alnSTToSeq*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnStruct with the alingment array
|    o 0 if had memory allocation error
\--------------------------------------------------------*/
struct alnStruct * dirMatrixToAlnST(
    struct seqStruct *refST,     /*Reference seq & length*/
    struct seqStruct *qryST,     /*Query seq & length*/
    struct scoresStruct *scoreST,
       /*Score to get alignment for*/
    /*Matrix to use in finding the alignment*/
    #if !defined BYTEMATRIX
       struct twoBitAry *dirMatrixST
    #else
       char *dirMatrixST
    #endif
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: getAlnAry
   '  - Builds an anlignment array for the input direction
   '    matrix
   '  o fun-02 sec-01:
   '    - VariAble declerations
   '  o fun-02 sec-02:
   '    - Assign memory and initalize variables
   '  o fun-02 sec-03:
   '    - Build the alignment array
   '  o fun-02 sec-04:
   '    - Find the starting and ending softmasked regions
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-02 Sec-01:
  ^  - Variable declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  long qryIndexL = 0;     /*Index of query base*/
  long refIndexL = 0;     /*Index of reference base*/
  unsigned long startIndexUL = 0;
  long lastRefMatchSnpL = 0;
  long lastQryMatchSnpL = 0;
  long lenRefL = 
      (long) refST->endAlnUL - (long) refST->offsetUL + 1;

  /* Find the starting query and reference base
  `  - (lenQryUL or lenRefL) + 1:
  `    o gives the length of each column in the matrix
  `  - bestScoreUI (/ or %) (lenRefL + 1)
  `    o gives the index of the base for the best score
  `  - -1:
  `    o converts the index 1 output to index 0
  */

  char *qrySeqStr =
      qryST->seqCStr + scoreST->qryEndUL;
    /*qryST->offsetUL
    + (indexUL / (lenRefL + 1))
    - 1;*/
  char *refSeqStr =
      refST->seqCStr + scoreST->refEndUL;
    /*refST->offsetUL
    + (indexUL % (lenRefL + 1))
    - 1;*/

  struct alnStruct *alnST = 0;

  /*Variables for searching through the direction matrix*/
  #if !defined BYTEMATRIX
     uint8_t bitElmUC = 0;
  #endif

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-02 Sec-02:
  ^  - Assign memory and initalize variables
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  alnST = calloc(1, sizeof(struct alnStruct));
  if(alnST == 0) return 0;
    /*alnST is already initilized for non pointer values*/

  alnST->refAlnStr =
     calloc(refST->lenSeqUL + 1, sizeof(char));

  if(alnST->refAlnStr == 0)
  { /*If I had a memory error*/
    freeAlnST(alnST, 1);
    alnST = 0;
    return 0;
  } /*If I had a memory error*/

  alnST->qryAlnStr =
     calloc(qryST->lenSeqUL + 1, sizeof(char));

  if(alnST->refAlnStr == 0)
  { /*If I had a memory error*/
    freeAlnST(alnST, 1);
    alnST = 0;
    return 0;
  } /*If I had a memory error*/

  /*Get the length of the reference and query*/
  alnST->refLenUL = refST->lenSeqUL;
  alnST->qryLenUL = qryST->lenSeqUL;

  /*Get the end of the query and reference alignment*/
  alnST->refEndAlnUL = refSeqStr - refST->seqCStr;
  alnST->qryEndAlnUL = qrySeqStr - qryST->seqCStr;

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-02 Sec-03:
  ^  - Build the alignment array
  ^  o fun-02 sec-5 sub-1:
  ^    - Find the best score cell
  ^  o fun-02 sec-5 sub-2:
  ^    - Find the best path
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*******************************************************\
  * Fun-02 Sec-5 Sub-1:
  *  - Find the cell with the best score
  \*******************************************************/

  startIndexUL = (1 + scoreST->qryEndUL);
     /* 1 + to account for qryEndUL indel column*/
     /* Not accounting for qryEndUL being index 0, because
     `  bestScoreST.refEndUL as the reference position in
     `  the row*/
  startIndexUL *= (1 + lenRefL);
     /*1 + to account for deletion row*/
  startIndexUL += 1 + scoreST->refEndUL;
     /* 1 + to account for refEndUL being index 0*/
  

  /*get direction for the best element*/
  #if !defined BYTEMATRIX
    twoBitMvXElmFromStart(dirMatrixST, startIndexUL);
    bitElmUC = getTwoBitElm(dirMatrixST);
  #else
     dirMatrixST += startIndexUL;
  #endif
    /*I have not figured out how to convert end
    ` coordinates to an index. For some odd reason they
    ` are off when not at the end.
    */

  /*Grab the index of at the end of the alignment*/
  qryIndexL = qrySeqStr - qryST->seqCStr;
  refIndexL = refSeqStr - refST->seqCStr;

  /*******************************************************\
  * Fun-02 Sec-5 Sub-2:
  *  - Find the best path
  \*******************************************************/

  #if !defined BYTEMATRIX
  while(bitElmUC != defMvStop)
  #else
  while(*dirMatrixST != defMvStop)
  #endif
  { /*While I have more bases in the alignment*/
    #if !defined BYTEMATRIX
    switch(bitElmUC)
    #else
    switch(*dirMatrixST)
    #endif
    { /*Switch: check if bases is gap, match, or snp*/
      case defMvStop: goto finishAlignment;
      case defMvIns:
      /*Case: insertion (defMvIns)*/
        *(alnST->qryAlnStr + qryIndexL) = defGapFlag;

        --qrySeqStr; /*insertion only query has a base*/
        --qryIndexL;
        ++(alnST->numInssUL);
        #if !defined BYTEMATRIX
           twoBitMvBackXElm(dirMatrixST, lenRefL + 1);
        #else
           dirMatrixST -= (lenRefL + 1);
        #endif
          /* lenRefL (index 1) cells per row; need + 1 to
          `  get to the cell above
          */
        break;
      /*Case: insertion (defMvIns)*/

      case defMvSnp:
      /*Case: match/snp (defMvSnp)*/

        lastRefMatchSnpL = (unsigned long) refIndexL;
        lastQryMatchSnpL = (unsigned long) qryIndexL;

        if(checkIfBasesMatch(qrySeqStr, refSeqStr))
        { /*If the bases were a match*/
          *(alnST->qryAlnStr + qryIndexL) = defMatchFlag;
          *(alnST->refAlnStr + refIndexL) = defMatchFlag;
          ++(alnST->numMatchesUL);
        } /*If the bases were a match*/

        else
        { /*Else was a SNP*/
          *(alnST->qryAlnStr + qryIndexL) = defSnpFlag;
          *(alnST->refAlnStr + refIndexL) = defSnpFlag;
          ++(alnST->numSnpsUL);
        } /*Else was a SNP*/

        #if !defined BYTEMATRIX
           twoBitMvBackXElm(dirMatrixST, lenRefL + 2);
        #else
           dirMatrixST -= (lenRefL + 2);
        #endif
          /*need + 2 to get to the next diagnol cell*/
        --qrySeqStr;
        --qryIndexL;

        --refSeqStr;
        --refIndexL;
        break;
      /*Case: match/snp (defMvSnp)*/

      case defMvDel:
      /*Case: deletion (defMvDel)*/
        *(alnST->refAlnStr + refIndexL) = defGapFlag;
        ++(alnST->numDelsUL);

        #if !defined BYTEMATRIX
           twoBitMvBackOneElm(dirMatrixST);
        #else
           --dirMatrixST;
        #endif
        --refSeqStr;
        --refIndexL;
        break;
      /*Case: deletion (defMvDel)*/
    } /*Switch: check if bases is gap, match, or snp*/

    /*Get the next direction to move*/
    ++(alnST->lenAlnUL);

    #if !defined BYTEMATRIX
       bitElmUC = getTwoBitElm(dirMatrixST);
    #endif
  } /*While I have more bases to add to the path*/

  finishAlignment:

  /*Get the ending position for the alignment
  ` I really just want to recored the last match/snp
  */
  alnST->refStartAlnUL = lastRefMatchSnpL;  
  alnST->qryStartAlnUL = lastQryMatchSnpL;  

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-02 Sec-04:
  ^  - Find the starting and ending softmasked regions
  ^  o fun-02 sec-04 sub-01:
  ^    - Add softmasking to the start of the alingment
  ^  o fun-02 sec-04 sub-02:
  ^    - Add softmasking to the end of the alingment
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*******************************************************\
  * Fun-02 Sec-04 Sub-01:
  *  - Add softmasking to the start of the alingment
  \*******************************************************/

   while(refIndexL >= 0)
   { /*While the reference has starting softmasked bases*/
     *(alnST->refAlnStr + refIndexL) = defSoftMaskFlag;
     --refIndexL;
   } /*While the reference has starting softmasked bases*/

   while(qryIndexL >= 0)
   { /*While the query has starting softmasked bases*/
     *(alnST->qryAlnStr + qryIndexL) = defSoftMaskFlag;
     --qryIndexL;
   } /*While the query has starting softmasked bases*/

  /*******************************************************\
  * Fun-02 Sec-04 Sub-02:
  *  - Add softmasking to the end of the alingment
  \*******************************************************/

   /*Calloc set all variables to 0 at the ends*/
   refIndexL = alnST->refLenUL - 1;
   while(*(alnST->refAlnStr + refIndexL) == 0)
   { /*While I have softmaskig at the references end*/
     *(alnST->refAlnStr + refIndexL) = defSoftMaskFlag;
     --refIndexL;
   } /*While I have softmaskig at the end*/

   qryIndexL = alnST->qryLenUL - 1;
   while(*(alnST->qryAlnStr + qryIndexL) == 0)
   { /*While I have softmaskig at the querys end*/
     *(alnST->qryAlnStr + qryIndexL) = defSoftMaskFlag;
     --qryIndexL;
   } /*While I have softmaskig at the querys end*/

   /*Add the flags marking the end of the alignment*/
   *(alnST->refAlnStr + alnST->refLenUL) = defEndAlnFlag;
   *(alnST->qryAlnStr + alnST->qryLenUL) = defEndAlnFlag;

   return alnST;
} /*dirMatrixToAlnST*/

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints out the input alingment
|  - Returns:
|    o 0 for success
|    o 1 If failed to get time for emboss format
|    o 64 for memory errors.
\--------------------------------------------------------*/
char printAln(
  FILE *outFILE,            /*File to print alingment to*/
  char *outStr,             /*print file name to header*/
  struct seqStruct *refST,  /*reference Id & sequence*/
  struct seqStruct *qryST,  /*query Id & sequence*/
  struct alnStruct *alnST,  /*Has alignment to print out*/
  long scoreL,              /*Score of the alignment*/
  struct alnSet *settings,  /*Settings used for alingment*/
  char *scoreMtxFileStr     /*printing out scoring matrix*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printAln
   '  - Prints out the alignment to an alignment file
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Print out the header for the alignment
   '  o fun-03 sec-03:
   '    - Allocate memory and copy read ids
   '  o fun-03 sec-04:
   '    - Final prep before printing
   '  o fun-03 sec-05:
   '    - Print out the alignment
   '  o fun-03 sec-06:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Variable declerations
   ^  o fun-03 sec-01 sub-01:
   ^    - Variables dealing with buffers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-01 Sub-01:
   *  - Variables dealing with output buffers
   \******************************************************/

   /*These buffers hold one line of the alignment*/
   char *refBuffStr = 0;
   char *qryBuffStr = 0;
   char *eqxBuffStr = 0;

   /*These are for adding new bases to the buffers*/
   char *tmpRefStr = 0;
   char *tmpQryStr = 0;
   char *tmpEqxStr = 0;

   /*These point to point in buffer to add new bases/
   ` numbers to
   */
   char *startRefStr = 0;
   char *startQryStr = 0;
   char *startEqxStr = 0;

   /*I am limitnig reference id and query id to 10 char*/
   char *endRefIdStr = 0;
   char *endQryIdStr = 0;
   char oldRefCapCharC = 0;
   char oldQryCapCharC = 0;

   /******************************************************\
   * Fun-03 Sec-01 Sub-02:
   *  - Other variables
   \******************************************************/

   char *refAlnStr = 0;
   char *qryAlnStr = 0;

   char *refSeqStr = 0;
   char *qrySeqStr = 0;

   unsigned short lineWrapUS = settings->lineWrapUS;
      /*This is here because I have a minimum line wrap*/
   unsigned short wrapUS = 0;
      /*This will be adjusted to number of bases per line*/

   unsigned long basesInBuffUL = 0;  /*Base priting out*/

   unsigned long refBaseUL = 0; /*Reference base on*/
   unsigned long qryBaseUL = 0; /*Query base on*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Print out the header for the alignment
   ^  o fun-03 sec-02 sub-01:
   ^    - Check line wrap to see if is 0 (no wrap)
   ^  o fun-03 sec-02 sub-02:
   ^    - Default format (expand cigar) header/prep
   ^  o fun-03 sec-02 sub-03:
   ^    - EMOBSS format header/prep
   ^  o fun-03 sec-02 sub-04:
   ^    - Clustal format prep
   ^  o fun-03 sec-02 sub-05:
   ^    - Fasta format prep
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-01:
   *  - Check to see if line wrap is 0 (no line wrap)
   \******************************************************/

   /*If the user wanted no line wrapping*/
   if(lineWrapUS == 0)
      lineWrapUS = refST->lenSeqUL + qryST->lenSeqUL + 128;
      /*128 to provide room for ids and numbers*/

   wrapUS = lineWrapUS;

   /******************************************************\
   * Fun-03 Sec-02 Sub-02:
   *  - Default format (expand cigar) header/prep
   \******************************************************/

   switch(settings->formatFlag)
   { /*Switch: Check with format type using*/
     case defExpandCig:
     /*Case: Printing out the expanded cigar header*/
        pExpandCigHead(
          outFILE,
          refST,
          qryST,
          scoreL,
          scoreMtxFileStr,
          settings,
          alnST
        ); /*Print out the epxnaded cigar header*/

        if(wrapUS < 42)
        { /*If I need to change user input to minimum*/
           wrapUS = 42;
           lineWrapUS = 42;
           /*My min line wrap (longest line in header is
           ` 42 characters)
           */
        } /*If I need to change user input to minimum*/

        wrapUS -= 5;
          /*Tags for each line: 3 chars + ':' + ' '*/
        if(settings->pBasePosBl) wrapUS -= 20;
           /*Print starting and ending base positions
           ` I am printing digits at both ends
           */
        break;
     /*Case: Printing out the expanded cigar header*/

     /****************************************************\
     * Fun-03 Sec-02 Sub-03:
     *  - EMOBSS format header/prep
     \****************************************************/

     case defEMBOSS:
     /*Case: Printing out the EMBOSS header*/
        endRefIdStr = capIdLen(refST, 10, &oldRefCapCharC);
        endQryIdStr = capIdLen(qryST, 10, &oldQryCapCharC);

        /*EMBOSS uses numbers for lines*/
        if(wrapUS < 42)
        { /*If I need to change user input to minimum*/
           wrapUS = 42;
           lineWrapUS = 42;
           /*My min line wrap for EMBOSS (10 bases/line)*/
        } /*If I need to change user input to minimum*/

        wrapUS -= 22;
            /*10 = 10 char + ':' + ' ' + 9 digits + ' '*/
        wrapUS -= 10; /*for ending ' ' + 9 digits*/

        if(
           pEMBOSSHead(
              outFILE,
              refST,
              qryST,
              scoreL,
              outStr,
              scoreMtxFileStr,
              settings,/*Settings for the alignment*/
              alnST,
              1       /*I want to include the file header*/
           )
        ) return 1; /*Could not get the date/time*/

        break;
     /*Case: Printing out the EMBOSS header*/

     /****************************************************\
     * Fun-03 Sec-02 Sub-04:
     *  - Clustal format prep
     \****************************************************/

     case defClustal:
     /*Case: Printing out in clustal format*/
        endRefIdStr = capIdLen(refST, 10, &oldRefCapCharC);
        endQryIdStr = capIdLen(qryST, 10, &oldQryCapCharC);

        if(wrapUS < 32)
        { /*If I need to increase the line wrap to min*/
           wrapUS = 32;
           lineWrapUS = 32; /*10 bases per line*/
        } /*If I need to increase the line wrap to min*/

        wrapUS -= 12; /*12 = 10 char + ':' + ' '*/
        if(settings->pBasePosBl) wrapUS -= 10;
          /*10 = ' ' + 9 digits*/
          /*Clustal format has digits at end*/

        /*Clustal format prints out 60 bases at most*/
        if(wrapUS > 60)
        { /*If the line wrap was to large*/ 
           lineWrapUS -= (wrapUS - 60);
           wrapUS = 60;
        } /*If the line wrap was to large*/ 

        break;
     /*Case: Printing out in clustal format*/

     /****************************************************\
     * Fun-03 Sec-02 Sub-05:
     *  - Fasta format prep
     \****************************************************/

     case defFasta:
     /*Case: Printing out a fasta file*/
        lineWrapUS = refST->lenSeqUL + qryST->lenSeqUL;

        if(wrapUS)
        { /*If the user wanted line wrapping*/
           if(wrapUS < 10) wrapUS = 10;

           lineWrapUS += 1 + (lineWrapUS / wrapUS);
             /*Account for new lines in line*/
        } /*If the user wanted line wrapping*/

        else wrapUS = lineWrapUS;
        break;
     /*Case: Printing out a fasta file*/
   } /*Switch: Check with format type using*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - Allocate memory and copy read ids
   ^  o fun-03 sec-03 sub-01:
   ^    - Allocate memory
   ^  o fun-03 sec-03 sub-02:
   ^    - Set up buffers (add "\n\0" to end + start)
   ^  o fun-03 sec-03 sub-03:
   ^    - Copy tags for expand cigar format
   ^  o fun-03 sec-03 sub-04:
   ^    - Copy read ids for EMBOSS format
   ^  o fun-03 sec-03 sub-05:
   ^    - Copy read ids for clustal format
   ^  o fun-03 sec-03 sub-06:
   ^    - When printing numbers; adding padding to end of 
   ^      the eqx line
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-03 Sub-01:
   *  - Allocate memory
   \******************************************************/

   refBuffStr = calloc(lineWrapUS + 2, sizeof(char));
   if(refBuffStr == 0) return 64; /*Memory error*/

   qryBuffStr = calloc(lineWrapUS + 2, sizeof(char));
   if(qryBuffStr == 0)
   { /*If had a memory error*/
      free(refBuffStr);
      refBuffStr = 0;
      return 64; /*Memory error*/
   } /*If had a memory error*/

   if(settings->formatFlag != defFasta)
   { /*If not printing a fasta file*/
      eqxBuffStr = calloc(lineWrapUS + 2, sizeof(char));
      if(eqxBuffStr == 0)
      { /*If had a memory error*/
         free(refBuffStr);
         free(qryBuffStr);

         refBuffStr = 0;
         qryBuffStr = 0;

         return 64; /*Memory error*/
      } /*If had a memory error*/
   } /*If not printing a fasta file*/

   /******************************************************\
   * Fun-03 Sec-03 Sub-02:
   *  - Set up buffers (add "\n\0" to end + start)
   \******************************************************/

   *(refBuffStr + lineWrapUS) = '\n';
   *(refBuffStr + lineWrapUS + 1) = '\0';

   *(qryBuffStr + lineWrapUS) = '\n';
   *(qryBuffStr + lineWrapUS + 1) = '\0';

   if(settings->formatFlag != defFasta)
   { /*If not printing a fasta file*/
      *(eqxBuffStr + lineWrapUS) = '\n';
      *(eqxBuffStr + lineWrapUS + 1) = '\0';
   } /*If not printing a fasta file*/

   startRefStr = refBuffStr;
   startQryStr = qryBuffStr;
   startEqxStr = eqxBuffStr;

   /******************************************************\
   * Fun-03 Sec-03 Sub-03:
   *  - Copy tags for expand cigar format
   \******************************************************/

   switch(settings->formatFlag)
   { /*Switch: Check if adding read ids to buffer*/
     case defExpandCig:
     /*Case: Adding Ref, Qry, and Eqx tags to buffers*/
        strcpy(startRefStr, "Ref: ");
        startRefStr += 5;

        strcpy(startQryStr, "Qry: ");
        startQryStr += 5;

        strcpy(startEqxStr, "Eqx: ");
        startEqxStr += 5;

        /*Check if printing out numbers*/
        if(settings->pBasePosBl)
        { /*If I need to add padding to the eqx line*/
           for(int eqxCntI = 0; eqxCntI < 10; ++eqxCntI)
           { /*Loop: add padding to eqx line*/
              *startEqxStr = ' ';
              ++startEqxStr;
           } /*Loop: add padding to eqx line*/
        } /*If I need to add padding to the eqx line*/

        break;
     /*Case: Adding Ref, Qry, and Eqx tags to buffers*/

     /****************************************************\
     * Fun-03 Sec-03 Sub-04:
     *  - Copy read ids for EMBOSS format
     \****************************************************/

     case defEMBOSS:
     /*Case: If printing alignment in EMBOSS format*/
        startRefStr =
           cpReadIdRPad(refST, startRefStr, ':', 12);
        startQryStr =
           cpReadIdRPad(qryST, startQryStr, ':', 12);

         for(int eqxCntI = 0; eqxCntI < 22; ++eqxCntI)
         { /*Loop: add padding to eqx line*/
            *startEqxStr = ' ';
            ++startEqxStr;
         } /*Loop: add padding to eqx line*/

        break;
     /*Case: If printing alignment in EMBOSS format*/

     /****************************************************\
     * Fun-03 Sec-03 Sub-05:
     *  - Copy read ids for clustal format
     \****************************************************/

     case defClustal:
     /*Case: Printing format in clustal format*/
         startRefStr =
            cpReadIdRPad(refST, startRefStr, ':', 12);
         startQryStr =
            cpReadIdRPad(qryST, startQryStr, ':', 12);

         for(int eqxCntI = 0; eqxCntI < 12; ++eqxCntI)
         { /*Loop: add padding to eqx line*/
            *startEqxStr = ' ';
            ++startEqxStr;
         } /*Loop: add padding to eqx line*/
     /*Case: Printing format in clustal format*/

     case defFasta: break; /*Handled at end*/
   } /*Switch: Check if adding read ids to buffer*/

   /******************************************************\
   * Fun-03 Sec-03 Sub-06:
   *  - When printing numbers; adding padding to end of the
   *    eqx line
   \******************************************************/

   /*Check if printing out numbers*/
   if(
      ( settings->pBasePosBl ||
        settings->formatFlag == defEMBOSS
      ) &&
      settings->formatFlag != defFasta /*No eqx line*/
   ){/*If I need to add padding to the end of the eqxLine*/
      tmpEqxStr = startEqxStr + wrapUS;

      for(int eqxCntI = 0; eqxCntI < 10; ++eqxCntI)
      { /*Loop: add padding to end of the eqx line*/
         *tmpEqxStr = ' ';
         ++tmpEqxStr;
      } /*Loop: add padding to end of the eqx line*/
   }/*If I need to add padding to the end of the eqx line*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-04:
   ^  - Final prep before printing
   ^  o fun-03 sec-04 sub-01:
   ^    - Find the position I am starting to print at
   ^  o fun-03 sec-04 sub-02:
   ^    - Add start base position (if asked for) to buffers
   ^  o fun-03 sec-04 sub-03:
   ^    - Print out the alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   ^  o fun-03 sec-04 sub-01:
   ^    - Find the position I am starting to print at
   \******************************************************/

   tmpRefStr = startRefStr;
   tmpQryStr = startQryStr;
   tmpEqxStr = startEqxStr;

   if(!settings->pFullAlnBl)
   { /*If only printing the aligned region*/
      refAlnStr = alnST->refAlnStr + alnST->refStartAlnUL;     
      qryAlnStr = alnST->qryAlnStr + alnST->qryStartAlnUL;     

      refSeqStr = refST->seqCStr + alnST->refStartAlnUL;
      qrySeqStr = qryST->seqCStr + alnST->qryStartAlnUL;

      refBaseUL = alnST->refStartAlnUL + 1;
      qryBaseUL = alnST->qryStartAlnUL + 1;
      /*Need a + 1 to convert index 0 to index 1*/
   } /*If only printing the aligned region*/

   else
   { /*Else I am printing the full alignment*/
      refAlnStr = alnST->refAlnStr;     
      qryAlnStr = alnST->qryAlnStr;     

      refSeqStr = refST->seqCStr;
      qrySeqStr = qryST->seqCStr;

      /*Not the best start, but it will do*/
      refBaseUL = 1;
      qryBaseUL = 1;
   } /*Else I am printing the full alignment*/

   /******************************************************\
   * Fun-03 Sec-04 Sub-02:
   *  - Add start base position (if asked for) to buffers
   \******************************************************/

   tmpRefStr = 
      addPosToBuff(tmpRefStr,refBaseUL,settings,1);
      /*1 Is for the first bases in buffer*/

   tmpQryStr = 
      addPosToBuff(tmpQryStr,qryBaseUL,settings,1);
      /*1 Is for the first bases in buffer*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-05:
   ^  - Print out the alignment
   ^  o fun-03 sec-05 sub-01:
   ^    - Check for insertions before other entries
   ^  o fun-03 sec-05 sub-02:
   ^    - Check query soft mask before reference checks
   ^  o fun-03 sec-05 sub-03:
   ^    - Check for reference softmasks entries
   ^  o fun-03 sec-05 sub-04:
   ^    - Check for deletions
   ^  o fun-03 sec-05 sub-05:
   ^    - Check for SNPs
   ^  o fun-03 sec-05 sub-06:
   ^    - Check for matches
   ^  o fun-03 sec-05 sub-07:
   ^    - Move to next base
   ^  o fun-03 sec-05 sub-08:
   ^    - Print out the buffer (when full)
   ^  o fun-03 sec-05 sub-09:
   ^    - Print out the fasta entry
   ^  o fun-03 sec-05 sub-10:
   ^    - Do the final non-fasta print
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-03 Sec-05 Sub-01:
    *  - Check for insertions before other entries
    \*****************************************************/

    while(
       *refAlnStr != defEndAlnFlag ||
       *qryAlnStr != defEndAlnFlag
    ){ /*Loop: Print out all bases*/
       /*Check if only printing out aligned portion*/
       if(!settings->pFullAlnBl)
         if(
           qryBaseUL > alnST->qryEndAlnUL + 1 &&
           refBaseUL > alnST->refEndAlnUL + 1
         ) goto printAlnFinalPrint;
         /*+1 to account for index 1*/

       if(*qryAlnStr == defGapFlag)
       { /*If the reference has a gap*/
          *tmpRefStr = '-';
          *tmpQryStr = *qrySeqStr;
          eqxAddGap(1, tmpEqxStr, settings->formatFlag);
            /*1 = insertion*/

          ++qryAlnStr;
          ++qrySeqStr;
          ++qryBaseUL;

          goto pAlnIncBuff;
       } /*If the reference has a gap*/

       /**************************************************\
       * Fun-03 Sec-05 Sub-02:
       *  - Check query soft mask before reference checks
       \**************************************************/

       if(*qryAlnStr == defSoftMaskFlag)
       { /*If have a soft masked region*/
          if(*refAlnStr == defSoftMaskFlag)
          { /*If there is also a soft masked ref base*/
             *tmpRefStr = *refSeqStr;
             ++refAlnStr;
             ++refSeqStr;
             ++refBaseUL;
          } /*If there is also a soft masked ref base*/

          else *tmpRefStr = '-';

          *tmpQryStr = *qrySeqStr;
          eqxAddSMask(tmpEqxStr, settings->formatFlag);

          ++qryAlnStr;
          ++qrySeqStr;
          ++qryBaseUL;

          goto pAlnIncBuff;
       } /*If have a soft masked region*/

       /**************************************************\
       * Fun-03 Sec-05 Sub-03:
       *  - Check for reference softmasks entries
       \**************************************************/

       switch(*refAlnStr)
       { /*Switch: Check if ref is gap,softmask,snp,match*/
          case defSoftMaskFlag:
          /*Case: Reference has a softmasked base*/
             /*Already check if query was sofmasked*/
             *tmpQryStr = '-';
             *tmpRefStr = *refSeqStr;

             eqxAddSMask(tmpEqxStr, settings->formatFlag);

             ++refAlnStr;
             ++refSeqStr;
             ++refBaseUL;

             break;
          /*Case: Reference has a softmasked base*/

          /***********************************************\
          * Fun-03 Sec-05 Sub-04:
          *  - Check for deletions
          \***********************************************/
         
          case defGapFlag:
          /*Case: deletion (reference maps to gap)*/
             *tmpQryStr = '-';
             *tmpRefStr = *refSeqStr;
             eqxAddGap(0, tmpEqxStr, settings->formatFlag);
                /*0 = del*/

             ++refAlnStr;
             ++refSeqStr;
             ++refBaseUL;

             break;
          /*Case: deletion (reference maps to gap)*/

          /***********************************************\
          * Fun-03 Sec-05 Sub-05:
          *  - Check for SNPs
          \***********************************************/

          case defSnpFlag:
          /*Case: Reference and query have snps*/
             *tmpQryStr = *qrySeqStr;
             *tmpRefStr = *refSeqStr;

             eqxAddSnp(tmpEqxStr, settings->formatFlag);

             ++refAlnStr;
             ++refSeqStr;
             ++refBaseUL;

             ++qryAlnStr;
             ++qrySeqStr;
             ++qryBaseUL;

             break;
          /*Case: Reference and query have snps*/

          /***********************************************\
          * Fun-03 Sec-05 Sub-06:
          *  - Check for matches
          \***********************************************/

          case defMatchFlag:
          /*Case: Reference and query have a match*/
             *tmpQryStr = *qrySeqStr;
             *tmpRefStr = *refSeqStr;

             eqxAddMatch(tmpEqxStr, settings->formatFlag);

             ++refAlnStr;
             ++refSeqStr;
             ++refBaseUL;

             ++qryAlnStr;
             ++qrySeqStr;
             ++qryBaseUL;

             break;
          /*Case: Reference and query have a match*/
       } /*Switch: Check if ref is gap,softmask,snp,match*/

       /**************************************************\
       * Fun-03 Sec-05 Sub-07:
       *  - Move to next base
       \**************************************************/

       pAlnIncBuff:

       /*Incurment buffers*/
       ++tmpQryStr;
       ++tmpRefStr;
       if(settings->formatFlag != defFasta) ++tmpEqxStr;
       ++basesInBuffUL;

       /**************************************************\
       * Fun-03 Sec-05 Sub-08:
       *  - Print out the buffer (when full)
       *  o fun-03 sec-05 sub-08 cat-01:
       *    - Add the ending base positions for printing
       *  o fun-03 sec-05 sub-08 cat-02:
       *    - Print out buffers and check if finished
       *  o fun-03 sec-05 sub-08 cat-03:
       *    - Reset buffers for the next line
       \**************************************************/

       /*+++++++++++++++++++++++++++++++++++++++++++++++++\
       + Fun-03 Sec-05 Sub-08 Cat-01:
       +  - Add the ending base positions for printing
       \+++++++++++++++++++++++++++++++++++++++++++++++++*/

       if(basesInBuffUL >= wrapUS)
       { /*If I need to print out the buffer*/
          basesInBuffUL = 0;

          if(settings->formatFlag == defFasta)
          { /*If working on a fasta file*/
             *tmpQryStr = '\n';
             *tmpRefStr = '\n';
             ++tmpQryStr;
             ++tmpRefStr;
             continue;
          } /*If working on a fasta file*/
  
          /*Add the ending base positions if needed*/
          tmpRefStr = 
            addPosToBuff(tmpRefStr,refBaseUL-1,settings,0);
             /*0 Is for the last base in buffer*/
             /*-1 to account for being on the next base*/

          tmpQryStr = 
            addPosToBuff(tmpQryStr,qryBaseUL-1,settings,0);
             /*0 Is for the last bases in buffer*/
             /*-1 to account for being on the next base*/

          /*++++++++++++++++++++++++++++++++++++++++++++++\
          + Fun-03 Sec-05 Sub-08 Cat-02:
          +  - Print out buffers and check if finished
          \++++++++++++++++++++++++++++++++++++++++++++++*/

          /*Print out my buffers*/
          fprintf(outFILE, "%s", refBuffStr);

          if(settings->formatFlag == defEMBOSS)
          { /*If printing in emboss format*/
             fprintf(outFILE, "%s", eqxBuffStr);
             fprintf(outFILE, "%s\n", qryBuffStr);
             /*\n to add a blank line between entries*/
          } /*If printing in emboss format*/
 
          else
          { /*Else (eqx line is last*/
             fprintf(outFILE, "%s", qryBuffStr);
             fprintf(outFILE, "%s\n", eqxBuffStr);
             /*\n to add a blank line between entries*/
          } /*Else (eqx line is last*/

          if(
             *refAlnStr == defEndAlnFlag &&
             *qryAlnStr == defEndAlnFlag
          ) goto pAlnReturn; /*finished*/
            /*Currently on next ref/query base*/

          /*Check if only printing out aligned portion*/
          if(!settings->pFullAlnBl)
            if(
              qryBaseUL > alnST->qryEndAlnUL + 1 &&
              refBaseUL > alnST->refEndAlnUL + 1
            ) goto pAlnReturn;
            /*+1 to account for being on the next base
            ` & comparing index 1 to index 0
            */

          /*++++++++++++++++++++++++++++++++++++++++++++++\
          + Fun-03 Sec-05 Sub-08 Cat-03:
          +  - Reset buffers for the next line
          \++++++++++++++++++++++++++++++++++++++++++++++*/

          /*Reset my buffers*/
          tmpRefStr = startRefStr;
          tmpQryStr = startQryStr;
          tmpEqxStr = startEqxStr;

          /*Add the starting base positions if needed*/
          /*This is already on the next base*/
          tmpRefStr = 
            addPosToBuff(tmpRefStr,refBaseUL,settings,1);

          tmpQryStr = 
            addPosToBuff(tmpQryStr,qryBaseUL,settings,1);
       } /*If I need to print out the buffer*/
    }  /*Loop: Print out all bases*/

   /******************************************************\
   * Fun-03 Sec-05 Sub-09:
   *  - Print out the fasta entry
   \******************************************************/

    printAlnFinalPrint:

    if(settings->formatFlag == defFasta)
    { /*If printing out fasta output*/
       /*Print out the reference header*/
       if(*(refST->idCStr + refST->lenIdUL - 1) == '\n')
          *(refST->idCStr + refST->lenIdUL - 1) = '\0';

       if(*refST->idCStr == '>')
          fprintf(outFILE,"%s %li",refST->idCStr,scoreL);
       else
          fprintf(outFILE,">%s %li",refST->idCStr,scoreL);

       if(*(refST->idCStr + refST->lenIdUL - 1) == '\0')
          *(refST->idCStr + refST->lenIdUL - 1) = '\n';

       if(!settings->pFullAlnBl)
          fprintf(
             outFILE,
             " %li %li\n",
             alnST->refStartAlnUL + 1,
             alnST->refEndAlnUL + 1
          ); /*If I only printed the ailgned region*/

       else fprintf(outFILE, " 1 %li\n", refST->lenSeqUL);

       /*Print out the reference sequence*/
       *tmpRefStr = '\0';
       fprintf(outFILE, "%s\n\n", refBuffStr);

       /*Print out the query header*/

       if(*(qryST->idCStr + qryST->lenIdUL - 1) == '\n')
          *(qryST->idCStr + qryST->lenIdUL - 1) = '\0';

       if(*qryST->idCStr == '>')
          fprintf(outFILE,"%s %li",qryST->idCStr,scoreL);
       else
          fprintf(outFILE,">%s %li",qryST->idCStr,scoreL);

       if(*(qryST->idCStr + qryST->lenIdUL - 1) == '\0')
          *(qryST->idCStr + qryST->lenIdUL - 1) = '\n';

       if(!settings->pFullAlnBl)
          fprintf(
             outFILE,
             " %li %li\n",
             alnST->qryStartAlnUL + 1,
             alnST->qryEndAlnUL + 1
          ); /*If I only printed the ailgned region*/

       else fprintf(outFILE, " 1 %li\n", qryST->lenSeqUL);

       /*Print out the query sequence*/
       *tmpQryStr = '\0';
       fprintf(outFILE, "%s\n", qryBuffStr);

       goto pAlnReturn;
    } /*If printing out fasta output*/

   /******************************************************\
   * Fun-03 Sec-05 Sub-10:
   *  - Do the final non-fasta print
   \******************************************************/

    /*Add the ending base positions if needed*/
    tmpRefStr = 
       addPosToBuff(tmpRefStr, refBaseUL - 1, settings, 0);
       /*0 Is for the last base in buffer*/
       /*-1 for being one base ahead*/

    tmpQryStr = 
       addPosToBuff(tmpQryStr, qryBaseUL - 1, settings, 0);
       /*0 Is for the last base in buffer*/
       /*-1 for being one base ahead*/

   *tmpRefStr = '\n';
   *tmpQryStr = '\n';
   *tmpEqxStr = '\n';

   *(tmpRefStr + 1) = '\0';
   *(tmpQryStr + 1) = '\0';
   *(tmpEqxStr + 1) = '\0';

   fprintf(outFILE, "%s", refBuffStr);

   if(settings->formatFlag == defEMBOSS)
      fprintf(outFILE,"%s%s",eqxBuffStr,qryBuffStr);
      /*\n to add a blank line between entries*/
 
   else fprintf(outFILE,"%s%s",qryBuffStr,eqxBuffStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-06:
   ^  - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   pAlnReturn:

   free(refBuffStr);
   free(qryBuffStr);
   free(eqxBuffStr);

   /*Revert ids back to full length*/
   if(endRefIdStr != 0) *endRefIdStr = oldRefCapCharC;
   if(endQryIdStr != 0) *endQryIdStr = oldQryCapCharC;

   return 0;
} /*printAln*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o All variables in alnST to be 0
| Note:
|  - This does not free alnAryUc, so only call this for
|    new alnST structures
\--------------------------------------------------------*/
void initAlnST(
  struct alnStruct *alnST // Strucutre to initialize
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-01: initAlnST
   '  - Initalize all values in alnST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   alnST->refAlnStr = 0;
   alnST->qryAlnStr = 0;
   alnST->refLenUL = 0;
   alnST->qryLenUL = 0;
   alnST->lenAlnUL = 0;

   alnST->refStartAlnUL = 0;
   alnST->refEndAlnUL = 0;
   alnST->qryStartAlnUL = 0;
   alnST->qryEndAlnUL = 0;

   alnST->numInssUL = 0;
   alnST->numDelsUL = 0;
   alnST->numSnpsUL = 0;
   alnST->numMatchesUL = 0;

   return;
} // initAlnST

/*--------------------------------------------------------\
| Output:
|  - Frees
|    o alnST, including all variables in in alnST
|    0 Or if heapBl is 0; frees all variables in alnST, but 
|      does not free alnST
\--------------------------------------------------------*/
void freeAlnST(
  struct alnStruct *alnST, // Strucutre to free
  char heapBl // 0: free variables in alnST, but keep alnST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-01: freeAlnST
   '  - Frees alnST and all variables in alnST
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(alnST->refAlnStr != 0) free(alnST->refAlnStr);
   if(alnST->qryAlnStr != 0) free(alnST->qryAlnStr);

   if(heapBl != 0) free(alnST);

   return;
} // initAlnST

/*--------------------------------------------------------\
| Output:
|  - Prints:
|    o The EMBOSS format headers to a file
|  - Returns:
|    o 1 if their was an error when getting the time
|    o 0 for no errors
\--------------------------------------------------------*/
char pEMBOSSHead(
  FILE *outFILE,          /*File to print the header to*/
  struct seqStruct *refST,/*Has reference sequence id*/
  struct seqStruct *qryST,/*Has query sequence id*/
  long scoreL,           /*Score of the alignment*/
  char *outStr,          /*Name of the output file*/
  char *scoreMtxFileStr, /*Name of file with score matrix*/ 
  struct alnSet *setting,/*Settings for the alignment*/
  struct alnStruct *alnST, /*Results from aligment*/
  char pFileHeaderBl
    /*1: Print the file header and entry header
    ` 0: Just print the entry header*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: pEMBOSSHead
   '  - Prints out the EMBOSS header to a file
   '  o fun-06 sec-01:
   '    - Variable declerations (& get time)
   '  o fun-06 sec-02:
   '    - Print out the file header
   '  o fun-06 sec-03:
   '    - Print out the entry header
   \*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-01:
   ^  - Variable declerations (& get time)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Get the date*/
   unsigned long lenAlnUL = 0;
   time_t timeST = time(NULL);
   struct tm *localTimeST = localtime(&timeST);
   char dateStr[64];
   size_t errT =
     strftime(
        dateStr,
        64,
        "%a %b %d %T %Y",
         localTimeST
     ); /*Format the date to my output format*/
     /*strftime formats the time
       %a is abbriviated day name
       %b/%h is abbriviated month name
       %d is the day
       %T = %H:%M:%T (Military time format)
       %Y is the full year
     */

   if((!errT) & pFileHeaderBl) return 1;
     /*If I could not get the time when printing the file
     ` header
     */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-02:
   ^  - Print out the file header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(pFileHeaderBl)
   { /*If printing out the file header*/
     fprintf(
       outFILE,
       "########################################\n"
     );

     if(setting->useNeedleBl == 1)
       fprintf(outFILE, "# Program: alnSeq-Needle\n");
     else if(setting->useWaterBl == 1)
       fprintf(outFILE, "# Program: alnSeq-Water\n");
     else if(setting->useHirschBl == 1)
       fprintf(outFILE, "# Program: alnSeq-Hirschberg\n");

     fprintf(outFILE, "# Rundate:  %s\n", dateStr);

     if(outStr == 0)
        fprintf(outFILE, "# Report_file: stdout\n");
     else fprintf(outFILE, "# Report_file: %s\n", outStr);

     fprintf(
       outFILE,
       "########################################\n"
     );
   } /*If printing out the file header*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-03:
   ^  - Print out the entry header
   ^  o fun-06 sec-03 sub-01:
   ^    - Print out settings
   ^  o fun-06 sec-03 sub-02:
   ^    - Print out the alignment stats
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Print out settings
   \******************************************************/

   fprintf(
      outFILE,
      "#=======================================\n"
   );

   fprintf(outFILE, "#\n# Aligned_sequences: 2\n");
   fprintf(outFILE, "# 1: %s\n", refST->idCStr + 1);
   fprintf(outFILE, "# 2: %s\n", qryST->idCStr + 1);
     /*+1 will get off the header marker (>)*/

   /*Check if using the default scoring matrix*/
   if(!scoreMtxFileStr)
      fprintf(outFILE, "# Matrix: %s\n", defMatrixNameStr);
   else fprintf(outFILE,"# Matrix: %s\n",scoreMtxFileStr);

   #if !defined NOGAPOPEN
      fprintf(
         outFILE,
         "# Gap_penalty: %i\n",
         setting->gapOpenI
      );
   #endif

   fprintf(
     outFILE,
     "# Extend_penalty: %i\n",
     setting->gapExtendI
   );

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Print out the alignment stats
   \******************************************************/

   lenAlnUL =
        alnST->numMatchesUL
      + alnST->numInssUL
      + alnST->numDelsUL;

   fprintf(
      outFILE,
      "#\n# Length: %lu\n",
      lenAlnUL
   );

   fprintf(
     outFILE,
     "# Identity:  %9li/%li (%.1f%%)\n",
     alnST->numMatchesUL,
     lenAlnUL,
     ((double) (100 * alnST->numMatchesUL) / lenAlnUL)
   );

   /*In may case Similarity (>51% of bases the same) is
   ` always equal to idenity.
   */
   fprintf(
     outFILE,
     "# Similarity:%9li/%li (%.1f%%)\n",
     alnST->numMatchesUL,
     lenAlnUL,
     ((double) (100 * alnST->numMatchesUL) / lenAlnUL)
   );

   fprintf(
     outFILE,
     "# Gaps:      %9li/%li (%.1f%%)\n",
     alnST->numInssUL + alnST->numDelsUL,
     lenAlnUL,
     ((double)
        (100*(alnST->numInssUL+alnST->numDelsUL))/lenAlnUL)
   );

   fprintf(outFILE, "# Score: %li\n", scoreL);

   fprintf(outFILE, "#\n#\n");
   fprintf(
     outFILE,
     "#=======================================\n\n"
   );

   return 0;
} /*pEMBOSSHead*/

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o The expanded cigar header to outFILE
\--------------------------------------------------------*/
void pExpandCigHead(
  FILE *outFILE,          /*File to print header to*/
  struct seqStruct *refST,/*Has reference sequence id*/
  struct seqStruct *qryST,/*Has query sequence id*/
  long scoreL,           /*Score of the alignment*/
  char *scoreMtxFileStr, /*Name of file with score matrix*/ 
  struct alnSet *setting,/*Settings for the alignment*/
  struct alnStruct *alnST/*Results from aligment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: pExpandCigHead
   '  - Prints out the expanded cigar header entry
   '  o fun-07 sec-01:
   '    - Print out the query and reference information
   '  o fun-07 sec-02:
   '    - Print out the settings
   '  o fun-07 sec-03:
   '    - Print out the alignment stats
   '  o fun-07 sec-04:
   '    - Print out the legend
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-01:
   ^  - Print out the query and reference information
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
     outFILE,
     "###########################################"
   );
   fprintf(outFILE, "\n");

   if(*(qryST->idCStr + qryST->lenIdUL - 1) == '\n')
      fprintf(outFILE, "# Qry = %s", qryST->idCStr + 1);
   else
      fprintf(outFILE, "# Qry = %s\n", qryST->idCStr + 1);

   fprintf(
     outFILE,
     "#    - Query aligned bases: %lu to %lu\n",
     alnST->qryStartAlnUL + 1, /*+1 to convert to index 1*/
     alnST->qryEndAlnUL + 1    /*+1 to convert to index 1*/
   );
     
   if(*(qryST->idCStr + qryST->lenIdUL - 1) == '\n')
      fprintf(outFILE, "# ref = %s", refST->idCStr + 1);
   else
      fprintf(outFILE, "# ref = %s\n", refST->idCStr + 1);

   fprintf(
     outFILE,
     "#    - Reference aligned bases: %lu to %lu\n",
     alnST->refStartAlnUL + 1, /*+1 to convert to index 1*/
     alnST->refEndAlnUL + 1    /*+1 to convert to index 1*/
   );

   fprintf(outFILE, "#\n"); /*Separator*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-02:
   ^  - Print out the settings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(setting->useNeedleBl == 1)
     fprintf(outFILE, "# Program: alnSeq-Needle\n");
   else if(setting->useWaterBl == 1)
     fprintf(outFILE, "# Program: alnSeq-Water\n");
   else if(setting->useHirschBl == 1)
     fprintf(outFILE, "# Program: alnSeq-Hirschberg\n");

   /*Check if using the default scoring matrix*/
   if(!scoreMtxFileStr)
     fprintf(outFILE, "# Matrix: %s\n", defMatrixNameStr);
   else fprintf(outFILE,"# Matrix: %s\n", scoreMtxFileStr);

   #if !defined NOGAPOPEN
      fprintf(
        outFILE,
        "# Gap open: %i\n",
         setting->gapOpenI
      );
   #endif

   fprintf(
      outFILE,
      "# Gap extend: %i\n",
      setting->gapExtendI
   );

   fprintf(outFILE, "#\n"); /*Separator*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-03:
   ^  - Print out the alignment stats
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      outFILE,
      "# Matches:    %lu\n",
      alnST->numMatchesUL
   );

   fprintf(
      outFILE,
      "# Mismatches: %lu\n",
      alnST->numSnpsUL
   );

   fprintf(
      outFILE,
      "# Insertions: %lu\n",
      alnST->numInssUL
   );

   fprintf(
      outFILE,
     "# Deletions:  %lu\n",
     alnST->numDelsUL
   );

   fprintf(
      outFILE,
      "# Total:      %lu\n",
      alnST->lenAlnUL
   );

   fprintf(outFILE,"# Alignment Score: %li\n", scoreL);

   fprintf(outFILE, "#\n"); /*Spacer*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-04:
   ^  - Print out the legend
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(outFILE, "# Eqx line legend\n");
   fprintf(outFILE, "#   - = is match\n");
   fprintf(outFILE, "#   - X is mismatch\n");
   fprintf(outFILE, "#   - I is insertion\n");
   fprintf(outFILE, "#   - D is deletion\n");
   fprintf(outFILE, "#   - S is soft mask\n");
   fprintf(
     outFILE,
     "###########################################\n\n"
   );

   return;
} /*pExpandCigHead*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o seqST->idCStr to end at the first white space or
|      at index maxIdLenI - 1 bases
|    o oldCharC to hold the character that was changed to
|      '\0'
|  - Returns:
|    o New end of seqST->idCStr (position changed to '\0')
\--------------------------------------------------------*/
char * capIdLen(
  struct seqStruct *seqST, /*Has id to limit length on*/
  int maxIdLenI,           /*Max id length (index 1)*/
  char *oldCharC    /*Will hold character changed to '\0'*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: capIdLen
   '  - Caps id length in seqST first white space or the
   '    character at maxIdLenI - 1.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   char *idIterStr = seqST->idCStr;

   /*Get off the header markers*/
   if(*idIterStr == '>' || *idIterStr == '@') ++idIterStr;

   /*Get off white space*/
   while(*idIterStr < 33) ++idIterStr;

   /*Find the first white space*/
   while(*idIterStr > 32) ++idIterStr;

   if(idIterStr - seqST->idCStr >= maxIdLenI)
      idIterStr = seqST->idCStr + maxIdLenI - 1;

   *oldCharC = *idIterStr; /*Save the old character*/
   *idIterStr = '\0';

   return idIterStr;
} /*capIdLen*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o buffStr to hold the new base position if the
|      settings support position printing.
|  - Returns:
|    o Pointer to chacter after added base position or 
|      input buffStr address if buffStr was not changed.
\--------------------------------------------------------*/
char * addPosToBuff(
   char *buffStr,          /*Buffer to add position to*/
   unsigned long basePosUL,/*Position to add to buffer*/
   struct alnSet *settings,
     /*Has boleans to tell if printing first base and
     ` format printing to
     */
   char startBl            /*1: is first base in buffer*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: addPosToBuff
   '  - Adds a base position to a buffer. This uses the
   '    printing format and position as a guide to
   '    determine if printing.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(settings->formatFlag)
   { /*Switch: Check if adding a base position*/
      case defExpandCig:
      /*Case: Using the expanded cigar format*/
         if(settings->pBasePosBl)
         { /*If adding the base position*/
            /*Check if printing at start or end*/
            if(startBl) sprintf(buffStr,"%9lu ",basePosUL);
            else sprintf(buffStr, " %9lu\n", basePosUL);
               /*sprintf adds '\0' to end*/

            return buffStr + 10;
         } /*If adding the base position*/

         return buffStr;
      /*Case: Using the expanded cigar format*/

      case defEMBOSS:
      /*Case: Using emboss format*/
         if(startBl) sprintf(buffStr,"%9lu ",basePosUL);
         else sprintf(buffStr, " %9lu\n", basePosUL);
            /*sprintf adds '\0' to end*/

         return buffStr + 10;
      /*Case: Using emboss format*/

      case defClustal:
         if(settings->pBasePosBl && !startBl)
         { /*If base position is at end (only clustal)*/
            sprintf(buffStr," %9lu\n", basePosUL);
               /*sprintf adds '\0' to end*/
            return buffStr + 10;
         } /*If base position is at end (only clustal)*/

         return buffStr;

      case defFasta: return buffStr; /*No position entry*/
   } /*Switch: Check if adding starting bases position*/

   return buffStr;
} /*addPosToBuff*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o eqxBuffStr to have the correct gap entry based on
|      formatFlag.
\--------------------------------------------------------*/
void eqxAddGap(
   char insBl,         /*1:insertion; 0: deletion)*/
   char *eqxBuffStr,   /*Buffer with eqx entry*/
   char formatFlag     /*Has format to add gap as*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: eqxAddGap
   '  - Adds a gap to an eqx buffer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(formatFlag)
   { /*Switch: Check if adding a base position*/
      case defExpandCig:
      /*Case: Using the expanded cigar format*/
         if(insBl) *eqxBuffStr = 'I';
         else *eqxBuffStr = 'D';
         return;
      /*Case: Using the expanded cigar format*/

      case defEMBOSS:
      case defClustal:
         *eqxBuffStr = ' ';
         return;

     case defFasta: return; /*No eqx entry*/
   } /*Switch: Check if adding starting bases position*/

   return;
} /*eqxAddGap*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o eqxBuffStr to have the correct soft mask entry based
|      on formatFlag.
\--------------------------------------------------------*/
void eqxAddSMask(
   char *eqxBuffStr,   /*Buffer with eqx entry*/
   char formatFlag     /*Has format to add soft mask as*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: eqxAddSMask
   '  - Adds an soft mask to an eqx buffer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(formatFlag)
   { /*Switch: Check if adding a base position*/
      case defExpandCig:
      /*Case: Using the expanded cigar format*/
         *eqxBuffStr = 'S';
         return;
      /*Case: Using the expanded cigar format*/

      case defEMBOSS:
      case defClustal:
         *eqxBuffStr = ' ';
         return;

     case defFasta: return; /*No eqx entry*/
   } /*Switch: Check if adding starting bases position*/

   return;
} /*eqxAddSMask*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o eqxBuffStr to have the correct SNP entry based on
|      formatFlag.
\--------------------------------------------------------*/
void eqxAddSnp(
   char *eqxBuffStr,   /*Buffer with eqx entry*/
   char formatFlag     /*Has format to add soft mask as*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: eqxAddSnp
   '  - Adds an SNP entry to an eqx buffer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(formatFlag)
   { /*Switch: Check if adding a base position*/
      case defExpandCig:
      /*Case: Using the expanded cigar format*/
         *eqxBuffStr = 'X';
         return;
      /*Case: Using the expanded cigar format*/

      case defEMBOSS:
      case defClustal:
         *eqxBuffStr = ' ';
         return;

     case defFasta: return; /*No eqx entry*/
   } /*Switch: Check if adding starting bases position*/

   return;
} /*eqxAddSnp*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o eqxBuffStr to have the correct match entry based on
|      formatFlag.
\--------------------------------------------------------*/
void eqxAddMatch(
   char *eqxBuffStr,   /*Buffer with eqx entry*/
   char formatFlag     /*Has format to add soft mask as*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: eqxAddMatch
   '  - Adds an match entry to an eqx buffer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(formatFlag)
   { /*Switch: Check if adding a base position*/
      case defExpandCig:
      /*Case: Using the expanded cigar format*/
         *eqxBuffStr = '=';
         return;
      /*Case: Using the expanded cigar format*/

      case defEMBOSS:
         *eqxBuffStr = '|';
         return;
      case defClustal:
         *eqxBuffStr = '*';
         return;

     case defFasta: return; /*No eqx entry*/
   } /*Switch: Check if adding starting bases position*/

   return;
} /*eqxAddMatch*/
