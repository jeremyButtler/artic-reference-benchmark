/*########################################################
# Name memWater
# Use:
#  o Holds functions doing a memory efficent Smith Waterman
#    pairwise alignments. These aligners return positions,
#    which can be then run through a global/local aligner
#    used to find the actual alignment.
# Includes:
#   - "waterman.h"
#   o "generalAlnFun.h"
#   o "alnStruct.h"
#   o "alnMatrixStruct.h"
#   o "twoBitArrays.h"
#   o "scoresST.h"
#   o "seqStruct.h"
#   o "alnSetStruct.h"
#   o "alnSeqDefaults.h"
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>
#   0 <string.h>
########################################################*/

#include "memWater.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' memWater SOF: Start Of Functions
' o fun-01 memWaterAln:
'   - Run a memory efficent Waterman Smith alignment on
'     input sequences
' o fun-02 mamWaterAltAln:
'   - Run a memory efficent Waterman Smith alignment that
'     returns alternative alignmetns
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: memWaterAln
| Call: memWaterAln(qryST, refST, settings);
| Use:
|   - Performs a memory efficent Smith Waterman alignment
|     on a pair of sequences
| Input;
|   - qryST:
|     o SeqStruct with the query sequence and index 0
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment.
|   - refST:
|     o SeqStruct with the reference sequence and index 0
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment
|   - settings:
|     o alnSet structure with the setttings, such as
|       gap open, gap extend, scoring matrix, and preffered
|       direction.
| Output:
|  - Returns:
|    o scoresStruct with the best score
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
struct scoresStruct * memWaterAln(
    struct seqStruct *qryST, /*query sequence*/
    struct seqStruct *refST, /*ref sequence*/
    struct alnSet *settings  /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: memWaterAln
   '  - Run a memory efficent Waterman Smith alignment on
   '    input sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in initial negatives for ref
   '  o fun-01 sec-04:
   '    - Fill the matrix with scores
   '  o fun-01 sec-05:
   '    - Set up for returing the matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-01: Variable declerations
   ^  o fun-01 sec-01 sub-01:
   ^    - Variables dealing with the query and reference
   ^      starting positions
   ^  o fun-01 sec-01 sub-02:
   ^    - Variables holding the scores (only two rows)
   ^  o fun-01 sec-01 sub-03:
   ^    - Directinol matrix variables
   ^  o fun-01 sec-01 sub-04:
   ^    - Variables for building returend alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-01 Sub-01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \******************************************************/

   /*Get start & end of the query and reference sequences*/
   char *refStartStr = refST->seqCStr + refST->offsetUL;
   char *refEndStr = refST->seqCStr + refST->endAlnUL;

   char *qryStartStr = qryST->seqCStr + qryST->offsetUL;
   char *qryEndStr = qryST->seqCStr + qryST->endAlnUL;

   char *qryIterStr = 0;
   char *refIterStr = 0;

   unsigned long lenRefUL =
       refST->endAlnUL - refST->offsetUL + 1;
     /*The + 1 is to account for index 0 of endAlnUL*/

   /*Set up counters for the query and reference base
   `  index
   */
   /******************************************************\
   * Fun-01 Sec-01 Sub-02:
   *  - Variables holding the scores (only two rows)
   \******************************************************/

   long insScoreL = 0;   /*Score for doing an insertion*/
   long snpScoreL = 0;   /*Score for doing an match/snp*/
   long delScoreL = 0;   /*Score for doing an deletion*/
   long nextSnpSL = 0;   /*Score for the next match/snp*/

   // Marks when to reset score buffer (every second row)
   long *scoreRowLP = 0; /*matrix to use in alignment*/
   long *scoreOnLP = 0;  /*Score I am working on*/

   /******************************************************\
   * Fun-01 Sec-01 Sub-03:
   *  - Directinol matrix variables
   \******************************************************/

   /*Direction matrix (one cell holds a single direction)*/
   #if defined TWOBITMSW
      struct twoBitAry *dirRow = 0;/*Holds directions*/
      unsigned char lastDirC = 0;
      unsigned char lastLastDirC = 0;
   #else
      char *dirRow = 0;  /*Holds directions*/
      char *firstDir = 0; /*Holds directions*/
      char lastDirC = 0;
      char lastLastDirC = 0;
   #endif

   struct scoresStruct *bestScoreST = 0;

   /*For recording the start position*/
   unsigned long *refStartUL = 0;
   unsigned long *qryStartUL = 0;
   unsigned long *refStartFirstIndexUL = 0;
   unsigned long *qryStartFirstIndexUL = 0;
   unsigned long lastRefStartUL = refST->offsetUL;
   unsigned long lastQryStartUL = qryST->offsetUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   ^  o fun-01 sec-02 sub-01:
   ^    - Allocate memory for the alignment
   ^  o fun-01 sec-02 sub-02:
   ^    - Allocate memory for alternative alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-02 Sub-01:
   *  - Allocate memory for the alignment
   \******************************************************/

   /*Make struct array for every base in reference*/
   bestScoreST = calloc(1, sizeof(struct scoresStruct));
   if(bestScoreST == 0) return 0;
   initScoresST(bestScoreST);

   scoreRowLP = calloc((lenRefUL + 1), sizeof(long));
   /*+ 1 is for the indel column*/
   if(scoreRowLP == 0)
   { /*If I had a memory error*/
     freeScoresST(bestScoreST, 0);
     return 0;
   } /*If I had a memory error*/

   #if defined TWOBITMSW
      dirRow = makeTwoBit(lenRefUL+1 , 0);
     /* Calls calloc and adds an extra element at end
     `  - lenRefUL + 1 accounts for insertion reference row
     */
   #else
      dirRow = calloc(lenRefUL + 1 ,sizeof(char));
      firstDir = dirRow;
   #endif

   if(dirRow == 0)
   { /*If I do not have a direction matrix for each cell*/
     freeScoresST(bestScoreST, 0);
     free(scoreRowLP);
     return 0;
   } /*If I do not have a direction matrix for each cell*/

   /******************************************************\
   * Fun-01 Sec-02 Sub-02:
   *  - Allocate memory for alternative alignments
   \******************************************************/

   refStartUL = calloc(lenRefUL,sizeof(unsigned long));
   refStartFirstIndexUL = refStartUL;
   if(refStartUL == 0)
   { /*If had a memory error*/
     freeScoresST(bestScoreST, 0);
      free(scoreRowLP);
      return 0;
   } /*If had a memory error*/

   /*One query position recoreded per refference position*/
   qryStartUL = calloc(lenRefUL, sizeof(unsigned long));
   qryStartFirstIndexUL = qryStartUL;
   if(qryStartUL == 0)
   { /*If had a memory error*/
      freeScoresST(bestScoreST, 0);
      free(scoreRowLP);
      free(refStartUL);
      return 0;
   } /*If had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Fill in the indel column in the indel row*/
   #if defined TWOBITMSW
      changeTwoBitElm(dirRow, defMvStop);
      twoBitMvToNextElm(dirRow);
   #else
      *dirRow = defMvStop;
      ++dirRow;
   #endif

   refIterStr = refStartStr;
   while(refIterStr <= refEndStr)
   { /*loop; till have initalized the first row*/
     #if defined TWOBITMSW
        changeTwoBitElm(dirRow, defMvStop);
        twoBitMvToNextElm(dirRow);
     #else
        *dirRow = defMvStop;
        ++dirRow;
     #endif

     ++scoreOnLP; /*Already set to 0 by calloc*/
     ++refIterStr; /*Move though the next base*/
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-01 sec-04 sub-01:
   ^    - Get the initial scores
   ^  o fun-01 sec-04 sub-02:
   ^    - Do the first move
   ^  o fun-01 sec-04 sub-03:
   ^    - Fill out the matrix
   ^  o fun-01 sec-04 sub-04:
   ^    - Find the next matches score
   ^  o fun-01 sec-04 sub-05:
   ^    - Find the best score for the last round
   ^  o fun-01 sec-04 sub-06:
   ^    - Find the score for the next deletion
   ^  o fun-01 sec-04 sub-07:
   ^    - Check if is an alternative base best score
   *  o fun-01 sec-04 sub-08:
   *    - Check if is an alternative base best score
   ^  o fun-01 sec-4 sub-09:
   ^    - Move to the next reference base
   ^  o fun-01 sec-04 sub-10:
   ^    - Find the scores for the next insertion
   ^  o fun-01 sec-04 sub-11:
   ^    - Find the best score for the last base
   ^  o fun-01 sec-04 sub-12:
   ^    - Is the last base in row an alternative alignment?
   ^  o fun-01 sec-04 sub-13:
   ^    - Move to the indel column
   ^  o fun-01 sec-04 sub-14:
   ^    - Find the scores for the first base in the row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Get the initial scores
   \******************************************************/

   #if defined TWOBITMSW
      twoBitMvXElmFromStart(dirRow, 0);
   #else
      dirRow = firstDir;
   #endif

   qryIterStr = qryStartStr;
   refIterStr = refStartStr;
   scoreOnLP = scoreRowLP;
   refStartUL = refStartFirstIndexUL;
   qryStartUL = qryStartFirstIndexUL;

   nextSnpSL =
        getBaseScore(qryIterStr, refIterStr, settings)
      + *scoreOnLP;

   /*These are always negative*/
   *scoreOnLP = 0;
   delScoreL = 0;
   insScoreL = 0;
   lastLastDirC = 0;
   ++scoreOnLP;

   /******************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Do the first move
   \******************************************************/

   #if defined TWOBITMSW
      changeTwoBitElm(dirRow, defMvStop);
      twoBitMvToNextElm(dirRow);
   #else
     *dirRow = defMvStop;
     ++dirRow;
   #endif

   /******************************************************\
   * Fun-01 Sec-04 Sub-03:
   *  - Fill out the matrix
   \******************************************************/

   /*Starting on the first sequence row*/
   while(qryIterStr <= qryEndStr)
   { /*loop; compare query base against all ref bases*/

     refIterStr = refStartStr + 1;

     /*First reference bases column*/
     while(refIterStr <= refEndStr)
     { /* loop; compare one query to one reference base*/

       /**************************************************\
       * Fun-01 Sec-04 Sub-04:
       *  - Find the next matches score
       \**************************************************/

       snpScoreL = nextSnpSL;

       nextSnpSL =
          getBaseScore(
              qryIterStr,
              refIterStr,
              settings
       ); /*Find the score for the two base pairs*/

       nextSnpSL += *scoreOnLP;

       /**************************************************\
       * Fun-01 Sec-04 Sub-06:
       *  - Find the best score for the last round
       \**************************************************/

       lastDirC = lastLastDirC;
       #if defined TWOBITMSW
          lastLastDirC = getTwoBitElm(dirRow);
          waterTwoBitMaxScore(
            dirRow,
            settings,
            &insScoreL,
            &snpScoreL,
            &delScoreL,
            scoreOnLP
          ); /*Update the scores*/

       #else
          lastLastDirC = *dirRow;
          waterByteMaxScore(
            dirRow,
            settings,
            &insScoreL,
            &snpScoreL,
            &delScoreL,
            scoreOnLP
          ); /*Update the scores*/
       #endif
       /**************************************************\
       * Fun-01 Sec-04 Sub-06:
       *  - Find the the next deletion score
       \**************************************************/

       /*Need to move the direction here, so I have
       ` the previous bases direction.
       */
       #if defined TWOBITMSW && !defined NOGAPOPEN
          indelScore(
             delScoreL,
             getTwoBitElm(dirRow),
             *scoreOnLP,
             settings
          );
       #elif !defined NOGAPOPEN
          indelScore(
             delScoreL,
             *dirRow,
             *scoreOnLP,
             settings
          );
       #else
          delScoreL = *scoreOnLP + settings->gapExtendI;
       #endif

       /**************************************************\
       * Fun-01 Sec-04 Sub-07:
       *  - Determine if is a best score (keep as primary)
       \**************************************************/

       updateStartPos(
          #if defined TWOBITMSW
             getTwoBitElm(dirRow),
          #else
             *dirRow,
          #endif
          lastDirC,
          lastRefStartUL,
          lastQryStartUL,
          refStartUL,
          qryStartUL,
          refIterStr - 1,
          refST->seqCStr,
          qryIterStr,
          qryST->seqCStr
       ); /*macro in water.h*/
          
       if(*scoreOnLP > bestScoreST->scoreL)
       { /*If this is the current best score*/
          bestScoreST->scoreL = *scoreOnLP;

          bestScoreST->refStartUL = *refStartUL;
          bestScoreST->qryStartUL = *qryStartUL;

          bestScoreST->refEndUL =
            refIterStr - refST->seqCStr - 1;
          bestScoreST->qryEndUL= qryIterStr-qryST->seqCStr;
       } /*If this is the current best score*/


       /***********************************************\
       * Fun-01 Sec-04 Sub-09:
       *  - Move to next reference base
       \***********************************************/
          
       ++refStartUL;
       ++qryStartUL;
       ++refIterStr; /*Move to next reference base*/
       ++scoreOnLP;  /*Move to the next score*/

       #if defined TWOBITMSW
          twoBitMvToNextElm(dirRow);
       #else
          ++dirRow;
       #endif

       /**************************************************\
       * Fun-01 Sec-04 Sub-10:
       *  - Find the the next insertion score
       \**************************************************/

       /*Find the next insertion score (Is one score
       ` ahead of the just filled score).
       */

       /*Get the new last insertion direction*/
       #if defined TWOBITMSW && !defined NOGAPOPEN
          indelScore(
             insScoreL,
             getTwoBitElm(dirRow),
             *scoreOnLP,
             settings
          );
       #elif !defined NOGAPOPEN
          indelScore(
             insScoreL,
             *dirRow,
             *scoreOnLP,
             settings
          );
       #else
          insScoreL = *scoreOnLP + settings->gapExtendI;
       #endif
     } /*loop; compare one query to one reference base*/

       /**************************************************\
       * Fun-01 Sec-04 Sub-11:
       *  - Find the best score for the last base
       \**************************************************/

       /*Find the best score for the last base. In this
       ` case, this score can only apply to indels. So,
       ` I need to move off it to avoid overwirting it
       */
       lastDirC = lastLastDirC;
       #if defined TWOBITMSW
          waterTwoBitMaxScore(
            dirRow,
            settings,
            &insScoreL,
            &nextSnpSL,
            &delScoreL,
            scoreOnLP
          ); /*Update the score and direction*/
       #else
          waterByteMaxScore(
            dirRow,
            settings,
            &insScoreL,
            &nextSnpSL,
            &delScoreL,
            scoreOnLP
          ); /*Update the score and direction*/
       #endif

     /****************************************************\
     * Fun-01 Sec-04 Sub-12:
     *  - Is the last base in row an alternative alignment?
     \****************************************************/

     updateStartPos(
        #if defined TWOBITMSW
           getTwoBitElm(dirRow),
        #else
           *dirRow,
        #endif
        lastDirC,
        lastRefStartUL,
        lastQryStartUL,
        refStartUL,
        qryStartUL,
        refIterStr - 1,
        refST->seqCStr,
        qryIterStr,
        qryST->seqCStr
     );

     if(*scoreOnLP > bestScoreST->scoreL)
     { /*If this is the current best score*/
        bestScoreST->scoreL = *scoreOnLP;

        bestScoreST->refStartUL = *refStartUL;
        bestScoreST->qryStartUL = *qryStartUL;

        bestScoreST->refEndUL= refIterStr-refST->seqCStr-1;
        bestScoreST->qryEndUL= qryIterStr-qryST->seqCStr;
     } /*If this is the current best score*/

     refStartUL = refStartFirstIndexUL;
     qryStartUL = qryStartFirstIndexUL;

     /**************************************************\
     *  Fun-01 Sec-04 Sub-13:
     *   - Move to the indel column
     \**************************************************/

      /*I need to move dirRow and insDir to the first
      ` first base in the next row (skip indel column).
      ` The score for the indel column is updated in the
      ` next subsection. I need to find the score for a
      ` deletion before the update.
      */
      #if defined TWOBITMSW
         twoBitMvXElmFromStart(dirRow, 0);
         changeTwoBitElm(dirRow, defMvIns);
         twoBitMvToNextElm(dirRow);
      #else
         dirRow = firstDir;
         *dirRow = defMvIns;
         ++dirRow;
      #endif

      /**************************************************\
      * Fun-01 Sec-04 Sub-14:
      *  - Find the scores for the first base in the row
      \**************************************************/

      /*Move to indel column and apply gap extension*/
      scoreOnLP = scoreRowLP;
      ++qryIterStr; /*Move to the next query base*/

      /*Find the next score for an snp/match*/
      nextSnpSL = 
           getBaseScore(qryIterStr,refStartStr,settings)
         + *scoreOnLP;

      /*Update the indel column and find next deletion*/
      lastLastDirC = 0;
      *scoreOnLP += settings->gapExtendI;
      delScoreL = *scoreOnLP + settings->gapExtendI;
      ++scoreOnLP; /*Move to the first base pair*/

     #if defined TWOBITMSW && !defined NOGAPOPEN
        indelScore(
           insScoreL,
           getTwoBitElm(dirRow),
           *scoreOnLP,
           settings
        );
     #elif !defined NOGAPOPEN
        indelScore(
           insScoreL,
           *dirRow,
           *scoreOnLP,
           settings
        );
     #else
        insScoreL = *scoreOnLP + settings->gapExtendI;
     #endif
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Move back to the lower right conor cell
   ` This is not needed, but is nice.
   */
   #if defined TWOBITMSW
      twoBitMvBackOneElm(dirRow);
      changeTwoBitElm(dirRow, defMvStop);
      twoBitMvBackOneElm(dirRow);
      freeTwoBit(dirRow, 0, 0); /*0 frees everything*/
   #else
      *(dirRow - 1) = defMvStop;
      free(firstDir);
   #endif

   free(scoreRowLP);
   free(qryStartUL);
   free(refStartUL);
   
   scoreRowLP = 0;
   qryStartUL = 0;
   refStartUL = 0;

   return bestScoreST;
} /*memWaterAln*/

/*--------------------------------------------------------\
| Name: memWaterAltAln
| Call: memWaterAltAln(qryST, refST, settings);
| Use:
|   - Performs a memory efficent Smith Waterman alignment
|     on a pair of sequences
| Input;
|   - qryST:
|     o SeqStruct with the query sequence and index 0
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment
|   - refST:
|     o SeqStruct with the reference sequence and
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment
|   - settings:
|     o alnSet structure with the setttings, such as
|       gap open, gap extend, min score, scoring matrix,
|       and preffered direction.
| Output:
|  - Returns:
|    o alnMatrixStruct with the best score and alternative
|      alignments (includes best score)
|      - This does not return a direction matrix
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
struct alnMatrixStruct * memWaterAltAln(
    struct seqStruct *qryST,   /*query sequence and data*/
    struct seqStruct *refST,   /*ref sequence and data*/
      /* both qryST and refST have the sequence,
      `  they also have the point to start the alignment
      `  seqST->offsetUL (index 0) and the point to end
      `  the alignment seqST->endAlnUL (index 0).
      */
    struct alnSet *settings /*Settings for the alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: memWaterAltAln
   '  - Run a memory efficent Waterman Smith alignment that
   '    returns alternative alignmetns
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

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-01: Variable declerations
   ^  o fun-02 sec-01 sub-01:
   ^    - Variables dealing with the query and reference
   ^      starting positions
   ^  o fun-02 sec-01 sub-02:
   ^    - Variables holding the scores (only two rows)
   ^  o fun-02 sec-01 sub-03:
   ^    - Directinol matrix variables
   ^  o fun-02 sec-01 sub-04:
   ^    - Variables for building returend alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-02 Sec-01 Sub-01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \******************************************************/

   /*Get start & end of the query and reference sequences*/
   char *refStartStr = refST->seqCStr + refST->offsetUL;
   char *refEndStr = refST->seqCStr + refST->endAlnUL;

   char *qryStartStr = qryST->seqCStr + qryST->offsetUL;
   char *qryEndStr = qryST->seqCStr + qryST->endAlnUL;

   char *qryIterStr = 0;
   char *refIterStr = 0;

   /*Find the length of the reference and query*/
   unsigned long lenQryUL =
       qryST->endAlnUL - qryST->offsetUL + 1;
     /*The + 1 is to account for index 0 of endAlnUL*/

   unsigned long lenRefUL =
       refST->endAlnUL - refST->offsetUL + 1;
     /*The + 1 is to account for index 0 of endAlnUL*/

   /*Set up counters for the query and reference base
   `  index
   */
   /******************************************************\
   * Fun-02 Sec-01 Sub-02:
   *  - Variables holding the scores (only two rows)
   \******************************************************/

   long insScoreL = 0;   /*Score for doing an insertion*/
   long snpScoreL = 0;   /*Score for doing an match/snp*/
   long delScoreL = 0;   /*Score for doing an deletion*/
   long nextSnpSL = 0;   /*Score for the next match/snp*/

   // Marks when to reset score buffer (every second row)
   long *scoreRowLP = 0; /*matrix to use in alignment*/
   long *scoreOnLP = 0;  /*Score I am working on*/

   /******************************************************\
   * Fun-02 Sec-01 Sub-03:
   *  - Directinol matrix variables
   \******************************************************/

   /*Direction matrix (one cell holds a single direction)*/
   #if defined TWOBITMSW
      struct twoBitAry *dirRow = 0;/*Holds directions*/
      unsigned char lastDirC = 0;
      unsigned char lastLastDirC = 0;
   #else
      char *dirRow = 0;  /*Holds directions*/
      char *firstDir = 0; /*Holds directions*/
      char lastDirC = 0;
      char lastLastDirC = 0;
   #endif

   /*The structure to return (has results)*/
   struct alnMatrixStruct *retMtxST = 0;

   /*For recording alternative alignments*/
   struct scoresStruct *qryBasesST = 0;
   struct scoresStruct *refBasesST = 0;
   struct scoresStruct *bestScoreST = 0;

   unsigned long *refStartUL = 0;
   unsigned long *qryStartUL = 0;
   unsigned long *refStartFirstIndexUL = 0;
   unsigned long *qryStartFirstIndexUL = 0;
   unsigned long lastRefStartUL = refST->offsetUL;
   unsigned long lastQryStartUL = qryST->offsetUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-02:
   ^  - Allocate memory for alignment
   ^  o fun-02 sec-02 sub-01:
   ^    - Allocate memory for the alignment
   ^  o fun-02 sec-02 sub-02:
   ^    - Allocate memory for alternative alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-01:
   *  - Allocate memory for the alignment
   \******************************************************/

   retMtxST = calloc(1, sizeof(struct alnMatrixStruct));
   if(retMtxST == 0) return 0;
   initAlnMatrixST(retMtxST);
   /*+ 1 is for the indel column*/

   scoreRowLP = calloc((lenRefUL + 1), sizeof(long));
   if(scoreRowLP == 0)
   { // If I had a memory error
     freeAlnMatrixST(retMtxST);
     return 0;
   } // If I had a memory error

   #if defined TWOBITMSW
      dirRow = makeTwoBit(lenRefUL+1 , 0);
     /* Calls calloc and adds an extra element at end
     `  - lenRefUL + 1 accounts for insertion reference row
     */
   #else
      dirRow = calloc(lenRefUL + 1 ,sizeof(char));
      firstDir = dirRow;
   #endif

   if(dirRow == 0)
   { /*If I do not have a direction matrix for each cell*/
     freeAlnMatrixST(retMtxST);
     free(scoreRowLP);
     return 0;
   } /*If I do not have a direction matrix for each cell*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-02:
   *  - Allocate memory for alternative alignments
   \******************************************************/

   refStartUL = calloc(lenRefUL,sizeof(unsigned long));
   refStartFirstIndexUL = refStartUL;
   if(refStartUL == 0)
   { /*If had a memory error*/
      freeAlnMatrixST(retMtxST);
      free(scoreRowLP);
      return 0;
   } /*If had a memory error*/

   /*One query position recoreded per refference position*/
   qryStartUL = calloc(lenRefUL,sizeof(unsigned long));
   qryStartFirstIndexUL = qryStartUL;
   if(qryStartUL == 0)
   { /*If had a memory error*/
      freeAlnMatrixST(retMtxST);
      free(scoreRowLP);
      free(refStartUL);
      return 0;
   } /*If had a memory error*/

   /*Make struct array for every base in reference*/
   retMtxST->refBasesST =
     calloc(lenRefUL, sizeof(struct scoresStruct));
   if(retMtxST->refBasesST == 0)
   { /*If had memory error*/
      freeAlnMatrixST(retMtxST);
      free(scoreRowLP);
      free(refStartUL);
      free(qryStartUL);
      return 0;
   } /*If had memory error*/

   /*Make struct array for every base in the query*/
   retMtxST->qryBasesST =
     calloc(lenQryUL, sizeof(struct scoresStruct));
   if(retMtxST->qryBasesST == 0)
   { /*If had memory error*/
     freeAlnMatrixST(retMtxST);
     free(scoreRowLP);
     free(refStartUL);
     free(qryStartUL);
     return 0;
   } /*If had memory error*/

   qryBasesST = retMtxST->qryBasesST;
   refBasesST = retMtxST->refBasesST;
   retMtxST->lenRefScoresUL = lenRefUL;
   retMtxST->lenQueryScoresUL = lenQryUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Fill in the indel column in the indel row*/
   #if defined TWOBITMSW
      changeTwoBitElm(dirRow, defMvStop);
      twoBitMvToNextElm(dirRow);
   #else
      *dirRow = defMvStop;
      ++dirRow;
   #endif

   refIterStr = refStartStr;
   while(refIterStr <= refEndStr)
   { /*loop; till have initalized the first row*/
     #if defined TWOBITMSW
        changeTwoBitElm(dirRow, defMvStop);
        twoBitMvToNextElm(dirRow);
     #else
        *dirRow = defMvStop;
        ++dirRow;
     #endif

     refBasesST->refStartUL = refIterStr - refST->seqCStr;
     refBasesST->qryStartUL = qryST->offsetUL;

     ++refBasesST;

     ++scoreOnLP; /*Already set to 0 by calloc*/
     ++refIterStr; /*Move though the next base*/
   } /*loop; till have initalized the first row*/

   qryIterStr = qryStartStr;
   while(qryIterStr <= qryEndStr)
   { /*Loop: Initalize the query scores*/
     qryBasesST->refStartUL = refIterStr - refST->seqCStr;
     qryBasesST->qryStartUL = qryST->offsetUL;
     ++qryBasesST;
     ++qryIterStr;
   } /*Loop: Initalize the query scores*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-02 sec-04 sub-01:
   ^    - Get the initial scores
   ^  o fun-02 sec-04 sub-02:
   ^    - Do the first move
   ^  o fun-02 sec-04 sub-03:
   ^    - Fill out the matrix
   ^  o fun-02 sec-04 sub-04:
   ^    - Find the next matches score
   ^  o fun-02 sec-04 sub-05:
   ^    - Find the best score for the last round
   ^  o fun-02 sec-04 sub-06:
   ^    - Find the score for the next deletion
   ^  o fun-02 sec-04 sub-07:
   ^    - Check if is an alternative base best score
   *  o fun-02 sec-04 sub-08:
   *    - Check if is an alternative base best score
   ^  o fun-02 sec-4 sub-09:
   ^    - Move to the next reference base
   ^  o fun-02 sec-04 sub-10:
   ^    - Find the scores for the next insertion
   ^  o fun-02 sec-04 sub-11:
   ^    - Find the best score for the last base
   ^  o fun-02 sec-04 sub-12:
   ^    - Is the last base in row an alternative alignment?
   ^  o fun-02 sec-04 sub-13:
   ^    - Move to the indel column
   ^  o fun-02 sec-04 sub-14:
   ^    - Find the scores for the first base in the row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-02 Sec-04 Sub-01:
   *  - Get the initial scores
   \******************************************************/

   #if defined TWOBITMSW
      twoBitMvXElmFromStart(dirRow, 0);
   #else
      dirRow = firstDir;
   #endif

   qryIterStr = qryStartStr;
   refIterStr = refStartStr;
   scoreOnLP = scoreRowLP;
   refStartUL = refStartFirstIndexUL;
   qryStartUL = qryStartFirstIndexUL;

   qryBasesST = retMtxST->qryBasesST;
   refBasesST = retMtxST->refBasesST;

   nextSnpSL =
        getBaseScore(qryIterStr, refIterStr, settings)
      + *scoreOnLP;

   /*These are always negative*/
   *scoreOnLP = 0;
   delScoreL = 0;
   insScoreL = 0;
   lastLastDirC = 0;
   ++scoreOnLP;

   /******************************************************\
   * Fun-02 Sec-04 Sub-01:
   *  - Do the first move
   \******************************************************/

   #if defined TWOBITMSW
      changeTwoBitElm(dirRow, defMvStop);
      twoBitMvToNextElm(dirRow);
   #else
     *dirRow = defMvStop;
     ++dirRow;
   #endif

   /******************************************************\
   * Fun-02 Sec-04 Sub-03:
   *  - Fill out the matrix
   \******************************************************/

   /*Starting on the first sequence row*/
   while(qryIterStr <= qryEndStr)
   { /*loop; compare query base against all ref bases*/

     refIterStr = refStartStr + 1;

     /*First reference bases column*/
     while(refIterStr <= refEndStr)
     { /* loop; compare one query to one reference base*/

       /**************************************************\
       * Fun-02 Sec-04 Sub-04:
       *  - Find the next matches score
       \**************************************************/

       snpScoreL = nextSnpSL;

       nextSnpSL =
          getBaseScore(
              qryIterStr,
              refIterStr,
              settings
       ); /*Find the score for the two base pairs*/

       nextSnpSL += *scoreOnLP;

       /**************************************************\
       * Fun-02 Sec-04 Sub-06:
       *  - Find the best score for the last round
       \**************************************************/

       lastDirC = lastLastDirC;
       #if defined TWOBITMSW
          lastLastDirC = getTwoBitElm(dirRow);
          waterTwoBitMaxScore(
            dirRow,
            settings,
            &insScoreL,
            &snpScoreL,
            &delScoreL,
            scoreOnLP
          ); /*Update the scores*/

       #else
          lastLastDirC = *dirRow;
          waterByteMaxScore(
            dirRow,
            settings,
            &insScoreL,
            &snpScoreL,
            &delScoreL,
            scoreOnLP
          ); /*Update the scores*/
       #endif

       /**************************************************\
       * Fun-02 Sec-04 Sub-06:
       *  - Find the the next deletion score
       \**************************************************/

       /*Need to move the direction here, so I have
       ` the previous bases direction.
       */
       #if defined TWOBITMSW && !defined NOGAPOPEN
          indelScore(
             delScoreL,
             getTwoBitElm(dirRow),
             *scoreOnLP,
             settings
          );
       #elif !defined NOGAPOPEN
          indelScore(
             delScoreL,
             *dirRow,
             *scoreOnLP,
             settings
          );
       #else
          delScoreL = *scoreOnLP + settings->gapExtendI;
       #endif

       /**************************************************\
       * Fun-02 Sec-04 Sub-07:
       *  - Determine if is a best score (keep as primary)
       \**************************************************/

          updateStartPos(
             #if defined TWOBITMSW
                getTwoBitElm(dirRow),
             #else
                *dirRow,
             #endif
             lastDirC,
             lastRefStartUL,
             lastQryStartUL,
             refStartUL,
             qryStartUL,
             refIterStr - 1,
             refST->seqCStr,
             qryIterStr,
             qryST->seqCStr
          );
          
       keepAltScore(
          *scoreOnLP,
          settings->minScoreL,
          qryBasesST,
          refBasesST,
          *refStartUL, /*Is in index 0*/
          refIterStr - refST->seqCStr - 1, 
          *qryStartUL,
          qryIterStr - qryST->seqCStr
       ); /*Macro in waterman.h*/

       /***********************************************\
       * Fun-02 Sec-04 Sub-09:
       *  - Move to next reference base
       \***********************************************/
          
       ++refBasesST;
       ++refStartUL;
       ++qryStartUL;
       ++refIterStr; /*Move to next reference base*/
       ++scoreOnLP;  /*Move to the next score*/

       #if defined TWOBITMSW
          twoBitMvToNextElm(dirRow);
       #else
          ++dirRow;
       #endif

       /**************************************************\
       * Fun-02 Sec-04 Sub-10:
       *  - Find the the next insertion score
       \**************************************************/

       /*Find the next insertion score (Is one score
       ` ahead of the just filled score).
       */

       /*Get the new last insertion direction*/
       #if defined TWOBITMSW && !defined NOGAPOPEN
          indelScore(
             insScoreL,
             getTwoBitElm(dirRow),
             *scoreOnLP,
             settings
          );
       #elif !defined NOGAPOPEN
          indelScore(
             insScoreL,
             *dirRow,
             *scoreOnLP,
             settings
          );
       #else
          insScoreL = *scoreOnLP + settings->gapExtendI;
       #endif
     } /*loop; compare one query to one reference base*/

       /**************************************************\
       * Fun-02 Sec-04 Sub-11:
       *  - Find the best score for the last base
       \**************************************************/

       /*Find the best score for the last base. In this
       ` case, this score can only apply to indels. So,
       ` I need to move off it to avoid overwirting it
       */
       lastDirC = lastLastDirC;
       #if defined TWOBITMSW
          waterTwoBitMaxScore(
            dirRow,
            settings,
            &insScoreL,
            &nextSnpSL,
            &delScoreL,
            scoreOnLP
          ); /*Update the score and direction*/
       #else
          waterByteMaxScore(
            dirRow,
            settings,
            &insScoreL,
            &nextSnpSL,
            &delScoreL,
            scoreOnLP
          ); /*Update the score and direction*/
       #endif

     /****************************************************\
     * Fun-02 Sec-04 Sub-12:
     *  - Is the last base in row an alternative alignment?
     \****************************************************/

     updateStartPos(
        #if defined TWOBITMSW
           getTwoBitElm(dirRow),
        #else
           *dirRow,
        #endif
        lastDirC,
        lastRefStartUL,
        lastQryStartUL,
        refStartUL,
        qryStartUL,
        refIterStr - 1,
        refST->seqCStr,
        qryIterStr,
        qryST->seqCStr
     );

     keepAltScore(
        *scoreOnLP,
        settings->minScoreL,
        qryBasesST,
        refBasesST,
        *refStartUL,
        refIterStr - refST->seqCStr - 1, 
        *qryStartUL,
        qryIterStr - qryST->seqCStr
     ); /*Macro in waterman.h*/

     ++qryBasesST;
     refBasesST = retMtxST->refBasesST;

     refStartUL = refStartFirstIndexUL;
     qryStartUL = qryStartFirstIndexUL;

     /**************************************************\
     *  Fun-02 Sec-04 Sub-13:
     *   - Move to the indel column
     \**************************************************/

      /*I need to move dirRow and insDir to the first
      ` first base in the next row (skip indel column).
      ` The score for the indel column is updated in the
      ` next subsection. I need to find the score for a
      ` deletion before the update.
      */
      #if defined TWOBITMSW
         twoBitMvXElmFromStart(dirRow, 0);
         changeTwoBitElm(dirRow, defMvIns);
         twoBitMvToNextElm(dirRow);
      #else
         dirRow = firstDir;
         *dirRow = defMvIns;
         ++dirRow;
      #endif

      /**************************************************\
      * Fun-02 Sec-04 Sub-14:
      *  - Find the scores for the first base in the row
      \**************************************************/

      /*Move to indel column and apply gap extension*/
      scoreOnLP = scoreRowLP;
      ++qryIterStr; /*Move to the next query base*/

      /*Find the next score for an snp/match*/
      nextSnpSL = 
           getBaseScore(qryIterStr,refStartStr,settings)
         + *scoreOnLP;

      /*Update the indel column and find next deletion*/
      lastLastDirC = 0;
      *scoreOnLP += settings->gapExtendI;
      delScoreL = *scoreOnLP + settings->gapExtendI;
      ++scoreOnLP; /*Move to the first base pair*/

     #if defined TWOBITMSW && !defined NOGAPOPEN
        indelScore(
           insScoreL,
           getTwoBitElm(dirRow),
           *scoreOnLP,
           settings
        );
     #elif !defined NOGAPOPEN
        indelScore(
           insScoreL,
           *dirRow,
           *scoreOnLP,
           settings
        );
     #else
        insScoreL = *scoreOnLP + settings->gapExtendI;
     #endif
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   ^  o fun-02 sec-05 sub-01:
   ^    - clean up
   ^  o fun-02 sec-05 sub-02:
   ^    - find the best score
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-02 Sec-05 Sub-01:
   *  - clean up
   \******************************************************/

   /*Move back to the lower right conor cell
   ` This is not needed, but is nice.
   */
   #if defined TWOBITMSW
      twoBitMvBackOneElm(dirRow);
      changeTwoBitElm(dirRow, defMvStop);
      twoBitMvBackOneElm(dirRow);
      freeTwoBit(dirRow, 0, 0); /*0 frees everything*/
   #else
      *(dirRow - 1) = defMvStop;
      free(firstDir);
   #endif

   free(scoreRowLP);
   free(qryStartUL);
   free(refStartUL);
   
   scoreRowLP = 0;
   qryStartUL = 0;
   refStartUL = 0;

   /******************************************************\
   * Fun-02 Sec-05 Sub-02:
   *  - Find the best score
   \******************************************************/

   bestScoreST = &retMtxST->bestScoreST;
   refBasesST = retMtxST->refBasesST;
   for(unsigned long iterUL=0; iterUL < lenRefUL; ++iterUL)
   { /*Loop: Find the highest reference score*/
      if(bestScoreST->scoreL < refBasesST->scoreL)
      { /*If I have a new best score*/
         bestScoreST->scoreL = refBasesST->scoreL;
         bestScoreST->refStartUL = refBasesST->refStartUL;
         bestScoreST->refEndUL = refBasesST->refEndUL;
         bestScoreST->qryStartUL = refBasesST->qryStartUL;
         bestScoreST->qryEndUL = refBasesST->qryEndUL;
      } /*If I have a new best score*/

      ++refBasesST;
   } /*Loop: Find the highest reference score*/

   qryBasesST = retMtxST->qryBasesST;
   for(unsigned long iterUL=0; iterUL < lenQryUL; ++iterUL)
   { /*Loop: Find the highest reference score*/
      if(bestScoreST->scoreL < qryBasesST->scoreL)
      { /*If I have a new best score*/
         bestScoreST->scoreL = qryBasesST->scoreL;
         bestScoreST->refStartUL = qryBasesST->refStartUL;
         bestScoreST->refEndUL = qryBasesST->refEndUL;
         bestScoreST->qryStartUL = qryBasesST->qryStartUL;
         bestScoreST->qryEndUL = qryBasesST->qryEndUL;
      } /*If I have a new best score*/

      ++qryBasesST;
   } /*Loop: Find the highest reference score*/

   return retMtxST;
} /*memWaterAltAln*/

