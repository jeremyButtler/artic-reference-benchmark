/*#########################################################
# Name waterman
# Use:
#  o Holds functions doing a Waterman-Smith pairwise
#    alignments. This version outputs a single alignment
# Libraries:
#   - "generalAlnFun.h"
#   - "alnStruct.h"
#   - "alnMatrixStruct.h"
#   o "twoBitArrays.h"
#   o "scoresST.h"
#   o "seqStruct.h"
#   o "alnSetStruct.h"
#   o "alnSeqDefaults.h"
# C Standard libraries Used:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>  // by alnSetStructure.h
#   - <string.h>
#########################################################*/

#include "waterman.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-01 WatermanSmithAln:
'     o Perform a Waterman Smith alignment on input
'       sequences
'  o fun-02 addBestBaseScore:
'    - Adds a score and index to the kept scores list
'  o fun-03 printMatrixCig:
'    - Prints out a cigar for an single path in a
'      direction matrix
'  o fun-04 printAltWaterAlns:
'    - Prints out the best aligment and the saved
'       alterantive alignments  (best alignment for each
'       base) to a file
'  - fun-05 updateDirScoreWaterSingle:
'     o Picks the best score and direction for the current
'       base pairs being compared in a Waterman-Smith
'       alignment
'    - Inlined function is in header at bottom
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnMatrixStruct with the direction matrix and scores
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
struct alnMatrixStruct * WatermanAln(
    struct seqStruct *qryST, /*query sequence and data*/
    struct seqStruct *refST,   /*ref sequence and data*/
      /* both qryST and refST have the sequence,
      `  they also have the point to start the alignment
      `  seqST->offsetUL (index 0) and the point to end
      `  the alignment seqST->endAlnUL (index 0).
      */
    struct alnSet *settings /*Settings for the alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: WatermanAln
   '  - Run a Waterman Smith alignment on input sequences
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
   #if !defined BYTEMATRIX && !defined NOGAPOPEN
      struct twoBitAry *dirMatrix = 0;/*Direction matrix*/
      struct twoBitAry insDir;     /*Direction above cell*/
   #elif !defined BYTEMATRIX
      struct twoBitAry *dirMatrix = 0;/*Direction matrix*/
   #elif !defined NOGAPOPEN 
      char *dirMatrix = 0;/*Direction matrix*/
      char *insDir;       /*Direction above cell*/
   #else
      char *dirMatrix = 0;    /*Direction matrix*/
   #endif

   /*The structure to return (has results)*/
   struct alnMatrixStruct *retMtxST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

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

   #if !defined BYTEMATRIX
      dirMatrix= makeTwoBit((lenRefUL+1) * (lenQryUL+1),0);
     /* Calls calloc and adds an extra element at end
     `  - lenRefUL + 1 accounts for insertion reference row
     `  - lenQeurI + 1 accounts for insertion query column
     */
   #else
      dirMatrix =
         calloc((lenRefUL+1) *(lenQryUL+1)+2,sizeof(char));
   #endif

   if(dirMatrix == 0)
   { /*If I do not have a direction matrix for each cell*/
     freeAlnMatrixST(retMtxST);
     free(scoreRowLP);
     return 0;
   } /*If I do not have a direction matrix for each cell*/

   retMtxST->dirMatrixST = dirMatrix;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Get the first indel position*/
   #if !defined BYTEMATRIX && !defined NOGAPOPEN
      cpTwoBitPos(dirMatrix, &insDir);
   #elif !defined NOGAPOPEN
      insDir = dirMatrix;
   #endif

   refIterStr = refStartStr - 1;
   while(refIterStr <= refEndStr)
   { /*loop; till have initalized the first row*/
     #if !defined BYTEMATRIX
        changeTwoBitElm(dirMatrix, defMvStop);
        twoBitMvToNextElm(dirMatrix);
     #else
        *dirMatrix = defMvStop;
        ++dirMatrix;
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
   *    - Does not exist, but exists in WatermanAltAln
   ^  o fun-01 sec-04 sub-09:
   ^    - Find the scores for the next insertion
   ^  o fun-01 sec-4 sub-10:
   ^    - Move to the next reference base
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

   qryIterStr = qryStartStr;
   refIterStr = refStartStr;
   scoreOnLP = scoreRowLP;

   nextSnpSL =
        getBaseScore(qryIterStr, refIterStr, settings)
      + *scoreOnLP;

   /*These are always negative*/
   *scoreOnLP = 0;
   delScoreL = 0;
   insScoreL = 0;
   ++scoreOnLP;

   /******************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Do the first move
   \******************************************************/

   #if !defined BYTEMATRIX && !defined NOGAPOPEN
      changeTwoBitElm(dirMatrix, defMvStop);
      twoBitMvToNextElm(dirMatrix);
      twoBitMvToNextElm(&insDir);
      twoBitMvToNextElm(&insDir);
   #elif !defined BYTEMATRIX
      changeTwoBitElm(dirMatrix, defMvStop);
      twoBitMvToNextElm(dirMatrix);
   #elif !defined NOGAPOPEN
     *dirMatrix = defMvStop;
     ++dirMatrix;
     insDir += 2;
   #else
     *dirMatrix = defMvStop;
     ++dirMatrix;
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

       #if !defined BYTEMATRIX
          waterTwoBitMaxScore(
            dirMatrix,
            settings,
            &insScoreL,
            &snpScoreL,
            &delScoreL,
            scoreOnLP
          ); /*Update the scores*/

       #else
          waterByteMaxScore(
            dirMatrix,
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
       #if !defined BYTEMATRIX && !defined NOGAPOPEN
          indelScore(
             delScoreL,
             getTwoBitElm(dirMatrix),
             *scoreOnLP,
             settings
          );
       #elif !defined NOGAPOPEN
          indelScore(
             delScoreL,
             *dirMatrix,
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

        /*This is faster than the branchless option.
        ` I am guessing to much is done in the if and 
        ` that the if if fired rarely.
        */
        if(retMtxST->bestScoreST.scoreL < *scoreOnLP)
        { /*if have a new best score*/
           retMtxST->bestScoreST.scoreL = *scoreOnLP;

           /*This will slow me down a bit, but makes life
           ` easier
           */
           retMtxST->bestScoreST.refEndUL =
              refIterStr - refST->seqCStr - 1;

           retMtxST->bestScoreST.qryEndUL =
              qryIterStr - qryST->seqCStr;
        } /*if have a new best score*/

       /**************************************************\
       * Fun-01 Sec-04 Sub-09:
       *  - Find the the next insertion score
       \**************************************************/

       /*Find the next insertion score (Is one score
       ` ahead of the just filled score).
       */

       ++scoreOnLP;

       #if !defined BYTEMATRIX && !defined NOGAPOPEN
          indelScore(
             insScoreL,
             getTwoBitElm(&insDir),
             *scoreOnLP,
             settings
          );
       #elif !defined NOGAPOPEN
          indelScore(
             insScoreL,
             *insDir,
             *scoreOnLP,
             settings
          );
       #else
          insScoreL = *scoreOnLP + settings->gapExtendI;
       #endif

       /***********************************************\
       * Fun-01 Sec-04 Sub-10:
       *  - Move to next reference base
       \***********************************************/
     
       /*Move to the next cell to score*/
       ++refIterStr; /*Move to next reference base*/

       /*Move to the next direction*/
       #if !defined BYTEMATRIX && !defined NOGAPOPEN
          twoBitMvToNextElm(dirMatrix);
          twoBitMvToNextElm(&insDir);
       #elif !defined BYTEMATRIX
          twoBitMvToNextElm(dirMatrix);
       #elif !defined NOGAPOPEN
          ++dirMatrix;
          ++insDir;
       #else
          ++dirMatrix;
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
       #if !defined BYTEMATRIX
          waterTwoBitMaxScore(
            dirMatrix,
            settings,
            &insScoreL,
            &nextSnpSL,
            &delScoreL,
            scoreOnLP
          ); /*Update the score and direction*/
       #else
          waterByteMaxScore(
            dirMatrix,
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

     /*This is faster than the branchless option.
     ` I am guessing to much is done in the if and 
     ` that the if if fired rarely.
     */
     if(retMtxST->bestScoreST.scoreL < *scoreOnLP)
     { /*if have a new best score*/
        retMtxST->bestScoreST.scoreL = *scoreOnLP;

        retMtxST->bestScoreST.refEndUL =
           refIterStr - refST->seqCStr - 1;

        retMtxST->bestScoreST.qryEndUL =
           qryIterStr - qryST->seqCStr;
     } /*if have a new best score*/

      /**************************************************\
      *  Fun-01 Sec-04 Sub-13:
      *   - Move to the indel column
      \**************************************************/

      /*I need to move dirMatrix and insDir to the first
      ` first base in the next row (skip indel column).
      ` The score for the indel column is updated in the
      ` next subsection. I need to find the score for a
      ` deletion before the update.
      */
      #if !defined BYTEMATRIX && !defined NOGAPOPEN
         twoBitMvToNextElm(dirMatrix);
         changeTwoBitElm(dirMatrix, defMvIns);
         twoBitMvToNextElm(dirMatrix);

         twoBitMvToNextElm(&insDir);
      #elif !defined BYTEMATRIX
         twoBitMvToNextElm(dirMatrix);
         changeTwoBitElm(dirMatrix, defMvIns);
         twoBitMvToNextElm(dirMatrix);
      #elif !defined NOGAPOPEN
         ++dirMatrix;
         *dirMatrix = defMvIns;
         ++dirMatrix;
         ++insDir;
      #else
         ++dirMatrix;
         *dirMatrix = defMvIns;
         ++dirMatrix;
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
      *scoreOnLP += settings->gapExtendI;
      delScoreL = *scoreOnLP + settings->gapExtendI;
      ++scoreOnLP; /*Move to the first base pair*/

     /*At this point insDir is on the first base*/
     #if !defined BYTEMATRIX && !defined NOGAPOPEN
        indelScore(
           insScoreL,
           getTwoBitElm(&insDir),
           *scoreOnLP,
           settings
        );

        twoBitMvToNextElm(&insDir);
     #elif !defined NOGAPOPEN
        indelScore(
           insScoreL,
           *insDir,
           *scoreOnLP,
           settings
        );

        ++insDir;
     #else
        insScoreL = *scoreOnLP + settings->gapExtendI;
     #endif
     /*At this piont insDir is on the second base*/
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Move back to the lower right conor cell
   ` This is not needed, but is nice.
   */
   #if !defined BYTEMATRIX
      twoBitMvBackOneElm(dirMatrix);
      changeTwoBitElm(dirMatrix, defMvStop);
      twoBitMvBackOneElm(dirMatrix);
   #else
      *(dirMatrix - 1) = defMvStop;
   #endif

   free(scoreRowLP);
   return retMtxST;
} /*WatermanAln*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnMatrixStruct with the direction matrix, best
|      score, and alternative alignments (includes best
|      score)
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
struct alnMatrixStruct * WatermanAltAln(
    struct seqStruct *qryST, /*query sequence and data*/
    struct seqStruct *refST,   /*ref sequence and data*/
      /* both qryST and refST have the sequence,
      `  they also have the point to start the alignment
      `  seqST->offsetUL (index 0) and the point to end
      `  the alignment seqST->endAlnUL (index 0).
      */
    struct alnSet *settings /*Settings for the alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: WatermanAltAln
   '  - Run a Waterman Smith alignment on input sequences
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
   #if !defined BYTEMATRIX
      struct twoBitAry *dirMatrix = 0;/*Direction matrix*/
      struct twoBitAry insDir;     /*Direction above cell*/
      unsigned char lastDirC = 0;
      unsigned char lastLastDirC = 0;
   #else
      char *dirMatrix = 0;/*Direction matrix*/
      char *insDir = 0;   /*Direction above cell*/
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

   #if !defined BYTEMATRIX
      dirMatrix= makeTwoBit((lenRefUL+1) * (lenQryUL+1),0);
     /* Calls calloc and adds an extra element at end
     `  - lenRefUL + 1 accounts for insertion reference row
     `  - lenQeurI + 1 accounts for insertion query column
     */
   #else
      dirMatrix =
         calloc((lenRefUL+1) *(lenQryUL+1)+2,sizeof(char));
   #endif

   if(dirMatrix == 0)
   { /*If I do not have a direction matrix for each cell*/
     freeAlnMatrixST(retMtxST);
     free(scoreRowLP);
     return 0;
   } /*If I do not have a direction matrix for each cell*/

   retMtxST->dirMatrixST = dirMatrix;

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

   /*Get the first indel position*/
   #if !defined BYTEMATRIX
      cpTwoBitPos(dirMatrix, &insDir);
   #else
      insDir = dirMatrix;
   #endif

   /*Fill in the indel column in the indel row*/
   #if !defined BYTEMATRIX
      changeTwoBitElm(dirMatrix, defMvStop);
      twoBitMvToNextElm(dirMatrix);
   #else
      *dirMatrix = defMvStop;
      ++dirMatrix;
   #endif

   refIterStr = refStartStr;
   while(refIterStr <= refEndStr)
   { /*loop; till have initalized the first row*/
     #if !defined BYTEMATRIX
        changeTwoBitElm(dirMatrix, defMvStop);
        twoBitMvToNextElm(dirMatrix);
     #else
        *dirMatrix = defMvStop;
        ++dirMatrix;
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
   ^  o fun-02 sec-04 sub-09:
   ^    - Find the scores for the next insertion
   ^  o fun-02 sec-4 sub-10:
   ^    - Move to the next reference base
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

   #if !defined BYTEMATRIX
      changeTwoBitElm(dirMatrix, defMvStop);
      twoBitMvToNextElm(dirMatrix);
      twoBitMvToNextElm(&insDir);
      twoBitMvToNextElm(&insDir);
   #else
     *dirMatrix = defMvStop;
     ++dirMatrix;
     insDir += 2;
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
       #if !defined BYTEMATRIX
          lastLastDirC = getTwoBitElm(&insDir);
          waterTwoBitMaxScore(
            dirMatrix,
            settings,
            &insScoreL,
            &snpScoreL,
            &delScoreL,
            scoreOnLP
          ); /*Update the scores*/

       #else
          lastLastDirC = *insDir;
          waterByteMaxScore(
            dirMatrix,
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
       #if !defined BYTEMATRIX && !defined NOGAPOPEN
          indelScore(
             delScoreL,
             getTwoBitElm(dirMatrix),
             *scoreOnLP,
             settings
          );
       #elif !defined NOGAPOPEN
          indelScore(
             delScoreL,
             *dirMatrix,
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
          #if !defined BYTEMATRIX
             getTwoBitElm(dirMatrix),
          #else
             *dirMatrix,
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
          
       ++refBasesST;
       ++refStartUL;
       ++qryStartUL;

       /**************************************************\
       * Fun-02 Sec-04 Sub-09:
       *  - Find the the next insertion score
       \**************************************************/

       /*Find the next insertion score (Is one score
       ` ahead of the just filled score).
       */

       ++scoreOnLP;
         /*Last insertion is now the last snp direction*/

       /*Get the new last insertion direction*/
       #if !defined BYTEMATRIX && !defined NOGAPOPEN
          indelScore(
             insScoreL,
             getTwoBitElm(&insDir),
             *scoreOnLP,
             settings
          );
       #elif !defined NOGAPOPEN
          indelScore(
             insScoreL,
             *insDir,
             *scoreOnLP,
             settings
          );
       #else
          insScoreL = *scoreOnLP + settings->gapExtendI;
       #endif

       /***********************************************\
       * Fun-02 Sec-04 Sub-10:
       *  - Move to next reference base
       \***********************************************/
     
       /*Move to the next cell to score*/
       ++refIterStr; /*Move to next reference base*/

       /*Move to the next direction*/
       #if !defined BYTEMATRIX
          twoBitMvToNextElm(dirMatrix);
          twoBitMvToNextElm(&insDir);
       #else
          ++dirMatrix;
          ++insDir;
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
       #if !defined BYTEMATRIX
          waterTwoBitMaxScore(
            dirMatrix,
            settings,
            &insScoreL,
            &nextSnpSL,
            &delScoreL,
            scoreOnLP
          ); /*Update the score and direction*/
       #else
          waterByteMaxScore(
            dirMatrix,
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
        #if !defined BYTEMATRIX
           getTwoBitElm(dirMatrix),
        #else
           *dirMatrix,
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

      /*I need to move dirMatrix and insDir to the first
      ` first base in the next row (skip indel column).
      ` The score for the indel column is updated in the
      ` next subsection. I need to find the score for a
      ` deletion before the update.
      */
      #if !defined BYTEMATRIX
         twoBitMvToNextElm(dirMatrix);
         changeTwoBitElm(dirMatrix, defMvIns);
         twoBitMvToNextElm(dirMatrix);

         twoBitMvToNextElm(&insDir);
      #else
         ++dirMatrix;
         *dirMatrix = defMvIns;
         ++dirMatrix;
         ++insDir;
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
      *scoreOnLP += settings->gapExtendI;
      delScoreL = *scoreOnLP + settings->gapExtendI;
      lastLastDirC = 0;
      ++scoreOnLP; /*Move to the first base pair*/


     /*At this point insDir is on the first base*/
     #if !defined BYTEMATRIX && !defined NOGAPOPEN
        indelScore(
           insScoreL,
           getTwoBitElm(&insDir),
           *scoreOnLP,
           settings
        );

        twoBitMvToNextElm(&insDir);
     #elif !defined NOGAPOPEN
        indelScore(
           insScoreL,
           *insDir,
           *scoreOnLP,
           settings
        );

        ++insDir;
     #elif !defined BYTEMATRIX
        insScoreL = *scoreOnLP + settings->gapExtendI;
        twoBitMvToNextElm(&insDir);
     #else
        insScoreL = *scoreOnLP + settings->gapExtendI;
        ++insDir;
     #endif
     /*At this piont insDir is on the second base*/
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
   #if !defined BYTEMATRIX
      twoBitMvBackOneElm(dirMatrix);
      changeTwoBitElm(dirMatrix, defMvStop);
      twoBitMvBackOneElm(dirMatrix);
   #else
      *(dirMatrix - 1) = defMvStop;
   #endif

   free(scoreRowLP);
   free(qryStartUL);
   free(refStartUL);

   scoreRowLP = 0;
   qryStartUL = 0;
   refStartUL = 0;

   /******************************************************\
   * Fun-02 Sec-05 Sub-01:
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
} /*WatermanAltAln*/

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints the path of the input index in dirST
\*-------------------------------------------------------*/
void printMatrixCig(
  FILE *outFILE,            // File to print to
  struct twoBitAry *dirST,  // Has index to print
  unsigned long lenRefUL,   // Length of the reference
  long scoreL               // Score of the path
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printMatrixCig
   '  - Prints out a cigar for an single path in a
   '    direction matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   uint8_t lastDirUC = defMvStop;
   char lastDirCigC = 0;
   uint32_t numDirUI = 0;
   uint8_t curDirUC = 0;

   struct twoBitAry matrixST;

   // As index 1
   unsigned long qryEndUL =
     (twoBitGetIndex(dirST) / (lenRefUL + 1));

   unsigned long refEndUL =
     (twoBitGetIndex(dirST) % (lenRefUL + 1));

   cpTwoBitPos(dirST, &matrixST);
   curDirUC = getTwoBitElm(&matrixST);

   // Print out the score, position of ending query base,
   // & position of ending reference base
   fprintf(
     outFILE,
     "%ld\t%lu\t%lu\t",
     scoreL,
     qryEndUL,
     refEndUL
   );
   
   goto initializePrint;
   
   while(curDirUC != defMvStop)
   { // While I have a cigar to build
     if(curDirUC != lastDirUC && lastDirUC != defMvStop)
     { // If I need to print out the last direction
       if(numDirUI > 1)
         fprintf(outFILE, "%u%c", numDirUI, lastDirCigC);

       else fprintf(outFILE, "%c", lastDirCigC);

       initializePrint:

       lastDirUC = defMvStop;
       numDirUI = 0;
       lastDirUC = curDirUC;

       switch(curDirUC)
       { // Switch: find the cigar symbol
         case defMvIns:
           lastDirCigC = 'I';
           break;

         case defMvSnp:
           lastDirCigC = 'X';
           break;

         case defMvDel:
           lastDirCigC = 'D';
           break;
       } // Switch: find the cigar symbol
     } // If I need to print out the last direction

     switch(curDirUC)
     { // Switch: Check which way to move
       case defMvIns:
         twoBitMvBackXElm(&matrixST, lenRefUL + 1);
         break;

       case defMvSnp:
         twoBitMvBackXElm(&matrixST, lenRefUL + 2);
         break;

       case defMvDel:
         twoBitMvBackOneElm(&matrixST);
         break;
     } // Switch: Check which way to move

     ++numDirUI;
     curDirUC = getTwoBitElm(&matrixST);
   } // While I have a cigar to build

   // Print out the last score
   if(numDirUI > 1)
     fprintf(outFILE, "%u%c", numDirUI, lastDirCigC);

   else fprintf(outFILE, "%c", lastDirCigC);

   qryEndUL = (twoBitGetIndex(&matrixST) /(lenRefUL+1));
   refEndUL = (twoBitGetIndex(&matrixST) %(lenRefUL+1));

   // Print the starting position of query and reference
   fprintf(outFILE, "\t%lu\t%lu\n", qryEndUL, refEndUL);

   return;
} // printMatrixCig

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints out the score and the position of the first
~      refence base, first query base, last reference base
|      and last query base to outFILE.
\--------------------------------------------------------*/
void printAltWaterAlns(
  struct alnMatrixStruct *alnMtxST,
     /*Has alternative alignments*/
  long minScoreL, /*Min score to keep an alternative*/
  FILE *outFILE   /*File to write alignments to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: printAltWaterAlns
   '  - Prints out all saved alternatives alignments
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Open the output file
   '  o fun-04 sec-03:
   '    o Print out the reference alignments
   '  o fun-04 sec-04:
   '    - Print out the query aligments
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-04 Sec-01:
  ^  - Variable declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  struct scoresStruct *scoreST = 0;
  char *startStr = "Alt:";
     /*To distingush alterantives from primary*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-04 Sec-03:
  ^  - Print out reference alternative alignments positions
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*Print out the alternative alignment header*/
  fprintf(outFILE, "%s score fristRefBase ", startStr);
  fprintf(outFILE, " firstQueryBase lastRefBase");
  fprintf(outFILE, " lastQueryBase\n");
  scoreST = alnMtxST->refBasesST;

  for(
    unsigned long ulRefBase = 0;
    ulRefBase < alnMtxST->lenRefScoresUL;
    ++ulRefBase
  ){ /*For all reference bases in the alignment*/
    /*Check if even need to print out any reference alns*/
    if(scoreST->scoreL < minScoreL)
      goto nextRefScore;

    /*Prints out the score, first reference and query base
    ` in the alignment, and the last reference and query
    ` base in the alignment
    */
    fprintf(
       outFILE,
       "%s %li %lu %lu %lu %lu\n",
       startStr,
       scoreST->scoreL,
       scoreST->refStartUL,
       scoreST->qryStartUL,
       scoreST->refEndUL,
       scoreST->qryEndUL
    );

     nextRefScore:
    ++scoreST;  /*Move to the next entry*/
  } /*For all reference bases in the alignment*/
     
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-04:
  ^  - Print out the query alignments
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  scoreST = alnMtxST->qryBasesST;

  for(
    unsigned long ulQryBase = 0;
    ulQryBase < alnMtxST->lenQueryScoresUL;
    ++ulQryBase
  ){ /*For all reference bases in the alignment*/

    /*Check if even need to print out any reference alns*/
    if(scoreST->scoreL < minScoreL)
      goto nextQueryScore;

    /*Prints out the score, first reference and query base
    ` in the alignment, and the last reference and query
    ` base in the alignment
    */
    fprintf(
       outFILE,
       "%s %li %lu %lu %lu %lu\n",
       startStr,
       scoreST->scoreL,
       scoreST->refStartUL,
       scoreST->qryStartUL,
       scoreST->refEndUL,
       scoreST->qryEndUL
    );

    nextQueryScore:
    ++scoreST;  /*Move to the next entry*/
  } /*For all reference bases in the alignment*/

  return;
} /*printAltWaterAlns*/
