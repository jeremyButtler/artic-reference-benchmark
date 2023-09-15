/*#########################################################
# Name alignmentsFun
# Use:
#  o Holds functions for doing a pairwise Needleman Wunsch
#    alignment
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
#   o <stdio.h>
#########################################################*/

#include "needleman.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-01 NeedleManWunschAln:
'    o Perform a Needleman-Wunsch alignment on the two
'      input sequences
'  - fun-02 updateDirAndScore:
'    o Picks best score and direction for the current base
'      pairs being compared in a Needleman Wunsch alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnMatrixStruct with the directoin matrix and best
|      score struct pointing to the connor left cell in
|      the direction matrix
\--------------------------------------------------------*/
struct alnMatrixStruct * NeedlemanAln(
    struct seqStruct *qryST, // query sequence and data
    struct seqStruct *refST,  // ref sequence and data
      /* both qryST and refST have the sequence,
      `  they also have the point to start the alignment
      `  seqST->offsetUL (index 0) and the point to end
      `  the alignment seqST->endAlnUL (index 0).
      */
    struct alnSet *settings // Settings for the alignment
    // *startI and *endI paramaters should be index 1
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: NeedlemanAln
   '  - Perform a Needleman-Wunsch alignment on the two
   '    input sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in the initial negatives for the reference
   '  o fun-01 sec-04:
   '    - Fill the matrix with scores
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Get start & end of the query and reference sequences*/
   char *refStartStr = refST->seqCStr + refST->offsetUL;
   char *refEndStr = refST->seqCStr + refST->endAlnUL;

   char *qryStartStr = qryST->seqCStr + qryST->offsetUL;
   char *qryEndStr = qryST->seqCStr + qryST->endAlnUL;

   char *qryIterStr = 0;
   char *refIterStr = 0;

   /*Find the length of the reference and query. The +1
   ` is to account for offsetUL being index 0
   */
   unsigned long lenQryUL =
       qryST->endAlnUL - qryST->offsetUL + 1;
   unsigned long lenRefUL =
       refST->endAlnUL - refST->offsetUL + 1;

   /*Scoring variables*/
   long insScoreL = 0;   /*Score for doing an insertion*/
   long snpScoreL = 0;   /*Score for doing an match/snp*/
   long delScoreL = 0;   /*Score for doing an deletion*/
   long nextSnpSL = 0;   /*Score for the next match/snp*/

   long *scoreRowLP = 0; /*matrix to use in alignment*/
   long *scoreOnLP = 0;  /*Score I am working on*/

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

   /*Structure to return*/
   struct alnMatrixStruct *retMtxST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retMtxST = calloc(1, sizeof(struct alnMatrixStruct));
   if(retMtxST == 0) return 0;

   initAlnMatrixST(retMtxST);
   scoreRowLP = calloc(lenRefUL + 1, sizeof(long));
      /*I need two rows to keep track of the scores (2x)
      ` - lenRefUL + 1 is to account for insertion column
      */

   if(scoreRowLP == 0)
   { /*If I had a memory error*/
     freeAlnMatrixST(retMtxST);
     return 0;
   } /*If I had a memory error*/

   #if !defined BYTEMATRIX
      dirMatrix=
         makeTwoBit((lenRefUL+1) * (lenQryUL+1) + 2,0);
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
   ^  - Fill in the initial negatives for the reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Build up the indels for the reference row*/
   scoreOnLP = scoreRowLP;
   *scoreOnLP = 0;         /*Top left cell starts at 0*/
   ++scoreOnLP;            /*Move to next cell (score)*/

   /*The insertion/match/snp score is on the previous row*/
   #if !defined BYTEMATRIX
      #if !defined NOGAPOPEN
         cpTwoBitPos(dirMatrix, &insDir);
      #endif

      changeTwoBitElm(dirMatrix, defMvStop);
      twoBitMvToNextElm(dirMatrix);
      changeTwoBitElm(dirMatrix, defMvDel);
      twoBitMvToNextElm(dirMatrix);
   #else
      #if !defined NOGAPOPEN
         insDir = dirMatrix;
      #endif

      *dirMatrix = defMvStop;
      ++dirMatrix;
      *dirMatrix = defMvDel;
      ++dirMatrix;
   #endif

   /*2nd score (first indel in matrix)*/
   #ifdef NOGAPOPEN
      *scoreOnLP = settings->gapExtendI;
   #else
      *scoreOnLP = settings->gapOpenI;
   #endif

   /*Set up scores for remaning cells in the first row*/
   ++scoreOnLP;      /*Move to the 2nd reference base*/
   refIterStr = refStartStr;

   while(refIterStr < refEndStr)
   { /*loop; till have initalized the first row*/
     *scoreOnLP = *(scoreOnLP - 1) + settings->gapExtendI;

     /*Move to the next cell (ref base)*/
     #if !defined BYTEMATRIX
        changeTwoBitElm(dirMatrix, defMvDel);
        twoBitMvToNextElm(dirMatrix);
     #else
        *dirMatrix = defMvDel;
        ++dirMatrix;
     #endif

     ++scoreOnLP;
     ++refIterStr;
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-01 sec-04 sub-01:
   ^    - Set up for filling the rest of the matrix
   ^  o fun-01 sec-04 sub-02:
   ^    - Fill in the first cell (indel column)
   ^  o fun-01 sec-04 sub-03:
   ^    - Get scores for insertion, deletion, match
   ^  o fun-01 sec-04 sub-04:
   ^    - Move to the next refernce/query base
   ^  o fun-01 sec-04 sub-05:
   ^    - Handle the first cell (indel col) in new row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   ^ Fun-01 Sec-04 Sub-01:
   ^  - Set up for filling the rest of the matrix
   \******************************************************/

   qryIterStr = qryStartStr;
   refIterStr = refStartStr;

   /*Add the insertion to the first score
   ` The first score on the second row is in the insertion
   ` column and is always the start of a new insertion.
   */
   scoreOnLP = scoreRowLP;

   nextSnpSL =
        getBaseScore(qryIterStr,refIterStr,settings)
      + *scoreOnLP;

   /*Fill in the current indel column for this row*/
   #ifndef NOGAPOPEN
      *scoreOnLP = settings->gapOpenI;
      /**scoreOnLP = settings->gapExtendI;
      ` This was the old version. However, this was a
      ` mistake and can result in slightly worse alignments
      ` (when alignments have many gaps).
      */
   #else
      *scoreOnLP = settings->gapExtendI;
   #endif

   delScoreL = *scoreOnLP + settings->gapExtendI;
   ++scoreOnLP;
   insScoreL = *scoreOnLP + settings->gapExtendI;

   #if !defined BYTEMATRIX && !defined NOGAPOPEN
      changeTwoBitElm(dirMatrix, defMvIns);
      twoBitMvToNextElm(dirMatrix);
      twoBitMvToNextElm(&insDir);
      twoBitMvToNextElm(&insDir);
   #elif !defined BYTEMATRIX
      changeTwoBitElm(dirMatrix, defMvIns);
      twoBitMvToNextElm(dirMatrix);
   #elif !defined NOGAPOPEN
     *dirMatrix = defMvIns;
     ++dirMatrix;
     insDir += 2;
   #else
     *dirMatrix = defMvIns;
     ++dirMatrix;
   #endif

   /******************************************************\
   * Fun-01 Sec-04 Sub-02:
   *  - Fill in the first cell (indel column)
   \******************************************************/

   /******************************************************\
   * Fun-01 Sec-04 Sub-03:
   *  - Get scores for insertion, deletion, match
   \******************************************************/

   /*Starting on the first sequence row*/
   while(qryIterStr <= qryEndStr)
   { /*loop; fill the direction matrix with socres*/

       /*I have already found the first bases match*/
       refIterStr = refStartStr + 1;

       /*Find scores/directions for all basepairs in row*/
       /*I am doing one of the end to avoid caculating
       ` extra scores for each base
       */
       while(refIterStr <= refEndStr)
       { /*loop; compare one query to all reference bases*/
           /*Find the score for a match/snp*/
           snpScoreL = nextSnpSL;

           /*Find the next score for an snp/match*/
           nextSnpSL = 
              getBaseScore(
                 qryIterStr,
                 refIterStr,
                 settings
           );
           nextSnpSL += *scoreOnLP;

           /*Find the best score*/
           #if !defined BYTEMATRIX
              twoBitMaxScore(
                dirMatrix,
                settings,
                &insScoreL,
                &snpScoreL,
                &delScoreL,
                scoreOnLP
              ); /*Update the score and direction*/
           #else
              charMaxScore(
                dirMatrix,
                settings,
                &insScoreL,
                &snpScoreL,
                &delScoreL,
                scoreOnLP
              ); /*Update the score and direction*/
           #endif

           /*Find the next deletion score*/
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
              delScoreL = *scoreOnLP+settings->gapExtendI;
           #endif

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
              insScoreL = *scoreOnLP +settings->gapExtendI;
           #endif

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

           /**********************************************\
           * Fun-01 Sec-04 Sub-04:
           *  - Move to the next refernce/query base
           \**********************************************/
       
           /*Move to the next reference base*/
           ++refIterStr;
       } /*loop; compare one query to all reference bases*/

       /**************************************************\
       * Fun-01 Sec-04 Sub-05:
       *  - Find the best score for the last base
       \**************************************************/

       /*Find the best score for the last base. In this
       ` case, this score can only apply to indels. So,
       ` I need to move off it to avoid overwirting it
       */
       #if !defined BYTEMATRIX
          twoBitMaxScore(
            dirMatrix,
            settings,
            &insScoreL,
            &nextSnpSL,
            &delScoreL,
            scoreOnLP
          ); /*Update the score and direction*/
       #else
          charMaxScore(
            dirMatrix,
            settings,
            &insScoreL,
            &nextSnpSL,
            &delScoreL,
            scoreOnLP
          ); /*Update the score and direction*/
       #endif

       /**************************************************\
       *  Fun-01 Sec-04 Sub-05:
       *   - Update the indel columns directoin
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
       * Fun-01 Sec-04 Sub-06:
       *  - Find the first bases scores
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
      #elif !defined BYTEMATRIX
         insScoreL = *scoreOnLP + settings->gapExtendI;
      #endif
      /*At this piont insDir is on the second base*/
   } /*loop; fill the direction matrix with socres*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Account for being two cells out of bounds*/
   #if !defined BYTEMATRIX
      twoBitMvBackOneElm(dirMatrix);
      changeTwoBitElm(dirMatrix, defMvStop);
      twoBitMvBackOneElm(dirMatrix);
   #else
      *(dirMatrix - 1) = defMvStop;
   #endif

   /*Set the best score to the cornor right cell*/
   retMtxST->bestScoreST.refEndUL =
      refIterStr - refST->seqCStr - 1;

   retMtxST->bestScoreST.qryEndUL =
      qryIterStr - qryST->seqCStr - 1;

   retMtxST->bestScoreST.scoreL = *(scoreRowLP + lenRefUL);

   /*Clean UP*/
   free(scoreRowLP);
   return retMtxST;
} /*NeeldeManWunschAln*/
