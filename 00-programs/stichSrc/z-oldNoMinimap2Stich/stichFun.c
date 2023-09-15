/*#########################################################
# Name: stichFun
# Use:
#  - Holds functions used in stiching amplicons into a
#    consensus.
# Libraries:
#  - "alnSeqSrc/memWater.h"
#  - "alnSeqSrc/hirschberg.h"
#  - "sitchAmpStruct.h"            (No .c file)
#  - "samFunSrc/trimSam.h"
#  o "samFunSrc/samEntryStruct.h"  (No .c file)
#  o "samFunSrc/cStrToNumberFun.h" (No .c file)
#    o Also in alnSeqSrc
#  o "alnSeqSrc/generalAlnFun.h"
#  o "alnSeqSrc/alnStruct.h"
#  o "alnSeqSrc/alnMatrixStruct.h"
#  o "alnSeqSrc/twoBitArrays.h"    (No .c file)
#  o "alnSeqSrc/scoresST.h"        (No .c file)
#  o "alnSeqSrc/seqStruct.h"
#  o "alnSeqSrc/alnSetStruct.h"
#  o "alnSeqSrc/alnSeqDefaults.h"
#  o "alnSeqSrc/dataTypeShortHand.h"         (No .c file)
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdint.h>
#  o <stdio.h>
#  o <string.h>
#########################################################*/

#include "stichFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' stichFun SOF: Start Of Functions
'  - Functions to stich amplicons into a consensus sequence
'  o fun-01: getAmpPos 
'    - Find the starting positions and ending postions of
'      each amplicon on the reference
'  o fun-02: stichAmpCon 
'    - Uses a reference to stiche togther amplicons into a
'      consensuses
'  o fun-03: alnAmpToRef 
'    - Aligns an amplicon sequence to a reference sequence.
'  o fun-04: hirschbergToSeq 
'    - This is used in alnAmpToRef. It is used to convert
'      the output from HirschbergFun to an query sequence
'      alignment.
'  o fun-05: stichAmps 
'    - Stiches amplions into a consensus
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: getAmpPos (Fun-01:)
| Use:
|   - Find the starting positions and ending postions of
|     each amplicon on the reference
| Input:
|  - refST:
|    o seqStruct with reference sequence
|  - qryST:
|    o seqStruct with query sequence
|  - ampFaFILE:
|    o Fasta file handle with amplicons to find positions
|  - numAmpsUL:
|    o This will hold the number of amplicons in ampFaFILE
|  - indexAryUL:
|    o This will hold the index for every amplicon sequence
|      in the file.
|  - settings:
|    o alnSet struct with the settings to use for aligment
| Output:
|   - Modifies:
|     o numAmpsUL to hold the number of amplicons
|     o indexAryUL to hold the index of every sequence.
|   - Returns:
|     o Pointer to array of scoresStructs with the
|       coordinates of each alignment
|     o 0 for memory errors
\--------------------------------------------------------*/
struct scoresStruct * getAmpPos(
   struct seqStruct *refST,   /*Reference sequence*/
   struct seqStruct *qryST,  /*Blank struct to work with*/
   FILE *ampFaFILE,          /*File with amplicons*/
   ulong *numAmpsUL, /*Will have number amplicons*/
   ulong **indexAryUL, /*Will hold file index for scores*/
   struct alnSet *settings
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: getAmpPos
   '  - Gets the position of each amplicon on the reference
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Memory allocation an initialize variables
   '  o fun-01 Sec-03:
   '     - Align each sequence to the reference
   '  o fun-01 sec-04:
   '    - Reset file pointer and sort the sequences by
   '      reference starting positions
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ulong curPosUL = ftell(ampFaFILE);
   ulong *indexIterUL = 0;
   struct scoresStruct *ampsST = 0;
   struct scoresStruct *tmpScoreST = 0;
   struct scoresStruct *ampIterST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Memory allocation an initialize variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *numAmpsUL = 0;

   /*Find the number of amplicons*/
   while(readFaSeq(ampFaFILE, qryST) & 1)
      ++(*numAmpsUL);

   /*Move back to start*/
   fseek(ampFaFILE, curPosUL, SEEK_SET);

   *indexAryUL = malloc(sizeof(ulong) * *numAmpsUL);
   if(*indexAryUL == 0) return 0;

   ampsST=malloc(sizeof(struct scoresStruct) * *numAmpsUL);

   if(ampsST == 0)
   { /*If: I had a memory error*/
      free(*indexAryUL);
      *indexAryUL = 0;
      return 0; /*Memory error*/
   } /*If: I had a memory error*/

   ampIterST = ampsST;
   indexIterUL = *indexAryUL;
   *indexIterUL = ftell(ampFaFILE);
   ++indexIterUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Align each sequence to the reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(readFaSeq(ampFaFILE, qryST) & 1)
   { /*Loop: Find the starting position for each amplicon*/
      *indexIterUL = ftell(ampFaFILE);
      qryST->endAlnUL = qryST->lenSeqUL;
      tmpScoreST = memWaterAln(qryST, refST, settings);

      if(tmpScoreST == 0)
      { /*If: I had a memory error*/
         free(*indexAryUL);
         freeScoresSTAry(ampsST, *numAmpsUL, 0);

         ampsST = 0;
         *indexAryUL = 0;

         return 0;
      } /*If: I had a memory error*/

      cpScoreST(ampIterST, tmpScoreST);
      ++ampIterST;
      freeScoresST(tmpScoreST, 0); /*No longer need*/
      ++indexIterUL;
   } /*Loop: Find the starting position for each amplicon*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Reset file pointer and sort the sequences by
   ^    reference starting positions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Make sure the last index (not sorted) is at the end
   ` of the file
   */
   fseek(ampFaFILE, curPosUL, SEEK_END);
   *indexIterUL = ftell(ampFaFILE);
   fseek(ampFaFILE, curPosUL, SEEK_SET);

   sortScoresStartLenIndex(
      &ampsST,
      0,
      *numAmpsUL - 1,
      *indexAryUL
   );

   return ampsST;
} /*getAmpPos*/

/*--------------------------------------------------------\
| Name: stichAmpCon (Fun-02:)
| Use:
|   - Uses a reference to stiche togther amplicons into a
|     consensuses
| Input:
|   - refST:
|     o seqStruct with the reference sequence and the
|       offset to start the alignment (offsetUL).
|   - ampST:
|     o seqStruct with the amplicon sequence.
|     o Primers should be removed.
|   - conAlnST:
|     o 0 or alnStruct having the alignment for the current
|       consensuses.
|     o This funciton will build this up, so start with 0.
|   - conHasPrioBl:
|     o This is for when bases overlap between the amplicon
|       and the consensus.
|     o 1: Always keep the consensuses bases
|     o 0: Always keep the amplicons bases
|   - settings:
|     o Settings for the alignment
| Output:
|  - Returns:
|    o Pointer to conAlnST or a alnStrcut if conAlnST is 0
|    o 0 for error
|  - Modifies:
|    o conSeqStr to hold the sitched together consensus
|    o conAlnST to have the alignment updated to inlucde
|      the amplicon
|    o refST->offsetUL to be 0.
|    o refST->endAlnUL to be refST->lenSeqUL.
|    o ampST->offsetUL to be 0.
|    o ampST->endAlnUL to be ampST->lenSeqUL.
\--------------------------------------------------------*/
char * stichAmpCon(
   FILE *ampFaFILE,            /*File with sequences*/
   struct scoresStruct *ampsAryST, /*Array of scores*/
   ulong *seqIndexAryUL,       /*Index for every seq*/
   struct seqStruct *refST,    /*Reference sequence*/
   struct seqStruct *ampST,    /*Amplicon sequence*/
   struct alnSet *alnSetST,     /*Settings for alnSeq*/ 
   struct stichSet *stichSetST  /*Settings for stich*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: stichAmpCon
   '  - Stiches together amplicon consensuses to to make
   '    an consensus genome
   '  o fun-02 sec-01:
   '     - Variable declerations
   '  o fun-02 sec-02:
   '    - Align each amplicon sequence
   '  o fun-02 sec-03:
   '    - Stich together each amplicon sequence
   '  o fun-02 sec-04:
   '    - Deal with disagreements in the amplicons
   '  o fun-02 sec-05:
   '    - Clean up and exit
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *alnSeqStr = 0;
   char *conSeqStr = 0;
   char *conIterStr = 0;

   ulong refEndUL = 0;
   ulong depthUL = 0;
   ulong insSupUL = 0;
   ulong delSupUL = 0;
   ulong snpSupUL = 0;
   ulong supportUL = 0;

   ulong *indexIterUL = seqIndexAryUL;
   ulong curPosUL = ftell(ampFaFILE);

   struct scoresStruct *ampScoreST = ampsAryST;
   struct stichAmpST *conST = 0;
   struct stichAmpST *stichIterST = 0;
   struct stichAmpST *altStichST = 0;
   struct stichAmpST *lastConBaseST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-02:
   ^  - Align each amplicon sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fseek(ampFaFILE, *indexIterUL, SEEK_SET);

   while(readFaSeq(ampFaFILE, ampST) & 1)
   { /*While I have sequences to stich together*/
      alnSeqStr =
         alnAmpToRef(
            refST,
            ampST,
            ampScoreST,
            alnSetST
         );

      if(alnSeqStr == 0)
      { /*If: I had a memory error*/
         freeStichAmpSTList(&conST);
         return 0;
      } /*If: I had a memory error*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun-02 Sec-03:
      ^  - Stich together each amplicon sequence
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      conST = 
         stichAmps(
            alnSeqStr,
            ampScoreST,
            conST,
            &refEndUL,    /*Last ref base in consensus*/
            stichSetST->maskC
         );

      if(conST == 0)
      { /*If: I had a memory error*/
         free(alnSeqStr);
         alnSeqStr = 0;
         /*conST has already been freeded*/
         return 0;
      } /*If: I had a memory error*/

      ++ampScoreST;
      ++indexIterUL;
      fseek(ampFaFILE, *indexIterUL, SEEK_SET);
      free(alnSeqStr);
      alnSeqStr = 0;
   } /*While I have sequences to stich together*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-04:
   ^  - Deal with disagreements in the amplicons
   *  o fun-02 sec-04 sub-01:
   *    - Set up for colapsing the consensus list
   *  o fun-02 sec-04 sub-02:
   *    - Find the amount of support for each error type
   *  o fun-02 sec-04 sub-03:
   *    - Handle cases were have values beneath users
   *      min depth (needs 100% agreement)
   *  o fun-02 sec-04 sub-04:
   *    - Handle cases with values at or above the
   *      users min read depth
   *  o fun-02 sec-04 sub-05:
   *    - I only have one base, go with that
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-02 Sec-04 Sub-01:
    *  - Set up for colapsing the consensus list
    \*****************************************************/

   --ampScoreST; /*Get back on the last score*/
   while(conST->lastBase != 0) conST = conST->lastBase;

   stichIterST = conST;

   conSeqStr =
      malloc(sizeof(char) * (ampScoreST->refEndUL << 1));

   conIterStr = conSeqStr;

   if(conSeqStr == 0)
   { /*If I had a memory error*/
      free(alnSeqStr);
      alnSeqStr = 0;
      freeStichAmpSTList(&conST);
      return 0;
   } /*If I had a memory error*/

   stichNextBase:

   while(stichIterST != 0)
   { /*While I have bases to convert to chars*/

      /***************************************************\
      * Fun-02 Sec-04 Sub-02:
      *  - Find the amount of support for each error type
      \***************************************************/

      lastConBaseST = stichIterST;
      altStichST = stichIterST;
      stichIterST = stichIterST->nextBase;
      depthUL = altStichST->depthUL; /*# amplicons*/

      delSupUL = 0;
      insSupUL = 0;
      snpSupUL = 0;

      if(altStichST)
      { /*If: I have alternative bases*/
         while(altStichST != 0)
         { /*Loop: Find the total support for position*/
            switch(altStichST->errC)
            { /*Switch: Check the error type*/
               case defMvStop: break; /*Never fires*/

               case defMvDel:
                  delSupUL += altStichST->supportUL;
                  break;

               case defMvIns:
                  insSupUL += altStichST->supportUL;
                  break;

               case defMvSnp:
                  snpSupUL += altStichST->supportUL;
                  break;
 
               /*This only happens for missing regios*/
               case defMvMask:
                  snpSupUL += altStichST->supportUL;
                  break;
            } /*Switch: Check the error type*/

            altStichST = altStichST->altBase;
         } /*Loop: Find the total support for position*/

         /************************************************\
         * Fun-02 Sec-04 Sub-03:
         *  - Handle cases were have values beneath users
         *    min depth (needs 100% agreement)
         *  o fun-02 sec-04 sub-03 cat-01:
         *    - Handle insertion cases for low dephts
         *  o fun-02 sec-04 sub-03 cat-02:
         *    - Handle deltion cases for low dephts
         *  o fun-02 sec-04 sub-03 cat-03:
         *    - Handle snp/match or snp/match with
         *      deletions cases for low depths
         \************************************************/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-02 Sec-04 Sub-03 Cat-01:
         +  - Handle insertion cases for low depths
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         altStichST = lastConBaseST;
 
         if(depthUL < stichSetST->minDepthUL)
         { /*If I am beneath the min depth limits*/
            if(altStichST->errC == defMvIns)
            { /*If: I had an insertion*/
               if(insSupUL < depthUL)
               { /*If: ins was not 100% supported*/
                  while(altStichST != 0)
                     freeStichAmpST(&altStichST);
                  goto stichNextBase;   
               } /*If: ins was not 100% supported*/

               else if(insSupUL > 0)
               { /*Else if: have 100% support for ins*/
                  /*Check if there are alterantive bases*/
                  if(altStichST->supportUL < depthUL)
                     *conIterStr = stichSetST->maskC | 32;
                  else *conIterStr=altStichST->baseC | 32;

                  ++conIterStr;

                  while(altStichST != 0)
                     freeStichAmpST(&altStichST);
                  goto stichNextBase;   
               } /*Else if: have 100% support for ins*/
            } /*If: I had an insertion*/


            /*++++++++++++++++++++++++++++++++++++++++++++\
            + Fun-02 Sec-04 Sub-03 Cat-02:
            +  - Handle deletion cases for low dephts
            \++++++++++++++++++++++++++++++++++++++++++++*/

            if(delSupUL > snpSupUL && delSupUL >= depthUL)
            { /*Else if: A Deletion is 100% supported*/ 
                while(altStichST != 0)
                  freeStichAmpST(&altStichST);
               goto stichNextBase;   
            } /*Else if: A Deletion is 100% supported*/ 

            /*++++++++++++++++++++++++++++++++++++++++++++\
            + Fun-02 Sec-04 Sub-03 Cat-03:
            +  - Handle snp/match or snp/match with
            +    deletions cases for low depths
            \++++++++++++++++++++++++++++++++++++++++++++*/

            else
            { /*Else: there is an snp or snp/del combo*/
               while(
                     altStichST->errC != defMvSnp
                  && altStichST->errC != defMvMask
               ){freeStichAmpST(&altStichST);}

               depthUL -= delSupUL;

               if(altStichST->supportUL >= depthUL)
                  *conIterStr = altStichST->baseC & (~32);
               else *conIterStr= stichSetST->maskC & (~32);

               ++conIterStr;

               while(altStichST != 0)
                  freeStichAmpST(&altStichST);
               goto stichNextBase;   
            } /*Else: there is an snp or snp/del combo*/
         } /*If I am beneath the min depth limits*/

         /************************************************\
         * Fun-02 Sec-04 Sub-04:
         *  - Handle cases with values at or above the
         *    users min read depth
         *  o fun-02 sec-04 sub-04 cat-01:
         *    - Handle insertion has to low support cases
         *  o fun-02 sec-04 sub-04 cat-02:
         *    - Handle insertion has enough support cases
         *  o fun-02 sec-04 sub-04 cat-03:
         *    - Handle deletion has enough support cases
         *  o fun-02 sec-04 sub-04 cat-04:
         *    - Handle snp/match cases
         \************************************************/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-02 Sec-04 Sub-04 Cat-01:
         +  - Handle insertion has to low support cases
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         supportUL = (insSupUL * 100) / depthUL;

         if(
               altStichST->errC == defMvIns
            && supportUL < stichSetST->minSupportUL
           )
         { /*If this insertion is not supported*/
            while(altStichST != 0)
               freeStichAmpST(&altStichST);

            goto stichNextBase;   
         } /*If this insertion is not supported*/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-02 Sec-04 Sub-04 Cat-02:
         +  - Handle insertion has enough support cases
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         else if(altStichST->errC == defMvIns)
         { /*Else if: This insertion is supported*/
            while(altStichST != 0)
            { /*Loop: Find the best position*/
               supportUL =
                    (altStichST->supportUL * 100)
                  / depthUL;

               if(supportUL >= stichSetST->minSupportUL)
               { /*If: alt base is the supported base*/
                  *conIterStr = altStichST->baseC | 32;
                  ++conIterStr;

                  while(altStichST != 0)
                     freeStichAmpST(&altStichST);

                  goto stichNextBase;   
               } /*If: alt base is the supported base*/

               freeStichAmpST(&altStichST);
            } /*Loop: Find the best position*/

            /*No base had enough support to keep a base*/
            *conIterStr = stichSetST->maskC | 32;
            ++conIterStr;
            goto stichNextBase;   
         } /*Else if: This insertion is supported*/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-02 Sec-04 Sub-04 Cat-03:
         +  - Handle deletion has enough support cases
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         supportUL = (delSupUL * 100) / depthUL;
         if(
                delSupUL > snpSupUL
             && supportUL >= stichSetST->minSupportUL
         ){ /*If this base is a deletion*/
            while(altStichST != 0)
               freeStichAmpST(&altStichST);

            goto stichNextBase;   
         } /*If this base is a deletion*/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-02 Sec-04 Sub-04 Cat-04:
         +  - Handle snp/match cases
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         while(altStichST != 0)
         { /*See if there is good support for position*/
            supportUL = 
              (altStichST->supportUL * 100)/depthUL;

            if(
                  altStichST->errC != defMvDel 
               && supportUL >= stichSetST->minSupportUL
            ){ /*If: alt base is the supported base*/
               *conIterStr = altStichST->baseC;
               ++conIterStr;

               while(altStichST != 0)
                  freeStichAmpST(&altStichST);

               goto stichNextBase;   
            } /*If: alt base is the supported base*/

            freeStichAmpST(&altStichST);
         } /*See if there is good support for position*/

         /*Else not support for base*/
         *conIterStr = stichSetST->maskC & (~32);
         ++conIterStr;
         goto stichNextBase;   
      } /*If: I have alternative bases*/

      /***************************************************\
      * Fun-02 Sec-04 Sub-05:
      *  - I only have one base, go with that
      \***************************************************/

      else
      { /*Else this is the only base*/
         if(altStichST->errC != defMvDel)
         { /*If this base is not a deletion*/
            *conIterStr = altStichST->baseC & (~32);
            ++conIterStr;
         } /*If this base is not a deletion*/

         freeStichAmpST(&altStichST);
      } /*Else this is the only base*/
   } /*While I have bases to convert to chars*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-05:
   ^  - Clean up and exit
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fseek(ampFaFILE, curPosUL, SEEK_SET);
   return conSeqStr;
} /*stichAmpCon*/

/*--------------------------------------------------------\
| Name: alnAmpToRef (Fun-03:)
| Use:
|  - Aligns an amplicon sequence to a reference sequence.
| Input:
|  - refST:
|    o Has the reference sequence
|  - qryST:
|    o Has the query sequence
|  - ampScoreST:
|    o Has the starting position of the alignments
|  - alnSetST:
|    o Has the settings for the alignment.
| Output:
|  - Modifies:
|    o ampScoreST->refStratUL to hold the first reference
|      base in the alignment.
|    o ampScoreST->qryStratUL to hold the first query base
|      in the alignment.
|  - Returns:
|    o c-string with the aligned amplicon sequence
|    o 0 for memory errors
\--------------------------------------------------------*/
char * alnAmpToRef(
   struct seqStruct *refST,      /*Has reference to align*/
   struct seqStruct *ampST,      /*Has amplicon to align*/
   struct scoresStruct *ampScoreST, /*Score for alignment*/
   struct alnSet *alnSetST       /*Settings for alnSeq*/ 
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: alnAmpToRef:
   '  - Aligns the ampicon to the reference using the
   '    cordinates in ampScoreST.
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Memory allocation (set up for Hirschberg)
   '  o fun-03 sec-03:
   '    - Run the hirschberg alignment
   '  o fun-03 sec-04:
   '    - Clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-03 Sec-01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   char *alnStr = 0;
   unsigned long lenRefUL =
     ampScoreST->refEndUL - ampScoreST->refStartUL + 1;
     /*+1 to convert to index 1 (values are index 0)*/
   unsigned long lenQryUL =
     ampScoreST->qryEndUL - ampScoreST->qryStartUL + 1;
     /*+ 1 to convert to index 1 (values are index 0)*/

   long *forwardScoreRowL = 0;
   long *reverseScoreRowL = 0;

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
   ^  Fun-03 Sec-02:
   ^   - Memory allocation (set up for Hirschberg)
   ^   o fun-03 sec-02 sub-01:
   ^     - Initalize the ouput alignment structure 
   ^   o fun-03 sec-02 sub-02:
   ^     - Initalize the scoring rows
   ^   o fun-03 sec-02 sub-03:
   ^     - Initalize the direction rows
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-01:
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
   * Fun-03 Sec-02 Sub-02:
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
   ^  Fun-03 Sec-03:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Sening in offset values, because alignment array is
   ` sized to the alignmnet region
   */
   HirschbergFun(
     refST->seqCStr + ampScoreST->refStartUL,
     0,                /*1st reference base to align*/
     lenRefUL,         /*Length of ref region to align*/
     ampST->seqCStr + ampScoreST->qryStartUL,
     0,                /*1st query base to align*/
     lenQryUL,         /*length of query target region*/
     forwardScoreRowL, /*For scoring*/
     reverseScoreRowL, /*For scoring*/
     refAln,      /*Holds the reference alignment*/
     qryAln,      /*Holds the query alignment*/
     dirRow,      /*Direction row for thread safe scoring*/
     alnSetST     /*Settings for the alignment*/
   );
     /*dirRow becomes a dummy variable for -DNOGAPOPEN*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-03 Sec-04:
   ^    - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   free(forwardScoreRowL);
   free(reverseScoreRowL);

   /*Convert hirschberFun output to aligned query*/
   alnStr=hirschbergToSeq(ampST,ampScoreST,refAln,qryAln);
 
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

   return alnStr; /*Is either 0 or the aligned sequence*/
} /*alnAmpToRef*/

/*--------------------------------------------------------\
| Name: hirschbergToSeq (Fun-04:)
| Use:
|  - This is used in alnAmpToRef. It is used to convert
|    the output from HirschbergFun to an query sequence
|    alignment.
| Input:
|  - refST:
|    o Has the reference sequence
|  - qryST:
|    o Has the query sequence
|  - ampScoreST:
|    o Has the starting position of the alignments
|  - refAln:
|    o The reference alignment from HirschberFun
|  - qryAln:
|    o The query alignment from HirschberFun
| Output:
|  - Modifies:
|    o ampScoreST->refStratUL to hold the first reference
|      base in the alignment.
|    o ampScoreST->qryStratUL to hold the first query base
|      in the alignment.
|  - Returns
|    o A c-string with the aligned query sequence.
|      Insertions are in lower case, with deletions as '-'.
\--------------------------------------------------------*/
char * hirschbergToSeq(
  struct seqStruct *qryST,         /*Query sequence*/
  struct scoresStruct *ampScoreST, /*Alignment positions*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *refAln,
     struct twoBitAry *qryAln
  #else
     char *refAln, /*has reference alignment*/
     char *qryAln  /*has query alignment*/
  #endif
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: hirschbergToSeq
   '  - Converts the twobit or char arries used in the
   '    Hirschberg to an aligned query sequence, with
   '    insertions being lower case. This called in
   '    alnAmpToRef
   '  o fun-04 sec-01:
   '     - Variable declerations
   '  o fun-04 sec-02:
   '     - Allocate memory
   '  o fun-04 sec-03:
   '     - Find the first aligned reference base
   '  o fun-04 sec-04:
   '     - Find the first aligned query base
   '  o fun-04 sec-05:
   '     - Build the alignment
   '  o fun-04 sec-06:
   '     - Clean up and return alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *qryStr = qryST->seqCStr + ampScoreST->qryStartUL;
   char *alnStr = 0;
   char *alnIterStr = 0;

   uint8_t bitUC = 0;
   ulong lenQryUL =
     ampScoreST->qryEndUL - ampScoreST->qryStartUL + 1;
     /*+ 1 to convert to index 1 (values are index 0)*/

   long refFirstAlnBaseL = ampScoreST->refStartUL;
   long qryFirstAlnBaseL = ampScoreST->qryStartUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-02:
   ^  - Allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   alnStr = malloc(sizeof(char) * (lenQryUL << 1));
   if(alnStr == 0) return 0; /*Memory error*/
   alnIterStr = alnStr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-03:
   ^  - Find the first aligned reference base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #ifdef HIRSCHTWOBIT
     twoBitMvXElmFromStart(refAln, 0);
     bitUC = getTwoBitElm(refAln);
   #else
     bitUC = (uint8_t) *refAln;
   #endif

   while(bitUC == defGapFlag)
   { /*Loop: Find the first aligned reference base*/
      #ifdef HIRSCHTWOBIT
         twoBitMvToNextElm(refAln);
         bitUC = getTwoBitElm(refAln);
      #else
         ++refAln;
         bitUC = (uint8_t) *refAln;
      #endif

      ++refFirstAlnBaseL;
   } /*Loop: Find the first aligned reference base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-04:
   ^  - Find the first aligned query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #ifdef HIRSCHTWOBIT
     twoBitMvXElmFromStart(qryAln, 0);
     bitUC = getTwoBitElm(qryAln);
   #else
     bitUC = (uint8_t) *qryAln;
   #endif

   while(bitUC == defGapFlag)
   { /*Loop: Find the first aligned query base*/
      #ifdef HIRSCHTWOBIT
         twoBitMvToNextElm(qryAln);
         bitUC = getTwoBitElm(qryAln);
      #else
         ++qryAln;
         bitUC = (uint8_t) *qryAln;
      #endif

      ++qryFirstAlnBaseL;
      ++qryStr;
   } /*Loop: Find the first aligned query base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-05:
   ^  - Build the alignment
   ^  o fun-04 sec-05 sub-01:
   ^    - Check if reference has a gap (deletion in query)
   ^  o fun-04 sec-05 sub-02:
   ^    - Check if positions is an snp/match
   ^  o fun-04 sec-05 sub-03:
   ^    - Check if positions is insertion (query only)
   ^  o fun-04 sec-05 sub-04:
   ^    - Get the next base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-04 Sec-05 Sub-01:
   *  - Check if reference has a gap (deletion in query)
   \******************************************************/

   while(bitUC != defEndAlnFlag)
   { /*Loop: Add query aligned bases to alnStr*/
      #ifdef HIRSCHTWOBIT
      if(getTwoBitElm(refAln) == defGapFlag)
      #else
      if(*refAln == defGapFlag)
      #endif
      { /*If this is a deletion*/
         *alnIterStr = '-';

         #ifdef HIRSCHTWOBIT
            twoBitMvToNextElm(refAln);
         #else
            ++refAln;
         #endif

         ++alnIterStr;
         continue;
      } /*If this is a deletion*/

      /***************************************************\
      * Fun-04 Sec-05 Sub-02:
      *  - Check if positions is an snp/match
      \***************************************************/

      switch(bitUC)
      { /*Switch: Check the error type*/
         case 0: break;

         case defSnpFlag:
         case defMatchFlag:
         /*Case: match/snp*/
            *alnIterStr = *qryStr & (~32); /*Upper case*/
            ++alnIterStr;
            ++qryStr;

            #ifdef HIRSCHTWOBIT
               twoBitMvToNextElm(qryAln);
               twoBitMvToNextElm(refAln);
            #else
               ++qryAln;
               ++refAln;
            #endif
            break;
         /*Case: match/snp*/

         /************************************************\
         * Fun-04 Sec-05 Sub-03:
         *  - Check if positions is insertion (query only)
         \************************************************/

         case defGapFlag:
         /*Case: Insertions*/
            *alnIterStr = *qryStr | 32; /*Lower case*/
            ++alnIterStr;
            ++qryStr;

            #ifdef HIRSCHTWOBIT
               twoBitMvToNextElm(qryAln);
            #else
               ++qryAln;
            #endif

            break;
         /*Case: inerstions*/
      } /*Switch: Check the error type*/

      /***************************************************\
      * Fun-04 Sec-05 Sub-04:
      *  - Get the next base
      \***************************************************/

      #ifdef HIRSCHTWOBIT
         bitUC = getTwoBitElm(qryAln);
      #else
         bitUC = (uint8_t) *qryAln;
      #endif
   } /*Loop: Add query aligned bases to alnIterStruct*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-06:
   ^  - Clean up and return alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Nothing to clean up*/

    if(ampScoreST->refStartUL != refFirstAlnBaseL)
       ampScoreST->refStartUL = refFirstAlnBaseL;

    if(ampScoreST->qryStartUL != qryFirstAlnBaseL)
       ampScoreST->qryStartUL = qryFirstAlnBaseL;

    *alnIterStr = '\0'; /*Mark end of alignment*/
    --alnIterStr;

    /*Indels at the end are likely a sign of low quality
    ` bases, which are best trimmed off. Also patterns of
    ` indel-Match-indel-Match at the end is also likely
    ` a sign of a low quality section
    ` Or they are just an artifact of a global alignment.
    ' This issue is likely an error in my Hirschberg.
    ' TODO: Seed why I need this.
    */
    while(
          *alnIterStr == '-'       /*deletion*/
       || *alnIterStr & 32         /*insertion*/
       || *(alnIterStr - 1) & 32   /*low quality end*/
       || *(alnIterStr - 1) == '-' /*low quality end*/
    ) { /*Loop: trim the indels off the end*/
       *alnIterStr = '\0';
       --alnIterStr;

       if(*alnIterStr == '-' || !(*alnIterStr & 32))
          --(ampScoreST->refEndUL);
    } /*Loop: trim the indels off the end*/

    return alnStr;
} /*hirschbergToSeq*/

/*--------------------------------------------------------\
| Name: stichAmps (Fun-05:)
| Use:
|  - Stiches amplions into a consensus
| Input:
|  - alnSeqStr:
|    o The aligned amplicon sequence to stich into the
|      consensus
|  - ampScoreST:
|    o scoresStruct with the first reference base alnSeqStr
|      starts on (ampScoreST->refStartUL).
|    o ampScoreST->refStartUL should be in index 0
|  - conST:
|    o stichAmpST on the last added base in the consensus.
|      If 0, this will start a new consensus.
|  - refEndUL:
|    o The position of the last added base in conST
|    o This should be index 0;
|  - maskC:
|    o Character to mask amplicon with
| Output:
|  - Modifies:
|    o refEndUL to have to the position last added base
|      added by stichAmpST.
|  - Returns:
|    o A pointer to the last added base in conST.
|      - conST is a double linked list, so you can get to
|        the first base from the last base or the last base
|        from the first base.
|    o 0 if had a memory error (the list in conST is
|      freeded)
|  - Frees:
|    o the list in conST if had a memory error
| Note:
|  - stichAmps assumes that the amplicons have been sorted
|    by starting position on reference. With the amplicon
|    mapping to the first reference base coming first.
|    o This can be done with sortScoresStartLen() or
|      sortScoresStartLenIndex() in alnSeqSrc/scoresST.h.
\--------------------------------------------------------*/
struct stichAmpST * stichAmps(
   char *alnSeqStr,          /*Amplicon sequence*/
   struct scoresStruct *ampScoreST, /*Score for amplicon*/
   struct stichAmpST *conST, /*Consensus sequence (list)*/
   ulong *refEndUL,          /*Last ref base in consensus*/
   char maskC                /*What to use for maksing*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: stichAmps
   '  o fun-05 sec-01:
   '    - Variable declerations
   '  o fun-05 sec-02:
   '    - Check if this is the first amplicon or have a gap
   '      between the consensus and next amplicon 
   '  o fun-05 sec-03:
   '    - Add bases to the overlap
   '  o fun-05 sec-04:
   '    - Add bases to the overlap
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ulong overlapUL = 0;
   ulong gapStartUL = 0;
   struct stichAmpST *altBaseST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-02:
   ^  - Check if this is the first amplicon or have a gap
   ^    between the consensus and next amplicon 
   ^  o fun-05 sec-02 sub-01:
   ^    - Check if first amplicon, if so make first struct
   ^  o fun-05 sec-02 sub-02:
   ^    - Add masking to missing bases at start
   ^  o fun-05 sec-02 sub-03:
   ^    - Add in the bases in the amplicon
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-05 Sec-02 Sub-01:
   *  - Check if this is first amplicon
   \******************************************************/

   if(conST == 0 || *refEndUL < ampScoreST->refStartUL)
   { /*If: This is the first amplicon*/
      if(conST == 0)
      { /*If: this is the first base*/
         conST = malloc(sizeof(struct stichAmpST));
         if(conST == 0) return 0;
         initStichAmpST(conST);
         gapStartUL = 0;
      } /*If: this is the first base*/

      else
      { /*Else: There is a gap between the con and amp*/
         conST->nextBase=malloc(sizeof(struct stichAmpST));

         if(conST->nextBase == 0)
         { /*If: I had a memory error*/
            freeStichAmpSTList(&conST);
            return 0;
         } /*If: I had a memory error*/

         initStichAmpST(conST->nextBase);
         conST->nextBase->lastBase = conST;
         conST = conST->nextBase;
         gapStartUL = *refEndUL + 1;
      } /*Else: There is a gap between the con and amp*/

      /***************************************************\
      * Fun-05 Sec-02 Sub-02:
      *  - Add masking to missing bases at start
      \***************************************************/

      for(
         ulong iMask = gapStartUL;
         iMask < ampScoreST->refStartUL;
         ++iMask
      ){ /*Loop: Add in masking*/
         conST->baseC = maskC;
         conST->errC = (defMvMask & (~32));
            /*& ~32 is to make sure upper case. Lower case
            ` marks an insertion
            */
         conST->supportUL = 0;
         conST->depthUL = 0;

         conST->nextBase=malloc(sizeof(struct stichAmpST));

         if(conST->nextBase == 0)
         { /*If: I had a memory error*/
            freeStichAmpSTList(&conST);
            return 0;
         } /*If: I had a memory error*/

         initStichAmpST(conST->nextBase);
         conST->nextBase->lastBase = conST;
         conST = conST->nextBase;
      } /*Loop: Add in maksing*/

      goto stichAmpAddBases;
   } /*If: This is the first amplicon*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-03:
   ^  - Add bases to the overlap
   ^  o fun-05 sec-03 sub-01:
   ^    - Find the start of the overlap on the consensus
   ^  o fun-05 sec-03 sub-02:
   ^    - Add in the overlap bases
   ^  o fun-05 sec-03 sub-03:
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-05 Sec-03 Sub-01:
   *  - Find the start of the overlap on the consensus
   \******************************************************/

   overlapUL = *refEndUL - ampScoreST->refStartUL;

   while(overlapUL > 0)
   { /*Loop: Find start of overlap on the consensus*/
      if(conST->lastBase == 0) break; /*At start of seq*/
      ++(conST->depthUL); /*Add in the overlap*/

      switch(conST->errC)
      { /*Switch: Check the error type*/
         case 0: break; /*Blank entry, should not happen*/
         case defMvIns: break;
         case defMvDel: --overlapUL; break;
         case defMvSnp: --overlapUL; break;
      } /*Switch: Check the error type*/

      conST = conST->lastBase;
   } /*Loop: Find start of overlap on the consensus*/

   /******************************************************\
   * Fun-05 Sec-03 Sub-02:
   *  - Add in the overlap bases
   *  o fun-05 sec-03 sub-02 cat-01:
   *    - Handle insertions in both consensus and amplicon
   *  o fun-05 sec-03 sub-02 cat-02:
   *    - Case: insertion present, but have a new
   *      alternative base.
   *  o fun-05 sec-03 sub-02 cat-03:
   *    - add a new insertion into the consensus
   *  o fun-05 sec-03 sub-02 cat-04:
   *    - Consensus has insertion, but not the amplicon
   *  o fun-05 sec-03 sub-02 cat-05:
   *    - Handle deletions/snps/matchs/masked bases
   \******************************************************/

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun-05 Sec-03 Sub-02 Cat-01:
   +  - Handle insertions in both consensus and amplicon
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   while(conST->nextBase != 0)
   { /*Loop: Add in the overlapping bases*/
      if(*alnSeqStr & 32 && *alnSeqStr != '-')
      { /*If: I have an insertion*/
         if(conST->errC == defMvIns)
         { /*If: the overlap supports this*/
            altBaseST = conST;

            while(altBaseST->baseC != *alnSeqStr)
            { /*Loop: Find the matching alternative base*/
               if(altBaseST->altBase == 0) break;
               altBaseST = altBaseST->altBase;
            } /*Loop: Find the matching alternative base*/

            /*++++++++++++++++++++++++++++++++++++++++++++\
            + Fun-05 Sec-03 Sub-02 Cat-02:
            +  - Case: insertion present, but have a new
            +    alternative base.
            \++++++++++++++++++++++++++++++++++++++++++++*/

            if(altBaseST == 0)
            { /*If: This is a new alternative base*/
               altBaseST->altBase =
                  malloc(sizeof(struct stichAmpST));

               if(altBaseST->altBase == 0)
               { /*If: I had a memory error*/
                  freeStichAmpSTList(&conST);
                  return 0;
               } /*If: I had a memory error*/

               initStichAmpST(altBaseST->altBase);
               altBaseST = altBaseST->altBase;               
               altBaseST->baseC = *alnSeqStr;
               altBaseST->errC = defMvIns;
            } /*If: This is a new alternative base*/

            conST = conST->nextBase;
            ++(altBaseST->supportUL);
            ++alnSeqStr;
            continue;
         } /*If: the overlap supports this*/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-05 Sec-03 Sub-02 Cat-03:
         +  - add a new insertion into the consensus
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         else
         { /*Else: The consnesus does not have an ins*/   
            do { /*Loop: Add in the new insertions*/
               altBaseST=malloc(sizeof(struct stichAmpST));

               if(altBaseST == 0)
               { /*If: I had a memory error*/
                  freeStichAmpSTList(&conST);
                  return 0;
               } /*If: I had a memory error*/

               initStichAmpST(altBaseST);

               altBaseST->lastBase = conST;
               altBaseST->nextBase = conST->nextBase;
               conST->nextBase = altBaseST;
               conST = altBaseST;

               conST->baseC = *alnSeqStr;
               conST->errC = defMvIns;
               conST->supportUL = 1;
               conST->depthUL = conST->lastBase->depthUL;
           
               ++alnSeqStr;
            } while(*alnSeqStr & 32 && *alnSeqStr != '-');

            conST = conST->nextBase;
            ++alnSeqStr;
            continue;
         } /*Else: The consnesus does not have an ins*/   
      } /*If: I have an insertion*/

      /*++++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun-05 Sec-03 Sub-02 Cat-04:
      +  - Consensus has insertion, but not the amplicon
      \++++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(conST->errC == defMvIns)
      { /*Else if: the consensus has an insertion*/
         while(conST->errC == defMvIns)
            conST = conST->nextBase;
         continue;
      } /*Else if: the consensus has an insertion*/

      /*++++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun-05 Sec-03 Sub-02 Cat-05:
      +  - Handle deletions/snps/matchs/masked bases
      \++++++++++++++++++++++++++++++++++++++++++++++++++*/

      altBaseST = conST;

      while(altBaseST->baseC != *alnSeqStr)
      { /*Loop: Find the matching alternative base*/
         if(altBaseST->altBase == 0) break;
         altBaseST = altBaseST->altBase;
      } /*Loop: Find the matching alternative base*/

      if(altBaseST == 0)
      { /*If: This is a new alternative base*/
         altBaseST->altBase =
            malloc(sizeof(struct stichAmpST));

         if(altBaseST->altBase == 0)
         { /*If: I had a memory error*/
            freeStichAmpSTList(&conST);
            return 0;
         } /*If: I had a memory error*/

         initStichAmpST(altBaseST->altBase);
         altBaseST = altBaseST->altBase;               
         altBaseST->baseC = *alnSeqStr;

         if(*alnSeqStr == '-') altBaseST->errC = defMvDel;
         else altBaseST->errC = defMvSnp;
      } /*If: This is a new alternative base*/

      ++(altBaseST->supportUL);
      ++alnSeqStr;
      conST = conST->nextBase;
   } /*Loop: Add in the overlapping bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-04:
   ^  - Add bases to the overlap
   ^  o fun-05 sec-03 sub-01:
   ^    - Add non-overlapping amplicon bases to consensus
   ^  o fun-05 sec-03 sub-02:
   ^    - Clean up the last base (empty) and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-05 Sec-02 Sub-01:
   *  - Add non-overlapping amplicon bases to consensus
   \******************************************************/

   stichAmpAddBases:

   while(*alnSeqStr != '\0')
   { /*Loop: Add bases to the consensus*/
      conST->baseC = *alnSeqStr;

      if(*alnSeqStr == '-')
         conST->errC = defMvDel;
      else if(*alnSeqStr & 32)
         conST->errC = defMvIns;
      else
         conST->errC = defMvSnp;

      ++alnSeqStr;

      conST->depthUL = 1;
      conST->supportUL = 1;

      /*Add a new base to the list*/
      conST->nextBase=malloc(sizeof(struct stichAmpST));

      if(conST->nextBase == 0)
      { /*If: I had a memory error*/
         freeStichAmpSTList(&conST);
         return 0;
      } /*If: I had a memory error*/

      initStichAmpST(conST->nextBase);
      conST->nextBase->lastBase = conST;
      conST = conST->nextBase;
   } /*Loop: Add bases to the consensus*/

   /***************************************************\
   * Fun-05 Sec-04 Sub-02:
   *  - Clean up the last base and return
   \***************************************************/

   /*I have on extra base that I need to free*/
   conST = conST->lastBase;
   freeStichAmpST(&(conST->nextBase));
   conST->nextBase = 0;
   *refEndUL = ampScoreST->refEndUL;

   return conST;
} /*stichAmps*/
