/*#########################################################
# Name: stichMinimapFun
# Use:
#  - Holds functions used in stiching amplicons into a
#    consensus. These functions use minimap2
# Libraries:
#  - "sitchAmpStruct.h"              (No .c file)
#  - "samFunSrc/trimSam.h"
#  - "samFunSrc/cStrFun.h"           (No .c file)
#  - "stichSetStruct.h"              (No .c file)
#  o "stichDefaults.h"               (No .c file)
#  o "samFunSrc/samEntryStruct.h"    (No .c file)
#  o "samFunSrc/cStrToNumberFun.h"   (No .c file)
#  o "alnSeqSrc/dataTypeShortHand.h" (No .c file)
#  o "alnSeqSrc/seqStruct.h"
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdint.h>
#  o <stdio.h>
#  o <string.h>
# Requires:
#  - Minimap2 be in your file path
#########################################################*/

#include "stichMinimapFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' stichMinimapFun SOF: Start Of Functions
'  - Holds the functions used to stich amplicons together.
'  - This variation uses minimap2 for alignment
'  o fun-01 getAmpPosMinimap:
'    - Gets the position of each amplicon on the reference
'  o fun-02 stichAmpConMinimap:
'    - Stiches together amplicon consensuses to to make
'      an consensus genome
'  o fun-03 stichAmpsMinimap:
'    - Stiches amplions into a consensus
'  o fun-04 samEntryToAlnSeq:
'    - This uses a samEntry struct to get an aligned
'  o fun-05 freeSamEntryAry:
'    - Frees an array of samEntry structures
'  o fun-06 stichAmpConToCStr:
'    - Convert a stichAmpST consensus to a merged, c-string
'      consensus
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: getAmpPosMinmap (Fun-01:)
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
|  - settings:
|    o stichSet struct with the number of threads to use
|  - outErrUC:
|    o Holds any errors that happened (see output)
| Output:
|   - Modifies:
|     o numAmpsUL to hold the number of amplicons
|     o indexAryUL to hold the index of every sequence.
|     o outErrUC to hold the error, if an error occured
|       - 0: No error
|       - 1: amplicon fasta file error
|       - 2: reference fasta file error
|       - 3: reference fasta has multiple sequences
|       - 4: Could not allocate memory
|       - 5: minimap2 errored out or memory error
|       - 6: No amplicons were kept (none mapped to ref)
|   - Returns:
|     o Pointer to array of samEntry structures
|     o 0 for memory errors
\--------------------------------------------------------*/
struct samEntry * getAmpPosMinimap(
   char *refFileStr,  /*Name of fasta file with reference*/
   char *ampFileStr,  /*Name of fasta file with amplicons*/
   ulong *numAmpsUL,  /*Will have number amplicons kept*/
   struct stichSet *settings, /*Number of threads to use*/
   uchar *outErrUC   /*Reports the error type*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: getAmpPosMinimap
   '  - Gets the position of each amplicon on the reference
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - File checks and memory allocation
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

   char minimap2CmdStr[2048];
   char threadsStr[8];
   char *tmpStr = 0;
   uchar errUC = 0;

   struct samEntry *ampsST = 0;
   struct samEntry *ampIterST = 0;
   struct seqStruct seqST;

   FILE *stdinFILE = 0; /*Minimap2 output*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Do checks and memory alloctions
   ^  o fun-01 sec-02 sub-01:
   ^    - Check amplicon fasta and get number amplicons
   ^  o fun-01 sec-02 sub-02:
   ^    - Check if the reference fasta is valid
   ^  o fun-01 sec-02 sub-03:
   ^    - Make sure only one sequence in reference file
   ^  o fun-01 sec-02 sub-04:
   ^    - Do memory allocations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-02 Sub-01:
   *  - Check amplicon fasta and get number amplicons
   \******************************************************/

   initSeqST(&seqST);
   *numAmpsUL = 0;

   stdinFILE = fopen(ampFileStr, "r");

   if(ampFileStr == 0) 
   { /*If: the amplicon file could not be opened*/
      *outErrUC = 1;
      return 0;
   } /*If: the amplicon file could not be opened*/

   errUC = 1;

   /*Find the number of amplicons*/
   while(errUC & 1)
   { /*Loop: Find number of amplicons & check file*/
      errUC = readFaSeq(stdinFILE, &seqST);

      if(errUC > 1)
      { /*If: this is an invalid fasta file*/
         *outErrUC = 1;
         fclose(stdinFILE);
         stdinFILE = 0;
         return 0;
      } /*If: this is an invalid fasta file*/

      ++(*numAmpsUL);
   } /*Loop: Find number of amplicons & check file*/

   fclose(stdinFILE);
   stdinFILE = 0;

   /******************************************************\
   * Fun-01 Sec-02 Sub-02:
   *  - Check if the reference fasta is valid
   \******************************************************/
 
   stdinFILE = fopen(refFileStr, "r");

   if(stdinFILE == 0)
   { /*If: I could not open the reference file*/
      *outErrUC = 2;
      return 0;
   } /*If: I could not open the reference file*/

   errUC = readFaSeq(stdinFILE, &seqST);

   if(errUC > 1)
   { /*If: this is an invalid fasta file*/
      fclose(stdinFILE);
      stdinFILE = 0;
      freeSeqST(&seqST, 0); /*Struct is on stack*/
      *outErrUC = 2;
      return 0;
   } /*If: this is an invalid fasta file*/

   /******************************************************\
   * Fun-01 Sec-02 Sub-03:
   *  - Make sure only one sequence in reference file
   \******************************************************/

   errUC = readFaSeq(stdinFILE, &seqST);

   /*At this point my checks are done*/
   fclose(stdinFILE);
   stdinFILE = 0;
   freeSeqST(&seqST, 0); /*Struct is on stack*/

   if(errUC)
   { /*If: this has multipe entries*/
      *outErrUC = 3;
      return 0;
   } /*If: this has multipe entries*/

   /******************************************************\
   * Fun-01 Sec-02 Sub-04:
   *  - Do memory allocations
   \******************************************************/

   ampsST = malloc(sizeof(struct samEntry) * *numAmpsUL);

   if(ampsST == 0)
   { /*If: I had a memory error*/
      *outErrUC = 4;
      return 0; /*Memory error*/
   } /*If: I had a memory error*/

   ampIterST = ampsST;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Align each sequence to the reference
   ^  o fun-01 sec-03 sub-01:
   ^    - Run minimap2 and check if have output
   ^  o fun-01 sec-03 sub-02:
   ^    - Move past the header to the first sequence
   ^  o fun-01 sec-03 sub-03:
   ^    - Read in all minimap2 entries
   ^  o fun-01 sec-03 sub-04:
   ^    - final error checks
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-03 Sub-01:
   *  - Run minimap2 and check if have output
   \******************************************************/


    tmpStr=cStrCpInvsDelm(minimap2CmdStr,stichMinimap2CMD);

    uCharToCStr(threadsStr, settings->threadsUC);
    tmpStr = cpParmAndArg(tmpStr, "-t", threadsStr);

    tmpStr = cpParmAndArg(tmpStr, refFileStr, ampFileStr);

   /*Run minimap2*/
   stdinFILE = popen(minimap2CmdStr, "r");
   if(stdinFILE == 0) goto stichAlnAmpsErr;

   initSamEntry(ampIterST);
   errUC = readSamLine(ampIterST, stdinFILE);

   /*Check if I have something to work with*/
   if(!(errUC & 1))
   { /*If: I had an error*/
      pclose(stdinFILE);

      stichAlnAmpsErr:

      while(ampIterST >= ampsST)
      { /*Loop: Free all internal variables*/
         freeStackSamEntry(ampIterST);
         --ampIterST;
      } /*Loop: Free all internal variables*/

      free(ampsST); /*Free the entire array*/
      ampsST = 0;

      *outErrUC = 5;
      return 0;
   } /*If: I had an error*/

   /******************************************************\
   * Fun-01 Sec-03 Sub-02:
   *  - Move past the header to the first sequence
   \******************************************************/

   while(*ampIterST->samEntryCStr == '@')
   { /*Loop: Get past the header*/
      initSamEntry(ampIterST);
      errUC = readSamLine(ampIterST, stdinFILE);

      if(!(errUC & 1))
      { /*If: I had an error or file has only headers*/
         pclose(stdinFILE);

         while(ampIterST >= ampsST)
         { /*Loop: Free all internal variables*/
            freeStackSamEntry(ampIterST);
            --ampIterST;
         } /*Loop: Free all internal variables*/

         free(ampsST); /*Free the entire array*/
         ampsST = 0;

         *outErrUC = 5;
         return 0;
      } /*If: I had an error or file has only headers*/
   } /*Loop: Get past the header*/

   /******************************************************\
   * Fun-01 Sec-03 Sub-03:
   *  - Read in all minimap2 entries
   \******************************************************/
   
   *numAmpsUL = 0;

   while(errUC & 1)
   { /*Loop: Find the starting position for each amplicon*/
      trimSamEntry(ampIterST);

      /*Check if have an alignment I can use*/
      if(ampIterST->flagUSht & (4 | 256 | 2048))
         errUC = readSamLine(ampIterST, stdinFILE);
         /*4=unmapped, 256=secondary, 2048=supplemental*/
     
      ++(*numAmpsUL); /*Counter for number amplicons kept*/
      ++ampIterST;
      initSamEntry(ampIterST);

      errUC = readSamLine(ampIterST, stdinFILE);
   } /*Loop: Find the starting position for each amplicon*/

   pclose(stdinFILE);

   /******************************************************\
   * Fun-01 Sec-03 Sub-04:
   *  - Final error checks
   \******************************************************/

   if(errUC > 2)
   { /*If: I had an error (only option is memory)*/
      while(ampIterST >= ampsST)
      { /*Loop: Free all internal variables*/
         freeStackSamEntry(ampIterST);
         --ampIterST;
      } /*Loop: Free all internal variables*/

      free(ampsST); /*Free the entire array*/
      ampsST = 0;

      *outErrUC = 5;
      return 0;
   } /*If: I had an error (only option is memory)*/

   if(*numAmpsUL < 1)
   { /*If: no amplicons were kept*/
      freeStackSamEntry(ampIterST);

      free(ampsST); /*Free the entire array*/
      ampsST = 0;

      *outErrUC = 6;
      return 0;
   } /*If: no amplicons were kept*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Reset file pointer and sort the sequences by
   ^    reference starting positions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sortSamEntriesStartLen(&ampsST, 0, *numAmpsUL - 1);

   *outErrUC = 0; /*No error*/
   return ampsST;
} /*getAmpPosMinmap*/

/*--------------------------------------------------------\
| Name: stichAmpConMinimap (Fun-02:)
| Use:
|  - Stiches together amplicon consensuses to to make
|    an stichAmpST scaffold with alternative bases
| Input:
|  - ampsAryST:
|    o array of samEntry structs wich contain the sam file
|      entry for all amplicons to stich together
|  - numAmpsUL:
|    o Number of amplicons in ampsAryST
|  - stichSetST:
|    o Settings for stiching amplicons together.
| Output:
|  - Returns:
|    o List of stichAmpST structs that have the
|      uncollapsed (with alternative bases) scaffold
|    o 0 for error
| Note:
|  - This function works with output from getAmpPosMinimap
\--------------------------------------------------------*/
struct stichAmpST * stichAmpConMinimap(
   struct samEntry *ampsAryST,  /*Array of amplicons*/
   ulong numAmpsUL,             /*Number amplicons*/
   struct stichSet *stichSetST  /*Settings for stich*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: stichAmpConMinimap
   '  - Stiches together amplicon consensuses to to make
   '    an stichAmpST scaffold with alternative bases
   '  o fun-02 sec-01:
   '     - Variable declerations
   '  o fun-02 sec-02:
   '    - Align each amplicon sequence
   '  o fun-02 sec-03:
   '    - Stich together each amplicon sequence
   '  o fun-02 sec-04:
   '    - Move to first base and return consensus
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *alnSeqStr = 0;
   ulong refEndUL = 0;
   struct stichAmpST *conST = 0;
   struct samEntry *ampIterST = ampsAryST;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-02:
   ^  - Align each amplicon sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(numAmpsUL > 0)
   { /*While I have sequences to stich together*/
      --numAmpsUL;

      alnSeqStr = samEntryToAlnSeq(ampIterST);

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
         stichAmpsMinimap(
            alnSeqStr,
            (ulong) ampIterST->posOnRefUInt - 1,
            conST,
            &refEndUL,    /*Last ref base in consensus*/
            stichSetST->maskC
      );
      /*(ulong) ampIterST->posOnRefUInt - 1 is to account
      ` for minimap2 being index 1
      */

      if(conST == 0)
      { /*If: I had a memory error*/
         free(alnSeqStr);
         alnSeqStr = 0;
         /*conST has already been freeded*/
         return 0;
      } /*If: I had a memory error*/

      ++ampIterST;
      free(alnSeqStr);
      alnSeqStr = 0;
   } /*While I have sequences to stich together*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-04:
   ^  - Move to first base and return consensus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(conST->lastBase != 0) conST = conST->lastBase;
   return conST;
} /*stichAmpConMinimap*/

/*--------------------------------------------------------\
| Name: stichAmpsMinimap (Fun-03:)
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
|  - stichAmpsMinimap assumes that the amplicons have been sorted
|    by starting position on reference. With the amplicon
|    mapping to the first reference base coming first.
|    o This can be done with sortScoresStartLen() or
|      sortScoresStartLenIndex() in alnSeqSrc/scoresST.h.
\--------------------------------------------------------*/
struct stichAmpST * stichAmpsMinimap(
   char *alnSeqStr,          /*Amplicon sequence*/
   ulong ampStartUL,        /*First ref base in amplicon*/
   struct stichAmpST *conST, /*Consensus sequence (list)*/
   ulong *refEndUL,          /*Last ref base in consensus*/
   char maskC                /*What to use for maksing*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: stichAmpsMinimap
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Check if this is the first amplicon or have a gap
   '      between the consensus and next amplicon 
   '  o fun-03 sec-03:
   '    - Add bases to the overlap
   '  o fun-03 sec-04:
   '    - Add bases to the overlap
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ulong overlapUL = 0;
   ulong gapStartUL = 0;
   ulong insDepthUL = 0; /*Depth of new ins in overlap*/
   struct stichAmpST *altBaseST = 0;
   struct stichAmpST *lastBaseST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Check if this is the first amplicon or have a gap
   ^    between the consensus and next amplicon 
   ^  o fun-03 sec-02 sub-01:
   ^    - Check if first amplicon, if so make first struct
   ^  o fun-03 sec-02 sub-02:
   ^    - Add masking to missing bases at start
   ^  o fun-03 sec-02 sub-03:
   ^    - Add in the bases in the amplicon
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-01:
   *  - Check if this is first amplicon
   \******************************************************/

   if(conST == 0 || *refEndUL < ampStartUL)
   { /*If: This is the first amplicon*/
      if(conST == 0)
      { /*If: this is the first base*/
         conST = malloc(sizeof(struct stichAmpST));
         if(conST == 0) return 0;

         initStichAmpST(conST);
         conST->baseC = maskC & (~32);

         gapStartUL = 1;
         *refEndUL = 0; /*This will be overcounted*/
      } /*If: this is the first base*/

      else
      { /*Else: There is a gap between the con and amp*/
         gapStartUL = *refEndUL + 1; /*avoid overcounting*/
      } /*Else: There is a gap between the con and amp*/

      /***************************************************\
      * Fun-03 Sec-02 Sub-02:
      *  - Add masking to missing bases at start
      \***************************************************/

      for(
         ulong iMask = gapStartUL;
         iMask < ampStartUL;
         ++iMask
      ){ /*Loop: Add in masking*/
         conST->nextBase=malloc(sizeof(struct stichAmpST));

         if(conST->nextBase == 0)
         { /*If: I had a memory error*/
            freeStichAmpSTList(&conST);
            return 0;
         } /*If: I had a memory error*/

         initStichAmpST(conST->nextBase);
         conST->nextBase->lastBase = conST;
         conST = conST->nextBase;

         conST->baseC = maskC & (~32);
            /*& ~32 is to make sure upper case. Lower case
            ` marks an insertion
            */
         ++(*refEndUL);
      } /*Loop: Add in maksing*/

      /*--(*refEndUL);*/
      goto stichAmpAddBases;
   } /*If: This is the first amplicon*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - Add bases to the overlap
   ^  o fun-03 sec-03 sub-01:
   ^    - Find the start of the overlap on the consensus
   ^  o fun-03 sec-03 sub-02:
   ^    - Add in the overlap bases
   ^  o fun-03 sec-03 sub-03:
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-03 Sub-01:
   *  - Find the start of the overlap on the consensus
   \******************************************************/

   overlapUL = *refEndUL - ampStartUL;
   lastBaseST = conST;

   while(overlapUL > 0)
   { /*Loop: Find start of overlap on the consensus*/
      if(conST->lastBase == 0) break; /*At start of seq*/
      ++(conST->depthUL); /*Add in the overlap*/

      if(conST->baseC < 64 + 32) --overlapUL;
      /*'-' < 64, snp/match is uppercase (64 to 90)
      ` So this only does not fire for insertions.
      `*/

      conST = conST->lastBase;
   } /*Loop: Find start of overlap on the consensus*/

   /******************************************************\
   * Fun-03 Sec-03 Sub-02:
   *  - Add in the overlap bases
   *  o fun-03 sec-03 sub-02 cat-01:
   *    - Handle insertions in both consensus and amplicon
   *  o fun-03 sec-03 sub-02 cat-02:
   *    - Case: insertion present, but have a new
   *      alternative base.
   *  o fun-03 sec-03 sub-02 cat-03:
   *    - add a new insertion into the consensus
   *  o fun-03 sec-03 sub-02 cat-04:
   *    - Consensus has insertion, but not the amplicon
   *  o fun-03 sec-03 sub-02 cat-05:
   *    - Handle deletions/snps/matchs/masked bases
   \******************************************************/

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun-03 Sec-03 Sub-02 Cat-01:
   +  - Handle insertions in both consensus and amplicon
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   while(conST != 0)
   { /*Loop: Add in the overlapping bases*/
      if(*alnSeqStr >= 64 + 32)
      { /*If: I have an insertion*/
         if(conST->baseC >= 64 + 32)
         { /*If: This isnertion is in the consensus*/
            altBaseST = conST;

            while(altBaseST->altBase != 0)
            { /*Loop: Find the matching alternative base*/
               if(altBaseST->baseC == (*alnSeqStr | 32))
                  break;
               altBaseST = altBaseST->altBase;
            } /*Loop: Find the matching alternative base*/

            /*++++++++++++++++++++++++++++++++++++++++++++\
            + Fun-03 Sec-03 Sub-02 Cat-02:
            +  - Case: insertion present, but have a new
            +    alternative base.
            \++++++++++++++++++++++++++++++++++++++++++++*/

            if(altBaseST->baseC != (*alnSeqStr | 32))
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
            } /*If: This is a new alternative base*/

            conST = conST->nextBase;
            ++(altBaseST->supportUL);
            ++alnSeqStr;
            continue;
         } /*If: This isnertion is in the consensus*/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-03 Sec-03 Sub-02 Cat-03:
         +  - add a new insertion into the consensus
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         else
         { /*Else: The consnesus does not have an ins*/   
            bitMax(
               insDepthUL,
               conST->depthUL,
               conST->lastBase->depthUL
            ); /*Get the minimum*/
              /* I figure that the overlap should have
              `  detected the insertions. So, it should
              `  be the max of the previous or current
              `  base. This is not percise.
              */
           
            do { /*Loop: Add in the new insertions*/
               altBaseST=malloc(sizeof(struct stichAmpST));

               if(altBaseST == 0)
               { /*If: I had a memory error*/
                  freeStichAmpSTList(&conST);
                  return 0;
               } /*If: I had a memory error*/

               initStichAmpST(altBaseST);

               /*Insert the insertion behind the current
               ` base. conST will be one base ahead of the
               ` insertion.
               */
               altBaseST->lastBase = conST->lastBase;
               altBaseST->nextBase = conST;
               conST->lastBase = altBaseST;

               altBaseST->baseC = *alnSeqStr;
               altBaseST->supportUL = 1;
               altBaseST->depthUL = insDepthUL;

               ++alnSeqStr;
            } while(*alnSeqStr & 32 && *alnSeqStr != '-');

            conST = conST->nextBase;
            ++alnSeqStr;
            continue;
         } /*Else: The consnesus does not have an ins*/   
      } /*If: I have an insertion*/

      /*++++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun-03 Sec-03 Sub-02 Cat-04:
      +  - Consensus has insertion, but not the amplicon
      \++++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(conST->baseC >= 64 + 32)
      { /*Else if: the consensus has an insertion*/
         while(conST->baseC >= 64 + 32)
            conST = conST->nextBase;
         continue;
      } /*Else if: the consensus has an insertion*/

      /*++++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun-03 Sec-03 Sub-02 Cat-05:
      +  - Handle deletions/snps/matchs/masked bases
      \++++++++++++++++++++++++++++++++++++++++++++++++++*/

      altBaseST = conST;

      while(altBaseST->altBase != 0)
      { /*Loop: Find the matching alternative base*/
         if(altBaseST->baseC == *alnSeqStr) break;
           /*Every base should be in uppercase. I want
           ` to avoid & ~32 because that could confuse
           ` deletions with actual bases.
           */
         altBaseST = altBaseST->altBase;
      } /*Loop: Find the matching alternative base*/

      if(altBaseST->baseC != *alnSeqStr)
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
         /*Deletions are set up as '-' and matchs/snps as
         ` uppercase already. So, this will detect both
         */
      } /*If: This is a new alternative base*/

      ++(altBaseST->supportUL);
      ++alnSeqStr;
      conST = conST->nextBase;
   } /*Loop: Add in the overlapping bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-04:
   ^  - Add bases to the overlap
   ^  o fun-03 sec-03 sub-01:
   ^    - Add non-overlapping amplicon bases to consensus
   ^  o fun-03 sec-03 sub-02:
   ^    - Clean up the last base (empty) and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-01:
   *  - Add non-overlapping amplicon bases to consensus
   \******************************************************/

   stichAmpAddBases:

   if(lastBaseST != 0) conST = lastBaseST;

   while(*alnSeqStr != '\0')
   { /*Loop: Add bases to the consensus*/
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

      conST->baseC = *alnSeqStr;
      conST->depthUL = 1;
      conST->supportUL = 1;

      if(*alnSeqStr == '-') ++(*refEndUL); /*Deletion*/
      else if(*alnSeqStr < 64 + 32) ++(*refEndUL); /*snp*/

      ++alnSeqStr;
   } /*Loop: Add bases to the consensus*/

   /***************************************************\
   * Fun-03 Sec-04 Sub-02:
   *  - Clean up the last base and return
   \***************************************************/

   /*Account for overcounting*/
   /*if(*refEndUL != gapStartUL)*/ --(*refEndUL);

   /*I have on extra base that I need to free*/
   conST = conST->lastBase;
   conST->nextBase = 0;

   return conST;
} /*stichAmpsMinimap*/

/*--------------------------------------------------------\
| Name: samEntryToAlnSeq (Fun-04:)
| Use:
|  - This uses samEntry struct to get an aligned sequence.
|  - This function assumes that the samEntry has a sequence
|    and cigar entry. This means the sequence is mapped and
|    is not a secondary or supplemental alignment
| Input:
|  - ampST:
|    o samEntry struct with sequence to convert
| Output:
|  - Returns
|    o A c-string with the aligned query sequence.
|      Insertions are in lower case, with deletions as '-'.
\--------------------------------------------------------*/
char * samEntryToAlnSeq(
  struct samEntry *ampST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: samEntryToAlnSeq
   '  - This uses a samEntry struct to get an aligned
   '    sequence.
   '  - This function assumes that the samEntry has a
   '    sequence and cigar entry. This means that the
   '    sequence is mapped and is not a secondary or
   '    supplemental alignment
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

   char *seqStr = ampST->seqCStr;
   char *cigStr = ampST->cigarCStr;
   
   char *alnStr = 0;
   char *alnIterStr = 0;

   uint32_t lenEntryUI = 0;
   ulong lenAlnUL=ampST->readLenUInt + ampST->numDelUInt+1;
     /*+ 1 to convert to index 1 (values are index 0)*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-02:
   ^  - Allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   alnStr = malloc(sizeof(char) * lenAlnUL);
   if(alnStr == 0) return 0; /*Memory error*/
   alnIterStr = alnStr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-03:
   ^  - Make sure there is not softmask at the start
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   cigStr = cStrToUInt(cigStr, &lenEntryUI);

   if(*cigStr == 'S')
   { /*If: The samEntry was not trimmed*/
      for(uint32_t uiCig = 0; uiCig < lenEntryUI; ++uiCig)
         ++seqStr;

      /*Find the next cigar entry*/
      lenEntryUI = 0;
      ++cigStr; /*Get off letter*/
      cigStr = cStrToUInt(cigStr, &lenEntryUI);
   } /*If: The samEntry was not trimmed*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-05:
   ^  - Build the alignment
   ^  o fun-04 sec-05 sub-01:
   ^    - Check if positions is an snp/match
   ^  o fun-04 sec-05 sub-02:
   ^    - Check if reference has a gap (deletion in query)
   ^  o fun-04 sec-05 sub-03:
   ^    - Check if positions is insertion (query only)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-04 Sec-05 Sub-01:
   *  - Check if positions is an snp/match
   \******************************************************/

   while(*cigStr > 32) /*Ends in a tab (no spaces)*/
   { /*Loop: Add query aligned bases to alnStr*/

      if(*cigStr == '=' || *cigStr == 'X' || *cigStr =='M')
      { /*If: I have a match/snp*/
         for(uint32_t uiCig = 0;uiCig < lenEntryUI;++uiCig)
         { /*Loop: Add in all the snp/match entries*/
            *alnIterStr = *seqStr & (~32); /*Upper case*/
            ++alnIterStr;
            ++seqStr;
         } /*Loop: Add in all the snp/match entries*/

         lenEntryUI = 0;
         ++cigStr; /*Get off letter*/
         cigStr = cStrToUInt(cigStr, &lenEntryUI);
         continue;
      } /*If: I have a match/snp*/

      /***************************************************\
      * Fun-04 Sec-05 Sub-02:
      *  - Check if the cigar entry is a deletion
      \****************************************************/

      else if(*cigStr == 'D')
      { /*Else If: this is a deletion*/
         for(uint32_t uiCig = 0;uiCig < lenEntryUI;++uiCig)
         { /*Loop: Add in all the deletion entries*/
            *alnIterStr = '-';
            ++alnIterStr;
         } /*Loop: Add in all the deletion entries*/

         /*Find the next cigar entry*/
         lenEntryUI = 0;
         ++cigStr; /*Get off letter*/
         cigStr = cStrToUInt(cigStr, &lenEntryUI);
         continue;
      } /*Else If: this is a deletion*/

      /***************************************************\
      * Fun-04 Sec-05 Sub-03:
      *  - Check if positions is insertion
      \***************************************************/

      if(*cigStr=='I')
      { /*Else If: I have an insertion*/
         for(uint32_t uiCig = 0;uiCig < lenEntryUI;++uiCig)
         { /*Loop: Add in all the insertion entries*/
            *alnIterStr = *seqStr | 32; /*Lower case*/
            ++alnIterStr;
            ++seqStr;
         } /*Loop: Add in all the insertion entries*/

         lenEntryUI = 0;
         ++cigStr; /*Get off letter*/
         cigStr = cStrToUInt(cigStr, &lenEntryUI);
         continue;
      } /*Else If: I have an insertion*/
   } /*Loop: Add query aligned bases to alnIterStruct*/

   if(lenEntryUI > 0)
   { /*If: I have matches at the end*/
      for(uint32_t uiCig = 0;uiCig < lenEntryUI;++uiCig)
      { /*Loop: Add in all the match entries*/
         *alnIterStr = *seqStr & (~32); /*Upper case*/
         ++alnIterStr;
         ++seqStr;
      } /*Loop: Add in all the match entries*/
   } /*If: I have matches at the end*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-06:
   ^  - Clean up and return alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *alnIterStr = '\0'; /*Mark end of alignment*/
    return alnStr;
} /*samEntryToAlnSeq*/

/*--------------------------------------------------------\
| Name: freeSamEntryAry (Fun-06:)
| Use:
|  - Frees an array of samEntry structures
| Input:
|  - samAryST:
|    o Pointer to samEntry array to free
|  - numAmpsUL:
|    o Number of amplicons in the array
| Output:
|  - Frees:
|    - samAryST and its internal heap variables
|  - Modifies:
|    - samAryST to point to 0
\--------------------------------------------------------*/
void freeSamEntryAry(
   struct samEntry **samAryST,
   ulong numAmpsUL
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: freeSamEntryAry
   '  - Frees an array of samEntry structures
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   struct samEntry *ampIterST = *samAryST;

   while(numAmpsUL > 0)
   { /*Loop: Free all internal variables*/
      --numAmpsUL;
      freeStackSamEntry(ampIterST);
      ++ampIterST;
   } /*Loop: Free all internal variables*/

   free(*samAryST);
   *samAryST = 0;
   return;
} /*freeSamEntryAry*/

/*--------------------------------------------------------\
| Name: stichAmpConToCStr (Fun-06:)
| Use:
|  - Converts a list stichAmpST structs (consensus with
|    alternative bases) to a c-string with no alternative
|    bases.
| Input:
|  - conST:
|    o stichAmpST list to merge into a single consensus
|  - settings:
|    o stichSet struct with settings to use to when
|      merging the consensus
| Output:
|  - Returns:
|    o C-string with the merged consensus
|    o 0 for memory errors
\--------------------------------------------------------*/
char * stichAmpConToCStr(
   struct stichAmpST *conST, /*Has consensus to make*/
   struct stichSet *settings  /*Settings for stich*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: stichAmpConToCStr
   '  - Convert a stichAmpST consensus to a merged, c-string
   '    consensus
   '  o fun-06 sec-01:
   '    - Set up for colapsing the consensus list
   '  o fun-06 sec-02:
   '    - Find the possible largest consensus length
   '  o fun-06 sec-03:
   '    - Allocate memory for the consensus
   '  o fun-06 sec-04:
   '    - Merge alternative bases into a single consensus
   '  o fun-06 sec-04:
   '    - Clean up and exit
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-01:
   ^  - Set up for colapsing the consensus list
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *conSeqStr = 0;
   char *conIterStr = 0;

   ulong DELETEUL = 0;
   ulong lenAlnUL = 0;
   ulong depthUL = 0;

   ulong insSupUL = 0;
   ulong delSupUL = 0;
   ulong snpSupUL = 0;
   ulong supportUL = 0;
   ulong bestSupUL = 0;

   struct stichAmpST *stichIterST = 0;
   struct stichAmpST *altStichST = 0;
   struct stichAmpST *lastConBaseST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-02:
   ^  - Find the possible largest consensus length
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   stichIterST = conST;

   if(conST->lastBase == 0)
   { /*If: I am on the first base*/
      lenAlnUL = 0;

      while(stichIterST->nextBase != 0)
      { /*Loop: Get the sequence length*/
         stichIterST = stichIterST->nextBase;
         ++lenAlnUL;
      } /*Loop: Get the sequence length*/

      stichIterST = conST;
   } /*If: I am on the first base*/

   else if(conST->nextBase == 0)
   { /*Else if: I am on the last base*/
      lenAlnUL = 0;

      while(stichIterST->lastBase != 0)
      { /*Loop: Move to the first base*/
         stichIterST = stichIterST->lastBase;
         ++lenAlnUL;
      } /*Loop: Move to the first base*/
   } /*Else if: I am on the last base*/

   else
   { /*Else: I am not at a starting base*/
      /*Find the last base*/
      while(stichIterST->nextBase != 0)
         stichIterST = stichIterST->nextBase;

      /*Find the sequence length*/
      lenAlnUL = 0;

      while(stichIterST->lastBase != 0)
      { /*Loop: Get the sequence length*/
         stichIterST = stichIterST->lastBase;
         ++lenAlnUL;
      } /*Loop: Get the sequence length*/
   } /*Else: I am not at a starting base*/

   /*At this point stichIterST is on the first base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-03:
   ^  - Allocate memory for the consensus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   conSeqStr = malloc(sizeof(char) * (lenAlnUL + 2));
     /*+2 is to account for being off
     ` Frist + 1 is for '\0' at end
     ` Second + 1 is for lenAlnUL being index 0
     */
   conIterStr = conSeqStr;

   if(conSeqStr == 0) return 0;

   *(conSeqStr + lenAlnUL) = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-04:
   ^  - Merge all alternative bases into a single consensus
   ^  o fun-06 sec-04 sub-01:
   ^    - Find the amount of support for each error type
   ^  o fun-06 sec-04 sub-02:
   ^    - Handle cases were have values beneath users
   ^      min depth (needs 100% agreement)
   ^  o fun-06 sec-04 sub-03:
   ^    - Handle cases with values at or above the
   ^      users min read depth
   ^  o fun-06 sec-04 sub-04:
   ^    - I only have one base, go with that
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-06 Sec-04 Sub-01:
   *  - Find the amount of support for each error type
   \*******************************************************/

   while(stichIterST != 0)
   { /*While I have bases to convert to chars*/
      lastConBaseST = stichIterST;
      altStichST = stichIterST;
      stichIterST = stichIterST->nextBase;
      ++DELETEUL;
      depthUL = lastConBaseST->depthUL; /*# amplicons*/
      bestSupUL = 0;
      delSupUL = 0;
      insSupUL = 0;
      snpSupUL = 0;

      if(lastConBaseST->altBase)
      { /*If: I have alternative bases*/
         while(altStichST != 0)
         { /*Loop: Find the total support for position*/
            if(altStichST->baseC >= 64 + 32)
               insSupUL += altStichST->supportUL;

            else if(altStichST->baseC == '-')
               delSupUL += altStichST->supportUL;

            else snpSupUL += altStichST->supportUL;
               /*snp/match/mask*/

            altStichST = altStichST->altBase;
         } /*Loop: Find the total support for position*/

         /************************************************\
         * Fun-06 Sec-04 Sub-03:
         *  - Handle cases were have values beneath users
         *    min depth (needs 100% agreement)
         *  o fun-06 sec-04 sub-03 cat-01:
         *    - Handle insertion cases for low dephts
         *  o fun-06 sec-04 sub-03 cat-02:
         *    - Handle deltion cases for low dephts
         *  o fun-06 sec-04 sub-03 cat-03:
         *    - Handle snp/match or snp/match with
         *      deletions cases for low depths
         \************************************************/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-06 Sec-04 Sub-03 Cat-01:
         +  - Handle insertion cases for low depths
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         altStichST = lastConBaseST;
 
         if(depthUL < settings->minDepthUL)
         { /*If I am beneath the min depth limits*/
            if(altStichST->baseC >= 64 + 32)
            { /*If: I had an insertion*/
               if(insSupUL < depthUL) continue;

               else if(insSupUL > 0)
               { /*Else if: have 100% support for ins*/
                  /*Check if there are alterantive bases*/
                  if(altStichST->supportUL < depthUL)
                     *conIterStr = settings->maskC | 32;

                  else *conIterStr=altStichST->baseC | 32;
                  /*This requires 100% support to not
                  ` mask, so there are no alternatives
                  */

                  ++conIterStr;
                  continue;
               } /*Else if: have 100% support for ins*/
            } /*If: I had an insertion*/

            /*++++++++++++++++++++++++++++++++++++++++++++\
            + Fun-06 Sec-04 Sub-03 Cat-02:
            +  - Handle deletion cases for low dephts
            \++++++++++++++++++++++++++++++++++++++++++++*/

            if(delSupUL > snpSupUL && delSupUL >= depthUL)
               continue; /*Deletion is supported*/

            /*++++++++++++++++++++++++++++++++++++++++++++\
            + Fun-06 Sec-04 Sub-03 Cat-03:
            +  - Handle snp/match or snp/match with
            +    deletions cases for low depths
            \++++++++++++++++++++++++++++++++++++++++++++*/

            else
            { /*Else: there is an snp or snp/del combo*/
               depthUL -= delSupUL;
               *conIterStr = 0;

               while(altStichST != 0)
               { /*Loop: Find best base*/
                  if(altStichST->baseC != '-')
                  { /*If: I have an snp/match*/
                     if(altStichST->supportUL >= depthUL)
                        *conIterStr =
                           altStichST->baseC & (~32);
                     /*This requires 100% support*/
                  } /*If: I have an snp/match*/

                  altStichST = altStichST->altBase;
               } /*Loop: Find best base*/

               /*If no base had enough support*/
               if(*conIterStr == 0)
                  *conIterStr= settings->maskC & (~32);

               ++conIterStr;
               continue;
            } /*Else: there is an snp or snp/del combo*/
         } /*If I am beneath the min depth limits*/

         /************************************************\
         * Fun-06 Sec-04 Sub-04:
         *  - Handle cases with values at or above the
         *    users min read depth
         *  o fun-06 sec-04 sub-04 cat-01:
         *    - Handle insertion has to low support cases
         *  o fun-06 sec-04 sub-04 cat-02:
         *    - Handle insertion has enough support cases
         *  o fun-06 sec-04 sub-04 cat-03:
         *    - Handle deletion has enough support cases
         *  o fun-06 sec-04 sub-04 cat-04:
         *    - Handle snp/match cases
         \************************************************/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-06 Sec-04 Sub-04 Cat-01:
         +  - Handle insertion has to low support cases
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         supportUL = (insSupUL * 100) / depthUL;

         if(
               altStichST->baseC >= 64 + 32
            && supportUL < settings->minSupportUL
           ) continue; /*insertion not supported*/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-06 Sec-04 Sub-04 Cat-02:
         +  - Handle insertion has enough support cases
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         else if(altStichST->baseC >= 64 + 32)
         { /*Else if: This insertion is supported*/
            *conIterStr = 0;

            while(altStichST != 0)
            { /*Loop: Find the best position*/
               supportUL =
                    (altStichST->supportUL * 100)
                  / depthUL;

               if(supportUL >= settings->minSupportUL)
               { /*If: alt base is the supported base*/
                  /*Checking if last bases was better*/
                  if(supportUL > bestSupUL)
                  { /*If: this is the best base*/
                     *conIterStr = altStichST->baseC | 32;
                     bestSupUL = supportUL;
                  } /*If: this is the best base*/

                  else if(
                     *conIterStr == (settings->maskC | 32)
                  ) { /*Else if: The best base was a mask*/
                     *conIterStr = altStichST->baseC | 32;
                     bestSupUL = supportUL;
                  } /*Else if: The best base was a mask*/
               } /*If: alt base is the supported base*/

               altStichST = altStichST->altBase;
            } /*Loop: Find the best position*/

            /*No base had enough support to keep a base*/
            if(*conIterStr == 0)
               *conIterStr = settings->maskC | 32;

            ++conIterStr;
            continue;
         } /*Else if: This insertion is supported*/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-06 Sec-04 Sub-04 Cat-03:
         +  - Handle deletion has enough support cases
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         supportUL = (delSupUL * 100) / depthUL;
         if(
                delSupUL > snpSupUL
             && supportUL >= settings->minSupportUL
         ) continue; /*If is a deletion*/

         /*+++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun-06 Sec-04 Sub-04 Cat-04:
         +  - Handle snp/match cases (only case left)
         \+++++++++++++++++++++++++++++++++++++++++++++++*/

         *conIterStr = 0;

         while(altStichST != 0)
         { /*See if there is good support for position*/
            supportUL = 
              (altStichST->supportUL * 100)/depthUL;

            if(
                  altStichST->baseC != '-'
               && supportUL >= settings->minSupportUL
            ){ /*If: alt base is the supported base*/
               if(supportUL > bestSupUL)
               { /*If: this is the best base*/
                  *conIterStr = altStichST->baseC & (~32);
                  bestSupUL = supportUL;
               } /*If: this is the best base*/

               else if(
                  *conIterStr == (settings->maskC & ~32)
               ) { /*Else if: The best base was a mask*/
                  *conIterStr = altStichST->baseC & ~(32);
                  bestSupUL = supportUL;
               } /*Else if: The best base was a mask*/
            } /*If: alt base is the supported base*/

            altStichST = altStichST->altBase;
         } /*See if there is good support for position*/

         /*Else not support for base*/
         if(*conIterStr == 0)
            *conIterStr = settings->maskC & (~32);

         ++conIterStr;
         continue;
      } /*If: I have alternative bases*/

      /***************************************************\
      * Fun-06 Sec-04 Sub-05:
      *  - I only have one base, go with that
      \***************************************************/

      else
      { /*Else this is the only base*/
         if(altStichST->baseC >= 64 + 32)
         { /*If: this was an insertion*/
            *conIterStr = altStichST->baseC | 32;
            ++conIterStr;
         } /*If: this was an insertion*/

         else if(altStichST->baseC > 64)
         { /*If this base is an snp/match or mask*/
            *conIterStr = altStichST->baseC & (~32);
            ++conIterStr;
         } /*If this base is an snp/match or mask*/
 
         continue;
      } /*Else this is the only base*/
   } /*While I have bases to convert to chars*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-05:
   ^  - Clean up and exit
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*No clean up to do*/
   *conIterStr = '\0';
   return conSeqStr;
} /*stichAmpConToCStr*/
