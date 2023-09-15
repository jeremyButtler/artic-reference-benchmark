/*#########################################################
# Name: stichFun
# Use:
#  - Holds functions used in stiching amplicons into a
#    consensus.
# Libraries:
#  - "alnSeqSrc/memWater.h"
#  - "alnSeqSrc/hirschberg.h"
#  - "sitchAmpStruct.h"            (No .c file)
#  - "stichDefaults.h"             (No .c file)
#  o "alnSeqSrc/generalAlnFun.h"
#  o "alnSeqSrc/alnStruct.h"
#  o "alnSeqSrc/alnMatrixStruct.h"
#  o "alnSeqSrc/twoBitArrays.h"    (No .c file)
#  o "alnSeqSrc/scoresST.h"        (No .c file)
#  o "alnSeqSrc/seqStruct.h"
#  o "alnSeqSrc/alnSetStruct.h"
#  o "alnSeqSrc/alnSeqDefaults.h"
#  o "dataTypeShortHand.h"         (No .c file)
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdint.h>
#  o <stdio.h>
#  o <string.h>
#########################################################*/

#ifndef STICHFUN_H
#define STICHFUN_H

#include "alnSeqSrc/memWater.h"
#include "alnSeqSrc/hirschberg.h"
#include "stichAmpStruct.h"
#include "stichDefaults.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' stichFun SOH: Start Of Header
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
| Struct-01: stichSet
|  - Holds settings for stiching together amplicons
\--------------------------------------------------------*/
typedef struct stichSet{
   char maskC;       /*Mask to apply to low support bases*/
   char useMinimap2Bl; /*use minimap2*/
   ulong minDepthUL; /*Min depth to do voting (not 100%)*/
   ulong minSupportUL; /*Min support to keep disagreement*/
}stichSet;

/*--------------------------------------------------------\
| Macro-01: initStichSet
| Use:
|  - Initializes a stichSet structure with default settings
| Input:
|  - stichSetPtr
|    o Pointer to a stichSet structuer to initialize
| Output:
|  - Modifies:
|    o All values in stichSetPtr to be default values
|      (see stichDefaults.h for defaults)
\--------------------------------------------------------*/
#define initStichSet(stichSetPtr){\
   (stichSetPtr)->maskC = defMaskBase; \
   (stichSetPtr)->minDepthUL = defMinStichDepth; \
   (stichSetPtr)->minSupportUL = defMinStichSup; \
   (stichSetPtr)->useMinimap2Bl = defMinimap2Bl; \
} /*initStichSet*/

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
|  - settings:
|    o alnSet struct with the settings to use for aligment
| Output:
|   - Modifies:
|     o numAmpsUL to hold the number of amplicons
|   - Returns:
|     o Pointer to array of scoresStructs with the
|       coordinates of each alignment
|     o 0 for memory errors
\--------------------------------------------------------*/
struct scoresStruct * getAmpPos(
   struct seqStruct *refST,  /*Reference sequence*/
   struct seqStruct *qryST,  /*Blank struct to work with*/
   FILE *ampFaFILE,          /*File with amplicons*/
   unsigned long *numAmpsUL, /*Will have number amplicons*/
   ulong **indexAryUL, /*Will hold file index for scores*/
   struct alnSet *settings   /*Alignment settings*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01: getAmpPos
   '  - Gets the position of each amplicon on the reference
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
   ulong *seqIndexAryUL,      /*Index for every seq*/
   struct seqStruct *refST,    /*Reference sequence*/
   struct seqStruct *ampST,    /*Amplicon sequence*/
   struct alnSet *alnSetST,     /*Settings for alnSeq*/ 
   struct stichSet *stichSetST  /*Settings for stich*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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

#endif
