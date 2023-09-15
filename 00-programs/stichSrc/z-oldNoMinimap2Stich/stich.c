/*#########################################################
# Libraries:
#  - stichFun.h
#########################################################*/

#include "stichFun.h"

int main(
){ /*main*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *refStr = "../../ref.fa";
   char *ampStr = "../../tmp.fasta";
   char *conStr = 0;
   char *prefixStr = "test-consensus";
   ulong numAmpsUL = 0;
   ulong *seqIndexAryUL = 0;

   FILE *refFILE = fopen(refStr, "r");
   FILE *ampFILE = fopen(ampStr, "r");

   struct seqStruct refST;
   struct seqStruct qryST;

   struct scoresStruct *ampsAryST = 0;
   struct alnSet alnSetST;
   struct stichSet stichSetST;

   initAlnSet(&alnSetST);
   initStichSet(&stichSetST);
   initSeqST(&refST);
   initSeqST(&qryST);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-0?:
   ^  - Get the positions for each amplicon
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   readFaSeq(refFILE, &refST);
   fclose(refFILE);
   refFILE = 0;
   refST.endAlnUL = refST.lenSeqUL - 1;

   ampsAryST =
     getAmpPos(
        &refST,
        &qryST,
        ampFILE,
        &numAmpsUL,
        &seqIndexAryUL,
        &alnSetST
   );/*Find starting & ending positions for each sequence*/
   /*getAmpPos will also sort the postions*/

   if(ampsAryST == 0)
   { /*If I had a memory error*/
      fclose(ampFILE);
      freeSeqST(&refST, 0); /*0 For on stack*/
      freeSeqST(&qryST, 0); /*0 For on stack*/
      fprintf(stderr, "Failed to align all amplicons\n");
      exit(-1);
   } /*If I had a memor error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-0?:
   ^  - Stich the amplicons together
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   conStr = 
      stichAmpCon(
         ampFILE,
         ampsAryST,
         seqIndexAryUL,
         &refST,
         &qryST,
         &alnSetST,
         &stichSetST
   ); /*Build the consensus*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-0?:
   ^  - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fclose(ampFILE);
   freeSeqST(&refST, 0); /*0 For on stack*/
   freeSeqST(&qryST, 0); /*0 For on stack*/
   free(seqIndexAryUL);
   freeScoresSTAry(ampsAryST, numAmpsUL, 0);/*0 = on heap*/

   seqIndexAryUL = 0;
   ampsAryST = 0;

   if(conStr != 0)
   { /*If: I bulit a consensus*/
      fprintf(stdout, ">%s\n%s\n", prefixStr, conStr);
      free(conStr);
   } /*If: I bulit a consensus*/

   else
   { /*Else: I failed to build a consensus*/
      fprintf(
         stderr,
         "Failed to stich amplicons together\n"
      );
      exit(-1);
   } /*Else: I failed to build a consensus*/

   exit(0);
} /*main*/
