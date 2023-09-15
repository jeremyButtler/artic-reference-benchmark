/*#########################################################
# Name: stich
# Use:
#  - stiching amplicons into a consensus. This program
#    use minimap2
# Libraries:
#  - "stichMinimapFun.h"
#  - "stichInputAndHelp.h"           (No .c file)
#  o "stichSetStruct.h"              (No .c file)
#  o "stichDefaults.h"               (No .c file)
#  o "sitchAmpStruct.h"              (No .c file)
#  o "samFunSrc/trimSam.h"
#  o "samFunSrc/samEntryStruct.h"    (No .c file)
#  o "samFunSrc/cStrToNumberFun.h"   (No .c file)
#  o "samFunSrc/dataTypeShortHand.h" (No .c file)
#  o "samFunSrc/seqStruct.h"
#  o "cStrFun.h"                     (No .c file)
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdint.h>
#  o <stdio.h>
#  o <string.h>
# Requires:
#  - Minimap2 be in your file path
#########################################################*/

#include "stichMinimapFun.h"
#include "stichInputAndHelp.h"

int main(
   int lenArgsI,    /*Number of arguments input*/
   char *argsStr[]  /*Arguments/paramaters input*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '  - Runs stich
   '  o main sec-01:
   '    - Variable declerations
   '  o main sec-02:
   '    - Get and check user input
   '  o main sec-03:
   '    - Check if input files are valid
   '  o main sec-04:
   '    - Get the positions for each amplicon
   '  o main sec-05:
   '    - Stich the amplicons together
   '  o main sec-06:
   '    - Collapse the scaffold list into a single scaffold
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char refStr[256];
   char ampStr[256];
   char outStr[256];

   char prefixStr[256];

   ushort inputErrUS = 0;
   uchar errUC = 0;
   char *errStr = 0;

   char *conSeqStr = 0;
   ulong numAmpsUL = 0;

   struct samEntry *ampsAryST = 0;
   struct stichSet stichSetST;
   struct stichAmpST *conST = 0;

   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-02:
   ^  - Get and check user input
   ^  o main sec-02 sub-01:
   ^    - Get user input and set defaults
   ^  o main sec-02 sub-02:
   ^    - Check input (1st check if wanted help/version)
   ^  o main sec-02 sub-02:
   ^    - Error happened, check what kind of error
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   ^ Main Sec-02 Sub-01:
   ^  - Get user input and set defaults
   \******************************************************/

   initStichSet(&stichSetST);
   strcpy(prefixStr, defPrefix);
   refStr[0] = '\0';
   ampStr[0] = '\0';
   outStr[0] = '\0';

   errStr =
      stichGetInput(
         lenArgsI,
         argsStr,
         refStr,
         ampStr,
         outStr,
         prefixStr,
         &stichSetST,
         &inputErrUS
   ); /*Get the user input*/

   /******************************************************\
   ^ Main Sec-02 Sub-02:
   ^  - Check input (1st check if user wanted help/version)
   \******************************************************/

   if(inputErrUS)
   { /*If: there was a probelm*/
      if(inputErrUS & 1)
      { /*If: the help message was requested*/
         pStichHelp(stdout);

         if(inputErrUS & 2)
            fprintf(stdout, "%i\n", defVersion);
         exit(0); 
      } /*If: the help message was requested*/

      if(inputErrUS & 2)
      { /*If: the version number was requested*/
         fprintf(stdout, "%i\n", defVersion);
         exit(0); 
      } /*If: the version number was requested*/

      /***************************************************\
      ^ Main Sec-02 Sub-03:
      ^  - Error happened, check what kind of error
      \***************************************************/

      pStichHelp(stderr); /*Printing help message*/

      if(inputErrUS & 4)
      { /*If: Invalid input for min depth*/
         if(inputErrUS & 32)
           fprintf(
              stderr,
              "-min-depth; no argument was provided\n"
           );
         else
           fprintf(
             stderr,
             "Input for -min-depth is non-numeric\n"
           );
      } /*If: Invalid input for min depth*/

      if(inputErrUS & 8)
      { /*If: Invalid input for min support*/
         if(inputErrUS & 32)
           fprintf(
              stderr,
              "-min-support; no argument was provided\n"
           );
         else
           fprintf(
             stderr,
             "Input for -min-support is non-numeric\n"
           );
      } /*If: Invalid input for min support*/

      if(inputErrUS & 8)
      { /*If: Invalid input for number of threads*/
         if(inputErrUS & 32)
           fprintf(
              stderr,
              "-threads or -t; no argument was provided\n"
           );
         else
           fprintf(
             stderr,
             "Input for -threads or -t is non-numeric\n"
           );
      } /*If: Invalid input for number of threads*/

      if(inputErrUS & 64)
      { /*If: invalid input was input*/
         fprintf(
            stderr,
            "%s is not an valid parameter\n",
            errStr
         );
      } /*If: invalid input was input*/

      exit(-1);
   } /*If: there was a probelm*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-03:
   ^  - Check if input files are valid
   ^  o main sec-03 sub-01:
   ^    - Check the reference sequence
   ^  o main sec-03 sub-01:
   ^    - Check the amplicon sequences
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Main Sec-03 Sub-01:
   *  - Check the reference sequence
   \******************************************************/

   if(refStr[0] == '\0')
   { /*If: no reference sequence was provided*/
      fprintf(
         stderr,
         "No reference sequence was provided (use -ref)\n"
      );
      exit(-1);
   } /*If: no reference sequence was provided*/

   outFILE = fopen(refStr, "r");

   if(outFILE == 0)
   { /*If: reference file could not be opened*/
      fprintf(
         stderr,
         "Could not open -ref %s\n",
         refStr
      );
      exit(-1);
   } /*If: reference file could not be opened*/

   fclose(outFILE);
   outFILE = 0;

   /******************************************************\
   * Main Sec-03 Sub-02:
   *  - Check the amplicon sequences
   \******************************************************/

   if(ampStr[0] == '\0')
   { /*If: no amplicon sequences were provided*/
      fprintf(
        stderr,
        "No amplicon sequences were provided (use -amps)\n"
      );
      exit(-1);
   } /*If: no amplicon sequences were provided*/

   outFILE = fopen(ampStr, "r");

   if(outFILE == 0)
   { /*If: amplicons file could not be opened*/
      fprintf(
         stderr,
         "Could not open -amps %s\n",
         ampStr
      );
      exit(-1);
   } /*If: amplicons file could not be opened*/

   fclose(outFILE);
   outFILE = 0;

   /******************************************************\
   * Main Sec-03 Sub-03:
   *  - Check the output file
   \******************************************************/

   if(outStr[0] != '\0')
   { /*If: the user provided an output file*/
      outFILE = fopen(outStr, "r");

      if(outFILE != 0)
      { /*If: the output file is invalid*/
         if(!stichSetST.overwriteBl)
         { /*If: I am not overwriting files*/
            fprintf(
               stderr,
               "File -out %s already exits\n",
               outStr
            );
            fprintf(
               stderr,
               "Printing consensus to stdout\n"
            );

            outStr[0] = '\0';
         } /*If: I am not overwriting files*/

         fclose(outFILE);
         outFILE = 0;
         goto stichOutFileOk;
      } /*If: the output file is invalid*/

      outFILE = fopen(outStr, "w");

      if(outFILE == 0)
      { /*If: the output file is invalid*/
         fprintf(
            stderr,
            "Unable to create -out %s\n",
            outStr
         );
         fprintf(
            stderr,
            "Printing consensus to stdout\n"
         );

         outStr[0] = '\0';
      } /*If: the output file is invalid*/

      else
      { /*Else: This is a valid output file*/
         fclose(outFILE);
         outFILE = 0;
      } /*Else: This is a valid output file*/
   } /*If: the user provided an output file*/

   stichOutFileOk:

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-04:
   ^  - Get the positions for each amplicon
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ampsAryST =
     getAmpPosMinimap(
        refStr,
        ampStr,
        &numAmpsUL,
        &stichSetST,
        &errUC
   );/*Find starting & ending positions for each sequence*/
   /*getAmpPos will also sort the postions*/

   if(ampsAryST == 0)
   { /*If I had a memory error*/
      switch(errUC)
      { /*Switch: Find out what happened*/
         case 0: break; /*Never fires, nothing went wrong*/
         case 1:
            fprintf(
               stderr,
               "-amp-fa %s could not be opend\n",
               ampStr
            ); /*Invalid amplicon file*/
            exit(-1);
         case 2:
            fprintf(
               stderr,
               "-ref-fa %s could not be opend\n",
               refStr
            ); /*Invalid reference file*/
            exit(-1);
         case 3:
            fprintf(
               stderr,
               "-ref-fa %s has more than one sequence\n",
               refStr
            ); /*Invalid reference file*/
            exit(-1);
         case 4:
            fprintf(
               stderr,
               "Memory error (likely ran out of memory)\n"
            ); /*Invalid reference file*/
            exit(-1);
         case 5:
            fprintf(
               stderr,
               "Minimap2 errored out or memory error\n"
            ); /*Invalid reference file*/
            exit(-1);
         case 6:
            fprintf(
               stderr,
               "No amplicons mapped to the reference\n"
            ); /*Invalid reference file*/
            exit(-1);
         default: break;
      } /*Switch: Find out what happened*/

      fprintf(
         stderr,
         "Something happed in amplicon alignment\n"
      );
      exit(-1);
   } /*If I had a memor error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-05:
   ^  - Stich the amplicons together
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Build the consensus*/
   conST = 
      stichAmpConMinimap(ampsAryST,numAmpsUL,&stichSetST);

   freeSamEntryAry(&ampsAryST, numAmpsUL);
   ampsAryST = 0;

   if(conST == 0)
   { /*If: there was an error*/
      fprintf(
       stderr,
       "Error in stiching step, likely ran out of memory\n"
      ); /*Invalid reference file*/
      exit(-1);
   } /*If: there was an error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-06:
   ^  - Collapse the scaffold list into a single scaffold
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   conSeqStr = stichAmpConToCStr(conST, &stichSetST);
   freeStichAmpSTList(&conST);

   if(conSeqStr == 0)
   { /*If: there was an error*/
      fprintf(
       stderr,
       "Ran out of memory when collapsing the scaffold\n"
      ); /*Invalid reference file*/
      exit(-1);
   } /*If: there was an error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-07:
   ^  - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(outStr[0] == '\0') outFILE = stdout;
   else outFILE = fopen(outStr, "w");

   if(outFILE == 0)
   { /*If: I could not open the outfile*/
      fprintf(
         stderr,
         "Could not open output file (-out %s)\n",
         outStr
      );
      fprintf(
         stderr,
         "Printing consensus to stdout\n"
      );
      outFILE = stdout;
   } /*If: I could not open the outfile*/
   
   fprintf(outFILE, ">%s\n%s\n", prefixStr, conSeqStr);
   free(conSeqStr);
   if(outFILE != stdout) fclose(outFILE);

   exit(0);
} /*main*/
