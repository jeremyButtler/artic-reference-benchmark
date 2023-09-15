/*#########################################################
# Name: stichInputAndHelp
# Use:
#  - Holds functions to check user input & print out the
#    help message
# Libraries:
#  - "samFunSrc/cStrFun.h"
#  - "samFunSrc/cStrToNumberFun.h"
#  - "stichSetStruct.h"
#  o "samFunSrc/dataTypeShortHand.h"
# C Standard Libraries:
#  - <stdio.h>
#  - <string.h>
#  o <stdint.h>
#########################################################*/

#include "stichInputAndHelp.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' stichInputAndHelp SOF: Start Of Functions
'  - Checks user input and prints out help message
'  o fun-01 stichGetInput:
'    - Gets the user input (calls stichCheckArg)
'  o fun-02 stichCheckArg:
'    - Checks a single paramter/argument combination
'  o fun-03 pStichHelp:
'    - Prints the help message for stich
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: stichGetInput (Fun-01:)
| Use:
|  - Gets user input for stich
| Input:
|  - lenArgsI:
|    o length of argsStr(number arguments/parameters input)
|  - argsStr:
|    o Arguments and parameters the usesr input
|  - refFileStr:
|    o Buffer to hold the path to the fasta file with the
|      reference 
|  - ampFileStr:
|    o Buffer to hold the path to the fasta file with the
|      amplicons 
|  - outFileStr:
|    o Buffer to hold the path to the output file (fasta)
|  - prefixStr:
|    o Buffer to hold the prefix to name the scaffold
|  - settings:
|    o stichSet struct tho hold settings
|  - errUS:
|    o unsigned short to hold any errors
| Output:
|  - Returns:
|    o 0 for entry error or no error
|    o char pointer to unkown input if invalid input
|  - Modifies:
|    o refFileStr to have the refereces file path
|    o ampFileStr to have the amplicons file path
|    o outFileStr to have the output file path
|    o prefixStr to have the prefix
|    o settings: To have user input settings
|    o errUC has flags set
|      - 0: No flags set, no errors
|      - Flag 1: help message requested
|      - Flag 2: Version number requested
|      - Flag 4: -min-depth was non-numeric
|      - Flag 8: -min-support was non-numeric
|      - Flag 16: No numeric input for -threads/-t
|      - Flag 32 + 4: blank input for -min-depth
|      - Flag 32 + 8: blank input for -min-support
|      - Flag 32 + 16: blank input for -threads/-t
|      - Flag 64: Invalid input
\--------------------------------------------------------*/
char * stichGetInput(
   int lenArgsI,      /*Number of arguments input*/
   char *argsStr[],   /*Arguments/paramaters input*/
   char *refFileStr, /*Will hold the path to fasta file*/
   char *ampFileStr, /*Will hold the path to fasta file*/
   char *outFileStr, /*Will hold the output file*/
   char *prefixStr,  /*Holds prefix for consensus*/
   struct stichSet *settings,/*Settings specific to stich*/
   ushort *errUS       /*Error message*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   | Fun-01 TOC: stichGetInput
   |  - Get user input into stich
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   char *parmStr = 0;
   char *argStr = 0;
   char *tmpStr = 0;

   for(int iArg = 1; iArg < lenArgsI; ++iArg)
   { /*Loop: Search through user input*/
      parmStr = *(argsStr + iArg);
      argStr = *(argsStr + iArg + 1);

      tmpStr = 
         stichCheckArg(
            parmStr,
            argStr,
            &iArg,
            refFileStr,
            ampFileStr,
            outFileStr,
            prefixStr,
            settings,
            errUS
      );
   } /*Loop: Search through user input*/

   if(tmpStr == 0) return tmpStr;
   return 0;
} /*stichGetInput*/

/*--------------------------------------------------------\
| Name: stichCheckArg (Fun-02:)
| Use:
|  - Checks a single parameter/argument combination for
|    stich
| Input:
|  - parmStr:
|    o Points to the parameter to check
|  - argStr:
|    o Points to the argument to check
|  - countI:
|    o This is the current parameter on. This will be
|      incurmented by 1 if the paramter takes an agrugment
|  - lenArgsI:
|    o length of argsStr(number arguments/parameters input)
|  - argsStr:
|    o Arguments and parameters the usesr input
|  - refFileStr:
|    o Buffer to hold the path to the fasta file with the
|      reference 
|  - ampFileStr:
|    o Buffer to hold the path to the fasta file with the
|      amplicons 
|  - outFileStr:
|    o Buffer to hold the path to the output file (fasta)
|  - prefixStr:
|    o Buffer to hold the prefix to name the scaffold
|  - settings:
|    o stichSet struct tho hold settings
|  - errUS:
|    o unsigned short to hold any errors
| Output:
|  - Returns:
|    o 0 for entry error or no error
|    o char pointer to unkown input if invalid input
|  - Modifies:
|    o countI is incurmented by 1 if parameter has argument
|    o refFileStr to have the refereces file path
|    o ampFileStr to have the amplicons file path
|    o outFileStr to have the output file path
|    o prefixStr to have the prefix
|    o settings: To have user input settings
|    o errUC has flags set
|      - 0: No flags set, no errors
|      - Flag 1: help message requested
|      - Flag 2: Version number requested
|      - Flag 4: -min-depth was non-numeric
|      - Flag 8: -min-support was non-numeric
|      - Flag 16: No numeric input for -threads/-t
|      - Flag 32 + 4: blank input for -min-depth
|      - Flag 32 + 8: blank input for -min-support
|      - Flag 32 + 16: blank input for -threads/-t
|      - Flag 64: Invalid input
\--------------------------------------------------------*/
char * stichCheckArg(
   char *parmStr,    /*Parameter to check*/
   char *argStr,     /*Argument to check*/
   int *countI,      /*Counter to incurment*/
   char *refFileStr, /*Will hold the path to fasta file*/
   char *ampFileStr, /*Will hold the path to fasta file*/
   char *outFileStr, /*Will hold the output file*/
   char *prefixStr,  /*Holds prefix for consensus*/
   struct stichSet *settings,/*Settings specific to stich*/
   ushort *errUS       /*Error message*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   | Fun-02 TOC: stichCheckArg
   |  - Checks a single arg/paramter combination
   |  o fun-02 sec-01:
   |    - Variable declerations
   |  o fun-02 sec-02:
   |    - Check variables holding files or prefix
   |  o fun-02 sec-03:
   |    - Check variables holding alignment method
   |  o fun-02 sec-04:
   |    - Check variables holding scaffold colapse settings
   |  o fun-02 sec-05:
   |    - help message/version request or unkown input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *tmpStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-02:
   ^  - Check variables holding files or prefix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Using cStrCpInvsDelm to avoid openbsd warnings about
   ` strcpy is almost alwasy misused.
   */
   if(strcmp(parmStr, "-ref") == 0)
   { /*else if: The user input the reference file*/
      cStrCpInvsDelm(refFileStr, argStr);
      ++(*countI);
   } /*else if: The user input the reference file*/

   else if(strcmp(parmStr, "-amps") == 0)
   { /*else if: The user input the amplicon file*/
      cStrCpInvsDelm(ampFileStr, argStr);
      ++(*countI);
   } /*else if: The user input the amplicon file*/

   else if(
         strcmp(parmStr, "-o") == 0
      || strcmp(parmStr, "-out") == 0
   ){ /*Else If: The user specified an output file*/
      cStrCpInvsDelm(outFileStr, argStr);
      ++(*countI);
   } /*Else If: The user specified an output file*/

   else if(strcmp(parmStr, "-overwrite") == 0)
      settings->overwriteBl = 1;

   else if(strcmp(parmStr, "-no-overwrite") == 0)
      settings->overwriteBl = 0;

   else if(
         strcmp(parmStr, "-p") == 0
      || strcmp(parmStr, "-prefix") == 0
   ){ /*Else if: this is the prefix for the consensus*/
     cStrCpInvsDelm(prefixStr, argStr);
     ++(*countI);
   } /*Else if: this is the prefix for the consensus*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-03:
   ^  - Check variables holding alignment method
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   else if(strcmp(parmStr, "-minimap") == 0)
      settings->useMinimap2Bl = 1;

   /*TODO: This is not used. Add in non-minimap2 support*/
   else if(strcmp(parmStr, "-no-minimap") == 0)
      settings->useMinimap2Bl = 0;

   else if(
         strcmp(parmStr, "-threads") == 0
      || strcmp(parmStr, "-t") == 0
   ){ /*Else if: changing number of threads*/
       tmpStr=cStrToUChar(argStr, &settings->threadsUC);

       if(*tmpStr > 32) *errUS |= 16;
       if(*tmpStr == *argStr) *errUS |= (32 + 16);

       ++(*countI);
   } /*Else if: changing number of threads*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-04:
   ^  - Check variables holding scaffold colapse settings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   else if(strcmp(parmStr, "-min-depth") == 0)
   { /*Else if: user wanted to change voting depth*/
       tmpStr =
          strToULBase10(argStr, &settings->minDepthUL);

       if(*tmpStr > 32) *errUS |= 4;
       if(*tmpStr == *argStr) *errUS |= (32 + 4);

       ++(*countI);
   } /*Else if: user wanted to change voting depth*/

   else if(strcmp(parmStr, "-min-support") == 0)
   { /*Else if: user wanted to change min vote*/
       tmpStr =
          strToULBase10(argStr,&settings->minSupportUL);

       if(*tmpStr > 32) *errUS |= 8;
       if(*tmpStr == *argStr) *errUS |= (32 + 8);

       ++(*countI);
   } /*Else if: user wanted to change min vote*/

   else if(strcmp(parmStr, "-mask") == 0)
   { /*Else if: user wanted to change the mask char*/
       settings->maskC = *argStr;
       ++(*countI);
   } /*Else if: user wanted to change the mask char*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-05:
   ^  - help message/version request or unkown input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   else if(
         strcmp(parmStr, "-h") == 0
      || strcmp(parmStr, "--h") == 0
      || strcmp(parmStr, "-help") == 0
      || strcmp(parmStr, "--help") == 0
      || strcmp(parmStr, "help") == 0
   ) *errUS |= 1;

   else if(
         strcmp(parmStr, "-v") == 0
      || strcmp(parmStr, "-V") == 0
      || strcmp(parmStr, "--v") == 0
      || strcmp(parmStr, "--V") == 0
      || strcmp(parmStr, "-version") == 0
      || strcmp(parmStr, "-Version") == 0
      || strcmp(parmStr, "--version") == 0
      || strcmp(parmStr, "--Version") == 0
      || strcmp(parmStr, "version") == 0
      || strcmp(parmStr, "Version") == 0
   ) *errUS |= 2;

   else
   { /*Else: non-valid input*/
      *errUS |= 64;
      return argStr;
   } /*Else: non-valid input*/

   return 0;
} /*stichCheckArg*/

/*--------------------------------------------------------\
| Name: pStichHelp (Fun-03:)
| Use:
|  - Prints the help message for stich
| Input:
|  - outFILE:
|    o file to print help message to
| Output:
|  - Prints the help message to outFILE
\--------------------------------------------------------*/
void pStichHelp(
   FILE *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: pStichHelp
   '  - Prints the help message for stich
   '  o fun-03 sec-01:
   '    - Out header
   '  o fun-03 sec-02:
   '    - variables holding files or prefix
   '  o fun-03 sec-03:
   '    - alignment method
   '  o fun-03 sec-04:
   '    - scaffold colapse settings
   '  o fun-03 sec-05:
   '    - help message/version request or unkown input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Out header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      outFILE,
      "stich -ref ref.fasta -amps amplicons.fasta"
   );
   fprintf(outFILE, " [options...]\n");

   fprintf(outFILE, "Use:\n");
   fprintf(
      outFILE,
      "   - Uses a reference to stich amplicons into an"
   );
   fprintf(outFILE, " scaffold\n");

   fprintf(outFILE, "Input:\n");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - variables holding files or prefix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(outFILE, "   -ref: [Required]\n");
   fprintf(
    outFILE,
    "     o Reference (as fasta) to align amplicons with\n"
   );

   fprintf(outFILE, "   -amps: [Required]\n");
   fprintf(
    outFILE,
    "     o Amplicons (as fasta) to build scaffold with\n"
   );

   fprintf(outFILE, "   -out: [stdout]\n");
   fprintf(outFILE, "     o File to write consnesus to\n");

   if(defStichOverwrite == 1)
      fprintf(outFILE, "   -overwrite: [Yes]\n");
   else
      fprintf(outFILE, "   -overwrite: [No]\n");
   fprintf(
      outFILE,
      "     o Write over output file if already exists\n"
   );
   fprintf(
      outFILE,
      "     o This can be disabled with -no-overwrite\n"
   );

   fprintf(outFILE, "   -prefix: [%s]\n", defPrefix);
   fprintf(outFILE, "     o Prefix to name consensus\n");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - alignment method
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /* TODO: Allow this option
   if(defMinimap2Bl)
      fprintf(outFILE, "   -minimap: [Yes]\n");
   else
      fprintf(outFILE, "   -minimap: [No]\n");
   fprintf(outFILE, "     o Use minimap2 for alignment\n");
   fprintf(
      outFILE,
      "     o This is changed to alnSeq with -no-minimap\n"
   ); */

   fprintf(outFILE, "   -threads: [%i]\n", defThreads);
   fprintf(outFILE, "     o Number of threads to use\n");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-04:
   ^  - scaffold colapse settings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      outFILE,
      "   -min-depth: [%i]\n",
      defMinStichDepth
   );
   fprintf(
     outFILE,
     "     o Minimum depth needed to vote on"
   );
   fprintf(outFILE,  " disagreements\n");

   fprintf(
      outFILE,
      "   -min-support: [%i]\n",
      defMinStichSup
   );
   fprintf(
     outFILE,
     "     o Minimum support (as %%) needed to decide"
   );
   fprintf(outFILE, " which\n");
   fprintf(
      outFILE, 
      "       disagreement should be kept/not masked.\n");

   fprintf(
    outFILE,
   "     o Voting is done when total support >= -min-depth"
   );
   fprintf(outFILE, "\n");


   fprintf(outFILE, "   -mask: [%c]\n", defMaskBase);
   fprintf(
     outFILE,
     "     o Mask to apply to unsupported disagreements\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-05:
   ^  - help message/version request or unkown input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(outFILE, "   -h:\n");
   fprintf(outFILE, "     o Print this help message\n");

   fprintf(outFILE, "   -v:\n");
   fprintf(outFILE, "     o Print the version number\n");

   return;
} /*pStichHelp*/
