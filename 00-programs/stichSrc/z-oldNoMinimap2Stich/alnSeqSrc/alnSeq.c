/*#########################################################
# Name: alignSeq
# Use:
#  - Runs a Needlman Wunsch or Smith Waterman alignment on
#    a pair of fasta files
# Includes:
#  - "hirschberg.h"
#  - "needleman.h"
#  - "memWater.h"
#  o "waterman.h"
#  o "generalAlnFun.h"
#  o "alnStruct.h"
#  o "alnMatrixStruct.h"
#  o "twoBitArrays.h"
#  o "scoresST.h"
#  o "seqStruct.h"
#  o "alnSetStruct.h"
#  o "alnSeqDefaults.h"
#  o "twoBitArrays.h"
# C standard libraries:
#  o <string.h>
#  o <stdlib.h>
#  o <stdio.h>
#  o <stdint.h>
#########################################################*/

#include "hirschberg.h"
#include "needleman.h"
#include "memWater.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOP: Start Of Program
'  - main: Run this program
'  - fun-01 checkInput;
'    o Gets user input
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output: Modifies: Each input variable to hold user input
|  - Returns:
|     o 0: If no errors
|     o Pointer to parameter that was an issue
|     o Pointer to "-score-matrix" for invalid scoring
|       matrix input
|  - Prints to stdout when the scoring file is invalid
\--------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,        // Number arguments user input
    char *argsCStr[],       // Array with user arguments
    char **refFileCStr,     // file name of reference file
    char **queryFileCStr,   // File name of the query file
    char **outFileCStr,     // Name of the output file
    char **scoreMtrxFileStr, /*Holds scoring matrix file*/
    struct alnSet *settings // Aligment settings
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: checkInput
   '  - Checks & extracts user input
   '    fun-01 sec-1: Variable declerations
   '    fun-01 sec-2: Look through user input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Prints:
|    o help message to outFILE.
\--------------------------------------------------------*/
void printHelpMesg(
   FILE *outFILE,
   char breifBl   /*Print a shorter help message*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printHelpMesg
   '  - Prints out the help message to outFILE
   '  o fun-03 sec-01:
   '    - Usage block
   '  o fun-03 sec-02:
   '    - Input block
   '  o fun-03 sec-03:
   '    - Output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: printCompilerSettings
| Call: printCompilerSettings(file, 1/0)
| Use:
|   - Prints out the compiler flags used
|   - Prints out what each possible compiler flag does
| Input:
|   - outFILE:
|     o File to print flags and flag descriptions to
|   - pDescBl:
|     o 1: print out descriptions for each flag (all)
|     o 0: Do not print out any flag descriptions
| Output:
|   - Prints flags and description to outFILE
\--------------------------------------------------------*/
void printCompileSettings(
   FILE *outFILE, /*Output file*/
   char pDescBl   /*not 0: print out flag descriptions*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC:
   '  - Print out compile flags used and what each compiler
   '    flag does.
   '  o fun-04 sec-01:
   '    - Print out the compiled settings
   '  o fun-04 sec-02:
   '    - Print out what each compiler flag does
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


int main(
    int lenArgsInt,
    char *argsCStr[]
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '  - Main function that drives the Smith Waterman or
   '    Needlman Wunsch alignment
   '  o main sec-01:
   '    - Variable declerations
   '  o main sec-02:
   '    - Read in user input and check input
   '  o main sec-03:
   '    - read in the reference sequence
   '  o main sec-04:
   '    - read in the query sequence
   '  o main sec-05:
   '    - Do the alingment
   '  o main sec-06:
   '    - Print out multi-alignments or find the single
   '      alignment array
   '  o main sec-07:
   '    - Print out the alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // User input
   char *refFileCStr = 0;
   char *queryFileCStr = 0;
   char *outFileCStr = 0;
   char *inputCStr = 0;
   char *scoreMtrxFileStr = 0;

   // Holds the reference sequence
   struct seqStruct refST;
   struct seqStruct queryST;

   // Caputures error type from functions
   unsigned char errUC = 0;
   long bestScoreL = 0;  // Allows me to free score matrix

   // For holding settings
   struct alnSet settings;

   // For holding alignment output
   struct alnMatrixStruct *alnMtrxST = 0;
   struct alnStruct *alnST = 0;
   struct scoresStruct *bestScoreST = 0; /*For memWater*/

   FILE *faFILE = 0;
   FILE *outFILE = 0; /*Print alternative aligments to*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-2:
   ^  - Read in user input and check input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   initAlnSet(&settings);

   inputCStr =
       checkInput(
           &lenArgsInt,
           argsCStr,
           &refFileCStr,
           &queryFileCStr,
           &outFileCStr,
           &scoreMtrxFileStr,
           &settings
   ); // Get the user input

   if(inputCStr != 0)
   { // If had problematic input
        if(strcmp(inputCStr, "-h") == 0 ||
           strcmp(inputCStr, "--h") == 0 ||
           strcmp(inputCStr, "-help") == 0 ||
           strcmp(inputCStr, "--help") == 0 ||
           strcmp(inputCStr, "help") == 0
        ) { /*If user wanted the help message*/
            printHelpMesg(stdout, 1);/*Breif help message*/
            exit(0);
        } /*If user wanted the help message*/

        if(strcmp(inputCStr, "-h-all") == 0 ||
           strcmp(inputCStr, "--h-all") == 0 ||
           strcmp(inputCStr, "-help-all") == 0 ||
           strcmp(inputCStr, "--help-all") == 0 ||
           strcmp(inputCStr, "help-all") == 0
        ) { /*If user wanted the help message*/
            printHelpMesg(stdout, 0); /*full help message*/
            exit(0);
        } /*If user wanted the help message*/


        if(strcmp(inputCStr, "-V") == 0 ||
           strcmp(inputCStr, "-v") == 0 ||
           strcmp(inputCStr, "--V") == 0 ||
           strcmp(inputCStr, "--v") == 0 ||
           strcmp(inputCStr, "--version") == 0 ||
           strcmp(inputCStr, "--Version") == 0 ||
           strcmp(inputCStr, "-version") == 0 ||
           strcmp(inputCStr, "-Version") == 0 ||
           strcmp(inputCStr, "Version") == 0 ||
           strcmp(inputCStr, "version") == 0
        ) { /*if the user wanted the version number*/
            fprintf(
                stdout,
                "alnSeq version: %u\n",
                defVersion
            ); /*Print out the closest thing to a version*/
            exit(0);
        } /*Else if the user wanted the version number*/

        else if(strcmp(inputCStr, "-flags") == 0)
        { /*If user wanted to print the compiler flags*/
           printCompileSettings(stdout, 1);
           exit(0);
        } /*If user wanted to print the compiler flags*/

        /*Just the flags (no descriptions)*/
        else if(strcmp(inputCStr, "-flags-only") == 0)
        { /*If user wanted to print the compiler flags*/
            fprintf(
                stdout,
                "alnSeq version: %u\n",
                defVersion
            ); /*user will likely want the version number*/

           printCompileSettings(stdout, 0);
           exit(0);
        } /*If user wanted to print the compiler flags*/

        /*If file with scoring matrix was invalid*/
        else if(strcmp(inputCStr, "-score-matrix") == 0)
            exit(1);

        else if(inputCStr != 0)
        { /*If user had invalid input*/
            printHelpMesg(stderr, 1); /*short help*/
            fprintf(stderr, "%s is invalid\n", inputCStr);
            exit(1); /*Let user know their was an error*/
        } /*If user had invalid input*/
   } // If had problematic input

   if(outFileCStr != 0)
   { /*If printing output to a file*/
        outFILE = fopen(outFileCStr, "w");

        if(outFILE == 0)
        { /*If an invalid output file*/
            printf(
              "Output (-out %s) file is invalid.\n",
              outFileCStr
            ); /*Let user know about the invalid file*/

            exit(-1);
        } /*If an invalid output file*/

        fclose(outFILE);
        outFILE = 0;
   } /*If printing output to a file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-3:
   ^  - read in the reference sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   faFILE = fopen(refFileCStr, "r");

   if(faFILE == 0) 
   { // If reference file could not be opened
       fprintf(
         stderr,
         "Reference (-ref %s) could not be opend\n",
         refFileCStr
       );

       exit(-1);
   } // If reference file could not be opened

   // Read in the reference sequence
   initSeqST(&refST);
   errUC = readFaSeq(faFILE, &refST);
   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { // Invalid fasta file
       freeSeqST(&refST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Reference (-ref %s) is not valid\n",
         refFileCStr
       );

       exit(-1);
   } // Invalid fasta file

   if(errUC & 64)
   { // Invalid fasta file
       freeSeqST(&refST, 0); // 0 to makr on the stack
       fprintf(stderr, "Memory allocation error\n");
       exit(-1);
   } // Invalid fasta file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-04:
   ^  - read in the query sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   faFILE = fopen(queryFileCStr, "r");

   if(faFILE == 0) 
   { // If reference file could not be opened
       freeSeqST(&refST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Query (-query %s) could not be opend\n",
         queryFileCStr
       );

       exit(-1);
   } // If reference file could not be opened


   // Read in the query sequence
   initSeqST(&queryST);
   errUC = readFaSeq(faFILE, &queryST);
   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { // Invalid fasta file
       freeSeqST(&refST, 0); // 0 to makr on the stack
       freeSeqST(&queryST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Query (-query %s) is not valid\n",
         refFileCStr
       );
       exit(-1);
   } // Invalid fasta file

   if(errUC & 64)
   { // Invalid fasta file
      freeSeqST(&refST, 0); // 0 to makr on the stack
      freeSeqST(&queryST, 0); // 0 to makr on the stack
      fprintf(stderr, "Memory allocation error\n");
      exit(-1);
   } // Invalid fasta file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-05:
   ^  - Do the alingment
   ^  o main sec-05 sub-01:
   ^    - Set up for the alignment
   ^  o main sec-05 sub-02:
   ^    - Check if doing a Needleman alignment
   ^  o main sec-05 sub-02:
   ^    - Check if doing an Waterman alignment
   ^  o main sec-05 sub-04:
   ^    - Check if doing an memory efficent Waterman
   ^      alignment that prints out alternative alignment
   ^      positions (Gets combined with a Hirschberg)
   ^  o main sec-05 sub-05:
   ^    - Check if doing an memory efficent Waterman
   ^      alignment that prints out a single alignment
   ^      (uses a Hirschberg to get the alignment)
   ^  o main sec-05 sub-06:
   ^    - Check if doing an Hirschberg alignment
   ^  o main sec-05 sub-07:
   ^    - Let user know this was an invalid alignment
   ^  o main sec-05 sub-08:
   ^    - Check if the alignment failed
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Main Sec-05 Sub-01:
   *  - Set up for the alignment
   \******************************************************/

   /*Right know these are hardoced in, but at some piont
   `  it might be nice to allow the user the manipulate.
   */

   if(outFileCStr == 0) outFILE = stdout;
   else outFILE = fopen(outFileCStr, "w");

   seqToLookupIndex(&refST);
   seqToLookupIndex(&queryST);

   queryST.endAlnUL = queryST.lenSeqUL - 1;
   queryST.offsetUL = 0;

   refST.endAlnUL = refST.lenSeqUL - 1;
   refST.offsetUL = 0;

   /******************************************************\
   * Main Sec-05 Sub-02:
   *  - Check if doing an Needleman alignment
   \******************************************************/

   if(settings.useNeedleBl != 0)
     alnMtrxST = NeedlemanAln(&queryST, &refST, &settings);

   /******************************************************\
   * Main Sec-05 Sub-03:
   *  - Check if doing an Waterman alignment
   \******************************************************/

   else if(settings.useWaterBl != 0)
   { /*else if doing a waterman alignment*/
     if(settings.refQueryScanBl)
     { /*IF keeping some alternative alignments*/
       alnMtrxST=WatermanAltAln(&queryST,&refST,&settings);

       if(alnMtrxST == 0) goto alignmentFailed;

       printAltWaterAlns(
          alnMtrxST,
          settings.minScoreL,
          outFILE /*Prefix to name files*/
       ); /*Print out alternative alignment positions*/
     } /*IF keeping some alternative alignments*/

     else alnMtrxST=WatermanAln(&queryST,&refST,&settings);
   } /*else if doing a waterman alignment*/

   /******************************************************\
   * Main Sec-05 Sub-04:
   *  - Check if doing an memory efficent Waterman
   *    alignment that prints out alternative alignment
   *    positions (Gets combined with a Hirschberg)
   \******************************************************/

   else if(settings.memWaterBl != 0)
   { /*Else if; doing a memory waterman alignment*/
      if(settings.refQueryScanBl)
      { /*If doing an query-ref scan*/
         alnMtrxST =
            memWaterAltAln(&queryST, &refST, &settings);

         if(alnMtrxST == 0) goto alignmentFailed;

         if(settings.justScoresBl) goto printAlts;

         refST.offsetUL =
            alnMtrxST->bestScoreST.refStartUL;
         queryST.offsetUL =
            alnMtrxST->bestScoreST.qryStartUL;

         refST.endAlnUL = alnMtrxST->bestScoreST.refEndUL;
         queryST.endAlnUL =alnMtrxST->bestScoreST.qryEndUL;

         alnST = Hirschberg(&refST, &queryST, &settings);

         if(alnST == 0)
         { /*If the hirschberg alignment failed*/
            freeAlnMatrixST(alnMtrxST);
            alnMtrxST = 0;
            goto alignmentFailed;
         } /*If the hirschberg alignment failed*/

         printAlts:

         printAltWaterAlns(
            alnMtrxST,
            settings.minScoreL,
            outFILE /*Prefix to name files*/
         ); /*Print out alternative alignment positions*/

         if(settings.justScoresBl) goto noAlnOutFree;

         bestScoreL = alnMtrxST->bestScoreST.scoreL;
         freeAlnMatrixST(alnMtrxST); /*No longer need*/
         alnMtrxST = 0;

         goto noDirMatrix;
      } /*If doing an query-ref scan*/

      /***************************************************\
      * Main Sec-05 Sub-05:
      *  - Check if doing an memory efficent Waterman
      *    alignment that prints out a single alignment
      *    (uses a Hirschberg to get the alignment)
      \***************************************************/

      else
      { /*Else I am just finding the best alignment*/
         bestScoreST =
            memWaterAln(&queryST,&refST,&settings);

         if(bestScoreST == 0) goto alignmentFailed;

         if(settings.justScoresBl)
         { /*If I am just printing out coordinates*/
            fprintf(
               outFILE,
               "%li\t%lu\t%lu\t%lu\t%lu\n",
               bestScoreST->scoreL,
               bestScoreST->refStartUL,
               bestScoreST->qryStartUL,
               bestScoreST->refEndUL,
               bestScoreST->qryEndUL
            );

            goto noAlnOutFree;
         } /*If I am just printing out coordinates*/

         refST.offsetUL = bestScoreST->refStartUL;
         queryST.offsetUL = bestScoreST->qryStartUL;
         refST.endAlnUL = bestScoreST->refEndUL;
         queryST.endAlnUL = bestScoreST->qryEndUL;

         alnST = Hirschberg(&refST, &queryST, &settings);

         /*Free the current uneeded variables*/
         bestScoreL = bestScoreST->scoreL;
         freeScoresST(bestScoreST, 0);
         bestScoreST = 0;

         if(alnST == 0) goto alignmentFailed;
         else goto noDirMatrix;
      } /*Else I am just finding the best alignment*/
   } /*Else if; doing a memory waterman alignment*/

   /******************************************************\
   * Main Sec-05 Sub-06:
   *  - Check if doing an Hirschberg alignment
   \******************************************************/

   else if(settings.useHirschBl != 0)
   { /*Else if doing an Hirschberg alignment*/
     alnST = Hirschberg(&refST, &queryST, &settings);
     if(alnST == 0) goto alignmentFailed;
     goto noDirMatrix;

     /* The Hirschberg returns an alignment structure,
     `  instead of an directional matrix.
     */
   } /*Else if doing an Hirschberg alignment*/

   /******************************************************\
   * Main Sec-05 Sub-07:
   *  - Let user know this was an invalid alignment
   \******************************************************/

   else
   { /*If no aignment was requested*/
      freeSeqST(&refST, 0);   /*0 to mark on the stack*/
      freeSeqST(&queryST, 0); /*0 to mark on the stack*/

      printHelpMesg(stderr, 1); /*short help*/
      fprintf(
         stderr,
         "Invalid or no aligment method input\n"
      ); /*Let user know the error*/

      exit(-1);
   } /*If no aignment was requested*/

   /******************************************************\
   * Main Sec-05 Sub-08:
   *  - Check if the alignment failed
   \******************************************************/

   if(alnMtrxST == 0)
   { /*If did not have enough memory*/
      alignmentFailed:

      freeSeqST(&refST, 0); // 0 to makr on the stack
      freeSeqST(&queryST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Memory allocation error in aligment step\n"
       );
       exit(-1);
   } /*If did not have enough memory*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-06:
   ^  - Print out multi-alignments or find the single
   ^    alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(outFileCStr != 0) outFILE = fopen(outFileCStr, "w");

   else
   { /*Else putting output to stdout*/
     outFILE = stdout;
     outFileCStr = "out";
   } /*Else putting output to stdout*/

   alnST =
      dirMatrixToAlnST(
         &refST,
         &queryST,
         &alnMtrxST->bestScoreST,
         alnMtrxST->dirMatrixST
   );
   bestScoreL = alnMtrxST->bestScoreST.scoreL;
   freeAlnMatrixST(alnMtrxST); // No longer need
   alnMtrxST = 0;

   if(alnST == 0)
   { // IF i falied to make an alignment array
      freeSeqST(&refST, 0); // 0 to makr on the stack
      freeSeqST(&queryST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Memory error when making alignment array\n"
       );

       exit(-1);
   } // IF i falied to make an alignment array

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-07:
   ^  - Print out the alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   noDirMatrix: /*For aligners that return an alnStruct*/

   lookupIndexToSeq(&refST);
   lookupIndexToSeq(&queryST);

   if(
      printAln(
         outFILE,
         outFileCStr,
         &refST,
         &queryST,
         alnST,
         bestScoreL,
         &settings,
         scoreMtrxFileStr /*File name for scoring matrix*/
      )
   ){ /*If could not print out the alignmnet*/
      fprintf(stderr, "Failed to print alignment\n");

      fclose(outFILE);
      freeAlnST(alnST, 1); /*NEED TO SET UP*/

      freeSeqST(&refST, 0);   /*0 to specify on stack*/
      freeSeqST(&queryST, 0); /* 0 to do a stack free*/

      exit(1);
   } /*If could not print out the alignmnet*/

   freeAlnST(alnST, 1);

   noAlnOutFree: /*When memWater just printing positions*/

   fclose(outFILE);
   freeSeqST(&refST, 0);   /*0 to specify on stack*/
   freeSeqST(&queryST, 0); /* 0 to do a stack free*/

   exit(0);
} /*main*/

/*--------------------------------------------------------\
| Output: Modifies: Each input variable to hold user input
|  - Returns:
|     o 0: If no errors
|     o Pointer to parameter that was an issue
|     o Pointer to "-score-matrix" for invalid scoring
|       matrix input
|  - Prints to stdout when the scoring file is invalid
\--------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,        /*Number arguments user input*/
    char *argsCStr[],       /*Array with user arguments*/
    char **refFileCStr,     /*file name of reference file*/
    char **queryFileCStr,   /*File name of the query file*/
    char **outFileCStr,     /*Name of the output file*/
    char **scoreMtrxFileStr,/*Holds scoring matrix file*/
    struct alnSet *settings /*Aligment settings*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: checkInput
   '  - Checks & extracts user input
   '    fun-01 sec-1: Variable declerations
   '    fun-01 sec-2: Look through user input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    char *singleArgCStr = 0;
    char *dummyStr = 0;
    unsigned long scoreFileErrUL = 0;
    FILE *scoreFILE = 0; /*For loading the scoring matrix*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-01 Sec-2: Look through user input
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    for(int intArg = 1; intArg < *lenArgsInt; intArg += 2)
    { /*loop through all user input arguments*/
        /*The 0 index holds the program name*/
        singleArgCStr = *(argsCStr +intArg + 1);/*arg*/
        tmpCStr = *(argsCStr + intArg);        /*Paramter*/

        if(strcmp(tmpCStr, "-ref") == 0)
            *refFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-query") == 0)
            *queryFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-out") == 0)
            *outFileCStr = singleArgCStr;

        #if !defined NOGAPOPEN
           else if(strcmp(tmpCStr, "-gapopen") == 0)
               settings->gapOpenI =
                   strtol(singleArgCStr, &tmpCStr, 10);
        #endif

        else if(strcmp(tmpCStr, "-gapextend") == 0)
          settings->gapExtendI =
                strtol(singleArgCStr, &tmpCStr, 10);

        else if(strcmp(tmpCStr, "-print-aligned") == 0)
        { /*Else if printing only the aligned region*/
           settings->pFullAlnBl = 0;
           --intArg;
        } /*Else if printing only the aligned region*/

        
        else if(strcmp(tmpCStr, "-print-unaligned") == 0)
        { /*Else if printing only the aligned region*/
           settings->pFullAlnBl = 1;
           --intArg;
        } /*Else if printing only the aligned region*/

        else if(strcmp(tmpCStr, "-print-positions") == 0)
        { /*Else if printing out the base positions*/
           settings->pBasePosBl = 1;
           --intArg;
        } /*Else if printing out the base positions*/

        else if(strcmp(tmpCStr, "-no-positions") == 0)
        { /*Else if not printing out the base positions*/
           settings->pBasePosBl = 0;
           --intArg;
        } /*Else if not printing out the base positions*/

        else if(strcmp(tmpCStr, "-line-wrap") == 0)
          cStrToUSht(singleArgCStr, &settings->lineWrapUS);

        else if(strcmp(tmpCStr, "-use-needle") == 0)
        { /*Else if using a needleman alignment*/
            settings->useNeedleBl = 1;
            settings->useWaterBl = 0;
            settings->useHirschBl = 0;
            settings->memWaterBl = 0;
            --intArg;
        } /*Else if using a needleman alignment*/

        else if(strcmp(tmpCStr, "-use-water") == 0)
        { /*Else if doing a waterman smith alignment*/
            settings->useNeedleBl = 0;
            settings->useWaterBl = 1;
            settings->useHirschBl = 0;
            settings->memWaterBl = 0;
            --intArg;
        } /*Else if doing a waterman smith alignment*/

        else if(strcmp(tmpCStr, "-use-hirschberg") == 0)
        { /*Else if doing a Hirshberg alignment*/
            settings->useNeedleBl = 0;
            settings->useWaterBl = 0;
            settings->useHirschBl = 1;
            settings->memWaterBl = 0;
            --intArg;
        } /*Else if doing a Hirshberg alignment*/

        else if(strcmp(tmpCStr, "-use-mem-water") == 0)
        { /*Else if I am doing a memory effecient water*/
            settings->useNeedleBl = 0;
            settings->useWaterBl = 0;
            settings->useHirschBl = 0;
            settings->memWaterBl = 1;
            --intArg;
        } /*Else if I am doing a memory effecient water*/

        else if(strcmp(tmpCStr, "-only-scores") == 0)
        { /*Else if only printing scores for memWater*/
           settings->justScoresBl = 1;
           --intArg;
        } /*Else if only printing scores for memWater*/

        else if(strcmp(tmpCStr, "-scores-and-aln") == 0)
        { /*Else if only printing alignments for memWater*/
           settings->justScoresBl = 0;
           --intArg;
        } /*Else if only printing alignments for memWater*/

        else if(strcmp(tmpCStr, "-format-expand-cig") == 0)
        { /*Else if using the expanded cigar format*/
           settings->formatFlag = defExpandCig;
           --intArg;
        } /*Else if using an expanded cigar format*/

        else if(strcmp(tmpCStr, "-format-emboss") == 0)
        { /*Else if using an EMBOSS like format*/
           settings->formatFlag = defEMBOSS;
           --intArg;
        } /*Else if using an EMBOSS like format*/

        else if(strcmp(tmpCStr, "-format-clustal") == 0)
        { /*Else if using clustal format*/
           settings->formatFlag = defClustal;
           --intArg;
        } /*Else if using clustal format*/

        else if(strcmp(tmpCStr, "-format-fasta") == 0)
        { /*Else if outputing fasta format*/
           settings->formatFlag = defFasta;
           --intArg;
        } /*Else if outputing fasta format*/

        else if(
          strcmp(tmpCStr, "-query-ref-scan") == 0
        )
        { /*Else if doing more than the best alignment*/
          settings->refQueryScanBl = 1;
          settings->useNeedleBl = 0;
          settings->useHirschBl = 0;

          if(!(settings->useWaterBl))
             settings->memWaterBl = 1;

          --intArg;
        } /*Else if doing more than the best alignment*/

        else if(strcmp(tmpCStr, "-min-score") == 0)
          settings->minScoreL =
             strtol(singleArgCStr, &dummyStr, 10);

        /*Only do these checks when the user has not
        ` called a direction flag during compile time
        */
        #if defined SNPINSDEL
        #elif defined SNPDELINS
        #elif defined INSSNPDEL 
        #elif defined INSDELSNP
        #elif defined DELSNPINS
        #elif defined DELINSSNP
        #else
           else if(strcmp(tmpCStr, "-match-ins-del") == 0)
           { /*Else if matches->insertions->deletions*/
               settings->bestDirC = defSnpInsDel;
               --intArg;
           } /*Else if matches->insertions->deletions*/

           else if(strcmp(tmpCStr, "-match-del-ins") == 0)
           { /*Else if matches->deletions->insertions*/
               settings->bestDirC = defSnpDelIns;
               --intArg;
           } /* Else if matches->deletions->insertions*/

           else if(strcmp(tmpCStr, "-ins-match-del") == 0)
           { /*Else if insertions->matches->deletions*/
               settings->bestDirC = defInsSnpDel;
               --intArg;
           } /*Else if insertions->matches->deletions*/

           else if(strcmp(tmpCStr, "-del-match-ins") == 0)
           { /*Else if deletions->matches->insertions*/
               settings->bestDirC = defDelSnpIns;
               --intArg;
           } /*Else if deletions->matches->insertions*/

           else if(strcmp(tmpCStr, "-ins-del-match") == 0)
           { /*Else if insertions->deletions->matches*/
               settings->bestDirC = defInsDelSnp;
               --intArg;
           } /*Else if insertions->deletions->matches*/

           else if(strcmp(tmpCStr, "-del-ins-match") == 0)
           { /*Else if deletions->insertions->matches*/
               settings->bestDirC = defDelInsSnp;
               --intArg;
           } /*Else if deletions->insertions->matches*/

        #endif

        else if(strcmp(tmpCStr, "-score-matrix") == 0)
        { /*else if the user supplied a scoring matrix*/
            *scoreMtrxFileStr = singleArgCStr;
            scoreFILE = fopen(singleArgCStr, "r");

            if(scoreFILE == 0)
            { /*If I could not open the scoring file*/
                return tmpCStr; /*So user knows invalid*/

                fprintf(stderr,
                  "-score-matrix %s is not an file\n",
                  singleArgCStr
                ); /*Print out the problem*/
            } /*If I could not open the scoring file*/

            scoreFileErrUL =
              readInScoreFile(settings, scoreFILE);

            if(scoreFileErrUL != 0)
            { /*If the scoring file had an errors*/
              fprintf(
                stderr,
                "Invalid line (%lu) in -score-matrix %s\n",
                scoreFileErrUL,
                singleArgCStr
              ); /*Print out the problem*/

                return tmpCStr;    /*invalid file*/
            } /*If the scoring file had an error*/
        } /*else if the user supplied a scoring matrix*/
            
        else return tmpCStr; /*Invalid parameter*/
    } /*loop through all user input arguments*/

    return 0; /*input is valid*/
} /*checkInput*/

/*--------------------------------------------------------\
| Output:
|  - Prints:
|    o help message to outFILE.
\--------------------------------------------------------*/
void printHelpMesg(
   FILE *outFILE,
   char breifBl   /*Print a shorter help message*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printHelpMesg
   '  - Prints out the help message to outFILE
   '  o fun-03 Sec-01:
   '    - Usage block
   '  o fun-03 Sec-02:
   '    - Input block
   '  o fun-03 sec-03:
   '    - Output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Usage block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      outFILE,
      "alnSeq -query query.fasta -ref ref.fasta"
   );
   fprintf(outFILE, " [options...]\n");

   fprintf(outFILE, "Use:\n");
   fprintf(
     outFILE,
     "  - Does a pairwise alignment on two input sequences"
   );
   fprintf(outFILE, "\n");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Input block
   ^  o fun-03 sec-02 sub-01:
   ^     - General input/output block (includes line wrap)
   ^  o fun-03 sec-02 sub-02:
   ^     - Alignment algorithim selection block
   ^  o fun-03 sec-02 sub-03:
   ^     - Alignment paramaters block
   ^  o fun-03 sec-02 sub-04:
   ^    - File output settings (non-format)
   ^  o fun-03 sec-02 sub-05:
   ^     - File output format block
   ^  o fun-03 sec-02 sub-06:
   ^     - Waterman specific paramters block
   ^  o fun-03 sec-02 sub-07:
   ^     - Selecting alignment direction block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-01:
   *  - General input/output block (includes line wrap)
   \******************************************************/

   fprintf(outFILE, "Input:\n");

   fprintf(outFILE, "  -query: [Required]\n");
   fprintf(
      outFILE,
      "    o Fasta file with query sequence.\n"
   );

   fprintf(outFILE, "  -ref: [Required]\n");
   fprintf(
      outFILE,
      "    o Fasta file with reference sequence.\n"
   );

   fprintf(outFILE, "  -out: [stdout]\n");
   fprintf(outFILE,"    o File to output alignment to.\n");

   /******************************************************\
   * Fun-03 Sec-02 Sub-02:
   *  - Alignment algorithim selection block
   \******************************************************/

   if(defUseNeedle)
      fprintf(outFILE, "  -use-needle: [Yes]\n");
   else fprintf(outFILE, "  -use-needle: [No]\n");

   fprintf(
      outFILE,
      "    o Do a Needleman Wunsch alignment.\n" 
   );

   if(defUseWater)
      fprintf(outFILE, "  -use-water: [Yes]\n");
   else fprintf(outFILE, "  -use-water: [No]\n");

   fprintf(
      outFILE,
      "    o Do a Waterman Smith alignment.\n"
   );

   if(defUseHirsch)
       fprintf(outFILE, "  -use-hirschberg: [Yes]\n");
   else fprintf(outFILE, "  -use-hirschberg: [No]\n");

   fprintf(outFILE, "    o Do a Hirschberg alignment.\n");

   if(defUseMemWater)
       fprintf(outFILE, "  -use-mem-water: [Yes]\n");
   else fprintf(outFILE, "  -use-mem-water: [No]\n");

   fprintf(
      outFILE,
      "    o Do a memory efficent Waterman alignment.\n"
   );

   fprintf(
      outFILE,
      "    o Slower than Waterman, but uses O(N) memory.\n"
   );

   fprintf(
     outFILE,
     "    o Uses a Waterman to find the starting and"
   );
   fprintf(outFILE, " ending\n      ");
   fprintf(
      outFILE,
     "positions and a Hirschberg to find the alignment.\n"
   );

   /******************************************************\
   * Fun-03 Sec-02 Sub-03:
   *  - Alignment paramaters block
   \******************************************************/

   #if !defined NOGAPOPEN
      fprintf(outFILE, "  -gapopen: [%i]\n", defGapOpen);
      fprintf(
         outFILE,
         "    o Cost of starting an indel (as integer). A "
      );
      fprintf(
         outFILE,
         "negative\n      value is a penalty.\n"
      );
   #endif

   fprintf(outFILE,"  -gapextend: [%i]\n",defGapExtend);
   fprintf(
     outFILE,
     "    o Cost of extending an indel one base (< 0 is"
   );
   fprintf(outFILE, " penalty).\n");

   fprintf(
      outFILE,
      "  -score-matrix: [%s]\n",
      defMatrixNameStr
   );
   fprintf(
      outFILE,
      "    o File with scoring matrix to use. It should"
   );
   fprintf(
      outFILE,
      " have one\n      line for each score to change.\n"
   );
   fprintf(
      outFILE,
      "    o (see scoring-matrix.txt for an example).\n"
   );

   if(breifBl) goto helpOutputBlock;

   /******************************************************\
   * Fun-03 Sec-02 Sub-04:
   *  - File output settings (non-format)
   \******************************************************/

   if(defJustScoresBl)
      fprintf(outFILE,"  -only-scores: [Yes]\n");
   else
      fprintf(outFILE,"  -only-scores: [No]\n");

   fprintf(
      outFILE,
      "    o Only print out the score, starting, and"
   );
   fprintf(outFILE, " ending\n");
   fprintf(
      outFILE,
      "      coordinates of an alingnment (mem-water only)"
   );
   fprintf(outFILE, "\n");
   fprintf(
      outFILE,
      "    o Disable this setting with -scores-and-aln\n"
   );

   if(defPAln) fprintf(outFILE,"  -print-aligned: [No]\n");
   else fprintf(outFILE,"  -print-aligned: [Yes]\n");

   fprintf(
      outFILE,
      "    o Only print the aligned regions of each"
   );
   fprintf(outFILE, " sequence.\n");

   if(!defPAln)
      fprintf(outFILE,"  -print-unaligned: [Yes]\n");
   else fprintf(outFILE,"  -print-unaligned: [No]\n");

   fprintf(
      outFILE,
      "    o Print out the entire reference and query"
   );
   fprintf(outFILE, " sequence.\n");

   if(defPPos)
      fprintf(outFILE, "  -print-positions: [Yes]\n");
   else fprintf(outFILE, "  -print-positions: [No]\n");

   fprintf(
      outFILE,
      "    o Print the starting and ending position for"
   );
   fprintf(outFILE, " each line\n");
   fprintf(outFILE, "      in the alignment.\n");
   fprintf(
      outFILE,
     "    o Clustal format only prints the ending position"
   );
   fprintf(outFILE, ".\n");

   if(!defPPos)
      fprintf(outFILE, "  -no-positions: [Yes]\n");
   else fprintf(outFILE, "  -no-positions: [No]\n");

   fprintf(
      outFILE,
      "    o Turns off -print-positions\n"
   );
   fprintf(
      outFILE,
      "    o EMBOSS format always prints positions.\n"
   );


   fprintf(outFILE, "  -line-wrap: [%i]\n", defLineWrap);
   fprintf(
      outFILE,
      "    o Maximum characters per line in output file.\n"
   );
   fprintf(outFILE, "    o Input 0 for no line wrap.\n");
   fprintf(
      outFILE,
      "    o Minimum line wrap is 10 for fasta (header not"
   );
   fprintf(outFILE, "\n");
   fprintf(
      outFILE,
      "      wrapped), 32 for clustal, & 42 for expanded"
   );
   fprintf(outFILE, " cigar\n      and EMBOSS.\n");

   /******************************************************\
   * Fun-03 Sec-02 Sub-05:
   *  - File output format block
   \******************************************************/

   if(defFormat == defExpandCig)
      fprintf(outFILE, "  -format-expand-cig: [Yes]\n");
   else fprintf(outFILE, "  -format-expand-cig: [No]\n");

   fprintf(
      outFILE,
      "    o Prints the reference sequence, then query"
   );
   fprintf(outFILE, " sequence,\n");
   fprintf(outFILE, "      and then the eqx line.\n");
   fprintf(outFILE, "    o eqx line format:\n");
   fprintf(outFILE, "      - I = Insertion\n");
   fprintf(outFILE, "      - D = Deletion\n");
   fprintf(outFILE, "      - X = mismatch\n");
   fprintf(outFILE, "      - = = match\n");
   fprintf(outFILE, "      - S = soft mask\n");
   fprintf(outFILE, "    o reference/query lines\n");
   fprintf(
      outFILE,
      "      - Without -print-positions\n"
   );

   fprintf(outFILE, "        - Ref: sequence\n");
   fprintf(outFILE, "        - Qry: sequence\n");
   fprintf(
      outFILE,
      "      - With -print-positions\n"
   );
   fprintf(
      outFILE,
      "        - Ref: start-base sequence end-base\n"
   );
   fprintf(
      outFILE,
      "        - Qry: start-base sequence end-base\n"
   );


   if(defFormat == defEMBOSS)
      fprintf(outFILE, "  -format-emboss: [Yes]\n");
   else fprintf(outFILE, "  -format-emboss: [No]\n");

   fprintf(
      outFILE,
      "    o Prints the reference sequence, then eqx line,"
   );
   fprintf(outFILE,"and then\n      the query sequence.");
   fprintf(outFILE, "\n    o eqx line format:\n");
   fprintf(outFILE, "      - | = Match\n");
   fprintf(outFILE, "      - space = SNP/gap\n");
   fprintf(outFILE, "    o reference/query lines\n");
   fprintf(
      outFILE,
      "      - ID: start-base sequence end-base\n"
   );

   if(defFormat == defClustal)
      fprintf(outFILE, "  -format-clustal: [Yes]\n");
   else fprintf(outFILE, "  -format-clustal: [No]\n");

   fprintf(
      outFILE,
      "    o Prints the reference sequence, then query"
   );
   fprintf(outFILE, " sequence,\n");
   fprintf(outFILE, "      and then the eqx line.\n");
   fprintf(outFILE, "    o eqx line format:\n");
   fprintf(outFILE, "      - * = Match\n");
   fprintf(outFILE, "      - space = SNP/gap\n");
   fprintf(outFILE, "    o reference/query lines\n");
   fprintf(
      outFILE,
      "      - Without -print-positions\n"
   );
   fprintf(outFILE, "        - ID sequence\n");
   fprintf(
      outFILE,
      "      - With -print-positions\n"
   );

   fprintf(outFILE, "        - ID sequence end-base\n");

   if(defFormat == defFasta)
      fprintf(outFILE, "  -format-fasta: [Yes]\n");
   else fprintf(outFILE, "  -format-fasta: [No]\n");

   fprintf(
      outFILE, "    o Save alignment as a fasta file.\n"
   );
   fprintf(
      outFILE, "    o >id score first-base last-base\n"
   );
   fprintf(
      outFILE, "    o first-base is first aligned base\n"
   );
   fprintf(
      outFILE, "    o last-base is last aligned base\n"
   );

   /******************************************************\
   * Fun-03 Sec-02 Sub-06:
   *  - Waterman specific paramters block
   \******************************************************/

   if(defQueryRefScan)
       fprintf(outFILE,"  -query-ref-scan: [Yes]\n");
   else fprintf(outFILE,"  -query-ref-scan: [No]\n");

   fprintf(
      outFILE,
      "    o Waterman alignments only.\n"
   );
   fprintf(
      outFILE,
      "    o prints out the best score for each"
   );
   fprintf(
      outFILE,
      " reference and\n      query base.\n"
   );
   fprintf(
      outFILE,
      "    o Uses -use-mem-water if -use-water not used.\n"
   );

   fprintf(outFILE, "  -min-score: [%i]\n", defMinScore);
   fprintf(outFILE, "    o Waterman alignments only.\n");
   fprintf(
      outFILE,
      "    o Minimum score needed to keep an alternative\n"
   );
   fprintf(outFILE, "      alignment.\n");

   /******************************************************\
   * Fun-03 Sec-02 Sub-07:
   *  - Selecting alignment direction block
   \******************************************************/
   #if defined SNPINSDEL
   #elif defined SNPDELINS
   #elif defined INSSNPDEL 
   #elif defined INSDELSNP
   #elif defined DELSNPINS
   #elif defined DELINSSNP
   #else

      if(defBestDir == defSnpInsDel)
         fprintf(outFILE, "  -match-ins-del: [Yes]\n");
      else  fprintf(outFILE, "  -match-ins-del: [No]\n");

      fprintf(
         outFILE,
         "    o For equal scores choose matches/SNPs over "
      );
      fprintf(outFILE, "insertions\n      and");
      fprintf(outFILE, " insertions over deletions.\n");

      if(defBestDir == defSnpDelIns)
         fprintf(outFILE, "  -match-del-ins: [Yes]\n");
      else  fprintf(outFILE, "  -match-del-ins: [No]\n");

      fprintf(
         outFILE,
         "    o For equal scores choose matches/SNPs over "
      );
      fprintf(outFILE, "deletions\n      and");
      fprintf(outFILE, " deletions over insertions.\n");

      if(defBestDir == defInsSnpDel)
         fprintf(outFILE, "  -ins-match-del: [Yes]\n");
      else  fprintf(outFILE, "  -ins-match-del: [No]\n");

      fprintf(
         outFILE,
         "    o For equal scores choose insertions over "
      );
      fprintf(outFILE, "matches/SNPs\n      and");
      fprintf(outFILE, " matches/SNPs over deletions.\n");

      if(defBestDir == defInsDelSnp)
         fprintf(outFILE, "  -ins-del-match: [Yes]\n");
      else  fprintf(outFILE, "  -ins-del-match: [No]\n");

      fprintf(
         outFILE,
         "    o For equal scores choose insertions over "
      );
      fprintf(outFILE, "deletions\n      and");
      fprintf(outFILE, " deletions over matches/SNPs.\n");

      if(defBestDir == defDelSnpIns)
         fprintf(outFILE, "  -del-match-ins: [Yes]\n");
      else  fprintf(outFILE, "  -del-match-ins: [No]\n");

      fprintf(
         outFILE,
         "    o For equal scores choose deletions over "
      );
      fprintf(outFILE, "matches/SNPs\n      and");
      fprintf(outFILE, " matches/SNPs over insertions.\n");

      if(defBestDir == defDelInsSnp)
         fprintf(outFILE, "  -del-ins-match: [Yes]\n");
      else  fprintf(outFILE, "  -del-ins-match: [No]\n");

      fprintf(
         outFILE,
         "    o For equal scores choose deletions over "
      );
      fprintf(outFILE, "insertions\n      and");
      fprintf(outFILE, " insertions over matches/SNPs.\n");
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - Output block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   helpOutputBlock:

   if(breifBl)
     fprintf(
       outFILE,
       "  -h-all:\n    o Print the entire help message\n"
      );

   fprintf(
      outFILE,
      "  -v:\n    o Print alnSeq version.\n"
   );

   fprintf(
      outFILE,
      "  -flags:\n"
   );

   fprintf(
      outFILE,
      "    o Prints complier flags with descriptions\n"
   );

   fprintf(
      outFILE,
      "  -flags-only:\n"
   );

   fprintf(
      outFILE,
      "    o Prints the flags used and version number\n"
   );

   /*Output block*/
   fprintf(outFILE, "Output:\n");
   fprintf(
      outFILE,
      "  - Alignment to stdout or to file provided by"
   );
   fprintf(outFILE, " -out.\n");

   return;
} /*printHelpMesg*/

/*--------------------------------------------------------\
| Name: printCompilerSettings
| Call: printCompilerSettings(file, 1/0)
| Use:
|   - Prints out the compiler flags used
|   - Prints out what each possible compiler flag does
| Input:
|   - outFILE:
|     o File to print flags and flag descriptions to
|   - pDescBl:
|     o 1: print out descriptions for each flag (all)
|     o 0: Do not print out any flag descriptions
| Output:
|   - Prints flags and description to outFILE
\--------------------------------------------------------*/
void printCompileSettings(
   FILE *outFILE, /*Output file*/
   char pDescBl   /*not 0: print out flag descriptions*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC:
   '  - Print out compile flags used and what each compiler
   '    flag does.
   '  o fun-04 sec-01:
   '    - Print out the compiled settings
   '  o fun-04 sec-02:
   '    - Print out what each compiler flag does
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-01:
   ^  - Print out the compiled settings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(outFILE, "Compiler flags:\n");

   #if defined BYTEMATRIX
      fprintf(outFILE, "   -DBYTEMATRIX\n");
   #endif

   #if defined HIRSCHTWOBIT
      fprintf(outFILE, "   -DHIRSCHTWOBIT\n");
   #endif

   #if defined TWOBITMSW
      fprintf(outFILE, "   -DTWOBITMSW\n");
   #endif

   #if defined NOGAPOPEN
      fprintf(outFILE, "   -DNOGAPOPEN\n");
   #endif

   #if defined SNPINSDEL
      fprintf(outFILE, "   -DSNPINSDEL\n");
   #elif defined SNPDELINS
      fprintf(outFILE, "   -DSNPDELINS\n");
   #elif defined INSSNPDEL 
      fprintf(outFILE, "   -DINSSNPDEL\n");
   #elif defined INSDELSNP
      fprintf(outFILE, "   -DINSDELSNP\n");
   #elif defined DELSNPINS
      fprintf(outFILE, "   -DDELSNPINS\n");
   #elif defined DELINSSNP
      fprintf(outFILE, "   -DDELINSSNP\n");
   #endif

   if(!pDescBl) return;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-02:
   ^  - Print out what each compiler flag does
   ^  o fun-04 sec-02 sub-01:
   ^    - Print out the header
   ^  o fun-04 sec-02 sub-02:
   ^    - Print out what the -DBYTEMATRIX flag does
   ^  o fun-04 sec-02 sub-03:
   ^    - Print out what the -DHIRSCHTWOBIT flag does
   ^  o fun-04 sec-02 sub-04:
   ^    - Print out what the -DTWOBITMSW flag does
   ^  o fun-04 sec-02 sub-05:
   ^    - Print out what the -DNOGAPOPEN flag does
   ^  o fun-04 sec-02 sub-06:
   ^    - Print out what the -DSNPINSDEL flag does
   ^  o fun-04 sec-02 sub-07:
   ^    - Print out what the -DSNPDELINS flag does
   ^  o fun-04 sec-02 sub-08:
   ^    - Print out what the -DINSSNPDEL flag does
   ^  o fun-04 sec-02 sub-09:
   ^    - Print out what the -DINSDELSNP flag does
   ^  o fun-04 sec-02 sub-10:
   ^    - Print out what the -DDELSNPINS flag does
   ^  o fun-04 sec-02 sub-11:
   ^    - Print out what the -DDELINSSNP flag does
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-04 Sec-02 Sub-01:
   *  - Print out the header
   \******************************************************/

   fprintf(outFILE, "What each flag does:\n");

   fprintf(
      outFILE,
      "   - These flags are used to make alternative"
   );
   fprintf(outFILE, " versions\n     of alnSeq.\n");

   fprintf(
      outFILE,
      "   - Make alternate alnSeq:"
   );
   fprintf(
      outFILE,
      " make CFLAGS=\"flag1 flag2 ...\"\n"
   );
   /******************************************************\
   * Fun-04 Sec-02 Sub-02:
   *  - Print out what the -DBYTEMATRIX flag does
   \******************************************************/

   fprintf(outFILE, "   -DBYTEMATRIX:\n");
   fprintf(
      outFILE,
      "     - Makes the Needleman and Water man alingments"
   );

   fprintf(outFILE, " use a\n");
   fprintf(
      outFILE,
      "       byte directional matrix instead of a two bit"
   );
   fprintf(outFILE, " matrix.\n");

    fprintf(
       outFILE,
       "     - This gives a faster Needleman and"
    );
    fprintf(
       outFILE,
       " Waterman\n       alignment, but also increases "
    );
    fprintf( outFILE, " memory usage by 4x.\n");

   /******************************************************\
   * Fun-04 Sec-02 Sub-03:
   *  - Print out what the -DHIRSCHTWOBIT flag does
   \******************************************************/

   fprintf(outFILE, "   -DHIRSCHTWOBIT:\n");
   fprintf(
      outFILE,
      "     - This makes the Hirschberg alignment use two"
   );

   fprintf(outFILE, " bit\n");
   fprintf(
      outFILE,
      "       arrays instead of byte arrays.\n"
   );

   fprintf(
      outFILE,
      "     - This gives a small reduction in memory, but"
   );
   fprintf(
      outFILE,
      " also\n       makes alignments take 2x more time.\n"
   );

   /******************************************************\
   * Fun-04 Sec-02 Sub-04:
   *  - Print out what the -DTWOBITMSW flag does
   \******************************************************/

   fprintf(outFILE, "   -DTWOBITMSW\n");
   fprintf(
     outFILE,
     "     - Memory Waterman alignment uses two-bit arrays"
   );

   fprintf(outFILE, "\n");
   fprintf(
      outFILE,
      "       instead of byte arrays.\n"
   );

   fprintf(
      outFILE,
      "     - This gives a small reduction in memory, but"
   );
   fprintf(
      outFILE,
      " also\n       makes alignments take 2x more time.\n"
   );

   /******************************************************\
   * Fun-04 Sec-02 Sub-05:
   *  - Print out what the -DNOGAPOPEN flag does
   \******************************************************/

   fprintf(outFILE, "   -DNOGAPOPEN:\n");
   fprintf(
      outFILE,
      "     - Removes the gap opening penalty.\n"
   );
   fprintf(
      outFILE,
      "     - This decreases the alignment time, but may"
   );
   fprintf(
      outFILE,
      " also\n       decrease the alignment quality.\n"
   );

   /******************************************************\
   * Fun-04 Sec-02 Sub-06:
   *  - Print out what the -DSNPINSDEL flag does
   \******************************************************/

   fprintf(outFILE, "   -DSNPINSDEL:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This makes"
   );
   fprintf(
      outFILE,
      " alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: SNPs/Matches, then insertions, then"
   );
   fprintf(outFILE, "\n       deletions.\n");

   /******************************************************\
   * Fun-04 Sec-02 Sub-07:
   *  - Print out what the -DNSNPDELINS flag does
   \******************************************************/

   fprintf(outFILE, "   -DSNPDELINS:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This makes"
   );
   fprintf(
      outFILE,
      " alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: SNPs/Matches, then deletinos, then"
   );
   fprintf(outFILE, "\n       insertions.\n");

   /******************************************************\
   * Fun-04 Sec-02 Sub-08:
   *  - Print out what the -DINSSNPDEL flag does
   \******************************************************/

   fprintf(outFILE, "   -DINSSNPDEL:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This makes"
   );
   fprintf(
      outFILE,
      " alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: insertions, then SNPs/Matches, then"
   );
   fprintf(outFILE, "\n       deletions.\n");

   /******************************************************\
   * Fun-04 Sec-02 Sub-09:
   *  - Print out what the -DINSDELSNP flag does
   \******************************************************/

   fprintf(outFILE, "   -DINSDELSNP:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This makes"
   );
   fprintf(
      outFILE,
      " alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: insertions, then deletions, then"
   );
   fprintf(outFILE, "\n       SNPs/Matches.\n");

   /******************************************************\
   * Fun-04 Sec-02 Sub-10:
   *  - Print out what the -DDELSNPINS flag does
   \******************************************************/

   fprintf(outFILE, "   -DDELSNPINS:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This makes"
   );
   fprintf(
      outFILE,
      " alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: deletions, then SNPs/Matches, then"
   );
   fprintf(outFILE, "\n       deletions.\n");

   /******************************************************\
   * Fun-04 Sec-02 Sub-11:
   *  - Print out what the -DDELINSSNP flag does
   \******************************************************/

   fprintf(outFILE, "   -DDELINSSNP:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This makes"
   );
   fprintf(
      outFILE,
      " alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: deletions, then insertions, then"
   );
   fprintf(outFILE, "\n       SNPs/Matches.\n");

   return;
} /*printCompileSettings*/
