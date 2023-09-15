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

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' stichInputAndHelp SOH: Start Of Header
'  - Checks user input and prints out help message
'  o fun-01 stichGetInput:
'    - Gets the user input (calls stichCheckArg)
'  o fun-02 stichCheckArg:
'    - Checks a single paramter/argument combination
'  o fun-03 pStichHelp:
'    - Prints the help message for stich
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef STICHINPUTANDHELP_H
#define STICHINPUTANDHELP_H

#include "samFunSrc/cStrFun.h"
#include "samFunSrc/cStrToNumberFun.h"
#include "stichSetStruct.h"
#include <stdio.h>
#include <string.h>

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   | Fun-01 TOC: stichGetInput
   |  - Get user input into stich
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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

#endif
