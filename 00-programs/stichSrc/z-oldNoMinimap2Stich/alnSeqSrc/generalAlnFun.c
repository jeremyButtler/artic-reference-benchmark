/*#########################################################
# Name generalAlnFun
# Use:
#  o Holds general functions used in my Needleman Wunsch
#    or Waterman Smith alignment.
# Libraries:
#   - "alnSetStruct.h"
#   o "alignmentSettings.h"
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>  // Used by alnSetStructure.h
#########################################################*/

#include "generalAlnFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
' o fun-01 twoBitMaxScore:
'   - Picks the best score and direction for the current
'     base pairs being compared in an alignment. This
'     function is set up for two bit arrays.
'   - Function is inlied (in generalAlnFun.h)
' o fun-02 charMaxScore:
'   - Picks the best score and direction for the current
'     base pairs being compared in an alignment. This
'     function is set up for charters
'   - Function is inlined (in generalAlnFun.h)
'  o fun-03 alnMaxScore:
'    - Picks the best score for the current base pairs
'      being compared in an alignment.
'   - Function is inlied (in generalAlnFun.h)
'  o fun-03 checkIfBasesMatch:
'    - Are two bases are same? (includes anonymous bases)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o 1: if bases were a match
|    o 0 if bases do not mach
\--------------------------------------------------------*/
char checkIfBasesMatch(
    char *queryBaseC,// Query base to compare to reference
    char *refBaseC   // Reference base to compare to query
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-1 Sub-1: checkIfBasesMatch
   '  - Are two bases are same? (includes anonymous bases)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   
   /*-64 is to account for inputing lookup indexes*/
   // The switch should default to a look up table &
   // will be more clear
   switch(*queryBaseC & defToUper)
   { // Switch: Check if bases are the same
       case ('A' - 64):
       case 'A':
       // Case: Query is an A
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('A' - 64):
               case ('W' - 64):
               case ('M' - 64):
               case ('R' - 64):
               case ('D' - 64):
               case ('H' - 64):
               case ('V' - 64):
               case ('N' - 64):
               case ('X' - 64):
               case 'A':
               case 'W':
               case 'M':
               case 'R':
               case 'D':
               case 'H':
               case 'V':
               case 'N':
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an A

       case ('T' - 64):
       case 'T':
       // Case: Query is an T
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('T' - 64):
               case ('U' - 64):
               case ('W' - 64):
               case ('K' - 64):
               case ('B' - 64):
               case ('Y' - 64):
               case ('D' - 64):
               case ('H' - 64):
               case ('N' - 64):
               case ('X' - 64):

               case 'T': return 1;
               case 'U': return 1;
               case 'W': return 1;
               case 'K': return 1;
               case 'B': return 1;
               case 'Y': return 1;
               case 'D': return 1;
               case 'H': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an T

       case ('U' - 64):
       case 'U':
       // Case: Query is an U
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('T' - 64):
               case ('U' - 64):
               case ('W' - 64):
               case ('K' - 64):
               case ('B' - 64):
               case ('Y' - 64):
               case ('D' - 64):
               case ('H' - 64):
               case ('N' - 64):
               case ('X' - 64):

               case 'T': return 1;
               case 'U': return 1;
               case 'W': return 1;
               case 'K': return 1;
               case 'B': return 1;
               case 'Y': return 1;
               case 'D': return 1;
               case 'H': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an U

       case ('C' - 64):
       case 'C':
       // Case: Query is an C
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('C' - 64):
               case ('S' - 64):
               case ('M' - 64):
               case ('Y' - 64):
               case ('B' - 64):
               case ('H' - 64):
               case ('V' - 64):
               case ('N' - 64):
               case ('X' - 64):

               case 'C': return 1;
               case 'S': return 1;
               case 'M': return 1;
               case 'Y': return 1;
               case 'B': return 1;
               case 'H': return 1;
               case 'V': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an C

       case ('G' - 64):
       case 'G':
       // Case: Query is an G
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('G' - 64):
               case ('S' - 64):
               case ('K' - 64):
               case ('R' - 64):
               case ('B' - 64):
               case ('D' - 64):
               case ('V' - 64):
               case ('N' - 64):
               case ('X' - 64):

               case 'G': return 1;
               case 'S': return 1;
               case 'K': return 1;
               case 'R': return 1;
               case 'B': return 1;
               case 'D': return 1;
               case 'V': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an G

       case ('W' - 64):
       case 'W':
       // Case: Query is an W
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('C' - 64):
               case ('G' - 64):
               case ('S' - 64):

               case 'C': return 0;
               case 'G': return 0;
               case 'S': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an W

       case ('S' - 64):
       case 'S':
       // Case: Query is an S
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('A' - 64):
               case ('T' - 64):
               case ('U' - 64):
               case ('W' - 64):

               case 'A': return 0;
               case 'T': return 0;
               case 'U': return 0;
               case 'W': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an S

       case ('M' - 64):
       case 'M':
       // Case: Query is an M
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('G' - 64):
               case ('T' - 64):
               case ('U' - 64):
               case ('K' - 64):

               case 'G': return 0;
               case 'T': return 0;
               case 'U': return 0;
               case 'K': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an M

       case ('K' - 64):
       case 'K':
       // Case: Query is an K
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('A' - 64):
               case ('C' - 64):
               case ('M' - 64):

               case 'A': return 0;
               case 'C': return 0;
               case 'M': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an K

       case ('R' - 64):
       case 'R':
       // Case: Query is an R
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('C' - 64):
               case ('T' - 64):
               case ('U' - 64):
               case ('Y' - 64):

               case 'C': return 0;
               case 'T': return 0;
               case 'U': return 0;
               case 'Y': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an R

       case ('Y' - 64):
       case 'Y':
       // Case: Query is an Y
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('A' - 64):
               case ('G' - 64):
               case ('R' - 64):

               case 'A': return 0;
               case 'G': return 0;
               case 'R': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an Y

       case ('B' - 64):
       case 'B':
       // Case: Query is an B
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('A' - 64):
               case 'A': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an B

       case ('D' - 64):
       case 'D':
       // Case: Query is an D
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('C' - 64):
               case 'C': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an C

       case ('H' - 64):
       case 'H':
       // Case: Query is an D
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('G' - 64):
               case 'G': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an H

       case ('V' - 64):
       case 'V':
       // Case: Query is an D
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case ('T' - 64):
               case ('U' - 64):
               case 'T': return 0;
               case 'U': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an V

       case ('N' - 64):
       case ('X' - 64):
       case 'N':
       case 'X': return 1;
   } // Switch: Check if bases are the same

   return 0; // not a base
} // checkIfBasesMatch
