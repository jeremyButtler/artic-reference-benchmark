/*#########################################################
# Name: alnStruct
# Use:
#  - Holds the alingment structure that stores the
#    alignment. This also includes the functions needed
#    to maintain and print out the alignment
# Libraries:
#  - "seqStruct.h"
#  - "scoresST.h"
#  - "alnSetStruct.h"
#  - "generalAlnFun.h"
#  - include "alnMatrixStruct.h"
#  o "twoBitArrays.h"
#  o "alnSeqDefaults.h"
# C Standard Libraries:
#  - <time.h>
#  - <string.h>
#  o <stdio.h>
#  o <stdlib.h>
#  o <stdint.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o fun-01 alnSTToSeq:
'    - Makes an alignment for an single sequence
'      (reference or query)
'  o fun-02 dirMatrixToAlnST:
'    - Builds an anlignment array for the input direction
'  o fun-03 printAln:
'    - Prints out an alignment
'  o fun-04 initAlnST:
'    - Initalize all values in alnST to 0
'  o fun-05 freeAlnST:
'    - Frees alnST and all variables in alnST
'  o fun-06 pEMBOSSHead:
'    - Prints out the EMBOSS header to a file
'  o fun-07 pExpandCigHead:
'    - Prints out the expanded cigar header entry
'  o fun-08 capIdLen:
'    - Caps id length in seqST first white space or the
'      character at maxIdLenI - 1.
'  o fun-09 addPosToBuff:
'    - Adds a base position to a buffer. This uses the
'      printing format and position as a guide to
'      determine if printing.
'  o fun-10 eqxAddGap
'    - Adds a gap to an eqx buffer
'  o fun-11 eqxAddSMask
'    - Adds an soft mask to an eqx buffer
'  o fun-12 eqxAddSnp
'    - Adds an SNP entry to an eqx buffer
'  o fun-13 eqxAddMatch
'    - Adds an match entry to an eqx buffer
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef ALNSTRUCT_H
#define ALNSTRUCT_H

#include <time.h> /*For the EMBOSS format*/
#include <string.h>

#include "seqStruct.h"
#include "scoresST.h"
#include "generalAlnFun.h"
#include "alnMatrixStruct.h"

#define defEndAlnFlag 0

/*Keep these flag values between 1 and 3. This is so these
` flags can be used in twoBitArrays. Avoid using 0, because
` calloc will initialize all values to 0.
*/
#define defGapFlag 1
#define defSnpFlag 2
#define defMatchFlag 3

/*Not used by two bit arrays, can be set to any value,
` except 0
*/
#define defSoftMaskFlag 4

/*--------------------------------------------------------\
| Struct-01: alnStruct
|  - Holds the alignment array, which tells if a match,
|    SNP, del, ins, or masking was at a postion
\--------------------------------------------------------*/
typedef struct alnStruct
{ // alnStruct
  char *refAlnStr;
    /*Tells If each reference base is a match,snp,or gap*/
  char *qryAlnStr;
    /*Tells if each query base is a match, snp, or gap*/

  /*Length of the reference and query*/
  unsigned long refLenUL;
  unsigned long qryLenUL;
  unsigned long lenAlnUL;/*Length of alingnment with gaps*/

  /*Starting & ending position of alignment (index 0)*/
  unsigned long refStartAlnUL;
  unsigned long refEndAlnUL;
  unsigned long qryStartAlnUL;
  unsigned long qryEndAlnUL;

  /*This is for reporting similarity or other stats*/
  unsigned long numInssUL;   /*Number of insertions*/
  unsigned long numDelsUL;   /*Number of deletions*/
  unsigned long numSnpsUL;  /*Number of snps in alignment*/
  unsigned long numMatchesUL;/*Number of matches*/
}alnStruct;

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o Heap alloacted C-string with alignment for the
|      input sequence
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
char * alnSTToSeq(
    struct seqStruct *seqST,/*Has sequence to work with*/
    char qryBl,             /*1: working on query; 0: ref*/
    struct alnStruct *alnST,/*Has alignment array*/
    char extAlnRegionBl     /*1: Only keep aligned region*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: alnSTToSeq
   '  - Makes an alignment for an single sequence
   '    (reference or query)
   '  o fun-01 Sec-01:
   '    - Variable declerations 
   '  o fun-01 Sec-02:
   '    - Allocate memory & identify if ref or query seq
   '  o fun-01 Sec-03:
   '    - Add softmasking to the start
   '  o fun-01 Sec-04:
   '    - Make the aligned sequence
   '  o fun-01 Sec-05:
   '    - Add soft masking to the end
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnStruct with the alingment array
|    o 0 if had memory allocation error
\--------------------------------------------------------*/
struct alnStruct * dirMatrixToAlnST(
    struct seqStruct *refST,     /*Reference seq & length*/
    struct seqStruct *qryST,     /*Query seq & length*/
    struct scoresStruct *scoreST,
    /*Matrix to use in finding the alignment*/
    #if !defined BYTEMATRIX
       struct twoBitAry *dirMatrixST
    #else
       char *dirMatrixST
    #endif
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: getAlnAry
   '  - Builds an anlignment array for the input direction
   '    matrix
   '  o fun-02 sec-01:
   '    - VariAble declerations
   '  o fun-02 sec-02:
   '    - Assign memory and initalize variables
   '  o fun-02 sec-03:
   '    - Build the alignment array
   '  o fun-02 sec-04:
   '    - Find the starting and ending softmasked regions
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints out the input alingment
|  - Returns:
|    o 0 for success
|    o 1 If failed to get time for emboss format
|    o 64 for memory errors.
\--------------------------------------------------------*/
char printAln(
  FILE *outFILE,            /*File to print alingment to*/
  char *outStr,             /*print file name to header*/
  struct seqStruct *refST,  /*reference Id & sequence*/
  struct seqStruct *qryST,  /*query Id & sequence*/
  struct alnStruct *alnST,  /*Has alignment to print out*/
  long scoreL,              /*Score of the alignment*/
  struct alnSet *settings,  /*Settings used for alingment*/
  char *scoreMtxFileStr     /*printing out scoring matrix*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printAln
   '  - Prints out the alignment to an alignment file
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Print out the header for the alignment
   '  o fun-03 sec-03:
   '    - Allocate memory and copy read ids
   '  o fun-03 sec-04:
   '    - Final prep before printing
   '  o fun-03 sec-05:
   '    - Print out the alignment
   '  o fun-03 sec-06:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o All variables in alnST to be 0
| Note:
|  - This does not free alnAryUc, so only call this for
|    new alnST structures
\--------------------------------------------------------*/
void initAlnST(
  struct alnStruct *alnST // Strucutre to initialize
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-01: initAlnST
   '  - Initalize all values in alnST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Frees
|    o alnST, including all variables in in alnST
|    0 Or if heapBl is 0; frees all variables in alnST, but 
|      does not free alnST
\--------------------------------------------------------*/
void freeAlnST(
  struct alnStruct *alnST, // Strucutre to free
  char heapBl // 0: free variables in alnST, but keep alnST
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-01: freeAlnST
   '  - Frees alnST and all variables in alnST
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Prints:
|    o The EMBOSS format headers to a file
|  - Returns:
|    o 1 if their was an error when getting the time
|    o 0 for no errors
\--------------------------------------------------------*/
char pEMBOSSHead(
  FILE *outFILE,          /*File to print the header to*/
  struct seqStruct *refST,/*Has reference sequence id*/
  struct seqStruct *qryST,/*Has query sequence id*/
  long scoreL,           /*Score of the alignment*/
  char *outStr,          /*Name of the output file*/
  char *scoreMtxFileStr, /*Name of file with score matrix*/ 
  struct alnSet *setting,/*Settings for the alignment*/
  struct alnStruct *alnST, /*Results from aligment*/
  char pFileHeaderBl
    /*1: Print the file header and entry header
    ` 0: Just print the entry header*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: pEMBOSSHead
   '  - Prints out the EMBOSS header to a file
   '  o fun-06 sec-01:
   '    - Variable declerations (& get time)
   '  o fun-06 sec-02:
   '    - Print out the file header
   '  o fun-06 sec-03:
   '    - Print out the entry header
   \*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o The expanded cigar header to outFILE
\--------------------------------------------------------*/
void pExpandCigHead(
  FILE *outFILE,          /*File to print header to*/
  struct seqStruct *refST,/*Has reference sequence id*/
  struct seqStruct *qryST,/*Has query sequence id*/
  long scoreL,           /*Score of the alignment*/
  char *scoreMtxFileStr, /*Name of file with score matrix*/ 
  struct alnSet *setting,/*Settings for the alignment*/
  struct alnStruct *alnST/*Results from aligment*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: pExpandCigHead
   '  - Prints out the expanded cigar header entry
   '  o fun-07 sec-01:
   '    - Print out the query and reference information
   '  o fun-07 sec-02:
   '    - Print out the settings
   '  o fun-07 sec-03:
   '    - Print out the alignment stats
   '  o fun-07 sec-04:
   '    - Print out the legend
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o seqST->idCStr to end at the first white space or
|      at index maxIdLenI - 1 bases
|    o oldCharC to hold the character that was changed to
|      '\0'
|  - Returns:
|    o New end of seqST->idCStr (position changed to '\0')
\--------------------------------------------------------*/
char * capIdLen(
  struct seqStruct *seqST, /*Has id to limit length on*/
  int maxIdLenI,           /*Max id length (index 1)*/
  char *oldCharC    /*Will hold character changed to '\0'*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: capIdLen
   '  - Caps id length in seqST first white space or the
   '    character at maxIdLenI - 1.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o buffStr to hold the new base position if the
|      settings support position printing.
|  - Returns:
|    o Pointer to chacter after added base position or 
|      input buffStr address if buffStr was not changed.
\--------------------------------------------------------*/
char * addPosToBuff(
   char *buffStr,          /*Buffer to add position to*/
   unsigned long basePosUL,/*Position to add to buffer*/
   struct alnSet *settings,
     /*Has boleans to tell if printing first base and
     ` format printing to
     */
   char startBl            /*1: is first base in buffer*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: addPosToBuff
   '  - Adds a base position to a buffer. This uses the
   '    printing format and position as a guide to
   '    determine if printing.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o eqxBuffStr to have the correct gap entry based on
|      formatFlag.
\--------------------------------------------------------*/
void eqxAddGap(
   char insBl,         /*1:insertion; 0: deletion)*/
   char *eqxBuffStr,   /*Buffer with eqx entry*/
   char formatFlag     /*Has format to add gap as*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: eqxAddGap
   '  - Adds a gap to an eqx buffer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o eqxBuffStr to have the correct soft mask entry based
|      on formatFlag.
\--------------------------------------------------------*/
void eqxAddSMask(
   char *eqxBuffStr,   /*Buffer with eqx entry*/
   char formatFlag     /*Has format to add soft mask as*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: eqxAddSMask
   '  - Adds an soft mask to an eqx buffer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o eqxBuffStr to have the correct SNP entry based on
|      formatFlag.
\--------------------------------------------------------*/
void eqxAddSnp(
   char *eqxBuffStr,   /*Buffer with eqx entry*/
   char formatFlag     /*Has format to add soft mask as*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: eqxAddSnp
   '  - Adds an SNP entry to an eqx buffer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o eqxBuffStr to have the correct match entry based on
|      formatFlag.
\--------------------------------------------------------*/
void eqxAddMatch(
   char *eqxBuffStr,   /*Buffer with eqx entry*/
   char formatFlag     /*Has format to add soft mask as*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: eqxAddMatch
   '  - Adds an match entry to an eqx buffer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
