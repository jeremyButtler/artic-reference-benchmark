/*#########################################################
# Name: seqStruct
# Use:
#  o Holds functions for reading in or manipulating
#    sequences
# Libraries:
# C standard libraries:
#  - <stdlib.h>
#  - <stdio.h>
#  - <stdint.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  - struct-01: seqStruct
'    - Holds sequence and length of a input sequence
'  - fun-01 reverseComplementSeq:
'     o Reverse complement a sequence
'  - fun-02 complementBase:
'     o Returns the complement of a base
'  - fun-03 readFqSeq:
'    o Reads a fastq sequence from a fastq file
'  - fun-04 readFaSeq:
'     o Grabs the next read in the fasta file
'  - fun-05 addLineToBuffSeqFun:
'     o Add characters from file to buffer, if needed 
'       resize. This will only read in till the end of the
'       line
'  - fun-06 reverseCStr;
'     o Reverse a c-string to be backwards (here for
'       Q-score entries)
'  o fun-07 freeSeqST:
'     - Frees the seqST strucuter
'  o fun-08 TOC: Sec-01: initSeqST
'    - Sets vlues in seqST to zero
'  o fun-09 addStartEndToSeqST:
'    - Sets the start and ending corrdinates of a region
'      of interest in a sequence
'  o fun-10 initSeqST:
'     - Sets values in seqST to blank values
'  o fun-11 cpReadIdRPad:
'     - Copies read id to a buffer and adds in endIdC to
'       the end. If needed, this function will add right
'       padding of spaces to the end.
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef SEQSTRUCT_H
#define SEQSTRUCT_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/*--------------------------------------------------------\
| Struct-01: seqStruct
|  - Holds sequence and length of a input sequence
\--------------------------------------------------------*/
typedef struct seqStruct
{ // refQueryStruct
  char *idCStr;          // Id of th sequence
  unsigned long  lenIdUL;      // Length of the sequence id
  unsigned long  lenIdBuffUL;  // Lenght of Id buffer

  char *seqCStr;          // Sequence
  unsigned long  lenSeqUL;      // Length of the sequence
  unsigned long  lenSeqBuffUL;  // Lenght of sequence buffer

  char *qCStr;           // q-score entry
  unsigned long  lenQUL;       // Length of the Q-score
  unsigned long  lenQBuffUL;   // Length of Q-score buffer

  unsigned long  offsetUL;     // Offset for an alignment
  unsigned long  endAlnUL;     // Marks end of alignment
}seqStruct;

/*--------------------------------------------------------\
| Output:
|  - Modfies
|    o seqCStr to be the reverse complement sequence
\--------------------------------------------------------*/
void reverseComplementSeq(
  struct seqStruct *seqST
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: reverseComplementSeq
   '  - Reverse complement a sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns
|     o complement of the input base (0 if invalid base)
\--------------------------------------------------------*/
char complementBase(
    const char *baseC
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-1 Sub-1: complementBase
   '  - Return the complement of a base
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o refStruct to hold the read in fastq entry & sets
|      its pointers
|  - Returns:
|     o 0: if EOF
|     o 1: if succeded
|     o 2: If file was not a fastq file
|     o 130: If file had an invalide entry
|       - This error uses to flags, first it uses 2 to
|         specify that it is not a fastq file
|       - 2nd it uses 128 to specifty that it is not an
|         blank file
|     o 64: If malloc failed to find memory
\--------------------------------------------------------*/
uint8_t readFqSeq(
  FILE *fqFILE,       // fastq file to grab sequence from
  struct seqStruct *seqST // Will hold one fastq entry
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: readRefFqSeq
   '  -  Grabs the next read in the fastq file
   '  o fun-03 sec-1:
   '    - Variable declarations
   '  o fun-03 sec-2:
   '    - Check if need to allocate memory for buffer
   '  o fun-03 sec-3:
   '    - Read in the first data
   '  o fun-03 sec-4:
   '    - If not at file end, see if have the full entry
   '  o fun-03 sec-5:
   '    -Read till end of file, check if valid fastq entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o seqST to hold one fasta entry
|  - Returns:
|     o 0: if EOF
|     o 1: if succeded
|     o 2: for an invalid fasta entry
|     o 64: If malloc failed to find memory
| Note:
|   - This will remove new lines from the sequence.
|   - This will only remove spaces or tabs at the end of
|     each sequence entry, so "atg atc \n" will got to
|     "atc atc".
\--------------------------------------------------------*/
uint8_t readFaSeq(
  FILE *faFILE,           // Fasta file with sequence
  struct seqStruct *seqST // Will hold one fastq entry
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: readFaSeq
   '  -  Grabs the next read in the fasta file
   '  o fun-04 sec-1:
   '    - Variable declarations
   '  o fun-04 sec-2:
   '    - Check if need to allocate memory for buffer
   '  o fun-04 sec-3:
   '    - Read in the sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|   - Modifies:
|     o buffCStr to hold the next line.
|       - buffCStr is resizied if it is to small to hold
|         the next line.
|       - buffCStr + lenBuffUL - 2 will be '\0' or '\n'
|       - buffCStr set to 0 for memory allocation errors
|     o curBuffUL: Has the number of chars in the buffer
|     o lenBuffUL: Has the buffer size
|     o inFILE: Points to next line in file
|   - Returns:
|     o 0 if was end of file (EOF)
|     o 1 if read in the next line
|     o 64 if had a memory allocation error
\--------------------------------------------------------*/
unsigned char addLineToBuffSeqFun(
    char **buffCStr,          // Buffer to add data to
    unsigned long  *lenBuffUL, // Size of the buffer
    unsigned long  *curBuffUL, // Number of chars in buffer
    unsigned long resBuffUL,  // How much to resize buff by
    FILE *inFILE              // File to grab data from
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: addLineToBuffSeqFun
   '  - Read line of characters into the buffer.If needed
   '    this will resize the buffer.
   '   o fun-05 sec-1:
   '     - variable declerations
   '   o fun-05 sec-2:
   '     - Check if need to resize the buffer
   '   o fun-05 sec-3:
   '     - Read in the next line in the buffer
   '   o fun-05 sec-4:
   '     - If at end of file, update read in lengths
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o inCStr to be backwards (end at start, start at end)
\--------------------------------------------------------*/
void reverseCStr(
  char *inCStr,       // C-string to refeverse
  unsigned long lenCStrUI//Length of input string (index 1)
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-1 Sub-1: reverseCStr
   '  - Reverse a c-string to be backwards
   '    (here for Q-score entries)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Frees
|    o seqST from memory
| Notes:
|  - You will have to set the pointer to seqST to 0
|  - Make sure you have a pointer to seqST->seqCStr, and
|    seqST->idCStr if you did set freeSeqBl or freeIdbl to
|    0. Otherwise you will loose your handle to the data
\--------------------------------------------------------*/
void freeSeqST(
  struct seqStruct *seqST, // Struct to free
  char heapBl     // 0: seqST on stack only free variables
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: Sec-01: freeSeqST
   '  - Frees the seqST strucuter
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o all values in seqST to be 0
\--------------------------------------------------------*/
void initSeqST(
  struct seqStruct *seqST // Struct to initialize
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: Sec-01: initSeqST
   '  - Sets vlues in seqST to zero
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o seqST to have the start and end corrdiantes of the
|      target region in the sequence
\--------------------------------------------------------*/
void addStartEndToSeqST(
  unsigned long  startTargetUI,
    // Start of region of intreset
  unsigned long  endTargetUI,
     // End of region of interest
  struct seqStruct *seqST // Struct to add corrdinates to
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: Sec-01: addStartEndToSeqST
   '  - Sets the start and ending corrdinates of a region
   '    of interest in a sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o Sets id, sequence, and Q-score entreis to start
|      with '\0' and the id, sequence, and Q-score lengths
|      to 0. This does not change the buffer lengths.
\--------------------------------------------------------*/
void blankSeqST(
  struct seqStruct *seqST // Struct to Blank
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: Sec-01: initSeqST
   '  - Sets values in seqST to blank values
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o buffStr to hold the sequence id + endIdC. If the
|      id + endIdC is shorter than padRI, the copied id is
|      padded with spaces on the right till it is the same
|      size as padRI.
|  - Returns:
|    o Pointer to end of copied id or padding if padding
|      was applied.
\--------------------------------------------------------*/
char * cpReadIdRPad(
   struct seqStruct *seqST, /*Has read id to copy*/
   char *buffStr,           /*Buffer to add read id to*/
   char endIdC,    /*Char to add to end of id (0 to skip)*/
   int padRI       /*Padding to add to right of id*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: cpReadIdRPad
   '  - Copies read id to a buffer and adds in endIdC to
   '    the end. If needed, this function will add right
   '    padding of spaces to the end.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
