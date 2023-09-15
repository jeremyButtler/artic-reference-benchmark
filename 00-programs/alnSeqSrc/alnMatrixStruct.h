/*#########################################################
# Name: alnMatrixStruct
# Use:
#  - Holds the alnMatrix struct and its functions
#  - alnMatrix is used to store the results of an alignment
# Libraries:
#  - "scresST.h"
#  - "twoBitArrays.h"
# C Standard Libraries:
#  o <stdint.h>
#  o <stdlib.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o struct-01: alnMatrixStruct
'    - Holds the direction matrix and best score(s) for a
'      single aligment
'  o fun-01 initAlnMatrixST:
'    - Sets all variables in matrixST to 0
'  o fun-02 freeAlnMatrixST
'    - Sets all variables in matrixST to 0
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef ALNMATRIXSTRUCT_H
#define ALNMATRIXSTRUCT_H

#include "scoresST.h"
#include "twoBitArrays.h"

/*--------------------------------------------------------\
| Struct-01: alnMatrixStruct
|  - Holds the direction matrix and best score(s) for a
|    single aligment
\--------------------------------------------------------*/
typedef struct alnMatrixStruct
{ /*alnStruct*/
  struct scoresStruct bestScoreST;

  #if !defined BYTEMATRIX
     struct twoBitAry *dirMatrixST;
  #else
     char *dirMatrixST;
  #endif

  struct scoresStruct *refBasesST;
  unsigned long lenRefScoresUL;

  struct scoresStruct *qryBasesST;
  unsigned long lenQueryScoresUL;
}alnMatrixStruct;

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o All variables in matrixST to be 0
\--------------------------------------------------------*/
void initAlnMatrixST(
  struct alnMatrixStruct *matrixST /*Struct to initialize*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-01: initAlnMatrixST
   '  - Sets all variables in matrixST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o Frees alnMatrix and all of its held variables
\--------------------------------------------------------*/
void freeAlnMatrixST(
  struct alnMatrixStruct *matrixST /*Struct to free*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-01: freeAlnMatrixST
   '  - Sets all variables in matrixST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
