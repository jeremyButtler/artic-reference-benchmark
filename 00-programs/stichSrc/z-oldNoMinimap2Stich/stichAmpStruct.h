/*########################################################\
# Name: stichStruct
# Use:
#  - Holds the stichAmpST structer and its supporting
#    functions.
# Libraries:
# C Standard Libraries:
#  - <stdlib.h>
\########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start of Header
'  o struct-01 stichAmpST
'    - Holds the consensus sequence for the stiched
'      together amplicons.
'  o macro-01: initStichAmpST
'    - Initializes a stichAmpST structure to 0's
|  o macro-02: makeStichAmpST
|    - Makes an stich amp structure
'  o macro-03: freeStichAmpST
'    - Frees an stich amp structure
'  o macro-04: freeStichAmpSTList
'    - Frees a list of stich amp structures
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef STICHAMPSTRUCT_H
#define STICHAMPSTRUCT_H

#include <stdlib.h>

#define defMvMask 4 /*Error type is masked*/
  /*Using this format to reflect alnSeq*/

/*--------------------------------------------------------\
| Struct-01: stichsT
|  - Holds the consensus sequence for the stiched together
|    amplicons.
|  - This is set up to be a jagged matrix. Each base pionts
|    to the next base, the previous base, and a list of
|    alternative bases for that position. Each alt base
|    should only be for that position. Otherwise you will
|    beak the free functions.
\--------------------------------------------------------*/
typedef struct stichAmpST
{ /*stichAmpST*/
   char baseC;   /*Consensus base at position*/
   char errC;    /*Error type*/
   ulong supportUL; /*How much support?*/
   ulong depthUL;  /*How many amplicons had this region*/

   struct stichAmpST *nextBase;
   struct stichAmpST *lastBase;
   struct stichAmpST *altBase;
}stichAmpST;


/*--------------------------------------------------------\
| Macro-01: initStichAmpST
| Use:
|  - Initializes a stichAmpST structure to 0's
| Call: initStichAmpST(baseST);
| Input:
|  - baseST:
|    o Pointer to stichAmpST structure to initialize
| Output:
|  - Modifies:
|    o baseST to have all values set to defaults.
\--------------------------------------------------------*/
static inline void initStichAmpST(
   struct stichAmpST *baseST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Macro-01: initStichAmpST
   '  - Initializes a stichAmpST structure to 0's
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   baseST->baseC = 0;
   baseST->errC = 0;
   baseST->supportUL = 0;
   baseST->depthUL = 0;

   baseST->nextBase = 0;
   baseST->lastBase = 0;
   baseST->altBase = 0;
} /*initStichAmpST*/

/*--------------------------------------------------------\
| Macro-02: makeStichAmpST
| Use:
|  - Makes an stich amp structure
| Call: stichPtr = makeStichAmpST();
| Input:
| Output:
|  - Returns:
|    o An pointer to an initalized stichAmpST structer
|    o 0 For memory errors
\--------------------------------------------------------*/
static inline struct stichAmpST *  makeStichAmpST(
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Macro-02: makeStichAmpST
   '  - Makes an initialized stichAmpST structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   struct stichAmpST *retST = 0;
   retST = malloc(sizeof(struct stichAmpST));

   if(retST == 0) return 0;

   initStichAmpST(retST);
   return retST;
} /*makeStichAmpST*/

/*--------------------------------------------------------\
| Macro-03: swapStichAmpST
| Use:
|  - Swaps the values, but not pointers in two stichAmp
|    structs
| Call: swapStichAmpST(oneST, twoST);
| Input:
|  - oneST:
|    o Pointer to first stichAmpST structure to swap
|  - twoST:
|    o Pointer to second stichAmpST structure to swap
| Output:
|  - Modifies:
|    o oneST to have the values of twoST
|    o twoST to have the values of oneST
\--------------------------------------------------------*/
static inline void swapStichAmpST(
   struct stichAmpST *oneST,
   struct stichAmpST *twoST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Macro-03: swapStichAmpST
   '  - Frees a stichAmpST structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   ulong swapUL = 0;
   char swapC = 0;

   swapC = oneST->baseC;
   oneST->baseC = twoST->baseC;
   twoST->baseC = swapC;

   swapC = oneST->errC;
   oneST->errC = twoST->errC;
   twoST->errC = swapC;

   swapUL = oneST->supportUL;
   oneST->supportUL = twoST->supportUL;
   twoST->supportUL = swapUL;

   return;
} /*swapStichAmpST*/

/*--------------------------------------------------------\
| Macro-04: freeStichAmpST
| Use:
|  - Frees an stich amp structure
|  - This will move the alternative base down if there is
|    one.
| Call: freeStichAmpST(&baseST);
| Input:
|  - baseST:
|    o Pointer to stichAmpST structure to free
| Output:
|  - Frees:
|    o baseST from memory
|  - Sets:
|    o baseST to 0 if there are no alternative bases
|    o Puts the alternative base in baseST's position if
|      there are altnerative bases.
\--------------------------------------------------------*/
static inline void freeStichAmpST(
   struct stichAmpST **baseST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Macro-04: freeStichAmpST
   '  - Frees a stichAmpST structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   struct stichAmpST *tmpST = 0;

   if((*baseST)->altBase != 0)
   { /*If I have alternative bases in the swap*/
      swapStichAmpST(*baseST, (*baseST)->altBase);
      tmpST = (*baseST)->altBase;
      (*baseST)->altBase = tmpST->altBase;
      free(tmpST);
   } /*If I have alternative bases in the swap*/

   else
   { /*Else: there are not alternative bases*/
      if((*baseST)->nextBase!=0 && (*baseST)->lastBase!=0)
      { /*If: this is in the middle of a list*/
         (*baseST)->nextBase->lastBase=(*baseST)->lastBase;
         (*baseST)->lastBase->nextBase=(*baseST)->nextBase;
      } /*If: this is in the middle of a list*/

      free(*baseST);
      *baseST = 0;
   } /*Else: there are not alternative bases*/

   return;
} /*freeStichAmpST*/

/*--------------------------------------------------------\
| Macro-05: freeStichAmpSTList
| Use:
|  - Frees a list of stich amp structures
| Call: freeStichAmpSTList(&baseST);
| Input:
|  - baseST:
|    o Pointer to stichAmpST list to free
| Output:
|  - Frees:
|    o Every stichAmpST structer in baseListST from memory
|  - Sets:
|    o baseListST to 0.
\--------------------------------------------------------*/
static inline void freeStichAmpSTList(
   struct stichAmpST **baseListST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Macro-05: freeStichAmpSTList
   '  - Frees a list of stichAmpST structures
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   struct stichAmpST *iterST = 0;
   struct stichAmpST *nextBaseST = 0;

   if(iterST ==  0) return;

   iterST = *baseListST;

   while(iterST != 0)
   { /*While I have a next bases to free*/
      nextBaseST = iterST->nextBase;

      /*Free all alternative bases*/
      while(iterST != 0) freeStichAmpST(&iterST);
      
      iterST = nextBaseST;
   } /*While I have a next bases to free*/

   iterST = *baseListST;

   while(iterST != 0)
   { /*While I have a next bases to free*/
      nextBaseST = iterST->lastBase;

      /*Free all alternative bases*/
      while(iterST != 0) freeStichAmpST(&iterST);
      
      iterST = nextBaseST;
   } /*While I have a next bases to free*/

   while(*baseListST != 0) freeStichAmpST(baseListST);

   return;
} /*freeStichAmpSTList*/

#endif
