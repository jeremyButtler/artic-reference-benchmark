/*#########################################################
# Name: stichSetStruct
# Use:
#  - Holds the structure to hold the settings for stich
# Libraries:
#  - "samFunSrc/dataTypeShortHand.h"   (No .c file)
#  - "stichDefaults.h"                 (No .c file)
# C Standard Libraries:
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' stichSetStruct SOH: Start Of Header
'  o struct-01: stichSet
'    - Holds settings for stiching together amplicons
'  o macro-01: initStichSet
'    - Initializes a stichSet struct with default settings
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef STICHSETSTRUCT_H
#define STICHSETSTRUCT_H

#include "samFunSrc/dataTypeShortHand.h"
#include "stichDefaults.h"

/*--------------------------------------------------------\
| Struct-01: stichSet
|  - Holds settings for stiching together amplicons
\--------------------------------------------------------*/
typedef struct stichSet{
   char maskC;       /*Mask to apply to low support bases*/
   char useMinimap2Bl; /*use minimap2*/
   uchar overwriteBl;  /*Overwrite output file*/
   uchar threadsUC;
     /*The miniamp2 boolean is not used currently, but will
     ` be used if I ever get the alnSeq mem-water working
     ` correctly
     */
   ulong minDepthUL; /*Min depth to do voting (not 100%)*/
   ulong minSupportUL; /*Min support to keep disagreement*/
}stichSet;

/*--------------------------------------------------------\
| Macro-01: initStichSet
| Use:
|  - Initializes a stichSet structure with default settings
| Input:
|  - stichSetPtr
|    o Pointer to a stichSet structuer to initialize
| Output:
|  - Modifies:
|    o All values in stichSetPtr to be default values
|      (see stichDefaults.h for defaults)
\--------------------------------------------------------*/
#define initStichSet(stichSetPtr){\
   (stichSetPtr)->maskC = defMaskBase; \
   (stichSetPtr)->minDepthUL = defMinStichDepth; \
   (stichSetPtr)->minSupportUL = defMinStichSup; \
   (stichSetPtr)->useMinimap2Bl = defMinimap2Bl; \
   (stichSetPtr)->threadsUC = defThreads; \
   (stichSetPtr)->overwriteBl = defStichOverwrite; \
} /*initStichSet*/

#endif
