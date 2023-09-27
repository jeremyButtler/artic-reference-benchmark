/*######################################################################
# Use:
#   o Holds functions to build consensuses from fastq files or work
#     with consensuses
# Includes:
#   - "trimSam.h"
#   - "fqAndFaFun.h"
#   - "readExtract.h"
#   - "dataTypeShortHand.h"           (No .c file)
#   o "defaultSettings.h"
#   o "cStrFun.h"                     (No .c file)
#   o "cStrToNumberFun.h"             (No .c file)
#   o "printError.h"
#   o "samEntryStruct.h"              (No .c file)
#   o "findCoInftBinTree.h"
#   o "findCoInftChecks.h
#   o "scoreReadsFun.h"
#   o "fqGetIdsSearchFq.h"
#   o "fqGetIdsFqFun.h"
#   o "fqGetIdsStructs.h"
#   o "fqGetIdsHash.h"
#   o "fqGetIdsAVLTree.h"
# C standard Libraries:
#   o <time.h>
#   o <string.h>
#   o <stdlib.h>
#   o <stdio.h>
#   o <stdint.h>
######################################################################*/

#ifndef BUILDCONFUN_H
#define BUILDCONFUN_H

#include "trimSam.h"
#include "readExtract.h"
#include "fqAndFaFun.h"
#include "dataTypeShortHand.h"

/*---------------------------------------------------------------------\
| Struct-1: majCon
| Use:
|    - Holds information for building a majority consensus.
\---------------------------------------------------------------------*/
typedef struct majConStruct
{ /*majConStruct*/
   unsigned char useMajConBl; /*1: use the majority consensus step*/ 
   unsigned char minBaseQUC; /*Min q-score needed to keep an SNP/match*/ 
   unsigned char minInsQUC;  /*Min q-score needed to keep an insertion*/ 
   float minReadsPercBaseFlt;/*Min % of supporting reads to keep base*/
   float minReadsPercInsFlt; /*Min % of supporting reads to keep ins*/
   unsigned long lenConUL;   /*Holds length of ouput consensus*/
}majConStruct;

/*---------------------------------------------------------------------\
| Struct-2: raconStruct
| Use:
|    - Holds settings for Racon
\---------------------------------------------------------------------*/
typedef struct raconStruct
{ /*raconStruct*/
   unsigned char useRaconBl;  /*1: use racon consensus step*/ 
   unsigned char rndsRaconUC;/*Number of rounds to polish with racon*/ 
   unsigned long lenConUL;   /*Holds length of ouput consensus*/
}raconStruct;

/*---------------------------------------------------------------------\
| Struct-3: medakaStruct
| Use:
|    - Holds settings for Medaka
\---------------------------------------------------------------------*/
typedef struct medakaStruct
{ /*medakaStruct*/
   unsigned char useMedakaBl; /*1: use Medaka consensus step*/ 
   char modelCStr[64];        /*Model to use with medaka*/ 
   char condaBl;              /*1: use conda install, else python env*/
   unsigned long lenConUL;   /*Holds length of ouput consensus*/
}medakaStruct;

/*--------------------------------------------------------\
| Struct-4: ivarSet
| Use:
|  - Holds settings for ivar
\--------------------------------------------------------*/
typedef struct ivarSet
{ /*ivarSet*/
   char minSupStr[8];    /*Min support for snps*/
   char minInsSupStr[8]; /*Min support for insertions*/
   char minDepthStr[16]; /*Min read depth*/
   char maskC;           /*Base to mask with*/
   char minQStr[4];      /*Min Q-score to keep a base*/
   ulong lenConUL;       /*Length of consensus*/
   char useIvarBl;       /*Lets me know if ivar is used*/
}ivarSet;

/*--------------------------------------------------------\
| Name: initIvarSet
| Macro-01:
| Use:
|  - Initializes the settings for an ivarSet structuer
| Input:
|  - ivarSetST:
|    o ivarSet strucuture (not pointer) to initialize
| Output:
|  - Initializes
|    o All variables in ivarSetST
\--------------------------------------------------------*/
#define initIvarSet(ivarSetST){ \
   strcpy((ivarSetST).minSupStr, defIvarMinSupStr); \
   strcpy((ivarSetST).minInsSupStr, defIvarMinInsSupStr);\
   strcpy((ivarSetST).minDepthStr, defIvarMinDepthStr); \
   (ivarSetST).maskC = defIvarMask; \
   strcpy((ivarSetST).minQStr, defIvarMinQStr); \
   (ivarSetST).useIvarBl = (defConMethod == defUseIvar);\
} /*initIvarSet*/

/*---------------------------------------------------------------------\
| Struct-5: conBuildStruct
| Use:
|    - Holds settings for consensus building
\---------------------------------------------------------------------*/
typedef struct conBuildStruct
{ /*condBuildStruct*/
    uchar methodAryUC[32];/*Order to run methods*/
    uchar lenMethodUC;   /*Number consensus methods input*/

    char useStatBl;
         /*1: Use a stats file instead of read median-Q
         ` for finding a good read to build with
         */
    uchar clustUC; /*Cluster number of consensus*/
    uint numRndsToPolishUI;
       /*This is used to determine the number of times
       ` to rebuild the consensus.
       */
    float minPercMappedReadsFlt;
       /*Min % of total reads that mapped. This is used to
       ` determine if a consensus should be built when
       ` the subsample depth is smaller than the reads.
       */
    ulong minReadsToBuildConUL; /*Min reads for consensus*/
    ulong maxReadsToBuildConUL; /*max read subsample size*/
    ulong numReadsForConUL;
        /*Number of reads kept to buld a consensus.*/

    uint minConLenUI; /*Min consensus length*/
    ulong lenConUL;   /*length of ouput consensus*/

    /*Settings for the consensus building step*/
    struct majConStruct majConSet;
    struct raconStruct raconSet;
    struct medakaStruct medakaSet;
    struct ivarSet ivarSetST;
}condBuildStruct;

/*---------------------------------------------------------------------\
| Struct-6: baseStruct
| Use:
|    - Holds the base, error type (if insertion or deletion), & a
|      counter for number of reads with the same base.
\---------------------------------------------------------------------*/
typedef struct baseStruct
{ /*baseStruct*/
    char baseChar;     /*What is the base (0 for deletion)*/
    char errTypeChar;  /*0 = match, 1 = SNP, 2 = insertion*/
    
   /*Recored number of reads that supported this base (good quality)*/
   unsigned long numSupReadsUL;

   struct baseStruct *altBase; /*Alternate options for the base*/
   struct baseStruct *nextBase; /*For linked lists*/
}baseStruct;

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o 1 if built a consensus
|     o 2 if could not open a file (refs, stats?, or fastq)
|     o 8 if stat file was specified, but no stat file provided
|     o 16 if could not build a consensus, but no other errors
|     o 64 for memory allocation error
|   - File:
|     o consensus named fastqFileName--cluster-0--consensus.fasta
\---------------------------------------------------------------------*/
unsigned char buildCon(
    struct readBin *conData,
        /*Has fastq and other files needed to build the consensus
          WARNING: The fastq & stats file in conData will be modified*/
    char *refCStr,           /*Path to fasta file with reference*/
    char *threadsCStr,       /*Number threads to use with system calls*/
    struct conBuildStruct *conSet,   /*settings for building consensus*/
    struct samEntry *samStruct,      /*Will hold sam file data*/
    struct samEntry *bestReadSam,    /*For read median Q extraction*/
    struct minAlnStats *minReadReadStats,
        /*Minimum stats to keep a read/read mapping*/
    struct minAlnStats *minReadConStats
        /*Minimum stats needed to keep a read/consensus mapping*/
); /*Builds a consensus using input fastq file & best read or referece*/

/*---------------------------------------------------------------------\
| Output:
|    o Creates:
|      - Fasta file with the consensus
|    o Modifies:
|      - clustOn->consensusCStr to have file name of created consensus 
|    o Returns:
|      - 1 for success
|      - 2 or 4 for file errors
|      - 16 for minimap2 error
|      - 64 for memory allocation error
\---------------------------------------------------------------------*/
unsigned char buildSingleCon(
    char *threadsCStr,       /*Number threads to use with system calls*/
    struct readBin *clustOn, /*fastq file & reads to build consensus*/
    struct samEntry *samStruct,/*For reading in sequences or sam file*/
    struct conBuildStruct *conSet
        /*Settings to use while building the consensus*/
); /*Builds a consensus using the best read & top reads in clustOn*/

/*---------------------------------------------------------------------\
| Output:
|    - Returns:
|        o 1 If succeded
|        o 2 if could not open or read in the reference sequence
|        o 4 if could not open the bestReads file
|        o 16 if had to few sequences map to the selected read
|        o 32 if minimap2 crashed
|        o 64 for memory allocation errors
\---------------------------------------------------------------------*/
unsigned char simpleMajCon(
    unsigned char *clustUC,    /*Cluster number to assign to consensus*/
    char *threadsCStr,         /*Number of threads to use with minimap*/
    struct readBin *binStruct,  /*Has best read & top reads files*/
    struct samEntry *samStruct, /*For reading in sam file entries*/
    struct majConStruct *settings
        /*Has settings for building the consensus*/
); /*Builds a majority consensus from the best reads & top read*/

/*---------------------------------------------------------------------\
| Output:
|    Uses: Racon to build a consensus (file name in bin->consensusCStr)
\---------------------------------------------------------------------*/
void buildConWithRacon(
    struct raconStruct *settings,   /*Settings to use with racon*/
    unsigned char *clustUC,      /*Cluster on*/
    char *threadsCStr,              /*Number threads to use with racon*/
    struct readBin *conBin          /*Bin working on*/
); /*Builds a consensus using racon*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        binTree->consensusCStr to hold the polished consensus
|    Returns:
|        - 1 If built a consensus
|        - 2 If input file does not exists
|        - 4 IF consensus not built
\---------------------------------------------------------------------*/
unsigned char medakaPolish(
    struct medakaStruct *settings, /*Settins for running medaka*/
    unsigned char *clustUC,        /*Cluster on*/
    char *threadsCStr,             /*Number threads to use*/
    struct readBin *conBin,        /*bin with consensus & top reads*/
    struct samEntry *samStruct     /*For converting fastq to fasta*/
); /*Polish a consensus with medaka using the best reads*/

/*---------------------------------------------------------------------\
| Output:
|    - Returns:
|        - readBin Struct with the closest consensus
|    - Modifies:
|        - closesConUint to hold the score of the most similar consensus
| Note:
|    - clusters with ->balUChar < 0 will be ignored
\---------------------------------------------------------------------*/
struct readBin * cmpCons(
    struct readBin *conBin,       /*Bin with consensus to compare*/
    struct readBin *conBinTree,   /*Tree of consensus to compare to*/
    struct samEntry *samStruct,   /*Struct to hold input from minimap2*/
    struct samEntry *refStruct,   /*Struct to hold input from minimap2*/
    struct minAlnStats *minStats, /*Min stats needed to keep a error*/
    char *threadsCStr            /*Number threads to use with Minimap2*/
); /*Compares a consenses to a another consensus. This will do a
    recursive call if conBinTree has children*/

/*---------------------------------------------------------------------\
| Output:
|   o Returns:
|     - The next base in the list (if altBase != 0, returns altBase)
|   o Frees:
|     - The baseToFree from memory
|   o Modifies:
|     - If their is an alternative base (baseToFree->altBase != 0)
|         o lastBase->nextBase is set to baseToFree->altBase
|         o bastToFree is set to bastToFree->altBase
|     - If their is not alternative base (bastToFree->altBase == 0)
|         o lastBase->nextBase is set to baseToFree->nextBase
|         o bastToFree is set to 0
| Note:
|    o This function assumes that the next base pointer (nextBase) is
|      is alwasy 0 (not set) for teh alternate base pointer (altBase==0)
\---------------------------------------------------------------------*/
struct baseStruct * freeBaseStruct(
    struct baseStruct **baseToFree, /*Insertion list to free*/
    struct baseStruct *lastBase     /*Base to assign pointers to*/
); /*Frees an base from a linked list of bases*/

/*---------------------------------------------------------------------\
| Output: Modifies: majConStruct to have default settings
\---------------------------------------------------------------------*/
void initMajConStruct(
    struct majConStruct *majConSettings
    /*struct to set to default values in defaultSettings.h*/
); /*Sets input structers variables to default settings*/

/*---------------------------------------------------------------------\
| Output: Modifies: raconStruct to have default settings
\---------------------------------------------------------------------*/
void initRaconStruct(
    struct raconStruct *raconSettings
    /*struct to set to default values in defaultSettings.h*/
); /*Sets input structers variables to default settings*/

/*---------------------------------------------------------------------\
| Output: Modifies: medakaStruct to have default settings
\---------------------------------------------------------------------*/
void initMedakaStruct(
    struct medakaStruct *medakaSettings
    /*struct to set to default values in defaultSettings.h*/
); /*Sets input structers variables to default settings*/

/*---------------------------------------------------------------------\
| Output: Modifies: medakaStruct to have default settings
\---------------------------------------------------------------------*/
void initConBuildStruct(
    struct conBuildStruct *consensusSettings
    /*struct to set to default values in defaultSettings.h*/
); /*Sets input structers variables to default settings*/

/*--------------------------------------------------------\
| Name: ivarCon (Fun-13:)
| Use:
|  - Builds a consensus using ivar
| Input:
|  - settings:
|    o ivarSet structure with settings for ivar
|  - clustUC:
|    o Number of the cluster
|  - threadsStr:
|    o Number of threads to use (as string)
|  - conST:
|    o Has names of the consensus, best read, and top reads
|      fastq files.
|    o best read is the read to use in making the consensus
|      - fastq or fasta
|    o top reads are aligned to best read and used to build
|      the consensus
|    o consensus is the name of the output consensus
|      - Avoid periods in names
|  - samST:
|    o Used to read in sam entries.
| Output:
|    - Returns:
|        o 1 If succeded
|        o 2 for pipe error with samtools/ivar
|        o 4 Failed to make a consensus
|        o 32 for pipe error with minimap2
|        o 64 memory error
\--------------------------------------------------------*/
unsigned char ivarCon(
    struct ivarSet *settings,/*Settings to run ivar with*/
    uchar *clustUC,          /*Cluster number to use*/
    char *threadsStr,        /*Number of threads*/
    struct readBin *conST,   /*Has reads/consensus to use*/
    struct samEntry *samST  /*For getting sam file input*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: ivarCon
   '  - Builds a consensus using ivar from the best reads & top read
   '  o fun-13 sec-01:
   '    - Variable declerations
   '  o fun-13 sec-02:
   '    - Make a temporary file for trimming
   '  o fun-13 sec-03:
   '    - Prepare command to call ivar
   '  o fun-13 sec-04:
   '    - Run ivar
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
