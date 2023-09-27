/*#########################################################
# Name: samEntryStruct
# Use:
#  - Holds structer to hold a sam file entry. This also
#    includes the functions needed to support this
#    structer.
# Libraries
#  - "cStrToNumberFun.h"
# C Standard Libraries:
#  - <string.h>
#  - <stdlib.h>
#  - <stdio.h>
#  o <sdtint.h>
# Note:
#  - End of file has some general sam file basics
#########################################################*/

#ifndef SAMENTRYSTRUCT_H
#define SAMENTRYSTRUCT_H

#include <stdlib.h>          /*memory allocation*/
#include <string.h>
#include <stdio.h>
#include "cStrToNumberFun.h"

#define Q_ADJUST 33 /*offest to get q-score of 0*/
#define MAX_Q_SCORE 94 /*highest possible Q-score*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' samEntryStruct TOC: Table Of Contents
'  o struct-01 samEntry:
'    - Holds a single samfile entry
'  o struct-02 readStat:
'    - Holds query id, reference id, & stats for a read
'  o fun-01 blankSamEntry:
'    - Sets all non-alloacted variables in samEntryST to 0
'  o fun-02 initSamEntry:
'    - Initalize a samEntry struct to 0's
'  o fun-03 blankReadStats:
'    - modifes samEntryST to have every non-pointer
'      variable set to 0
'  o fun-04 deepCpSamReadStats:
'    - Copies the contents and buffer of one samEntry
'      struct to another samEntry struct.
'  o fun-05 cpSamEntry:
'    - Copies address of q-score & sequence pointers from
'      an old samEntry to an new SamEntry
'  o fun-06 freeStackSamEntry:
'    - Frees heap allocations in a stack allocated
'      samEntry struct
'  o fun-07 freeHeapSamEntry:
'    - Frees a samEntry structer (and sets to null)
'  o fun-08 processSamEntry:
'    - Sets Q-score, cigar, & sequence pointers and 
'      stats in a sam entry sotred in a samEntry struct.
'  o fun-09 readNextPartOfLine:
'    - Reads in next part of line when fgets did not
'      get a full line
'  o fun-10 findSamCig:
'    - Finds the cigar in same entry & also finds the
'      number of bases in query from cigar
'  o fun-11 readSamLine:
'    - Reads in sam entry into a samEntry struct
'  o fun-12 printSamEntry:
'    - Prints the sam file entry to a file. This does not
'      print any extra stats that were found.
'  o fun-13 samToFq:
'    - Prints sam entry as a fastq entry to a fastq file
'  o fun-14 printStatsHeader:
'    - Prints the stat file header made using a samEntry
'      struct
'  o fun-15 printSamStats
'    - Prints stats in samEntry structure to a tsv file
'  o fun-16 swapSamEntries:
'    - Swaps the values in two sam entries. This includes
'      histograms, so will take a bit of time.
'  o fun-17 sortScoresStartLen:
'    - Sorts an array of samEntry structs by reference
'      position (start) and aligned refence length with
'      shell short. Order is least to greatest.
'  o fun-18 readInConFa:
'    - Reads in a reference sequence from a fasta file
'      into a samEntry struct.
'    - I would recomend using seqStruct for this. This
'      was made before seqStruct & is kept here to make
'      sure nothing breaks
'  o fun-19 blankReadStat:
'    - Blanks (sets to 0) all values in a readStat struct
'  o fun-20 cpReadStat:
'    - Copies values from one readStat structure to
'      another read stat structuer
'  o fun-21 readStatsFileLine
'    - Reads in one line of read stats from a read stats
'      file into a readStat structure
'  o Fun-22 TOC: printReadStat
'    - Prints the stats in a read stat struct to an file
'  o fun-23 samEntryToReadStat:
'    - Copies stats in a samEntry structure to a
'      readStat structure
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Struct-01: samEntry
|  - Holds a single samfile entry
\--------------------------------------------------------*/
typedef struct samEntry
{ /*samEntry*/
    char *samEntryCStr;/*Holds c-string (entire sam entry)*/

    char *queryCStr;/*points to query id in samEntryCStr*/
    char *refCStr;  /*points to reference in samEntryCStr*/
    char *cigarCStr;/*points to cigar in samEntryCStr*/
    char *seqCStr;  /*points to seqeuence in samEntryCStr*/
    char *qCStr;    /*points to q-score in samEntryCStr*/

    uint8_t mapqUChar; /*Holds mapping quality*/

    uint16_t flagUSht; /*Holds the flag in samEntryCStr*/

    float medianQFlt;     /*holds median read Q-score*/
    float meanQFlt;       /*holds mean read Q-score*/
 
    /*The aligned Q-scores ignore low quality bases*/
    float medianAligQFlt; /*holds aligned median Q-score*/
    float meanAligQFlt;   /*holds mean aligned Q-score*/

    uint32_t posOnRefUInt;  /*First reference base*/
    uint32_t unTrimReadLenUInt;/*length of untrimmed read*/
    uint32_t readLenUInt;     /*Holds read length*/
    uint32_t readAligLenUInt; /*Aligned read length*/
    uint32_t numMatchUInt;    /*Holds number of matches*/
    uint32_t numKeptMatchUInt;
       /*number matches with Q-score > min Q*/
    uint32_t numKeptSNPUInt;  /*number of kept mismatches*/
    uint32_t numSNPUInt;     /*total number of mismatches*/
    uint32_t numKeptDelUInt; /*number of kept deletions*/
    uint32_t numDelUInt;     /*number of deletions*/
    uint32_t numKeptInsUInt; /*number of kept insertions*/
    uint32_t numInsUInt;     /*number of insertions*/

    uint32_t seqQHistUInt[MAX_Q_SCORE]; 
       /*Histogram of base Q-scores*/
    uint32_t seqQAlnHistUInt[MAX_Q_SCORE];
        /*Histogram of kept base Q-score (Qscore > minQ)*/


    uint64_t totalQScoreULng;
       /*Q-score of all bases added together*/
    uint64_t totalAlnQScoreULng;
       /*Q-score of kept bases added together*/

    unsigned long lenBuffULng;
       /*# bytes allocated to samEntryCStr*/
}samEntry;

/*--------------------------------------------------------\
| Struct-02: readStat
|  - Holds query id, reference id, & stats for a read
\--------------------------------------------------------*/
typedef struct readStat
{ /*readStat*/
    char queryIdCStr[128];
    char refIdCStr[128];

    unsigned char mapqUChar;

    float medianQFlt;     /*holds median Q-score*/
    float meanQFlt;       /*holds mean Q-score*/

    float meanAligQFlt;   /*holds mean aligned Q-score*/
    float medianAligQFlt; /*holds aligned median Q-score*/

    uint32_t readLenUInt;    /*read length*/
    uint32_t readAligLenUInt;/*Aligned read length*/
    uint32_t numMatchUInt;   /*number of matches*/
    uint32_t numKeptMatchUInt;
       /*number matches with Q-score > min Q*/
    uint32_t numKeptSNPUInt; /*number of kept mismatches*/
    uint32_t numSNPUInt;     /*number of mismatches*/
    uint32_t numKeptDelUInt; /*number of kept deletions*/
    uint32_t numDelUInt;     /*number of deletions*/
    uint32_t numKeptInsUInt; /*number of kept insertions*/
    uint32_t numInsUInt;     /*number of insertions*/
}readStat;

/*--------------------------------------------------------\
| Name: blankSamEntry (Fun-01:)
| Use:
|  - Sets all non-alloacted variables in samEntryST to 0
| Input:
|  - samEntryST:
|    - samEntry structure to blank/initialize
| Output:
|  - Modifies:
|    o Sets every variable but samEntryCStr to 0
|    o the first character in samEntryCStr to be '\0'
\--------------------------------------------------------*/
static inline void blankSamEntry(
    struct samEntry *samEntryST /*samEntry to blank*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: blankSamEntry 
   '  - Sets all non-alloacted variables in
   '    samEntryST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    samEntryST->queryCStr = 0;
    samEntryST->refCStr = 0;
    samEntryST->cigarCStr = 0;
    samEntryST->seqCStr = 0;
    samEntryST->qCStr = 0; 

    /*Make sure the c-strings are blanked*/
    if(samEntryST->samEntryCStr != 0)
        *samEntryST->samEntryCStr = '\0';

    samEntryST->flagUSht = 0;
    samEntryST->mapqUChar = 0;

    samEntryST->medianQFlt = 0;
    samEntryST->medianAligQFlt = 0;
    samEntryST->meanQFlt = 0;
    samEntryST->meanAligQFlt = 0;

    samEntryST->posOnRefUInt = 0;

    samEntryST->unTrimReadLenUInt = 0;
    samEntryST->readLenUInt = 0;
    samEntryST->readAligLenUInt = 0;

    samEntryST->numMatchUInt = 0;
    samEntryST->numKeptMatchUInt = 0;
    samEntryST->numKeptSNPUInt = 0;
    samEntryST->numSNPUInt = 0;
    samEntryST->numKeptInsUInt = 0;
    samEntryST->numInsUInt = 0;
    samEntryST->numKeptDelUInt = 0;
    samEntryST->numDelUInt = 0;

    samEntryST->totalQScoreULng = 0;
    samEntryST->totalAlnQScoreULng = 0;

    /*Reset the q-score histogram*/
    for(
       uint16_t uShtCnt = 0;
       uShtCnt < MAX_Q_SCORE;
       ++uShtCnt
    ) { /*Set all values in histograms to 0*/
        samEntryST->seqQHistUInt[uShtCnt] = 0;
        samEntryST->seqQAlnHistUInt[uShtCnt] = 0;
    } /*Set all values in histograms to 0*/

    return;
} /*blankSamEntry*/

/*--------------------------------------------------------\
| Name: initSamEntry (Fun-03:)
| Use:
|  - Sets all variables in a samEntry struct to 0 or NULL
|  - This should only be used after malloc or decleration.
| Input:
|  - samEntryST:
|    o pointer to samEntry struct to initialize
| Output:
|  - Modifies:
|    - samEntry to have all variables set to 0
| Warning:
|  - This sets samEntryCStr to 0 without freeing.
\--------------------------------------------------------*/
static inline void initSamEntry(
    struct samEntry *samEntryST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: initSamEntry
   '  - Initalize a samEntry struct to 0's
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    samEntryST->samEntryCStr = 0;
    samEntryST->lenBuffULng = 0;
    blankSamEntry(samEntryST);

    return;
} /*initSamEntry*/

/*--------------------------------------------------------\
| Name: blankReadStats (Fun-03:)
| Use:
|  - modifes samEntryST to have every non-pointer variable
|    set to 0
| Input:
|  - samEntryST:
|    o samEntry struct to blank
| Output:
|  - Modifies:
|    o Sets every non-pointer variable in blankReadStats
|      to 0. The only exeptions are flagUSht, mapqUChar,
|      readLenUInt, unTrimReadLenUInt
\--------------------------------------------------------*/
static inline void blankReadStats(
    struct samEntry *samEntryST /*Structure to blank*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: blankReadStats
   '  - modifes samEntryST to have every non-pointer
   '    variable set to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /*These are taken from the sam file*/
    /*samEntryST->flagUSht = 0;
    samEntryST->mapqUChar = 0;*/

    samEntryST->medianQFlt = 0;
    samEntryST->medianAligQFlt = 0;
    samEntryST->meanQFlt = 0;
    samEntryST->meanAligQFlt = 0;

    /* found when adding read entry to a samEnty struct*/
    /*samEntryST->readLenUInt = 0;
    samEntryST->unTrimReadLenUInt = 0;*/
    samEntryST->readAligLenUInt = 0;

    samEntryST->numMatchUInt = 0;
    samEntryST->numKeptMatchUInt = 0;
    samEntryST->numKeptSNPUInt = 0;
    samEntryST->numSNPUInt = 0;
    samEntryST->numKeptInsUInt = 0;
    samEntryST->numInsUInt = 0;
    samEntryST->numKeptDelUInt = 0;
    samEntryST->numDelUInt = 0;

    samEntryST->totalQScoreULng = 0;
    samEntryST->totalAlnQScoreULng = 0;

    /*Reset the q-score histogram*/
    for(
       uint16_t uShtCnt = 0;
       uShtCnt < MAX_Q_SCORE;
       ++uShtCnt
    ){ /*Set all values in histograms to 0*/
        samEntryST->seqQHistUInt[uShtCnt] = 0;
        samEntryST->seqQAlnHistUInt[uShtCnt] = 0;
    } /*Set all values in histograms to 0*/

    return;
} /*blankReadStats*/

/*--------------------------------------------------------\
| Name: deepCpSamReadStats (Fun-04:)
| Use:
|  - Copies the contents of one samEntry struct to another
|    struct.
|  - This is a deep copy, so buffer contents are copied.
|  - This does not copy alignment specific stats, such as
|    cigars and number of matches.
| Input:
|   - oldSamEntry:
|     o samEntry structer to copy
|   - newSamEntry:
|     o sam entry to copy values to
| Output: 
|  - Modifies:
|    o newSamEntry to hold variables in oldSamEntry
|  - Returns:
|    o 1 if suceeded
|    o 64 if failed (memory allowcation error)
\--------------------------------------------------------*/
static inline uint8_t deepCpSamReadStats(
    struct samEntry *oldSamEntry, /*sam entry to copy*/
    struct samEntry *newSamEntry  /*sam entry to copy to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: deepCpSamReadStats
   '  - Copies read stats and sam file line from on
   '    samEntry to another.
   '  - This ignores alignment specific items, such as
   '    cigars and number of matches.
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Check if need to resize arrays in newSamEntry
   '  o fun-04 sec-03:
   '    - Copy sam file line (entry) and stats for entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-04 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *lineCStr = oldSamEntry->samEntryCStr;
    char *cpLineCStr = newSamEntry->samEntryCStr;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-04 Sec-02:
    ^  - Check if need to resize arrays in newSamEntry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(oldSamEntry->lenBuffULng > newSamEntry->lenBuffULng)
    { /*If: need to resize the new array*/
        if(newSamEntry->samEntryCStr != 0)
            free(newSamEntry->samEntryCStr);

        newSamEntry->samEntryCStr =
           malloc(sizeof(char) * oldSamEntry->lenBuffULng);

        if(newSamEntry->samEntryCStr == 0)
            return 64;

        cpLineCStr = newSamEntry->samEntryCStr;
    } /*If: need to resize the new array*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-04 Sec-03:
    ^  - Copy the sam file line (entry) and stats for entry
    ^  o fun-04 sec-02 sub-01:
    ^    - Copy c-string with sam entry
    ^  o fun-04 sec-02 sub-02:
    ^    - Copy stats in from old entry to new
    ^    - Also set pointers to proper locations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-04 Sec-02 Sub-01:
    *  - Copy c-string with sam file line to new buffer
    \*****************************************************/

    while(*lineCStr != '\0') /*Avoid overwriting entry*/
    { /*Loop: copy the samfile line*/
        *cpLineCStr = *lineCStr; /*On a new sequence line*/
        cpLineCStr++;
        lineCStr++;
    } /*Loop: copy the samfile line*/

    *(cpLineCStr + 1) = '\0'; /*Make into a c-string*/

    /*****************************************************\
    * Fun-04 Sec-02 Sub-02:
    *  - Copy shared query stats in old entry to new
    *  - Also set pointers to proper locations
    \*****************************************************/

    /*Set pointers to sequence & Q-score line in copied
    ` sam entry
    */
    newSamEntry->seqCStr =
        newSamEntry->samEntryCStr +
        (oldSamEntry->seqCStr - oldSamEntry->queryCStr);

    (*newSamEntry).qCStr =
        newSamEntry->samEntryCStr +
        (oldSamEntry->qCStr - oldSamEntry->queryCStr);

    /*Record shared stats between oldSame querys*/
    newSamEntry->unTrimReadLenUInt = 0;
    newSamEntry->readLenUInt = oldSamEntry->readLenUInt;
    newSamEntry->meanQFlt = oldSamEntry->meanQFlt;
    newSamEntry->medianQFlt = oldSamEntry->medianQFlt;
    newSamEntry->mapqUChar = oldSamEntry->mapqUChar;
    newSamEntry->posOnRefUInt = 0;

    return 1;
} /*deepCpSamReadStats*/

/*--------------------------------------------------------\
| Name: cpSamEntry (Fun-05:)
| Use:
|  - Does a shallow copy of oldSamEntry to newSamEntry.
|  - Pointers are copied, but not the buffer.
| Input:
|   - oldSamEntry:
|     o samEntry structer to copy
|   - newSamEntry:
|     o sam entry to copy values to
| Output: 
|  - Modifies:
|    o newSamEntry to have some oldSamEntry values
|      (q-score, sequence, length's, mean/medianQ)
| Warning:
|  - The pointer addresses are copied instead of buffer
|    adressed. So, do not free oldSamEntry or newSamEntry
|    until you are done with both structs.
\--------------------------------------------------------*/
static inline void cpSamEntry(
    samEntry *oldSamEntry,/*Shallow copy from*/
    samEntry *newSamEntry/*Holds shallow copied variables*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: cpSamEntry
   '  - Copies address of q-score & sequence pointers from
   '    an old samEntry to an new SamEntry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*Set pointers to sequence and Q-score*/
    newSamEntry->seqCStr = oldSamEntry->seqCStr;
    newSamEntry->qCStr = oldSamEntry->qCStr;

    /*Record shared stats between oldSame querys*/
    newSamEntry->unTrimReadLenUInt = 0;
    newSamEntry->readLenUInt = oldSamEntry->readLenUInt;
    newSamEntry->meanQFlt = oldSamEntry->meanQFlt;
    newSamEntry->medianQFlt = oldSamEntry->medianQFlt;
    newSamEntry->mapqUChar = oldSamEntry->mapqUChar;
    newSamEntry->posOnRefUInt = 0;

    return;
} /*cpSamEntry*/

/*--------------------------------------------------------\
| Name: freeStackSamEntry (Fun-06:)
| Use:
|  - Frees all variables in samEntry, but not samEntry
| Input:
|  - samEntryST
|    o Sam entry struct to free the interal memory of
| Output:
|  - Frees:
|    o allocated memory in samEntryST
\--------------------------------------------------------*/
static inline void freeStackSamEntry(
    struct samEntry *samEntryST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: freeStackSamEntry
   '  - Frees heap allocations in a stack allocated
   '    samEntry struct
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(samEntryST->samEntryCStr != 0)
    { /*If I have a sam enty to free*/
        if(samEntryST->queryCStr==samEntryST->samEntryCStr)
        { /*If: I need to set values to 0*/
            samEntryST->queryCStr = 0;
            samEntryST->refCStr = 0;
            samEntryST->cigarCStr = 0;
            samEntryST->seqCStr = 0;
            samEntryST->qCStr = 0;
        } /*If: I need to set values to 0*/

        free(samEntryST->samEntryCStr);
        samEntryST->samEntryCStr = 0;
    } /*If I have a sam enty to free*/

    return;
} /*freeStackSamEntry*/

/*--------------------------------------------------------\
| Name: freeHeapSamEntry (Fun-07:)
| Use:
|  - Frees all a samEntry struct
| Input:
|  - samEntryST
|    o Sam entry struct to free
| Output:
|  - Frees:
|    o samEntryST and its variables with allocated memory
\--------------------------------------------------------*/
static inline void freeHeapSamEntry(
    struct samEntry **samEntryST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: freeHeapSamEntry
   '  - Frees a samEntry structer (and sets to null)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    freeStackSamEntry(*samEntryST);
    free(*samEntryST);
    *samEntryST = 0;

    return;
} /*freeStackSamEntry*/

/*--------------------------------------------------------\
| Name: processSamEntry (Fun-08:)
| Use:
|  - Extracts the data from the sam entry from the buffer
|    stored in samEntryST
| Input:
|  - samEntryST:
|    o Hs the sam entry with buffer to extract data from
| Output:
|  - Modifies:
|    o read id, reference id, q-score, cigar, & sequence
|      pointers to piont to their entires in the sam entry.
|      This also grabs the flag & mapq.
\--------------------------------------------------------*/
static inline void processSamEntry(
    struct samEntry *samEntryST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: processSamEntry
   '  - Sets Q-score, cigar, & sequence pointers and 
   '    stats in a sam entry sotred in a samEntry struct.
   '  o fun-08 sec-01:
   '    - Variable declerations
   '  o fun-08 sec-02:
   '    - Check if their is a sam entry
   '  o fun-08 sec-03:
   '    - Extract data from the sam entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-08 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *samIterUC = samEntryST->samEntryCStr;
    unsigned int cigEntryUI = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-08 Sec-02:
    ^  - Check if their is a sam entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samIterUC == 0 || *samIterUC == '@')
    { /*If their is nothing to do or is header entry*/
        samEntryST->queryCStr = 0;
        samEntryST->refCStr = 0;
        samEntryST->cigarCStr = 0;
        samEntryST->seqCStr = 0;
        samEntryST->mapqUChar = 0;
        samEntryST->qCStr = 0;
        samEntryST->flagUSht = 0;

        return;
    } /*If their is nothing to do or is header entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-08 Sec-03:
    ^  - Extract data from the sam entry
    ^  o fun-08 sec-03 sub-01:
    ^    - Find non-sequence & q-score entries
    ^  o fun-08 sec-03 sub-02:
    ^    - Find sequence length from the cigar entry
    ^  o fun-08 sec-03 sub-03:
    ^    - Find sequence & q-score entries
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-08 Sec-03 Sub-01:
    *  - Find non-sequence & q-score entries
    \*****************************************************/

    /*Store query id address*/
    samEntryST->queryCStr = samIterUC;

    /*Get the flag*/
    while(*samIterUC > 32) samIterUC++;
    samIterUC++;             /*Move of tab after query id*/
    samIterUC=cStrToUSht(samIterUC, &samEntryST->flagUSht);

    /*Find the reference id*/
    samIterUC++;           /*Move off tab after flag*/
    samEntryST->refCStr = samIterUC; 

    while(*samIterUC > 32) samIterUC++;
    samIterUC++;   /*Move of tab after reference id*/

    /*Find the starting position on the reference*/
    samIterUC =
       cStrToUInt(samIterUC, &samEntryST->posOnRefUInt);
    samIterUC++;   /*Move of tab after position*/

    /*Find the mapq entry*/
    samIterUC =
       cStrToUChar(samIterUC, &samEntryST->mapqUChar);

    /*Find the cigar entry*/
    samIterUC++;    /*Move of tab after mapq entry*/
    samEntryST->cigarCStr = samIterUC;

    /*****************************************************\
    * Fun-08 Sec-03 Sub-02:
    *  - Find sequence length from the cigar entry
    \*****************************************************/

     samEntryST->unTrimReadLenUInt = 0;

     while(*samIterUC > 32)
     { /*Loop: Get the stats from the cigar entry*/
        samIterUC = cStrToUInt(samIterUC, &cigEntryUI);

        switch(*samIterUC)
        { /*Switch: Count number of cigar entries*/
           case '=':
           case 'M':
           case '\t':
              samEntryST->unTrimReadLenUInt += cigEntryUI;
              samEntryST->numMatchUInt += cigEntryUI;
              break;
           case 'X':
              samEntryST->unTrimReadLenUInt += cigEntryUI;
              samEntryST->numSNPUInt += cigEntryUI;
              break;
           case 'I':
              samEntryST->unTrimReadLenUInt += cigEntryUI;
              samEntryST->numInsUInt += cigEntryUI;
              break;
           case 'D':
              samEntryST->numDelUInt += cigEntryUI;
              break;
           case 'S':
              samEntryST->unTrimReadLenUInt += cigEntryUI;
              break;
           default: break;
        } /*Switch: Count number of cigar entries*/

        cigEntryUI = 0;
        ++samIterUC; /*Get off entry character type*/
     } /*Loop: Get the stats from the cigar entry*/

     samEntryST->readLenUInt=samEntryST->unTrimReadLenUInt;
     ++samIterUC; /*Get off the tab*/

    /*****************************************************\
    * Fun-08 Sec-03 Sub-02:
    *  - Find sequence & q-score entries
    \*****************************************************/


    /*Find the sequence entry*/
    for(uint8_t uCharCnt = 0; uCharCnt < 3; ++uCharCnt)
    { /*Loop past RNEXT, PNEXT, & TLEN*/
        while(*samIterUC > 32) samIterUC++;

        samIterUC++;      /*Move of tab after last entry*/
    } /*Loop past RNEXT, PNEXT, & TLEN*/

    /*Recored sequence address*/
    samEntryST->seqCStr = samIterUC;

    /*Check if have a sequence entry*/
    if(*samIterUC != '*')
       samIterUC += samEntryST->readLenUInt + 1;
    else samIterUC += 2;

    samEntryST->qCStr = samIterUC; /*Set q-score pionter*/

    return;
} /*processSamEntry*/

/*--------------------------------------------------------\
| Name: readNextPartOfLine (Fun-09:)
| Use:
|  - Reads in the next part of the line.
|  - Marks the samStruct->lenBuffULng - 2 position to be
|    '\0\' or '\n' if a complete line was read in.
| Input:
|  - inFILE:
|    o File to read from
|  - samStruct:
|    o samEntry struct to copy buffer to
|  - buffChangeULng:
|    o Amount to increase buffer by when full
| Output:
|  - Returns:
|    o 0 for eof or memory error (no data read in)
|    o pointer to start of buffer with characters read in
|  - Sets:
|    o the second to last entry (samStruct->lenBuffULng-2)
|      in samStruct->samEntryCStr to '\0'.
\--------------------------------------------------------*/
static inline char * readNextPartOfLine(
    FILE *inFILE,
    struct samEntry *samStruct,
    uint64_t buffChangeULng 
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: readNextPartOfLine
   '  - Reads in next part of line when fgets did not
   '    get a full line
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char *tmpUCStr = 0;

    samStruct->samEntryCStr =
        realloc(
            samStruct->samEntryCStr, 
              sizeof(uint8_t)
            * (samStruct->lenBuffULng
            + buffChangeULng)
        ); /*Rellocate memory*/

    tmpUCStr =
       samStruct->samEntryCStr + samStruct->lenBuffULng -1;

    samStruct->lenBuffULng += buffChangeULng;

    while(*tmpUCStr == '\0') /*Make sure on first null*/
    { /*Loop: till not on null*/
        ++buffChangeULng;
        --tmpUCStr;
    } /*Loop: till not on null*/

    ++tmpUCStr;
      /*Not adjusting buff change because fgets reads in
      ` one off
      */

    if(samStruct->samEntryCStr == 0)
        return 0;

    /*Set my marker for non-entire line read in*/
    *(  samStruct->samEntryCStr
      + samStruct->lenBuffULng
      - 2
    ) = '\0';

    return fgets(tmpUCStr, buffChangeULng, inFILE);
} /*readNextPartOfLine*/

/*--------------------------------------------------------\
| Name findSamCig (Fun-10:)
| Use:
|  - Finds the cigar entry in a sam entry and used the
|    cigar to find the sequence length
| Input:
|  - inFILE:
|    o File that has the remaing part of the sam file entry
|  - samStruct:
|    o samEntry struct with buffer to hold cigar and its
|      related stats
| Output:
|  - Modifies:
|    o samStruct->samEntryCStr, samStruct->lenBuffULng, &
|      samStruct->unTrimedLenULng
|    o the cigar pointer in samStruct to point to the cigar
|  - Returns:
|    o 1 if succeded
|    o 4 if end of file or memory allocation error
| Note:
|  - This requires that samStruct->samEntryCStr has read
|    in at least part of the sam entry
\--------------------------------------------------------*/
static inline uint8_t findSamCig(
    FILE *inFILE, /*File to search for cigar in*/
    struct samEntry *samStruct
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: findSamCig
   '  - Finds the cigar in same entry & also finds the
   '    number of bases in query from cigar
   '  o fun-10 sec-01:
   '    - Variable declerations
   '  o fun-10 sec-02:
   '    - Set query pointer & read to start of flag
   '  o fun-10 sec-03:
   '    - Read in flag & move to start of reference id
   '  o fun-10 sec-04:
   '    - Set reference id pointer & move to position
   '  o fun-10 sec-05:
   '    - Read position & move to start of mapq
   '  o fun-10 sec-06:
   '    - Read in mapq & move to start of cigar
   '  o fun-10 sec-07:
   '    - Read in cigar
   '  o fun-10 sec-08:
   '    - Find length of sequence in the cigar
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-10 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *oldBufAddUChar = 0;/*Used to refind pointers*/
    char *samIterUChar = samStruct->samEntryCStr;
    uint32_t cigEntryUInt = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-10 Sec-02:
    ^  - Set query pointer & read to start of flag
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samStruct->queryCStr = samIterUChar;

    while(*samIterUChar != '\t')
    { /*Loop: till off the query id*/
        if(*samIterUChar == '\0')
        { /*If: I did not get a full line*/
            samIterUChar =
               readNextPartOfLine(
                 inFILE,
                 samStruct,
                 samStruct->lenBuffULng/*Double buff size*/
            ); /*Read in next section of line*/

            /*Check if I have an incomplete sam aligment*/
            if(samIterUChar == 0) return 4; 

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
        } /*If: I did not get a full line*/

        ++samIterUChar;
    } /*Loop: till off the query id*/

    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-10 Sec-03:
    ^  - Read in flag & move to start of reference id
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get how many characters mapq is in the sam entry*/
    cigEntryUInt = samIterUChar - samStruct->samEntryCStr;

    /*Need to make sure entire flag is read in*/ 
    for(
       uint32_t intIter = 0;
       *(samIterUChar+intIter) != '\t';
       ++intIter
    ) { /*Loop: till at end of number*/
        if(*(samIterUChar + intIter) == '\0')
        { /*If: at end of buffer*/
            samIterUChar =
               readNextPartOfLine(
                 inFILE,
                 samStruct,
                 samStruct->lenBuffULng/*Double buff size*/
            ); /*Read in next section of line*/

            /*check if is an incomplete sam aligment*/
            if(samIterUChar == 0) return 4;

            /*Make sure pionter is set start of flag*/
            samIterUChar =
               samStruct->samEntryCStr + cigEntryUInt;
            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
        } /*If: at end of buffer*/
    } /*Loop: till at end of number*/

    samIterUChar =
       cStrToUSht(samIterUChar, &samStruct->flagUSht);
    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-10 Sec-04:
    ^  - Set reference id pointer & move to position
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samStruct->refCStr = samIterUChar;

    while(*samIterUChar != '\t')
    { /*While on the reference id*/
        if(*samIterUChar == '\0')
        { /*If at end of buffer*/
            oldBufAddUChar = samStruct->samEntryCStr;

            samIterUChar =
              readNextPartOfLine(
                inFILE,
                samStruct,
                samStruct->lenBuffULng /*Double buff size*/
            ); /*Read in next section of line*/

            /*Check if is an incomplete sam aligment*/
            if(samIterUChar == 0) return 4;

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
    
            samStruct->refCStr =
                samStruct->samEntryCStr +
                (samStruct->refCStr - oldBufAddUChar);
        } /*If at end of buffer*/

        ++samIterUChar;
    } /*While on the reference id*/

    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-10 Sec-05:
    ^  - Get off position & move to start of mapq
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    for(
       uint32_t intIter = 0;
       *(samIterUChar+intIter) != '\t';
       ++intIter
    ){ /*Loop: till at start of of mapq*/
        if(*samIterUChar == '\0')
        { /*If: at end of buffer*/
            oldBufAddUChar = samStruct->samEntryCStr;

            samIterUChar =
              readNextPartOfLine(
                inFILE,
                samStruct,
                samStruct->lenBuffULng /*Double buff size*/
            ); /*Read in next section of line*/

            /*Check if is an incomplete sam aligment*/
            if(samIterUChar == 0) return 4;

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
    
            samStruct->refCStr =
                samStruct->samEntryCStr +
                (samStruct->refCStr - oldBufAddUChar);
        } /*If: at end of buffer*/
    } /*Loop: till at start of of mapq*/

    /*Record the position entry*/
    samIterUChar =
       cStrToUInt(samIterUChar, &samStruct->posOnRefUInt);
    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-10 Sec-06:
    ^  - Read in mapq & move to start of cigar
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get how many characters mapq is in the sam entry*/
    cigEntryUInt = samIterUChar - samStruct->samEntryCStr;

    /*Need to make sure entire mapq is read in*/ 
    for(
       uint32_t intIter = 0;
       *(samIterUChar+intIter) != '\t';
       ++intIter
    ){ /*Loop: till at end of mapq (is number)*/
        if(*(samIterUChar + intIter) == '\0')
        { /*If: at end of buffer*/
            oldBufAddUChar = samStruct->samEntryCStr;

            samIterUChar =
              readNextPartOfLine(
                inFILE,
                samStruct,
                samStruct->lenBuffULng /*Double buff size*/
            ); /*Read in next section of line*/

            /*Check if is an incomplete sam aligment*/
            if(samIterUChar == 0) return 4;

            /*Make sure pionter points to start of mapq*/
            samIterUChar =
               samStruct->samEntryCStr + cigEntryUInt;

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
    
            samStruct->refCStr =
                samStruct->samEntryCStr +
                (samStruct->refCStr - oldBufAddUChar);
        } /*If: at end of buffer*/
    } /*Loop: till at end of mapq (is number)*/

    samIterUChar =
       cStrToUChar(samIterUChar, &(samStruct->mapqUChar));
    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-10 Sec-07:
    ^  - Read in cigar (make sure I have a full entry)
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samStruct->cigarCStr = samIterUChar;

    while(*samIterUChar != '\t')
    { /*Loop: Till at end of cigar*/
        if(*samIterUChar == '\0')
        { /*If: at end of buffer*/
            oldBufAddUChar = samStruct->samEntryCStr;

            samIterUChar =
              readNextPartOfLine(
                inFILE,
                samStruct,
                samStruct->lenBuffULng /*Double buff size*/
            ); /*Read in next section of line*/

            /*Check if is an incomplete sam aligment*/
            if(samIterUChar == 0) return 4;

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
    
            samStruct->refCStr =
                samStruct->samEntryCStr +
                (samStruct->refCStr - oldBufAddUChar);

            samStruct->cigarCStr =
                samStruct->samEntryCStr +
                (samStruct->cigarCStr - oldBufAddUChar);
        } /*If: at end of buffer*/

        ++samIterUChar;
    } /*Loop: Till at end of cigar*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-10 Sec-08:
    ^  - Find length of sequence in the cigar
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samIterUChar = samStruct->cigarCStr;

    while(*samIterUChar != '\t')
    { /*Loop: though all cigar entries to read in*/
        samIterUChar =
           cStrToUInt(samIterUChar, &cigEntryUInt);

        switch(*samIterUChar)
        { /*Switch: Count number of cigar entries*/
           case '=':
           case 'M':
           case '\t':
              samStruct->unTrimReadLenUInt += cigEntryUInt;
              samStruct->numMatchUInt += cigEntryUInt;
              break;
           case 'X':
              samStruct->unTrimReadLenUInt += cigEntryUInt;
              samStruct->numSNPUInt += cigEntryUInt;
              break;
           case 'I':
              samStruct->unTrimReadLenUInt += cigEntryUInt;
              samStruct->numInsUInt += cigEntryUInt;
              break;
           case 'D':
              samStruct->numDelUInt += cigEntryUInt;
              break;
           case 'S':
              samStruct->unTrimReadLenUInt += cigEntryUInt;
              break;
           default: break;
        } /*Switch: Count number of cigar entries*/

        cigEntryUInt = 0;
        ++samIterUChar; /*Get off entry character type*/
    } /*Loop: though all cigar entries to read in*/

    samStruct->readLenUInt = samStruct->unTrimReadLenUInt;
    return 1;
} /*findSamCig*/

/*--------------------------------------------------------\
| Name: readSamLine (Fun-11:)
| Use:
|  - Reads in a sam file entry into a samEntry struct
| Input:
|  - samStruct:
|    o Pointer to samEntry struct to store the sam entry in
|  - inFILE:
|    o File with sam entry to get
| Output:
|  - Returns:
|    o 1 if succeded
|    o 2 if end of file
|    o 64 if memory allocation error
|  - Modifies:
|    o All pointers in samStruct & variables holding read
|      Lenths
\--------------------------------------------------------*/
static inline uint8_t readSamLine(
    struct samEntry *samStruct,
    FILE *inFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: readSamLine
   '  - Reads in sam entry into a samEntry struct
   '  o fun-11 sec-01:
   '    - Variable declerations
   '  o fun-11 sec-02:
   '    - Check if have buffer to read into/add buffer
   '  o fun-11 sec-03:
   '    - Check if can read line in one go
   '  o fun-11 sec-04:
   '    - Check if is header & need to read in more
   '  o fun-11 sec-05:
   '    - If not read in cigar & find sequence length
   '  o fun-11 sec-06:
   '    - Check if have full sequence from cigar read in
   '  o fun-11 sec-07:
   '    - If not find sequence length
   '  o fun-11 sec-08:
   '    - Use sequence length to get sequence & q-score
   '  o fun-11 sec-09:
   '    - If not, read rest of cigar by adding 1kb chunks
   '  o fun-11 sec-10:
   '    - Set sequence & q-score pointers
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *oldBufAddUChar = 0;
    char *samIterUChar = 0; /*fget error or end of file*/
    char checkC = 0;        /*Used for single char checks*/
    uint64_t newBuffSizeULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-02:
    ^  - Check if have buffer to read into & prepare buffer
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samStruct->samEntryCStr == 0)
    { /*If: I have no buffer*/
        samStruct->lenBuffULng = 1024;
        samStruct->samEntryCStr=malloc(sizeof(char)*1024);
    } /*If: I have no buffer*/

    /*Set second to last character in buffer to null.
    ' This tells me if I read in an entire line
    */
    *(
         samStruct->samEntryCStr
       + samStruct->lenBuffULng
       - 2
    ) = '\0';

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-03:
    ^  - Check if can read line in one go
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samIterUChar =
        fgets(
            samStruct->samEntryCStr,
            samStruct->lenBuffULng,
            inFILE
    ); /*See if can read in a single line*/

    if(samIterUChar == 0)
        return 2;           /*At end of the file*/

    checkC =
      *(samStruct->samEntryCStr +samStruct->lenBuffULng-2);

    if(checkC =='\0' || checkC == '\n')
    { /*If: I grabed the entire line*/
        processSamEntry(samStruct);
        return 1;
    } /*If: I grabed the entire line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-04:
    ^  - Check if is header & need to read in more
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*samStruct->samEntryCStr == '@')
    { /*If: this is a header line*/
        newBuffSizeULng = 1024; /*Expand buffer by*/
            /*I should be close to the full entry, so 
            ` expand the buffer by in small chunks only
            */
    
        while(!(checkC == '\0' || checkC == '\n'))
        { /*Loop: Get the next part of the line*/
            samIterUChar =
               readNextPartOfLine(
                  inFILE,
                  samStruct,
                  newBuffSizeULng
            ); /*Get the next or remnaing part of line*/

            if(samIterUChar == 0) return 64;

            checkC =
              *(
                   samStruct->samEntryCStr
                 + samStruct->lenBuffULng
                 - 2
            ); /*Find the ending character*/
        } /*Loop: Get the next part of the line*/

        return 1;
    } /*If: this is a header line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-05:
    ^  - Read in cigar & get sequence length from cigar
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Read in cigar (findSamCig) & check for memory error*/
    if(findSamCig(inFILE, samStruct) == 4)
        if(samStruct->samEntryCStr == 0) return 64;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-06:
    ^  - Check if have full sequence from cigar read in
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if have full sequence & set pionters*/
    checkC =
      *(samStruct->samEntryCStr +samStruct->lenBuffULng-2);

    if(checkC =='\0' || checkC == '\n')
    { /*If: I now have the entire line*/
        samIterUChar = samStruct->cigarCStr;

        /*I have to get past the cigar, RNEXT, PNEXT, &
        ` TLEN entires
        */
        for(uint32_t intTab = 0; intTab < 4; ++intTab)
        { /*Loop: till on sequence entry*/
            while(*samIterUChar != '\t')
                ++samIterUChar;
            ++samIterUChar; /*get off the tab*/
        } /*Loop: till on sequence entry*/

        samStruct->seqCStr = samIterUChar;
        samStruct->qCStr =
           samIterUChar + samStruct->readLenUInt + 1;
            /*reaLenUint to move past sequence; +1 to get
            ` off the tab
            */
        return 1;
    } /*If: I now have the entire line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-07:
    ^  - If not find sequence length
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samIterUChar = samStruct->cigarCStr;

    while(*samIterUChar != '\t')
    { /*Loop: till at end of cigar*/
        ++samIterUChar;
        ++newBuffSizeULng;
        /*Find length of the cigar (avoid guessing)*/
    } /*Loop: till at end of cigar*/

    /*MD entries will be bit longer than cigars, but I can
    ` double the cigar to get a rough size. This way I
    ` only have to read in once or twice
    */
    newBuffSizeULng = newBuffSizeULng << 1;

    /*Figure out buffer size with sequence & q-score*/
    newBuffSizeULng +=         /*newBuff has Cigar length*/
          (samStruct->readLenUInt << 1)  /*seq + q-score*/
        + samStruct->cigarCStr
        - samStruct->samEntryCStr /*pre cigar*/
        + 2048;                            /*Extra buffer*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-07:
    ^  - Use sequence length to get sequence & q-score
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(newBuffSizeULng > samStruct->lenBuffULng)
    { /*If: new buffer size with the sequence is bigger*/
        newBuffSizeULng -= samStruct->lenBuffULng;
        oldBufAddUChar = samStruct->samEntryCStr;

        samIterUChar =
           readNextPartOfLine(
              inFILE,
              samStruct,
              newBuffSizeULng
        ); /*See if I can finish the line*/

        if(samStruct->samEntryCStr == 0) return 64;

        /*Readjust previous pointers*/
        samStruct->queryCStr = samStruct->samEntryCStr;
    
        samStruct->refCStr =
            samStruct->samEntryCStr +
            (samStruct->refCStr - oldBufAddUChar);

        samStruct->cigarCStr =
            samStruct->samEntryCStr +
            (samStruct->cigarCStr - oldBufAddUChar);
    } /*If: new buffer size with the sequence is bigger*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-09:
    ^  - If not, read rest of cigar by adding 1kb chunks
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    newBuffSizeULng = 1024;
    /*I Should be close to the full line length.
    ` So I am going to read in small chunks at a time
    */
    checkC =
      *(samStruct->samEntryCStr +samStruct->lenBuffULng-2);

    while(!(checkC == '\0' || checkC == '\n'))
    { /*Loop: Finish reading in the sam entry*/
        oldBufAddUChar = samStruct->samEntryCStr;
        samIterUChar =
            readNextPartOfLine(
               inFILE,
               samStruct,
               newBuffSizeULng
        ); /*Read in the next part of the line*/

        if(samStruct->samEntryCStr == 0) return 64;

        checkC =
          *(
               samStruct->samEntryCStr
             + samStruct->lenBuffULng
             - 2
        ); /*Find the ending character*/

        /*Readjust previous pointers*/
        samStruct->queryCStr = samStruct->samEntryCStr;
    
        samStruct->refCStr =
            samStruct->samEntryCStr +
            (samStruct->refCStr - oldBufAddUChar);

        samStruct->cigarCStr =
            samStruct->samEntryCStr +
            (samStruct->cigarCStr - oldBufAddUChar);
    } /*Loop: Finish reading in the sam entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-10:
    ^  - Set sequence & q-score pointers
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*At this point I have the entire entry, but I still
    `  need to set some pointers
    */
    samIterUChar = samStruct->cigarCStr;

    /*TODO: SEQ IS ON Q-SCORE, FIX go though cigar entry
    ` I need to go through the cigar, RNEXT, PNEXT, & TLEN
    ` entires. This should not be the case and needs to be
    ` fixed.
    */
    for(uint32_t intTab = 0; intTab < 4; ++intTab)
    { /*Loop: Go to seq entry start*/
        while(*samIterUChar != '\t')
            ++samIterUChar;
        ++samIterUChar; /*get off the tab*/
    } /*Loop: Go to seq entry start*/

    samStruct->seqCStr = samIterUChar;

    if(*samStruct->seqCStr != '*')
    { /*If: I have a sequence entry*/
        if(*samStruct->cigarCStr != '*')
            samStruct->qCStr = samIterUChar + samStruct->readLenUInt +1;
        else
        { /*Else: I have no cigar entry*/
            while(*samIterUChar > 32)
            { /*Loop: Find lenght of sequence*/
                ++samStruct->unTrimReadLenUInt;
                ++samIterUChar;
            } /*Loop: Find lenght of sequence*/

            samStruct->readLenUInt =
               samStruct->unTrimReadLenUInt;
            ++samIterUChar;  /*Get off tab after sequence*/
            samStruct->qCStr = samIterUChar;
        } /*Else: I have no cigar entry*/
    } /*If: I have a sequence entry*/

    else samStruct->qCStr = samIterUChar +2;/*No sequence*/
   
    return 1;
} /*readSamLine*/

/*--------------------------------------------------------\
| Name: printSamEntry (Fun-12:)
| Use:
|  - Prints the sam file entry to a file. This does not
|    print any extra stats that were found.
| Input:
|  - samST
|    o samEntry struct with sam entry to print out
|  - outFILE:
|    o File to print the sam entry to
| Output:
|  - Prints:
|    o Sam file entry in samST to outFILE.
\--------------------------------------------------------*/
static inline void printSamEntry(
    struct samEntry *samST, /*sam entry to print out*/
    FILE *outFILE           /*File to print sam entry to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: printSamEntry
   '  - Prints the sam file entry to a file. This does not
   '    print any extra stats that were found.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*
    ` printing out query name since it is the start of the
    ` sam entry. This allows me to store the sam entry in a
    ` different samEntry structure, but still have a handle
    ` on the entire entry with the query. I do this when
    ` their might be multiple entries for a single read.
    ` In this case only the first entry has the sequence.
    */
    if(samST->queryCStr != 0)
       fputs(samST->queryCStr, outFILE);
       /*fprintf(outFILE, "%s", samST->queryCStr);*/
    else
       fputs(samST->samEntryCStr, outFILE);
       /*fprintf(outFILE, "%s", samST->samEntryCStr);*/
} /*printSamEntry*/


/*--------------------------------------------------------\
| Name: samToFq (Fun-13:)
| Use:
|  - Prints the sam entry as a fastq entry to a fastq file
| Input:
|  - samST:
|    o samEntry structure with fastq entry to print out
|  - outFILE:
|    o Fastq file to print the new fastq entry to
| Output:
|  - Prints:
|    o fastq entry from samST to outFILE
\--------------------------------------------------------*/
static inline void samToFq(
    struct samEntry *samST, /*sam entry to print as fastq*/
    FILE *outFILE           /*to print fastq entry to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: samToFq
   '  - Prints sam entry as a fastq entry to a fastq file
   '  o fun-13 sec-01:
   '    - Variable declerations
   '  o fun-13 sec-02:
   '    - Find length of query id
   '  o fun-13 sec-03:
   '    - Print sam entry as fastq entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-13 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpQueryCStr = samST->queryCStr;
    uint16_t lenQueryUSht = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-13 Sec-02:
    ^  - Find length of query id
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*tmpQueryCStr > 32)
    { /*wile not at end of query id*/
        tmpQueryCStr++; /*move to end of query name*/
        lenQueryUSht++; /*Find length of query id*/
    } /*wile not at end of query id*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-13 Sec-03:
    ^  - Print sam entry as fastq entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check If I need to print out the @ for the header*/
    if(*samST->queryCStr != '@')
        fwrite("@", sizeof(char), 1, outFILE);

    fwrite(
       samST->queryCStr,
       sizeof(char),
       lenQueryUSht,
       outFILE
    ); /*print out the query id*/

    fwrite("\n", sizeof(char), 1, outFILE);

    fwrite(
        samST->seqCStr,
        sizeof(char),
        samST->readLenUInt,
        outFILE
    ); /*Write sequence entry*/

    /*Print out the spacer entry*/
    fwrite("\n+\n", sizeof(char), 3, outFILE);

    fwrite(
        samST->qCStr,
        sizeof(char),
        samST->readLenUInt,
        outFILE
    ); /*Write Q-score entry*/

    /*Write new line for next entry*/
    fwrite("\n", sizeof(char), 1, outFILE);

    return;
} /*samToFq*/

/*--------------------------------------------------------\
| Name: printStatHeader (Fun-14:)
| Use:
|  - Prints the header for the stats file
| Input:
|  - statFILE:
|    o File to print the header to
| Output:
|  - Prints:
|    o stat file header to the input file
\--------------------------------------------------------*/
static inline void printStatHeader(
    FILE *statFILE /*File to print the hader to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-14 TOC: printStatsHeader
   '  - Prints the stat file header made using a samEntry
   '    struct
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(statFILE == 0)
      return;

    fprintf(statFILE, "Read\tRef\tMAPQ\treadLength");
    fprintf(statFILE, "\talignedLength\tmatches"); 
    fprintf(statFILE, "\tkeptMatches\tmismatches"); 
    fprintf(statFILE, "\tinsertions\tdeletions\tmedianQ"); 
    fprintf(statFILE, "\tmeanQ\talignedMedianQ");
    fprintf(statFILE, "\talignedMeanQ\tTotalMismatches");
    fprintf(statFILE, "\tTotalInsertions");
    fprintf(statFILE, "\tTotalDeletions\n");

    return;
} /*printStatHeader*/

/*--------------------------------------------------------\
| Name: printSamStats (Fun-15:)
| Use:
|  - Prints out the stats in a samEntry struct to a file
| Input:
|  - samEntryST:
|    o samEntry struct to print stats for
|  - pHeadBl:
|    o 1: Print the header for the stats tsv file
|    o 0: Do not print the header
|  - outFILE:
|    o TSV (tab deliminated) file to print stats to
| Output:
|  - Prints:
|    o line with stats from samEntryStruct to file (as tsv)
|  - Modifies:
|    o printHeaderChar to be 0 if set to 1
\--------------------------------------------------------*/
static inline void printSamStats(
    struct samEntry *samEntryST, /*Has stats to print out*/
    uint8_t *pHeadBl,            /*1: print header*/
    FILE *outFILE                /*File to output to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-15 TOC: printSamStats
   '  - Prints stats in samEntry structure to a tsv file
   '  o fun-15 sec-01:
   '    - Variable declerations
   '  o fun-15 sec-02:
   '    - Print out the header
   '  o fun-15 sec-03:
   '    - Convert query & reference ids to c-strings
   '  o fun-15 sec-04:
   '    - Print out stats from the samEntry structer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-15 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpQueryCStr = (*samEntryST).queryCStr;
    char *tmpRefCStr = (*samEntryST).refCStr;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-15 Sec-02:
    ^  - Print out the header
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*pHeadBl == 1)
    { /*If printing out the header entry*/
        printStatHeader(outFILE);
        *pHeadBl = 0;
    } /*If printing out the header entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-15 Sec-03:
    ^  - Convert query & reference ids to c-strings
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*tmpQueryCStr != '\t')
        tmpQueryCStr++; /*move to end of query name*/
    *tmpQueryCStr = '\0'; /*turn into c-string*/

    while(*tmpRefCStr != '\t')
        tmpRefCStr++; /*move to end of reference name*/
    *tmpRefCStr = '\0'; /*turn into c-string*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-15 Sec-04:
    ^  - Print out stats from the samEntry structer
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Print out the entry stats*/
    fprintf(
        outFILE,
        "%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%f\t%f",
        samEntryST->queryCStr, 
        samEntryST->refCStr,
        samEntryST->mapqUChar,
        samEntryST->readLenUInt,
        samEntryST->readAligLenUInt,
        samEntryST->numMatchUInt,
        samEntryST->numKeptMatchUInt,
        samEntryST->numKeptSNPUInt,
        samEntryST->numKeptInsUInt,
        samEntryST->numKeptDelUInt,
        samEntryST->medianQFlt,
        samEntryST->meanQFlt
    ); /*1st printf: print out stats*/

    fprintf(
        outFILE,
        "\t%f\t%f\t%u\t%u\t%u\n",
        samEntryST->medianAligQFlt,
        samEntryST->meanAligQFlt,
        samEntryST->numSNPUInt,
        samEntryST->numInsUInt,
        samEntryST->numDelUInt
    ); /*2nd fprintf: print out stats*/

    *tmpQueryCStr = '\t'; /*turn back into single string*/
    *tmpRefCStr = '\t'; /*turn back into singel string*/

    return;
} /*printSamStats*/

/*--------------------------------------------------------\
| Name: swapSamEntries (Fun-16:)
| Use:
|  - Swaps the values in two samEntry structs
| Input:
|  - samOneST:
|    o First samEntry structure to swap values
|  - samTwoST:
|    o Second samEntry structure to swap values
| Output:
|  - Modifies:
|    o samOnST to have samTwoST's values
|    o samTwoST to have samOneST's values
| Note:
|  - This does swap the histograms. So, there is a bit of
|    work in this swap.
\--------------------------------------------------------*/
static inline void swapSamEntries(
   struct samEntry *samOneST, /*First entry to swap*/
   struct samEntry *samTwoST  /*Second entry to swap*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-16 TOC: swapSamEntries
   '  - Swaps the values in two sam entries
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   char *swapStr = 0;
   uint8_t swapUC = 0;
   uint16_t swapUS = 0;
   uint32_t swapUI = 0;
   uint64_t swapUL = 0;
   float swapFlt = 0;

   swapStr = samOneST->samEntryCStr;
   samOneST->samEntryCStr = samTwoST->samEntryCStr;
   samTwoST->samEntryCStr = swapStr;

   swapStr = samOneST->queryCStr;
   samOneST->queryCStr = samTwoST->queryCStr;
   samTwoST->queryCStr = swapStr;

   swapStr = samOneST->refCStr;
   samOneST->refCStr = samTwoST->refCStr;
   samTwoST->refCStr = swapStr;

   swapStr = samOneST->cigarCStr;
   samOneST->cigarCStr = samTwoST->cigarCStr;
   samTwoST->cigarCStr = swapStr;

   swapStr = samOneST->seqCStr;
   samOneST->seqCStr = samTwoST->seqCStr;
   samTwoST->seqCStr = swapStr;

   swapStr = samOneST->qCStr;
   samOneST->qCStr = samTwoST->qCStr;
   samTwoST->qCStr = swapStr;

   swapUC = samOneST->mapqUChar;
   samOneST->mapqUChar = samTwoST->mapqUChar;
   samTwoST->mapqUChar = swapUC;

   swapUS = samOneST->flagUSht;
   samOneST->flagUSht = samTwoST->flagUSht;
   samTwoST->flagUSht = swapUS;

   swapFlt = samOneST->medianQFlt;
   samOneST->medianQFlt = samTwoST->medianQFlt;
   samTwoST->medianQFlt = swapFlt;

   swapFlt = samOneST->meanQFlt;
   samOneST->meanQFlt = samTwoST->meanQFlt;
   samTwoST->meanQFlt = swapFlt;

   swapFlt = samOneST->medianAligQFlt;
   samOneST->medianAligQFlt = samTwoST->medianAligQFlt;
   samTwoST->medianAligQFlt = swapFlt;

   swapFlt = samOneST->meanAligQFlt;
   samOneST->meanAligQFlt = samTwoST->meanAligQFlt;
   samTwoST->meanAligQFlt = swapFlt;

   swapUI = samOneST->posOnRefUInt;
   samOneST->posOnRefUInt = samTwoST->posOnRefUInt;
   samTwoST->posOnRefUInt = swapUI;

   swapUI = samOneST->unTrimReadLenUInt;

   samOneST->unTrimReadLenUInt =
      samTwoST->unTrimReadLenUInt;

   samTwoST->unTrimReadLenUInt = swapUI;

   swapUI = samOneST->readLenUInt;
   samOneST->readLenUInt = samTwoST->readLenUInt;
   samTwoST->readLenUInt = swapUI;

   swapUI = samOneST->readAligLenUInt;
   samOneST->readAligLenUInt = samTwoST->readAligLenUInt;
   samTwoST->readAligLenUInt = swapUI;

   swapUI = samOneST->numMatchUInt;
   samOneST->numMatchUInt = samTwoST->numMatchUInt;
   samTwoST->numMatchUInt = swapUI;

   swapUI = samOneST->numKeptMatchUInt;
   samOneST->numKeptMatchUInt = samTwoST->numKeptMatchUInt;
   samTwoST->numKeptMatchUInt = swapUI;

   swapUI = samOneST->numKeptSNPUInt;
   samOneST->numKeptSNPUInt = samTwoST->numKeptSNPUInt;
   samTwoST->numKeptSNPUInt = swapUI;

   swapUI = samOneST->numSNPUInt;
   samOneST->numSNPUInt = samTwoST->numSNPUInt;
   samTwoST->numSNPUInt = swapUI;

   swapUI = samOneST->numKeptDelUInt;
   samOneST->numKeptDelUInt = samTwoST->numKeptDelUInt;
   samTwoST->numKeptDelUInt = swapUI;

   swapUI = samOneST->numDelUInt;
   samOneST->numDelUInt = samTwoST->numDelUInt;
   samTwoST->numDelUInt = swapUI;

   swapUI = samOneST->numKeptInsUInt;
   samOneST->numKeptInsUInt = samTwoST->numKeptInsUInt;
   samTwoST->numKeptInsUInt = swapUI;

   swapUI = samOneST->numInsUInt;
   samOneST->numInsUInt = samTwoST->numInsUInt;
   samTwoST->numInsUInt = swapUI;

   swapUI = samOneST->numInsUInt;
   samOneST->numInsUInt = samTwoST->numInsUInt;
   samTwoST->numInsUInt = swapUL;

   swapUL = samOneST->totalAlnQScoreULng;
   samOneST->totalAlnQScoreULng =
      samTwoST->totalAlnQScoreULng;
   samTwoST->totalAlnQScoreULng = swapUL;

   swapUL = samOneST->lenBuffULng;
   samOneST->lenBuffULng =
      samTwoST->lenBuffULng;
   samTwoST->lenBuffULng = swapUL;
  
   for(uint32_t iQ = 0; iQ < MAX_Q_SCORE; ++iQ)
   { /*Loop: Swap Q-score histogram entries*/
     swapUI = samOneST->seqQHistUInt[iQ];

     samOneST->seqQHistUInt[iQ] =
        samTwoST->seqQHistUInt[iQ];

     samTwoST->seqQHistUInt[iQ] = swapUI;

     /*Deal with the aligned histogram*/
     swapUI = samOneST->seqQAlnHistUInt[iQ];

     samOneST->seqQAlnHistUInt[iQ] =
        samTwoST->seqQAlnHistUInt[iQ];

     samTwoST->seqQAlnHistUInt[iQ] = swapUI;
   } /*Loop: Swap Q-score histogram entries*/
} /*swapSamEntries*/

/*--------------------------------------------------------\
| Name: sortSamEntriesStartLen (Fun-17:)
| Use:
|  - Sorts an array of samEntry structs by reference
|    position and length
| Input:
|  - scoresST:
|    o Array of positions (in scoresStructs) to sort
|    o Each struct has the reference start and end
|  - firstElmUL:
|    o 1st element to sort in scoresST
|  - lastElmUL:
|    o last element to sort in scoresST
| Output:
|  - Modifies
|    o sorts the scoresST array by reference postion first
|      and alignment length second.
\--------------------------------------------------------*/
static inline void sortSamEntriesStartLen(
   struct samEntry **samST,  /*sam entries to sort*/
   unsigned long firstElmUL, /*1st element to sort*/
   unsigned long endElmUL    /*last element to sort*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-17 TOC: sortSamEntriesStartLen
   '  - Sorts an array of samEntry structs by reference
   '    position (start) and aligned refence length with
   '    shell short. Order is least to greatest.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-17 sec-02:
   '    - Find the max search value
   '  o fun-17 sec-03:
   '    - Sort the samEntry struts array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*Number of elements to sort*/
  unsigned long numElmUL = endElmUL - firstElmUL;
  unsigned long alnLenUL = 0;
  unsigned long swapAlnLenUL = 0;

  /*Number of sorting rounds*/
  unsigned long subUL = 0;
  unsigned long nextElmUL = 0;
  unsigned long lastElmUL = 0;
    /*Array of sizes for the sub arrays in shell sort*/
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-17 Sec-02:
  ^  - Find the max search value
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] + 1*/
  subUL = 1; /*Initialzie first array*/
  while(subUL < numElmUL) subUL = (3 * subUL) + 1;
  subUL = (subUL - 1) / 3; /*Account for overshooting*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-17 Sec-03:
  ^  - Sort the samEntry struts array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  while(subUL > 0)
  { /*loop trhough all sub arrays sort the subarrays*/
    for(unsigned long indexUL=0; indexUL<=subUL; ++indexUL)
    { /*For each element in the subarray*/
      for(unsigned long elmUL = indexUL;
          elmUL + subUL <= endElmUL;
          elmUL += subUL
      ) { /*Loop; swap each nth element of the subarray*/
        nextElmUL = elmUL + subUL;

        alnLenUL =
             (*samST + elmUL)->numMatchUInt
           + (*samST + elmUL)->numSNPUInt
           + (*samST + elmUL)->numDelUInt;
        swapAlnLenUL =
             (*samST + nextElmUL)->numMatchUInt
           + (*samST + nextElmUL)->numSNPUInt
           + (*samST + nextElmUL)->numDelUInt;

        if(
          (/*Check if the later position has a lower start*/ 
              (*samST + elmUL)->posOnRefUInt
            > (*samST + nextElmUL)->posOnRefUInt
          ) || 
          ((  /*If both scores have an equal start*/
                (*samST + elmUL)->posOnRefUInt
             == (*samST + nextElmUL)->posOnRefUInt
           ) && (alnLenUL > swapAlnLenUL)
          )
        ) { /*If I need to swap an element*/
          swapSamEntries(
            (*samST + elmUL),
            (*samST + nextElmUL)
          ); /*Swap scores around*/

          nextElmUL = elmUL;
          lastElmUL = elmUL;

          while(lastElmUL >= subUL)
          { /*loop; move swapped element back*/
             lastElmUL -= subUL;

             alnLenUL =
                  (*samST + nextElmUL)->numMatchUInt
                + (*samST + nextElmUL)->numSNPUInt
                + (*samST + nextElmUL)->numDelUInt;
             swapAlnLenUL =
                  (*samST + lastElmUL)->numMatchUInt
                + (*samST + lastElmUL)->numSNPUInt
                + (*samST + lastElmUL)->numDelUInt;

             if(
               (/*Check if later position has lower start*/ 
                   (*samST + nextElmUL)->posOnRefUInt
                 > (*samST + lastElmUL)->posOnRefUInt
               ) || 
               ((  /*If both scores have an equal start*/
                     (*samST + nextElmUL)->posOnRefUInt
                  == (*samST + lastElmUL)->posOnRefUInt
                ) && (alnLenUL > swapAlnLenUL)
               )
             ) break; /*If this element is positioned*/

             swapSamEntries(
               (*samST + nextElmUL),
               (*samST + lastElmUL)
             ); /*Swap scores around*/

             nextElmUL = lastElmUL;
          } /*loop; move swapped element back*/
        } /*If I need to swap elements*/
      } /*Loop; swap each nth element of the subarray*/
    } /*For each element in the subarray*/

    subUL = (subUL - 1) / 3; /*Move to the next round*/
  } /*loop through all sub arrays to sort the subarrays*/

  return; /*Finshed sorting the array*/
} /*sortSamEntriesStartLen*/

/*--------------------------------------------------------\
| Name: readInConFa (Fun-18:)
| Use:
|  - Reads in a reference sequence from a fasta file into
|    a samEntry struct.
|  - I would recomend using seqStruct for this. This was
|    set up before seqStruct, and is kept here to make sure
|    nothing breaks
| Input:
|  - conFaToReadCStr:
|    o Name of fasta file to get read from
|  - refST:
|    o samEntry struct to hold the reference sequence
| Output:
|  - Modifies
|    o refST: To hold the sequence (no header)
|  - Returns
|    o 1 if succeeded
|    o 2 if file does not exist
|    o 4 invalid file
|    o 64 memory allocation error
| Note:
|  - Fasta file should only have one sequence
\--------------------------------------------------------*/
static inline unsigned char readInConFa(
    char *conFaToReadCStr, /*Name of fasta file with ref*/
    struct samEntry *refST /*Sam struct to hold ref*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-18 TOC: readInConFa
   '  - Reads in a reference sequence from a fasta file
   '    into a samEntry struct.
   '  - I would recomend using seqStruct for this. This
   '    was made before seqStruct, & is kept here to make
   '    sure nothing breaks
   '  o fun-18 sec-01:
   '     - variable declerations
   '  o fun-18 sec-02:
   '     - Check if file exists and read in header
   '  o fun-18 sec-03:
   '     - Read in the sequence
   '  o fun-18 sec-04:
   '     - Set up q-score null entry, get length, & return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-18 Sec-01:
    ^  - variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    char *delStr = 0;       /*for removing white space*/
    char checkC = 0;        /*Used for single char checks*/
    unsigned int incBuffUInt = 1000;
    FILE *faFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-18 Sec-02:
    ^  - Check if file exists and move past header
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open the reference file for reading*/
    faFILE = fopen(conFaToReadCStr, "r");

    if(faFILE == 0)
        return 2;

    blankSamEntry(refST);

    /*Set up null endings for lines*/
    *(refST->samEntryCStr + refST->lenBuffULng - 1) = '\0';
    *(refST->samEntryCStr + refST->lenBuffULng - 2) = '\0';

    /*Read in the header*/
    fgets(
       refST->samEntryCStr,
       refST->lenBuffULng,
       faFILE
    );

    /*Find the ending character*/
    checkC=*(refST->samEntryCStr + refST->lenBuffULng - 2);

    while(!(checkC !='\0' || checkC !='\n'))
    { /*Loop: read in the header*/
        fgets(
           refST->samEntryCStr,
           refST->lenBuffULng,
           faFILE
        ); /*Read in the next part of the header*/

        /*Check if I got file input*/
        if(*refST->samEntryCStr == '\0') return 4;

        /*Resetup markers*/
        *(refST->samEntryCStr + refST->lenBuffULng-1)='\0';
        *(refST->samEntryCStr + refST->lenBuffULng-2)='\0';

       checkC =
          *(refST->samEntryCStr + refST->lenBuffULng - 2);
    } /*Loop: read in the header*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-18 Sec-03:
    ^  - Read in the sequence
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr = refST->samEntryCStr;
    *(refST->samEntryCStr + refST->lenBuffULng - 2) = '\0';

    /*Read in the first part of the sequence*/
    tmpCStr =
       fgets(
          refST->samEntryCStr,
          refST->lenBuffULng,
          faFILE
    ); /*Get the first part of the sequence*/

    while(tmpCStr != 0 && *tmpCStr != '>')
    { /*While on the sequence line*/
        refST->samEntryCStr =
            realloc(
                refST->samEntryCStr,
                refST->lenBuffULng + incBuffUInt
            );

        if(refST->samEntryCStr == 0)
        { /*memory allocation error*/
            fclose(faFILE);
            fclose(faFILE);
            return 64;
        } /*Memory allocation error*/

        /*Set pointer to new buffer*/
        tmpCStr =
             refST->samEntryCStr
           + refST->lenBuffULng
           - 1;

        /*Update buff size*/
        refST->lenBuffULng += incBuffUInt - 1;

        /*Read in next part of reference sequence*/
        tmpCStr = fgets(tmpCStr, incBuffUInt, faFILE);
    } /*While on the sequence line*/

    fclose(faFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-18 Sec-3:
    ^  - Set up q-score null entry, get length, and return
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    refST->seqCStr = refST->samEntryCStr;

    /*Find the end of the sequence*/
    tmpCStr = refST->seqCStr;
    delStr = tmpCStr;
    refST->readLenUInt = 0; /*Reset for new length*/

    while(*tmpCStr != '\0')
    { /*Find the end of the sequence*/
        if(*tmpCStr > 32)
        { /*If: this is not white space*/
           *delStr = *tmpCStr;
           ++delStr;
           ++refST->readLenUInt;
        } /*If: this is not white space*/

        ++tmpCStr;
    } /*Find the end of the sequence*/

    *delStr = '\0';

    /*Set up Q-score entry, so score reads knows is empty*/
    refST->qCStr = tmpCStr;
    strcpy(refST->qCStr, "\t*\t");

    return 1;
} /*readInConFa*/

/*--------------------------------------------------------\
| Name: blankReadStat (Fun-19:)
| Use:
|  - Blanks (sets to 0) all values in a readStat struct
| Input:
|  - readToBlank
|    o readStat structer to blank
| Output:
|  - Modifies:
|    o readToBlank to have all stats set to 0 & all
|      c-strings have the frist character be '\0'
\--------------------------------------------------------*/
static inline void blankReadStat(
    struct readStat *readToBlank
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-19 TOC: blankReadStat
   '  - Blanks (sets to 0) all values in a readStat struct
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*Blank Ids and mapping quality*/
    readToBlank->mapqUChar = 0;
    readToBlank->queryIdCStr[0] = '\0';
    readToBlank->refIdCStr[0] = '\0';

    /*Blank the Q-score stats*/
    readToBlank->medianQFlt = 0;
    readToBlank->medianAligQFlt = 0;
    readToBlank->meanQFlt = 0;
    readToBlank->meanAligQFlt = 0;

    /*Blank the read length stats*/
    readToBlank->readLenUInt = 0;
    readToBlank->readAligLenUInt = 0;

    /*Blank the total error stats*/
    readToBlank->numMatchUInt = 0;
    readToBlank->numSNPUInt = 0;
    readToBlank->numDelUInt = 0;
    readToBlank->numInsUInt = 0;

    /*Blank the kept error stats*/
    readToBlank->numKeptMatchUInt = 0;
    readToBlank->numKeptSNPUInt = 0;
    readToBlank->numKeptDelUInt = 0;
    readToBlank->numKeptInsUInt = 0;

    return;
} /*blankReadStat*/

/*--------------------------------------------------------\
| Name: cpReadStat (Fun-20:)
| Use:
|  - Copies values from one readStat structure to another
|    read stat structuer
| Input:
|  - newReadStat:
|    o readStat structure to copy values to
|  - oldReadStat:
|    o readStat strucuter to copy values from
| Output:
|  - Modifies:
|    o newReadList to have the same values as oldReadList
\--------------------------------------------------------*/
static inline void cpReadStat(
    struct readStat *newReadStat, /*copy stats to*/
    struct readStat *oldReadStat  /*copy stats from*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-20 TOC: cpReadStat
   '  - Copies values from one readStat structure to
   '    another read stat structuer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*Blank Ids and mapping quality*/
    newReadStat->mapqUChar = oldReadStat->mapqUChar;

    strcpy(
       newReadStat->queryIdCStr,
       oldReadStat->queryIdCStr
    );

    strcpy(newReadStat->refIdCStr, oldReadStat->refIdCStr);

    /*Blank the Q-score stats*/
    newReadStat->medianQFlt = oldReadStat->medianQFlt;

    newReadStat->medianAligQFlt =
       oldReadStat->medianAligQFlt;

    newReadStat->meanQFlt = oldReadStat->meanQFlt;
    newReadStat->meanAligQFlt = oldReadStat->meanAligQFlt;

    /*Blank the read length stats*/
    newReadStat->readLenUInt = oldReadStat->readLenUInt;

    newReadStat->readAligLenUInt =
       oldReadStat->readAligLenUInt;

    /*Blank the total error stats*/
    newReadStat->numMatchUInt = oldReadStat->numMatchUInt;
    newReadStat->numSNPUInt = oldReadStat->numSNPUInt;
    newReadStat->numDelUInt = oldReadStat->numDelUInt;
    newReadStat->numInsUInt = oldReadStat->numInsUInt;

    /*Blank the kept error stats*/
    newReadStat->numKeptMatchUInt =
      oldReadStat->numKeptMatchUInt;

    newReadStat->numKeptSNPUInt =
      oldReadStat->numKeptSNPUInt;

    newReadStat->numKeptDelUInt =
      oldReadStat->numKeptDelUInt;

    newReadStat->numKeptInsUInt =
      oldReadStat->numKeptInsUInt;

    return;
} /*cpReadStat*/

/*--------------------------------------------------------\
| Name: readStatsFileLine (Fun-21:)
| Use:
|  - Reads in one line of read stats from a read stats
|    file into a readStat structure
| Input:
|  - statsFILE:
|    o File to get stats from
|  - onHeadBl:
|    o 1: I am on the header of the file
|    o 0: I am not on the header of the file
|  - readST:
|    o readStat structure to hold the stats fro statsFILE
| Output:
|  - Modifies:
|    o Puts next line of stats from statsFILE into
|      readST
|    o If onHeadBl is one, it is set to 0
|    o statsFILE to point to the next entry
\--------------------------------------------------------*/
static inline uint8_t readStatsFileLine(
    FILE *statsFILE,        /*File with stats to read*/
    uint8_t *onHeadBl,      /*1: on header, skip one line*/
    struct readStat *readST /*Holds stats from statsFILE*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-21 TOC: readStatsFileLine
   '  - Reads in one line of read stats from a read stats
   '    file into a readStat structure
   '  o fun-21 sec-01:
   '    - variable declerations
   '  o fun-21 sec-02:
   '    - Initalize and read in first line
   '  o fun-21 sec-03:
   '    - If need to, read past header
   '  o fun-21 sec-04:
   '    - Copy read and reference id's
   '  o fun-21 sec-05:
   '    - Get mapq and read lengths
   '  o fun-21 sec-06:
   '    - Get number of matches and mismatches
   '  o fun-21 sec-07:
   '    - Get number of kept insertions and deletions
   '  o fun-21 sec-08:
   '    - Get Q-scores
   '  o fun-21 sec-09:
   '    - Get aligned Q-scores
   '  o fun-21 sec-10:
   '    - Get number of SNPs, insertions, & deletions
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-01:
    ^  - variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint16_t lenBuffUSht = 1024;
       /*This is larger than one full line*/

    char buffCStr[lenBuffUSht];
    char *tmpCStr = 0;   /*String to copy from*/
    char *cpTmpCStr = 0; /*String to copy to*/
    char *eofCStr = 0;   /*Tells me if at end of file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-02:
    ^  - Initalize and read in first line
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    blankReadStat(readST); /*Make sure no leftover data*/ 

    if(statsFILE == 0)
        return 8;

    eofCStr =
        fgets(
            buffCStr,
            lenBuffUSht,
            statsFILE
    ); /*Read in the line*/

    /*Check if end of file*/
    if(eofCStr == 0) return 2;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-03:
    ^  - If need to, read past header
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*onHeadBl == 1)
    { /*If on the header, move past*/
        /*Read in the next line*/
        eofCStr = fgets(buffCStr, lenBuffUSht, statsFILE);

        if(eofCStr == 0)
            return 2;   /*End of file, so only a header*/

        *onHeadBl = 0; /*No longer on the header*/
    } /*If on the header, move past*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-04:
    ^  - Copy read and reference id's
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr = buffCStr;
    cpTmpCStr = readST->queryIdCStr;

    while(*tmpCStr > 31)
    { /*While have a read id to copy*/
        *cpTmpCStr = *tmpCStr;
        ++cpTmpCStr;
        ++tmpCStr;

        if(*tmpCStr == '\0')
            return 16; /*Early end of file or line*/
    } /*While have a read id to copy*/

    *cpTmpCStr = '\0'; /*make sure a c-string*/
    ++tmpCStr; /*Get off the tab*/

    cpTmpCStr = readST->refIdCStr;

    while(*tmpCStr > 31)
    { /*While have a reference id to copy*/
        *cpTmpCStr = *tmpCStr;
        ++cpTmpCStr;
        ++tmpCStr;

        if(*tmpCStr == '\0')
            return 16; /*Early end of file or line*/
    } /*While have a reference id to copy*/

    *cpTmpCStr = '\0'; /*make sure a c-string*/
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-05:
    ^  - Copy mapq and read lengths
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the mapq (mapping quality)*/
    tmpCStr = cStrToUChar(tmpCStr, &readST->mapqUChar);
    ++tmpCStr; /*Get off the tab*/

    /*Get read length*/
    tmpCStr = cStrToUInt(tmpCStr, &readST->readLenUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the algined read length*/
    tmpCStr = cStrToUInt(tmpCStr,&readST->readAligLenUInt);
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-06:
    ^  - Get number of matches and mismatches
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the number of matches*/
    tmpCStr = cStrToUInt(tmpCStr, &readST->numMatchUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of kept matches*/
    tmpCStr=cStrToUInt(tmpCStr, &readST->numKeptMatchUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of kept SNPs*/
    tmpCStr = cStrToUInt(tmpCStr, &readST->numKeptSNPUInt);
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-07:
    ^  - Get number of insertions and deletions
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the number of kept insertions*/
    tmpCStr = cStrToUInt(tmpCStr, &readST->numKeptInsUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of kept deletions*/
    tmpCStr = cStrToUInt(tmpCStr, &readST->numKeptDelUInt);
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-08:
    ^  - Get mean and median aligned Q-scores
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the mean Q-score*/
    sscanf(tmpCStr, "%f", &readST->meanQFlt);
    while(*tmpCStr > 32) ++tmpCStr;
    ++tmpCStr; /*Get off the tab*/

    /*Get the median Q-score*/
    sscanf(tmpCStr, "%f", &readST->medianQFlt);
    while(*tmpCStr > 32) ++tmpCStr;
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-09:
    ^  - Get median and mean aligned Q-scores
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the mean aligned Q-score*/
    sscanf(tmpCStr, "%f", &readST->meanAligQFlt);
    while(*tmpCStr > 32) ++tmpCStr;
    ++tmpCStr; /*Get off the tab*/

    /*Get the median aligned Q-score*/
    sscanf(tmpCStr, "%f", &readST->medianAligQFlt);
    while(*tmpCStr > 32) ++tmpCStr;
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-21 Sec-10:
    ^  - Get total number of SNPs, insertions, & deletions
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the number of SNPs*/
    tmpCStr = cStrToUInt(tmpCStr, &readST->numSNPUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of insertions*/
    tmpCStr = cStrToUInt(tmpCStr, &readST->numInsUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of deletions*/
    tmpCStr = cStrToUInt(tmpCStr, &readST->numDelUInt);
    ++tmpCStr; /*Get off the tab*/

    return 1; /*Sucess*/
} /*readStatsFileLine*/

/*--------------------------------------------------------\
| Name: printReadStat (Fun-22:)
| Use:
|  - Prints the stats in a read stat struct to an file
| Input:
|  - readToPrint:
|    o readStat structure with read stats to print out
|  - outFILE:
|    o File to print stats to
| Output:
|  - Prints:
|    o Prints out variables in readToPrint structure to
|      outFILE
\--------------------------------------------------------*/
static inline void printReadStat(
    struct readStat *readToPrint,/*Read stats to print*/
    FILE *outFILE                /*File to print stats to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-22 TOC: printReadStat
   '  - Prints the stats in a read stat struct to an file
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(readToPrint == 0)
        return;

    if(outFILE == 0)
        return;

    /*Print out the entry stats*/
    fprintf(
        outFILE,
        "%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%f\t%f",
        readToPrint->queryIdCStr, 
        readToPrint->refIdCStr,
        readToPrint->mapqUChar,
        readToPrint->readLenUInt,
        readToPrint->readAligLenUInt,
        readToPrint->numMatchUInt,
        readToPrint->numKeptMatchUInt,
        readToPrint->numKeptSNPUInt,
        readToPrint->numKeptInsUInt,
        readToPrint->numKeptDelUInt,
        readToPrint->medianQFlt,
        readToPrint->meanQFlt
    ); /*1st printf: print out stats*/

    fprintf(
        outFILE,
        "\t%f\t%f\t%u\t%u\t%u\n",
        readToPrint->medianAligQFlt,
        readToPrint->meanAligQFlt,
        readToPrint->numSNPUInt,
        readToPrint->numInsUInt,
        readToPrint->numDelUInt
    ); /*2nd fprintf: print out stats*/

    return;
} /*printReadStat*/

/*--------------------------------------------------------\
| Name: samEntryToReadStat (Fun-23:)
| Use:
|  - Copies stats in a samEntry structure to a readStat
|    structure
| Input:
|  - newBin:
|    o readStat structs to hold the stats in samStruct
|  - samStruct:
|    o samEntry structure with stats to copy
| Output:
|  - Modifies:
|    o newBin to have stats in samStruct
\--------------------------------------------------------*/
static inline void samEntryToReadStat(
    struct readStat *newBin,   /*Copy stats to*/
    struct samEntry *samStruct /*Copy stats from*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-23 TOC: samEntryToReadStat
   '  - Copies stats in a samEntry structure to a
   '    readStat structure
   '  o fun-23 sec-02:
   '    - Copy stats
   '  o fun-23 sec-01:
   '    - Variable declerations
   '  o fun-23 sec-03:
   '    - Copy the query id
   '  o fun-23 sec-04:
   '    - Copy the reference id
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-23 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    char *tmpCpCStr = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-23 Sec-02:
    ^  - Copy stats
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    newBin->mapqUChar = samStruct->mapqUChar;

    newBin->medianQFlt = samStruct->medianQFlt;
    newBin->medianAligQFlt = samStruct->medianAligQFlt;
    newBin->meanQFlt = samStruct->meanQFlt;
    newBin->meanAligQFlt = samStruct->meanAligQFlt;

    
    newBin->readLenUInt = samStruct->readLenUInt;
    newBin->readAligLenUInt = samStruct->readAligLenUInt;

    newBin->numMatchUInt = samStruct->numMatchUInt;
    newBin->numSNPUInt = samStruct->numSNPUInt;
    newBin->numKeptSNPUInt = samStruct->numKeptSNPUInt;
    newBin->numDelUInt = samStruct->numDelUInt;
    newBin->numKeptDelUInt = samStruct->numKeptSNPUInt;
    newBin->numInsUInt = samStruct->numInsUInt;
    newBin->numKeptInsUInt = samStruct->numKeptSNPUInt;

    tmpCpCStr = newBin->queryIdCStr;
    tmpCStr = samStruct->queryCStr;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-23 Sec-03:
    ^  - Copy the query id
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*tmpCStr > 16)
    { /*While have the query id to copy over*/
        *tmpCpCStr = *tmpCStr;
        ++tmpCStr;
        ++tmpCpCStr;
    } /*While have the query id to copy over*/

    *tmpCpCStr = '\0'; /*Mark end of c-string*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-23 Sec-04:
    ^  - Copy the reference id
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCpCStr = newBin->refIdCStr;
    tmpCStr = samStruct->refCStr;

    while(*tmpCStr > 16)
    { /*While have the reference id to copy over*/
        *tmpCpCStr = *tmpCStr;
        ++tmpCStr;
        ++tmpCpCStr;
    } /*While have the reference id to copy over*/

    *tmpCpCStr = '\0'; /*Mark end of c-string*/
    return;
} /*samStructToReadBin*/

#endif

/*
Sam file table for first 11 columns (all sam files have)
| Col | Field |  Type  |        Brief description              |
|:---:|:-----:|:------:|:-------------------------------------:|
|  1  | QNAME | String |       Query template NAME             |
|  2  | FLAG  |  Int   |          bitwise FLAG                 |
|  3  | RNAME | String |     Reference sequence NAME           |
|  4  |  POS  |  Int   |  1-based leftmost mapping POSition    |
|  5  | MAPQ  |  Int   |          MAPping Quality              |
|  6  | CIGAR | String |            CIGAR string               |
|  7  | RNEXT | String | Reference name of the mate/next read  |
|  8  | PNEXT |  Int   |   Position of the mate/next read      |
|  9  | TLEN  |  Int   |      observed Template LENgth         |
| 10  |  SEQ  | String |          segment SEQuence             |
| 11  | QUAL  | String | ASCII of Phred-scaled base Quality+33 |
*/

/*eqx cigar entry
`   matches are '=' or blank (no symbol) at end or if only
`     matches
`   mimsatches are 'X'
`   deletions are 'D'
`   insertions are 'I'
`   soft masks are 'S'
`   Numbers come before entry (=, X, D, I, or S) & show
`     number of times the entry is repeated
`   Everything else is a hard mask (was removed)

`   EX: 10S5=1X701 (10 soft masked bases, 5 matches,
`       1 SNP, 701 matches)
*/

/* Sam file flag values tables
| Bit  | FLAG  |                        Description                                 |
|:----:|:-----:|:------------------------------------------------------------------:|
| 1    | 0x1   | template having multiple segments in sequencing                    |
| 2    | 0x2   | each segment properly aligned according to the aligner             |
| 4    | 0x4   | segment unmapped                                                   |
| 8    | 0x8   | next segment in the template unmapped                              |
| 16   | 0x10  | SEQ being reverse complemented                                     |
| 32   | 0x20  | SEQ of the next segment in the template being reverse complemented |
| 64   | 0x40  | the first segment in the template                                  |
| 128  | 0x80  | the last segment in the template                                   |
| 256  | 0x100 | secondary alignment                                                |
| 512  | 0x200 | not passing filters, such as platform/vendor quality controls      |
| 1024 | 0x400 | PCR or optical duplicate                                           |
| 2048 | 0x800 | supplementary alignment                                            |
*/
