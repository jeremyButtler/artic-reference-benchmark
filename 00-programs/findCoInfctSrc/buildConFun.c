/*######################################################################
# Use:
#   o Holds functions invovled in building consensus or working with
#     consensuses
######################################################################*/

#include "buildConFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF TOC: Start of Functions
'  o fun-1 buildCon:
'    - Builds a consensus from a fastq file using a
'      supplied reference or the best read (highest
'      medain Q) in the fastq file.
'  o fun-3 buildSingleCon:
'    - Builds a consensus using: the selected consensus
'      methods (ivar, medaka, racon, majCon)
'  o fun-4 simpleMajCon:
'    - Builds a majority consensus (ignores insertions)
'  o fun-5 buildConWithRacon:
'    - Buids a consensus using racon
'  o fun-6 medakaPolish:
'    - Polish a consensus with medaka using the best reads
'  o fun-7 cmpCons:
'    - Compares two consensus(recursive call if tree input)
'  o fun-9 freeBaseStruct:
'    - Frees a base struct from a linked list.
'    - The freeded base is set to the alternate base if
'      their is an alternate base
'    - The freed base is set to 0 if no alternative base is
'      present.
'  o fun-9 initMajConStruct:
'    - Set default settings for struct holding majority
'      consensus settings
'  o fun-10 initRaconStruct:
'    - Set default settings for Racon settings
'  o fun-10 initMedakaStruct:
'    - Set default settings for Racon settings
'  o fun-11 initMedakaStruct:
'    - Set default settings for Medaka settings
'  o fun-12 initConBuildStruct:
'    - Set default settings for struct holdoing consensus
'      building settings.
'  o fun-13 ivarCon:
'    - Builds a consensus using ivar
'  o macro-01 initIvarSet:
'    - Initializes the settings for an ivarSet structuer
'    - in buildConFun.h only
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "buildConFun.h"

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
) /*Builds a consensus using input fastq file & best read or referece*/
{ /*buildCon*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: buildCon
    '   fun-1 sec-1: Variable declerations
    '   fun-1 sec-2: If using a reference set up reference for polishing
    '   fun-1 sec-3: Make sure I know the number of reads in fastq file
    '   fun-1 sec-4: If not using consensus, find a good quality read
    '   fun-1 sec-5: Polish the read or reference
    '   fun-1 sec-6: Copy the best read back into the fastq file
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    char polishBl = 0;      /*Marks if need to build initial con*/
    char trueBl = 1;        /*So I can pass TRUE (1) to my functions*/
    char falseBl = 0;       /*So I can pass a FALSE (0) to functions*/
    char tmpBuffCStr[1024];
    uchar errUC = 0;
    ulong totalReadsUL = 0; /*Number reads in fastq file*/
    ulong mappedReadsUL =0;
      /*Reads that mapped and passed the quality checks,
      ` but may not be in the subsample
      */

    float percMappedFlt = 0; /*% of reads mapped*/
    struct samEntry *zeroSam = 0; /*holds the reference (0 to ignore)*/

    FILE *fqFILE = 0;
    FILE *bestReadFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-2: If using a reference set up reference for polishing
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refCStr != 0)
    { /*If have a reference to work with*/
        tmpCStr =
            cStrCpInvsDelm(conData->consensusCStr, conData->fqPathCStr);

        /*Get to end of .fastq or .fq*/
        if(*(tmpCStr - 2) == 't') tmpCStr -= 6; /*.fastq*/
        else tmpCStr -= 3;                      /*.fq*/

        /*Set up the best read name*/
        strcpy(conData->bestReadCStr, conData->consensusCStr);
        fqFILE = fopen(refCStr, "r");

        /*Need to make sure that I do not delete the orginal reference*/
        if(fqFILE == 0)
            return 2;

        bestReadFILE = fopen(conData->consensusCStr, "w");

        if(bestReadFILE == 0)
        { /*If have no write permision*/
            fclose(fqFILE);
            return 2;
        } /*If have no write permision*/

        while(fgets(tmpBuffCStr, 1024, fqFILE))
            fprintf(bestReadFILE, "%s", tmpBuffCStr);

        fclose(fqFILE);
        fclose(bestReadFILE);

        polishBl = 1; /*Polish the reference up*/
    } /*If have a reference to work with*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-3: Make sure I know the number of reads in fastq file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/


    if(conData->numReadsULng == 0)/*Find number of reads in fastq file*/
        conData->numReadsULng = getNumReadsInFq(conData->fqPathCStr);

    if(conData->numReadsULng == 0)
        return 2;

    totalReadsUL = conData->numReadsULng;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-5: Polish the read or reference
    ^   fun-1 sec-5 sub-1: Build the initial consensus if needed
    ^   fun-1 sec-5 sub-2: Polish the consensus
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-1 Sec-5 Sub-1: Build the initial consensus if needed
    \******************************************************************/

    while(!(errUC & 1))
    { /*While have not built a consensus*/

        if(conData->numReadsULng < conSet->minReadsToBuildConUL)
            return 16; /*To few reads to build a consensus*/

        if(!(polishBl & 1))
        { /*if need to find another read*/
            if(conSet->useStatBl & 1)
            { /*If using a stats file to extract reads*/
                if(conData->statPathCStr[0] == '\0')
                    return 8; /*If no stat file to extract reads with*/

                errUC = 128;

                while(errUC & 128)
                    errUC = extractBestRead(conData);
                   /*Want to loop till I get a read extracted. This is
                     needed, since the bin to consensus step does not
                     keep the stats file up to date*/
            } /*If using a stats file to extract reads*/
            
            else
                /*Extracts best read by medain Q-score*/
                errUC =
                    fqGetBestReadByMedQ(conData, samStruct,bestReadSam);

            if(!(errUC & 1))
            { /*If had an error*/
                remove(conData->bestReadCStr);
                return errUC;/*4, no Fq file, 8 write error, 64 memory*/
            } /*If had an error*/

            errUC = 
                findBestXReads(
                    &conSet->maxReadsToBuildConUL,
                    &conSet->numReadsForConUL, /*# of reads extracted*/
                    &mappedReadsUL,
                    threadsCStr,      /*Number threads for minimap2*/
                    &falseBl,        /*Do not using mapping quality*/
                    minReadReadStats, /*Min stats to keep reads*/
                    samStruct,  /*Struct to use for reading sam file*/
                    zeroSam,    /*Do not use reference in scoring*/
                    conData,
                    0,          /*Use the reference*/
                    1           /*Make a name using the input fastq*/
            );  /*Extract top reads that mapped to selected best read*/

            percMappedFlt=mappedReadsUL/(totalReadsUL*1.0);

            if(
                     conSet->numReadsForConUL 
                   < conSet->minReadsToBuildConUL
                 ||
                     percMappedFlt
                   < conSet->minPercMappedReadsFlt
            ){ /*If I did not extract enough reads*/
                errUC = 16;
                continue; /*If not enough reads to build consensus*/
            } /*If I did not extract enough reads*/

            /*Build consensus using best read & top reads in conData*/
            errUC =buildSingleCon(threadsCStr,conData,samStruct,conSet);

            if(errUC & 32)
            { /*If the cosensus is to small*/
                errUC = 16;
                continue;
            } /*If the cosensus is to small*/
        } /*if need to find another read*/

        /**************************************************************\
        * Fun-1 Sec-5 Sub-2: Polish the consensus
        \**************************************************************/

        /*Make sure bes read is fasta file, this will also save the
          best read from being wiped*/
        tmpCStr = conData->bestReadCStr;

        while(*tmpCStr != '\0')
            ++tmpCStr;
        --tmpCStr;
        *tmpCStr = 'a'; /*change fastq to fasta*/

        for(uint32_t rndUI = 0;rndUI< conSet->numRndsToPolishUI;++rndUI)
        { /*Loop till have done all the users requested polishing*/
            /*Set up the consnesus as the next best read*/
            fqFILE = fopen(conData->bestReadCStr, "r");

            if(fqFILE != 0 && !(polishBl & 1))
            { /*If need to remove the best read file*/
                fclose(fqFILE);
                fqFILE = 0;
                remove(conData->bestReadCStr);
            } /*If need to remove the best read file*/

            else if(fqFILE != 0)
                fclose(fqFILE); /*Make sure the file is closes*/

            rename(conData->consensusCStr, conData->bestReadCStr);
            conData->consensusCStr[0] = '\0'; /*Remove old name*/

            /*Extract the reads for the next rebuild*/
            errUC = 
                findBestXReads(
                    &conSet->maxReadsToBuildConUL,
                    &conSet->numReadsForConUL, /*# of reads extracted*/
                    &mappedReadsUL,
                    threadsCStr,     /*Number threads for minimap2*/
                    &trueBl,        /*will use mapq for reads here*/
                    minReadConStats,/*Min stats to keep reads*/
                    samStruct,  /*Struct to use for reading sam file*/
                    zeroSam,    /*Do not use reference in scoring*/
                    conData,
                    0,          /*Use the reference*/
                    1           /*Make a name using the input fastq*/
            );  /*Extract top reads that mapped to selected best read*/

            percMappedFlt=mappedReadsUL/(totalReadsUL*1.0);

            if(
                     conSet->numReadsForConUL 
                   < conSet->minReadsToBuildConUL
                 ||
                     percMappedFlt
                   < conSet->minPercMappedReadsFlt
            ){ /*If need to get a new best read*/
                polishBl = 0; /*do best read if reference fails*/
                errUC = 16; /*Could not build a consensus*/
                break; /*Get a new best read*/
            } /*If need to get a new best read*/

            if(errUC & 64)
                return 64; /*Memory allocation error*/

            /*Re-build consensus using consensus & top reads*/
            errUC =buildSingleCon(threadsCStr,conData,samStruct,conSet);

            if(errUC & 32)
            { /*If the cosensus is to small*/
                errUC = 16;
                break;
            } /*If the cosensus is to small*/

            if(errUC & 16)
                break; /*minimap2 failed, Get a new best read*/

            if(errUC & 64)
            { /*If I had a memory allocation error*/
                if(polishBl == 0)
                    remove(conData->bestReadCStr);

                return 64; /*Memory allocation error*/
            } /*If I had a memory allocation error*/
        } /*Loop till have done all the users requested polishing*/

        if(!(errUC & 1))    /*If need to do another round*/
            polishBl = 0; /*do best read if reference fails*/

        remove(conData->bestReadCStr); /*Remove temporary consensus*/
        *tmpCStr = 'q'; /*change best read back to fastq*/
    } /*While have not built a consensus*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-6: Copy the best read back into the fastq file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(!(errUC & 1))
    { /*If I could not build a consensus*/
        if(polishBl == 0)
            remove(conData->bestReadCStr);

        return errUC;  /*Failed to build a consensus*/
    } /*If I could not build a consensus*/

    if(polishBl == 0)
    { /*If had a best read file*/
        bestReadFILE = fopen(conData->bestReadCStr, "r");
        fqFILE = fopen(conData->fqPathCStr, "a");

        while(fgets(tmpBuffCStr, 1024, bestReadFILE))
            fprintf(fqFILE, "%s", tmpBuffCStr);

        ++conData->numReadsULng; /*Account for adding back best read*/

        fclose(fqFILE);
        fclose(bestReadFILE);
        remove(conData->bestReadCStr);
    } /*If had a best read file*/

    return 1;
} /*buildCon*/

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
|      - 32 if the consensus was to small
|      - 64 for memory allocation error
\---------------------------------------------------------------------*/
unsigned char buildSingleCon(
    char *threadsCStr,       /*Number threads to use with system calls*/
    struct readBin *clustOn, /*fastq file & reads to build consensus*/
    struct samEntry *samStruct,/*For reading in sequences or sam file*/
    struct conBuildStruct *conSet
        /*Settings to use while building the consensus*/
) /*Builds a consensus using the best read & top reads in clustOn*/
{ /*buildSingleCon*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-3 TOC: buildCon
    '  o fun-3 sec-1:
    '    - Variable declerations
    '  o Fun-3 Sec-02:
    '    - Buid a consensus
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-1:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t errUC = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-02:
    ^  - Buid a consensus
    ^  o fun-3 sec-02 sub-01:
    ^    - Check which methods using to build a consensus
    ^  o fun-3 sec-02 sub-02:
    ^    - Build the consensus with majority consensus
    ^  o fun-3 sec-02 sub-03:
    ^    - Build the consensus with racon
    ^  o fun-3 sec-2 sub-04:
    ^    - Build the consensus with medaka
    ^  o fun-3 sec-2 sub-05:
    ^    - Build consensus with ivar
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-3 Sec-2 Sub-01:
    *  - Check which methods using to build a consensus
    \*****************************************************/

    for(uchar iRnd = 0; iRnd < conSet->lenMethodUC; ++iRnd)
    { /*Loop: Build consensus according to user input*/
       switch(conSet->methodAryUC[iRnd])
       { /*Switch: Check wich consensus method I am using*/
          case defNoCon: goto builtOneCon;

          /***********************************************\
          * Fun-3 Sec-2 Sub-02:
          *  - Build majority consensus
          \***********************************************/

          case defUseMajCon:
          /*Case: Using the majority consensus method*/
              errUC =
                simpleMajCon(
                    &conSet->clustUC,
                    threadsCStr,
                    clustOn,
                    samStruct,
                    &conSet->majConSet
              ); /*Build the consensus*/

              if(!(errUC & 1))
                  return errUC;

              conSet->lenConUL =conSet->majConSet.lenConUL;

              if(conSet->lenConUL < conSet->minConLenUI)
                  return 32; /*The consensus was to short*/
             break;
          /*Case: Using the majority consensus method*/
    
          /***********************************************\
          * Fun-3 Sec-02 Sub-03:
          *  - Build consensus with racon
          \***********************************************/

          case defUseRacon:
          /*Case: using the racon method*/
              buildConWithRacon(
                  &conSet->raconSet,
                  &conSet->clustUC,
                  threadsCStr,
                  clustOn
              ); /*Builds a consensus using racon*/

              conSet->lenConUL = conSet->raconSet.lenConUL;

              if(conSet->lenConUL < conSet->minConLenUI)
                  return 32; /*The consensus was to short*/
              break;
          /*Case: using the racon method*/

          /***********************************************\
          * Fun-3 Sec-02 Sub-04:
          *  - Build consensus with medaka
          \***********************************************/

           case defUseMedaka:
           /*Case: Use the medaka consensus method*/
              errUC =
                  medakaPolish(
                      &conSet->medakaSet,
                      &conSet->clustUC,
                      threadsCStr,
                      clustOn,
                      samStruct
              ); /*Build the consensus with medaka*/

              if(!(errUC & 1)) return errUC;

              conSet->lenConUL= conSet->medakaSet.lenConUL;

              if(conSet->lenConUL < conSet->minConLenUI)
                  return 32; /*The consensus was to short*/
              break;
           /*Case: Use the medaka consensus method*/

          /***********************************************\
          * Fun-3 Sec-02 Sub-05:
          *  - Build consensus with ivar
          \***********************************************/

           case defUseIvar:
           /*Case: Use the ivar consensus method*/
              errUC =
                  ivarCon(
                      &conSet->ivarSetST,
                      &conSet->clustUC,
                      threadsCStr,
                      clustOn,
                      samStruct
              ); /*Build the consensus with medaka*/

              if(!(errUC & 1)) return errUC;

              conSet->lenConUL=conSet->ivarSetST.lenConUL;

              if(conSet->lenConUL < conSet->minConLenUI)
                  return 32; /*The consensus was to short*/
              break;
           /*Case: Use the ivar consensus method*/
       } /*Switch: Check wich consensus method I am using*/
    } /*Loop: Build consensus according to user input*/

    builtOneCon:

    return 1;
} /*buildSingleCon*/

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
    struct readBin *binStruct, /*Has best read & top reads files*/
    struct samEntry *samStruct,/*For reading in sam file entries*/
    struct majConStruct *settings
        /*Has settings for building the consensus*/
) /*Builds a majority consensus from the best reads & top read*/
{ /*simpleMajCon*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-4 TOC: simpleMajCon
    '    fun-4 sec-1: Variable declerations
    '    fun-4 sec-2: Check if the bestRead and topReads files exist
    '        - Also reads in the reference sequence
    '    fun-4 sec-3: Prepare the minimap2 command
    '    fun-4 sec-4: Make consensus array & initalize with reference
    '    fun-4 sec-5: Map reads to the reference read
    '    fun-4 sec-7: Print out cosensus & do clean up
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char qEntryBl = 0;           /*Marks if reference has Q-core entry*/
    char minimap2CmdCStr[2048];  /*Holds minimap2 command to run*/
    char *tmpCStr = 0;           /*Temp ptr for c-string manipulations*/
    char *cigCStr = 0;           /*Reading the cigar entry*/
    char *seqCStr = 0;           /*Manipulating/reading the sequence*/
    char *qCStr = 0;             /*Manipulating/reading q-score entry*/

    unsigned char qScoreUChar = 0; /*Holds the Q-score for a base*/
    uint8_t errUChar = 0;        /*Holds error messages from functions*/

    uint32_t cigEntryUInt = 0;  /*Holds number of bases in cigar entry*/

    unsigned long minNumBasesUL = 0; /*Min read support to keep a base*/
    unsigned long minInsUL = 0;   /*Min read support to keep insertion*/
    unsigned long numSupReadsUL = 0; /*NumberOfReads supporting a base*/
    unsigned long numSeqUL = 0;      /*Number of mapped sequences*/
    unsigned long numMisSeqUL = 0;   /*Number of mapped sequences*/

    struct baseStruct *headBase = 0; /*Head of the list of bases*/
    struct baseStruct *incBase = 0;  /*First base at a position*/
    struct baseStruct *lastBase = 0; /*Base before tmpBase*/
    struct baseStruct *tmpBase = 0;  /*Base position working on*/

    FILE *stdinFILE = 0;        /*For reading and writing files*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-2: Check if the bestRead and topReads files exist
    ^    - Also reads in the reference sequence
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    stdinFILE = fopen(binStruct->bestReadCStr, "r");

    if(stdinFILE == 0)
        return 2;

    blankSamEntry(samStruct); /*Make sure start with blank*/
    tmpCStr = binStruct->bestReadCStr;

    while(*tmpCStr != '\0')
        ++tmpCStr;

    if(*(tmpCStr - 1) == 'q')
    { /*If is a fastq file*/
        errUChar = readRefFqSeq(stdinFILE, samStruct, 1); /*fastq file*/
        qEntryBl = 1;
    } /*If is a fastq file*/

    else if(*(tmpCStr - 1) == 'a')
    { /*Else the best read is a fasta file*/
        errUChar = readInConFa(binStruct->bestReadCStr, samStruct);
        qEntryBl = 0;
    } /*Else the best read is a fasta file*/

    else
    { /*Else reference is the consensus*/
        errUChar = readInConFa(binStruct->consensusCStr, samStruct);
        qEntryBl = 0;
    } /*Else reference is the consensus*/

    if(errUChar & 64)
    { /*if had a memory allocation error*/
        fclose(stdinFILE);
        return 64; /*Memory allocation error*/
    } /*if had a memory allocation error*/

    if(!(errUChar & 1))
    { /*If had an issue with the consensus*/
        fclose(stdinFILE);
        return 2; /*issue with consensus*/
    } /*If had an issue with the consensus*/

    /*Make sure have doube the sequence + q-score size for insertions*/
    if(samStruct->readLenUInt * 2 > samStruct->lenBuffULng)
    { /*If need to increase the size of the buffer*/
        samStruct->samEntryCStr =
            realloc(
                samStruct->samEntryCStr,
                sizeof(char) * samStruct->readLenUInt * 2
        ); /*Resize the structs buffer*/

        if(samStruct->samEntryCStr == 0)
        { /*If had a memory allocation wrror*/
            fclose(stdinFILE);
            return 64;
        } /*If had a memory allocation wrror*/
    } /*If need to increase the size of the buffer*/

    /*Make the consensus array*/
    fclose(stdinFILE);

    if(errUChar & 64)
        return 64;
    if(!(errUChar & 1))
        return 4; /*Invalide reference*/

    /*Check if can open the top reads file*/
    stdinFILE = fopen(binStruct->topReadsCStr, "r");

     if(stdinFILE == 0)
        return 4;
     fclose(stdinFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-3: Prepare the minimap2 command
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr =
        cStrCpInvsDelm(binStruct->consensusCStr, binStruct->fqPathCStr);

    /*Get to end of .fastq or .fq*/
    if(*(tmpCStr - 2) == 't') tmpCStr -= 6; /*.fastq*/
    else tmpCStr -= 3;                      /*.fq*/

    tmpCStr = cStrCpInvsDelm(tmpCStr, "--clust-");
    tmpCStr = uCharToCStr(tmpCStr, *clustUC); /*Add cluster number*/
    tmpCStr=cStrCpInvsDelm(tmpCStr, "--con.fasta");

    /*Prepare the minimap2 command*/
    tmpCStr = cStrCpInvsDelm(minimap2CmdCStr, minimap2CMD);
    tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
    cpParmAndArg(
        tmpCStr,
        binStruct->bestReadCStr,
        binStruct->topReadsCStr
    );

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-4: Make the consensus array & initalize with reference
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    seqCStr = samStruct->seqCStr;
    qCStr = samStruct->qCStr;
    tmpBase = 0;

    /*Initalize base struct using the reference read*/

    while(*seqCStr > 16)
    { /*While have bases to read in*/
        if(*seqCStr < 33)
            continue;

        tmpBase = malloc(sizeof(struct baseStruct));

        if(tmpBase == 0)
        { /*If had a memory allocation error*/
            lastBase = 0;       /*So function knows no previous bases*/

            /*Free all the made bases*/
            while(headBase != 0)
                headBase = freeBaseStruct(&headBase, lastBase);

            return 64;
        } /*If had a memory allocation error*/

        if(headBase == 0)
            headBase = tmpBase;           /*If is the first base*/
        else
            lastBase->nextBase = tmpBase; /*Set up list*/

        tmpBase->nextBase = 0; /*Only matches at this point*/
        tmpBase->altBase = 0;  /*No inserts at this point*/
        tmpBase->baseChar = *seqCStr;    /*Record the base*/
        tmpBase->errTypeChar = 0;        /*match = 0, snp=1, ins=2*/

        qScoreUChar = *qCStr - Q_ADJUST; /*Get the Q-score*/

        /*qEntryBl ensures if only files if their was a Q-score entry*/
        if(qEntryBl && qScoreUChar < settings->minBaseQUC)
        { /*If base is to low of quality to keep, make a blank struct*/
            tmpBase->errTypeChar = 0;  /*match=0, snp=1, ins=2*/
            tmpBase->numSupReadsUL = 0;
        } /*If base is to low of quality to keep, make a blank struct*/

        else
            tmpBase->numSupReadsUL = 1;        /*match = 0, snp=1, ins=2*/

        lastBase = tmpBase;
        ++seqCStr;
        ++qCStr;
    } /*While have bases to read in*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-5: Map reads to the reference read
    ^    fun-4 sec-5 sub-1: run miniamp2 and read in first line
    ^    fun-4 sec-5 sub-2: Read in sequence and position on first base
    ^    fun-4 sec-5 sub-3: Add matches & SNPs to consensus array
    ^    Fun-4 Sec-5 sub-4: Get base matching cigar (match/snp or ins)
    ^    fun-4 sec-5 sub-5: Check if should keep base or discard
    ^    fun-4 sec-5 sub-6: Find matching base at position
    ^    fun-4 sec-5 sub-7: Make new base (no matching) or update count
    ^    fun-4 sec-5 sub-8: Make sure on an insertion
    ^    fun-4 sec-5 sub-9: Check if should keep ins
    ^    fun-4 sec-5 sub-10: Find matching ins at pos
    ^    fun-4 sec-5 sub-11: check if need to make ins
    ^    fun-4 sec-5 sub-12: Ingnore deletions
    ^    fun-4 sec-5 sub-13: Ingnore soft masking
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-4 Sec-5 Sub-1: run miniamp2 and read in first line
    \******************************************************************/

    stdinFILE = popen(minimap2CmdCStr, "r");
    blankSamEntry(samStruct); /*Make sure start with blank*/

    /*Read in a single sam file line to check if valid (header)*/
    errUChar = readSamLine(samStruct, stdinFILE);

    if(!(errUChar & 1))
    { /*If an error occured*/
        pclose(stdinFILE);
        return 32;
    } /*If an error occured*/

    /******************************************************************\
    * Fun-4 Sec-5 Sub-2: Read in sequence and position on first base
    \******************************************************************/

    while(errUChar & 1)
    { /*While have alignments to read in from the sam file*/
        if(*samStruct->samEntryCStr == '@')
        { /*If on a header entry, read in next entry*/
            blankSamEntry(samStruct);
            errUChar = readSamLine(samStruct, stdinFILE);
            continue;
        } /*If on a header entry, read in next entry*/

        cigCStr = samStruct->cigarCStr;
        seqCStr = samStruct->seqCStr;
        qCStr = samStruct->qCStr;

        if(*seqCStr == '*' || (*qCStr == '*' && *(qCStr + 1) == '\t'))
        { /*If no entry to check*/
            blankSamEntry(samStruct);
            errUChar = readSamLine(samStruct, stdinFILE);
            continue;
        } /*If no entry to check*/

        if(samStruct->flagUSht & 4)
        { /*If was an unampped read*/
            blankSamEntry(samStruct);
            errUChar = readSamLine(samStruct, stdinFILE);
            ++numMisSeqUL;
            ++numSeqUL;
            continue;
        } /*If was an unampped read*/

        if(samStruct->flagUSht & 2048)
        { /*If was an supplemental read*/
            blankSamEntry(samStruct);
            errUChar = readSamLine(samStruct, stdinFILE);
            continue; /*Do not bother with supplemental alignments*/
            /*My logic is that this is PCR data, so any supplemental
              alignment is a repeat region, from rolling circle
              replication (RCR), or is a chimera. For repeats I am just
              introducing noise, RCR messes me up to much, and chimeras
              I want to ignore.*/
        } /*If was an supplemental read*/

        /*Get first base in the list*/
        ++numSeqUL; /*Count the number of non-supplemental reads*/
        incBase = headBase;
        cigEntryUInt = samStruct->posOnRefUInt - 1;
            /*-1 for 1 index for posOnRef, but 0 index fo incBase*/

        while(cigEntryUInt > 0)
        { /*While not on the starting base*/
            incBase = incBase->nextBase;
            --cigEntryUInt;
        } /*While not on the starting base*/

        while(*cigCStr != '\t')
        { /*While not at the end of the sam alignment sequence*/
            /*Get the cigar entry*/
            readCigEntry(&cigCStr, &cigEntryUInt);

        /**************************************************************\
        * Fun-4 Sec-5 Sub-3: Add matches & SNPs to consensus array
        \**************************************************************/

            switch(*cigCStr)
            { /*switch: check the error type & add bases to consensus*/
                case 'X':              /*snp, similar loop to match*/
                case '=':              /*Match, similar loop to snp*/
                case '\t':             /*Match at end of cigar*/
                /*Switch: Match, end of cigar match, snp, or insertion*/
                    while(cigEntryUInt > 0)
                    { /*While have bases to add to the bases list*/
                        --cigEntryUInt;

                        /**********************************************\
                        * Fun-4 Sec-5 Sub-4: Find first non-insertion
                        \**********************************************/

                        while(incBase->errTypeChar & 2)
                        { /*Get off the insertion*/
                            lastBase = incBase;
                            incBase = incBase->nextBase;
                        } /*Get off the insertion*/

                        /**********************************************\
                        * Fun-4 Sec-5 Sub-5: Check if should keep base
                        \**********************************************/

                        tmpBase = incBase; /*Base position on*/
                        qScoreUChar = *qCStr - Q_ADJUST; /*get Q-score*/

                        if(qScoreUChar < settings->minBaseQUC)
                        { /*If the base is low quality, skip*/
                            /*Move to the next base in the list*/
                            lastBase = incBase;
                            incBase = incBase->nextBase;

                            /*Move to the next base in the sam entry*/
                            ++seqCStr;
                            ++qCStr;
                            continue; /*No support for position*/
                        } /*If the base is low quality, skip*/

                        /**********************************************\
                        * Fun-4 Sec-5 Sub-6: Find alt matching base
                        \**********************************************/

                        while(tmpBase->baseChar != *seqCStr)
                        { /*While the bases are not equal*/
                            if(tmpBase->baseChar == 0)
                                break; /*Empty base*/

                            lastBase = tmpBase;
                            tmpBase = tmpBase->altBase;

                            if(tmpBase == 0)
                                break; /*Need to create a new base*/
                        } /*While the bases are not equal*/

                        /**********************************************\
                        * Fun-4 Sec-5 Sub-7: Make or update alt base
                        \**********************************************/

                        if(tmpBase == 0)
                        { /*If need to create a new base*/
                            lastBase->altBase =
                                malloc(sizeof(struct baseStruct));

                            tmpBase = lastBase->altBase;
                            tmpBase->nextBase = 0; /*Is an alterantive*/
                            tmpBase->altBase = 0;  /*No alternatives*/

                            tmpBase->numSupReadsUL = 0;
                            tmpBase->baseChar = *seqCStr;
                            tmpBase->errTypeChar = 0;
                        } /*If need to create a new base*/

                        else if(tmpBase->baseChar == 0)
                            tmpBase->baseChar = *seqCStr;

                        ++tmpBase->numSupReadsUL;

                        /*Move to the next base*/
                        lastBase = incBase;
                        incBase = incBase->nextBase;
                        ++seqCStr;
                        ++qCStr;
                    } /*While have bases to add to the bases list*/
            
                    break;
                /*Switch: Match, end of cigar match, snp, or insertion*/


                case 'I':   /*Is an insertion*/
                /*Switch: for insertions*/

                    while(cigEntryUInt > 0)
                    { /*While have bases to add to the bases list*/
                        --cigEntryUInt;

                        /**********************************************\
                        * Fun-4 Sec-5 Sub-8: Make sure on an insertion
                        \**********************************************/

                        /*Insertions at the ends will likely be treated
                          as softmasks by minimap2. This means that
                          incBase will != 0 for insertions
                        */
                        if(!(incBase->errTypeChar & 2))
                        { /*If next base is not an insertion*/
                            tmpBase=malloc(sizeof(struct baseStruct));

                            if(lastBase->errTypeChar & 2)
                            { /*If the last base was an ins*/ 
                                tmpBase->nextBase = incBase->nextBase;
                                incBase->nextBase = tmpBase;
                                incBase = tmpBase;
                            } /*If the last base was an ins*/ 

                            else
                            { /*else the last base was not an ins*/ 
                                tmpBase->nextBase = incBase;
                                lastBase->nextBase = tmpBase;
                                incBase = tmpBase;
                            } /*else the last base was not an ins*/ 

                            /*Set defualts for the new base*/
                            incBase->altBase = 0;
                            incBase->errTypeChar = 2;
                            incBase->baseChar = 0;
                            incBase->numSupReadsUL = 0;
                        } /*If next base is not an insertion*/

                        /**********************************************\
                        * Fun-4 Sec-5 Sub-9: Check if should keep ins
                        \**********************************************/
    
                        tmpBase = incBase; /*Base position on*/
                        qScoreUChar = *qCStr - Q_ADJUST; /*get Q-score*/
    
                        if(qScoreUChar < settings->minInsQUC)
                        { /*If the base is low quality, skip*/
                            /*Move to the next base in the list*/
    
                            /*This avoids scattered insertions*/
                            if(
                                cigEntryUInt > 0 &&  /*More insertions*/
                                !(incBase->nextBase->errTypeChar & 2)
                                    /*^^Next base is not an insertion*/
                            ) { /*If have other insertions to process*/
                                /*Make a new base for next insertion*/
                                tmpBase =
                                    malloc(sizeof(struct baseStruct));
    
                                tmpBase->nextBase = incBase->nextBase;
                                incBase->nextBase = tmpBase;
                                lastBase = incBase;
                                incBase = tmpBase;
    
                                /*Set the defaulst for the new base*/
                                incBase->altBase = 0;
                                incBase->errTypeChar = 2;
                                incBase->baseChar = 0;
                                incBase->numSupReadsUL = 0;
                            } /*If have other insertions to process*/
    
                            else
                            { /*Else is safe to move to the next base*/
                                lastBase = incBase;
                                incBase = incBase->nextBase;
                            } /*Else is safe to move to the next base*/
    
                            /*Move to the next base in the sam entry*/
                            ++seqCStr;
                            ++qCStr;
                            continue;
                        } /*If the base is low quality, skip*/
    
                        /**********************************************\
                        * Fun-4 Sec-5 Sub-10: Find matching ins at pos
                        \**********************************************/
    
                        while(tmpBase->baseChar != *seqCStr)
                        { /*While the bases are not equal*/
                            if(tmpBase->baseChar == 0)
                                break; /*Empty base*/
    
                            lastBase = tmpBase;
                            tmpBase = tmpBase->altBase;
    
                            if(tmpBase == 0)
                                 break; /*Need to create a new base*/
                        } /*While the bases are not equal*/
    
                        /**********************************************\
                        * Fun-4 Sec-5 Sub-11: check if need to make ins
                        \**********************************************/
    
                        if(tmpBase == 0)
                        { /*If need to create a new base*/
                            lastBase->altBase =
                                    malloc(sizeof(struct baseStruct));
    
                            tmpBase = lastBase->altBase;
                            tmpBase->nextBase = 0; /*Is an alterantive*/
                            tmpBase->altBase = 0;  /*No alternatives*/
    
                            tmpBase->errTypeChar = 2; /*Insertion*/
                            tmpBase->baseChar = *seqCStr;
                            tmpBase->numSupReadsUL = 0;
                        } /*If need to create a new base*/
    
                        else if(tmpBase->baseChar == 0)
                            tmpBase->baseChar = *seqCStr;
    
                        ++tmpBase->numSupReadsUL;
    
                        /*This avoids scattered deletions*/
                        if(
                            cigEntryUInt > 0 &&  /*More insertions*/
                            !(incBase->nextBase->errTypeChar & 2)
                                    /*^^Next base is not an insertion*/
                        ) { /*If have other insertions to process*/
                            /*Make a new base for next insertion*/
                            tmpBase =
                                malloc(sizeof(struct baseStruct));
    
                            tmpBase->nextBase = incBase->nextBase;
                            incBase->nextBase = tmpBase;
                            lastBase = incBase;
                            incBase = tmpBase;
    
                            /*Set the defaulst for the new base*/
                            incBase->altBase = 0;
                            incBase->errTypeChar = 2;
                            incBase->baseChar = 0;
                            incBase->numSupReadsUL = 0;
                        } /*If have other insertions to process*/
    
                        else /*Can move to the next base*/
                        { /*Else is safe to move to the next base*/
                            lastBase = incBase;
                            incBase = incBase->nextBase;
                        } /*Else is safe to move to the next base*/
    
                        ++seqCStr;
                        ++qCStr;
                    } /*While have bases to add to the bases list*/
                
                    break;
                /*Switch: for insertions*/

                /******************************************************\
                * Fun-4 Sec-5 Sub-12: Ignore deletions
                \******************************************************/

                case 'D':
                    /*Deletions will be marked as empty and will be 
                      caught if their are to few bases, so can ignore.*/

                    while(cigEntryUInt > 0)
                    { /*While have bases to ignore*/
                        --cigEntryUInt;

                        /*Move past other sequences insertions*/
                        while(incBase->errTypeChar & 2)
                            incBase = incBase->nextBase; 

                        lastBase = incBase;
                        incBase = incBase->nextBase;
                    } /*While have bases to ignore*/

                    break;

                /******************************************************\
                * Fun-4 Sec-5 Sub-13: Ignore soft masking
                \******************************************************/

                case 'S':
                /*Switch: ingnore soft maskes 'S'*/
                    while(cigEntryUInt > 0)
                    { /*While have bases to ignore*/
                        --cigEntryUInt;
                        ++seqCStr;
                        ++qCStr;
                    } /*While have bases to ignore*/

                    break;
                /*Switch: ingnore soft maskes 'S'*/
            } /*switch: check the error type & add bases to consensus*/
        } /*While not at the end of the sam alignment sequence*/

        blankSamEntry(samStruct);
        errUChar = readSamLine(samStruct, stdinFILE);
    } /*While have alignments to read in from the sam file*/

    pclose(stdinFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-6: Merge bases into a single majority consensus
    ^    fun-4 sec-6 sub-1: Set up for deciding bases to keep
    ^    fun-4 sec-6 sub-2: Find best alternate base for each position
    ^    fun-4 sec-6 sub-3: Check if should keep base
    ^    fun-4 sec-6 sub-4: Add the base & Q-score to the sequence
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-4 Sec-6 Sub-1: Set up for deciding bases to keep
    \******************************************************************/

    /*Minimum number of bases needed keep an SNP, match, or insterion*/
    minNumBasesUL = numSeqUL * settings->minReadsPercBaseFlt;
    minInsUL = numSeqUL * settings->minReadsPercInsFlt;
    settings->lenConUL = 0;

    seqCStr = samStruct->samEntryCStr;
    samStruct->readLenUInt = 0;
    lastBase = 0;                     /*For freeing the list*/

    /******************************************************************\
    * Fun-4 Sec-6 Sub-2: Find the best alternate base for each position
    \******************************************************************/

    while(headBase != 0)
    { /*For all positions in the consensus, remove non-majority bases*/
        incBase = headBase;         /*Move to 1st alternative base*/
        tmpBase = headBase->altBase;
        numSupReadsUL = incBase->numSupReadsUL;
       
        while(tmpBase != 0)
        { /*While have bases to remove*/
            numSupReadsUL += tmpBase->numSupReadsUL;

            if(incBase->numSupReadsUL < tmpBase->numSupReadsUL) 
                headBase = freeBaseStruct(&incBase, lastBase);
                /*Freeing puts tmpBase were incBase is*/

            else
            { /*If removing the current base in the list*/
                /*Need to make sure freeBaseStruct does nothing funny*/
                incBase->altBase = tmpBase->altBase;
                tmpBase->altBase = 0;
                freeBaseStruct(&tmpBase, lastBase);
            } /*If removing the current base in the list*/

            tmpBase = incBase->altBase; /*Move to next alternate base*/
        } /*While have bases to remove*/

        /**************************************************************\
        * Fun-4 Sec-6 Sub-3: Check if should keep base
        \**************************************************************/

        if(headBase->errTypeChar & 2)
        { /*If is an insertion*/
            /*Find the % of reads/supplemental alignments needed to
              support this position*/

            if(numSupReadsUL < minInsUL)
            { /*If have to few insertions to have support*/
                headBase = freeBaseStruct(&headBase, lastBase);
                continue;                /*Move on to the next base*/       
            } /*If have to few insertions to have support*/
        } /*If is an insertion*/

        else
        { /*If is a match or SNP*/
            /*Find the % of reads/supplemental alignments needed to
              support this position*/

            if(numSupReadsUL<minNumBasesUL)
            { /*If have to few SNPs or matches to have support*/
                headBase = freeBaseStruct(&headBase, lastBase);
                continue;                /*Move on to the next base*/       
            } /*If have to few SNPs or matches to have support*/
        } /*If is a match or SNP*/

        /**************************************************************\
        * Fun-4 Sec-6 Sub-4: Add the base & Q-score to the sequence
        \**************************************************************/

        /*Put the base into the sequence*/
        *seqCStr = headBase->baseChar; /*Set the base*/
        ++seqCStr; /*Move to next sequence entry*/
        ++settings->lenConUL; /*count base in the consensus length*/

        /*Free base and move to the next base*/
        headBase = freeBaseStruct(&headBase, lastBase);
    } /*For all positions in the consensus, remove non-majority bases*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-7: Print out cosensus & do clean up
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(numMisSeqUL > minNumBasesUL)
        return 16;                        /*If had to few mapped reads*/

    *seqCStr = '\0';

    /*Reset the sequence & q-score line pointers (for clairity)*/
    seqCStr = samStruct->samEntryCStr;

    /*Write consensus as fasta file*/
    stdinFILE = fopen(binStruct->consensusCStr, "w");
    fprintf(stdinFILE, ">%s\n%s\n", binStruct->consensusCStr, seqCStr);
    fclose(stdinFILE);

    return 1;
} /*simpleMajCon*/

/*---------------------------------------------------------------------\
| Output:
|    Uses: Racon to build a consensus (file name in bin->consensusCStr)
\---------------------------------------------------------------------*/
void buildConWithRacon(
    struct raconStruct *settings,   /*Settings to use with racon*/
    unsigned char *clustUC,         /*Cluster on*/
    char *threadsCStr,              /*Number threads to use with racon*/
    struct readBin *conBin          /*Bin working on*/
) /*Builds a consensus using racon*/
{ /*buildConWithRacon*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ' Fun-5 TOC: buildConWithRacon
    '    fun-5 sec-1: Variable declerations
    '    fun-5 sec-2: Set up the consensus name
    '    fun-5 sec-3: Build the consensus using racon
    '    fun-5 sec-4: Rename the consensus if needed and find the length
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-5 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned short lenBuffUS = 2048;
    unsigned short lenTmpBuffUS = 2048;
    char
        *refFileCStr = conBin->bestReadCStr,
        minimap2CmdCStr[lenBuffUS], /*Holds command to run minimap2*/
        raconCmdCStr[lenBuffUS],    /*Holds command to run racon*/
        tmpConCStr[lenTmpBuffUS],   /*Hold the consensus name*/
        *tmpCStr = 0,
        *tmpFileCStr = 0,    /*For swapping consensus file names*/
        *tmpFastaFileCStr = "tmp-2023-12-01-con-1563577017123412.fasta",
        *tmpSamFileCStr = "tmp-2023-12-01-map-15635770171234122342.sam";

     unsigned long numBytesUL = 0; /*Number bytes read in by fread*/
     FILE *conFILE = 0; /*For counting the consensus length*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-5 Sec-2: Set up the consensus name
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if need to build the consensus name*/
    if(conBin->consensusCStr[0] == '\0')
    { /*If using the best read for the first round of racon*/
        tmpCStr =
            cStrCpInvsDelm(conBin->consensusCStr, conBin->fqPathCStr);

        /*Get to end of .fastq or .fq*/
        if(*(tmpCStr - 2) == 't') tmpCStr -= 6; /*.fastq*/
        else tmpCStr -= 3;                      /*.fq*/

        tmpCStr = cStrCpInvsDelm(tmpCStr, "--clust-");

        /*Copy cluster number into file name; tmpCStr will point to end*/
        tmpCStr = uCharToCStr(tmpCStr, *clustUC);
        tmpCStr=cStrCpInvsDelm(tmpCStr, "--con.fasta");

        tmpFileCStr = conBin->consensusCStr;
        refFileCStr = conBin->bestReadCStr;
    } /*If using the best read for the first round of racon*/

    else
    { /*I am using a consensus*/
        tmpFileCStr = tmpFastaFileCStr;
        cStrCpInvsDelm(tmpConCStr, conBin->consensusCStr);
        refFileCStr = conBin->consensusCStr;
    } /*I am using a consensus*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-5 Sec-3: Build the consensus using racon
    ^    fun-5 sec-3 sub-1: Build sam file with minimap2
    ^    fun-5 sec-3 sub-1: Build consensus with racon
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-5 Sec-3 Sub-1: Build sam file with minimap2
    \******************************************************************/

    for(uint8_t UCRnd = 0; UCRnd < settings->rndsRaconUC; ++UCRnd)
    { /*Loop until have met user requirment for racon*/

        /*Build and run the minimap2 command*/
        tmpCStr = cStrCpInvsDelm(minimap2CmdCStr, minimap2CMD);
        tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
        tmpCStr =
            cpParmAndArg(tmpCStr, refFileCStr, conBin->topReadsCStr);
        cpParmAndArg(tmpCStr, ">", tmpSamFileCStr);

        system(minimap2CmdCStr);

        /**************************************************************\
        * Fun-5 Sec-3 Sub-1: Build consensus with racon
        \**************************************************************/

        tmpCStr = cStrCpInvsDelm(raconCmdCStr, raconCMD);
        tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
        tmpCStr = 
            cpParmAndArg(tmpCStr, conBin->topReadsCStr, tmpSamFileCStr);

        tmpCStr = cpParmAndArg(tmpCStr, refFileCStr, "> ");
        tmpCStr = cStrCpInvsDelm(tmpCStr, tmpFileCStr);

        system(raconCmdCStr); /*Run racon*/

        if(tmpFileCStr == conBin->consensusCStr)
        { /*If named the last consensus after the final name*/
            tmpFileCStr = tmpFastaFileCStr;
            refFileCStr = conBin->consensusCStr;
        } /*If named the last consensus after the final name*/

        else
        { /*Else named the last consensus after the temporary name*/
            tmpFileCStr = conBin->consensusCStr;
            refFileCStr = tmpFastaFileCStr;
        } /*Else named the last consensus after the temporary name*/
    } /*Loop until have met user requirment for racon*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-5 Sec-4: Rename consensus if needed and find the length
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    remove(tmpSamFileCStr); /*No longer need this*/

    /*If next round would have output to the consensus file*/
    if(tmpFileCStr != conBin->consensusCStr)
    { /*If need to rename the final consensus*/
        remove(conBin->consensusCStr);
        rename(tmpFileCStr, conBin->consensusCStr);
    } /*If need to rename the final consensus*/

    else                    /*Last consensus saved to consensus file*/
        remove(tmpFastaFileCStr); /*No longer need*/

    settings->lenConUL = 0; /*Reset for counting*/
    conFILE = fopen(conBin->consensusCStr, "r");

    if(conFILE == 0)
        return; /*No consensus made*/

    minimap2CmdCStr[lenBuffUS - 1] = '\0';
    minimap2CmdCStr[lenBuffUS - 2] = '\0';

    while(fgets(minimap2CmdCStr, lenBuffUS, conFILE))
    { /*While have a header to read in*/
        if(minimap2CmdCStr[lenBuffUS - 2] == '\0' ||
           minimap2CmdCStr[lenBuffUS - 2] == '\n'
        ) break; /*If at the end of the line*/
    } /*While have a header to read in*/

    numBytesUL = fread(minimap2CmdCStr,sizeof(char),lenBuffUS, conFILE);
    while(numBytesUL != 0)
    { /*While have a sequence to read in*/
        /*Account for reading in a full buffer of data*/
        settings->lenConUL += numBytesUL;

        numBytesUL = /*Read in the next line*/
            fread(minimap2CmdCStr, sizeof(char), lenBuffUS, conFILE);
    } /*While have a sequence to read in*/

    /*Consensuses output by racon have only two lines, one for the
      header and another for the sequence. So I can get away with a 
      very simple read function*/

    fclose(conFILE);
    return;
} /*buildConWithRacon*/

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
) /*Polish a consensus with medaka using the best reads*/
{ /*medakaPolish*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-6 TOC: Sec-1 Sub-1: medakaPolish
    '    fun-6 sec-1: Variable declerations
    '    fun-6 sec-2: See if the best reads and consensus files exist
    '    fun-6 sec-3: Build the command to run medaka
    '    fun-6 sec-4: Run medaka & clean up extra files
    '    fun-6 sec-5: rename consensus & delete directory medaka made
    '    fun-6 sec-6: Find the consensus length
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-1: Variable declerations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    unsigned short lenBuffUS = 2048;
    char medakaConPathCStr[512];  /*Path to consesnsus medaka made*/
    char medDirCStr[256];         /*For temporary c-string building*/
    char medakaCmdCStr[lenBuffUS];/*Holds command to run medaka*/
    char *tmpCStr = medakaCmdCStr;/*Temporary, for manipulating cStrs*/
    unsigned char errUC = 0;      /*Holds errors*/

    unsigned long fileLenULng = 0; /*See if files have something*/
    unsigned long numBytesUL = 0;  /*Hold in bytes read by fread*/

    FILE *testFILE = 0; /*Test if files exist*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-2: See if the consensus or best read file exist
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Open file*/
    testFILE = fopen(conBin->consensusCStr, "r");

    if(testFILE == 0)
    { /*Check if using the consensus*/
        testFILE = fopen(conBin->bestReadCStr, "r");

        if(testFILE == 0)
            return 2;     /*Nothing to build a consensus with*/

        fclose(testFILE); /*Was only for checking*/

        tmpCStr =
            cStrCpInvsDelm(conBin->consensusCStr, conBin->fqPathCStr);

        /*Get to end of .fastq or .fq*/
        if(*(tmpCStr - 2) == 't') tmpCStr -= 6; /*.fastq*/
        else tmpCStr -= 3;                      /*.fq*/

        tmpCStr = cStrCpInvsDelm(tmpCStr, "--clust-");

        /*Copy cluster number into file name*/
        tmpCStr = uCharToCStr(tmpCStr, *clustUC);
        tmpCStr=cStrCpInvsDelm(tmpCStr, "--con.fasta");

        /*Determine if best read is a fasta file or fastq*/
        tmpCStr = conBin->bestReadCStr;

        while(*tmpCStr != '\0')
            ++tmpCStr;

        errUC = 1; /*Assume no error*/

        if(*(tmpCStr - 1) == 'q')
        { /*If is a fastq file, so is a read*/
            errUC =
                fqToFa(
                    conBin->bestReadCStr,
                    conBin->consensusCStr,
                    samStruct
            ); /*convert fastq best read to fasta read*/
        } /*If is a fastq file, so is a read*/

        else   /*Is a fasta file (likely a consensus*/
            rename(conBin->bestReadCStr, conBin->consensusCStr);

        if(errUC & 64)
            return 64; /*Memory allocation error*/
        if(!(errUC & 1))
            return 2;  /*File read/write error*/

        testFILE = fopen(conBin->consensusCStr, "r");
    } /*Check if using the consensus*/

    /*Finde the length of the file*/
    fseek(testFILE, 0L, SEEK_END);   /*Find end of file*/
    fileLenULng = ftell(testFILE);   /*Find offset (length) of end*/
    fclose(testFILE);                /*No longer need the file open*/

    if(fileLenULng < 10)
        return 2;    /*Nothing or little in the file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-3: Build the command to run medaka
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Make temporary director for medaka*/
    tmpCStr = cStrCpInvsDelm(medDirCStr, conBin->fqPathCStr);

    /*Get to end of .fastq or .fq*/
    if(*(tmpCStr - 2) == 't') tmpCStr -= 6; /*.fastq*/
    else tmpCStr -= 3;                      /*.fq*/

    tmpCStr = cStrCpInvsDelm(tmpCStr, "--medaka");

    /*Copy the path to the consensus medaka will build*/
    tmpCStr = cStrCpInvsDelm(medakaConPathCStr, medDirCStr);
    tmpCStr = cStrCpInvsDelm(tmpCStr, "/consensus.fasta");

    /*Copy the enviroment activate command*/
    if(settings->condaBl & 1)
        tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medCondCMD);
    else
        tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medakaCMD);

    *tmpCStr = ' ';
    ++tmpCStr;

    tmpCStr = cStrCpInvsDelm(tmpCStr, "medaka_consensus");
    tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
    tmpCStr = cpParmAndArg(tmpCStr, "-i", conBin->topReadsCStr);
    tmpCStr = cpParmAndArg(tmpCStr, "-d", conBin->consensusCStr);
    tmpCStr = cpParmAndArg(tmpCStr, "-m", settings->modelCStr);
    tmpCStr = cpParmAndArg(tmpCStr, "-o", medDirCStr);

    /*Copy the deactivate command in*/
    if(settings->condaBl & 1)
        tmpCStr = cStrCpInvsDelm(tmpCStr, medCondCMDEnd);
    else
        tmpCStr = cStrCpInvsDelm(tmpCStr, medakaCMDEnd);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-4: Run medaka & clean up extra files
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    system(medakaCmdCStr); /*Run medaka command*/

    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medDirCStr);
    cStrCpInvsDelm(tmpCStr, "/calls_to_draft.bam");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medDirCStr);
    cStrCpInvsDelm(tmpCStr, "/calls_to_draft.bam.bai");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medDirCStr);
    cStrCpInvsDelm(tmpCStr,"/consensus.fasta.gaps_in_draft_coords.bed");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medDirCStr);
    cStrCpInvsDelm(tmpCStr,"/consensus_probs.hdf");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    /*Remove the temporary mapping files made by medaka*/
    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, conBin->consensusCStr);
    tmpCStr = cStrCpInvsDelm(tmpCStr, ".fai");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    /*Remove the temporary mapping files made by medaka*/
    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, conBin->consensusCStr);
    tmpCStr = cStrCpInvsDelm(tmpCStr, ".map-ont.mmi");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-5: rename consensus & delete directory medaka made
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    testFILE = fopen(medakaConPathCStr, "r");

    if(testFILE == 0)
        return 4;    /*Medaka did not build a consensus*/

    /*See if the test file has anything in it*/
    fseek(testFILE, 0L, SEEK_END);   /*Find end of file*/
    fileLenULng = ftell(testFILE);   /*Find offset (length) of end*/
    fclose(testFILE);                /*No longer need the file open*/

    if(fileLenULng < 10)
        return 4;    /*Nothing or little in the file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-6: Find the consensus length
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    remove(conBin->consensusCStr); /*remove the old consensus*/
    rename(medakaConPathCStr, conBin->consensusCStr);
    remove(medDirCStr); /*Delete directory made by medaka*/

    settings->lenConUL = 0; /*Reset for counting*/
    testFILE = fopen(conBin->consensusCStr, "r");

    if(testFILE == 0)
        return 4; /*No consensus made*/

    /*Add in markers to mark the end of headers*/
    medakaCmdCStr[lenBuffUS - 1] = '\0';
    medakaCmdCStr[lenBuffUS - 2] = '\0';

    while(fgets(medakaCmdCStr, lenBuffUS, testFILE))
    { /*While have a header to read in*/
        if(medakaCmdCStr[lenBuffUS - 2] == '\0' ||
           medakaCmdCStr[lenBuffUS - 2] == '\n'
        ) break; /*If at the end of the line*/
    } /*While have a header to read in*/

    numBytesUL = fread(medakaCmdCStr,sizeof(char),lenBuffUS,testFILE);

    while(numBytesUL != 0)
    { /*While have a sequence to read in*/
        /*Account for reading in a full buffer of data*/
        settings->lenConUL += numBytesUL;

        numBytesUL = /*Read in the next line*/
            fread(medakaCmdCStr, sizeof(char), lenBuffUS, testFILE);
    } /*While have a sequence to read in*/

    /*Consensuses output by medaka have only two lines, one for the
      header and another for the sequence. So I can get away with a 
      very simple read function*/

    return 1;
} /*medakaPolish*/

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
) /*Compares a consenses to a another consensus. This will do a
    recursive call if conBinTree has children*/
{ /*cmpCons*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-7 TOC: cmpCons
    '     fun-7 sec-1: Variable declerations
    '     fun-7 sec-2: Check if have valid user input
    '     fun-7 sec-3: Run minimap2 and check first line
    '     fun-7 sec-5: Score the consensus to the reference
    '     fun-7 sec-6: Compare this score to the other consensuses
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-7 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t errUChar = 0; /*For holding error returns from functions*/
    uint8_t zeroUChar = 0;

    char *tmpCStr = 0;
    char minimap2CmdCStr[2048];

    uint32_t
        incBuffUInt = 10000; /*Amount to increase buff size each time*/

    FILE
        *stdinFILE = 0; /*File to see if input files are valid*/

    struct readBin
        *refBin = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-7 Sec-2: Check if have valid user input
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(conBin == 0)
        return 0;

    if(conBin->leftChild != 0)
        return 0; /*This is a bin, not a cluster*/

    if(conBin->balUChar < 0)
        return 0;

    if(conBinTree == 0)
        return 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-4: Read in target consensus sequence
    ^    fun-7 sec-4 sub-1: Read past the header
    ^    fun-7 sec-4 sub-2: Read in the sequence
    ^    fun-7 sec-4 sub-3: Count the sequence length
    ^ Note:
    ^    - I should change this to the consensus I am checking & then
    ^      have a marker that tells this function to use the input
    ^      reference
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-7 Sec-4 Sub-1: Read past the header
    \******************************************************************/

    blankSamEntry(refStruct);

    /*Set up null endings for lines*/
    *(refStruct->samEntryCStr + refStruct->lenBuffULng - 1) = '\0';
    *(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) = '\0';

    /*Open the reference file for reading*/
    stdinFILE = fopen(conBin->consensusCStr, "r");

    /*Read in the header, I know it worked, because of minimap2*/
    fgets(refStruct->samEntryCStr, refStruct->lenBuffULng, stdinFILE);

    while(
     !(*(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) !='\0' ||
       *(refStruct->samEntryCStr + refStruct->lenBuffULng -2)!='\n')
    ) { /*While have a header to read in*/
        /*Read in the next part of the header*/
        fgets(refStruct->samEntryCStr, refStruct->lenBuffULng, stdinFILE);

        /*Resetup markers*/
        *(refStruct->samEntryCStr + refStruct->lenBuffULng - 1) = '\0';
        *(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) = '\0';
    } /*While have a header to read in*/

    /******************************************************************\
    * Fun-7 Sec-4 Sub-2: Read in the sequnence
    \******************************************************************/

    tmpCStr = refStruct->samEntryCStr;
    *(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) = '\0';

    /*Read in the first part of the sequence*/
    fgets(refStruct->samEntryCStr, refStruct->lenBuffULng, stdinFILE);

    while(*(refStruct->samEntryCStr+refStruct->lenBuffULng-2) > 16)
    { /*While on the sequence line*/
        refStruct->samEntryCStr =
            realloc(
                refStruct->samEntryCStr,
                refStruct->lenBuffULng + incBuffUInt
            );

        if(refStruct->samEntryCStr == 0)
        { /*memory allocation error*/
            fclose(stdinFILE);
            return 0;
        } /*Memory allocation error*/

        /*Set pointer to new buffer*/
        tmpCStr = refStruct->samEntryCStr + refStruct->lenBuffULng;
        refStruct->lenBuffULng += incBuffUInt; /*Update buff size*/

        /*Reset new line markers*/
        *(refStruct->samEntryCStr + refStruct->lenBuffULng-2) ='\0';

        /*Read in next part of reference sequence*/
        fgets(tmpCStr, refStruct->lenBuffULng, stdinFILE);
    } /*While on the sequence line*/

    fclose(stdinFILE);

    /******************************************************************\
    * Fun-7 Sec-4 Sub-3: Count the sequence length
    \******************************************************************/

    refStruct->seqCStr = refStruct->samEntryCStr;

    /*scoreAln assumes the Q-score entry is not null, so I am setting
      it to something to avoid errors. This is ok since I am telling
      scoreAln that their is not Q-score for the reference
    */
    refStruct->qCStr = refStruct->seqCStr;
    tmpCStr = refStruct->seqCStr;
    refStruct->readLenUInt = 0;

    while(*tmpCStr != '\0')
    { /*While have bases to count*/
        ++refStruct->readLenUInt;
        ++tmpCStr;
    } /*While have bases to count*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-7 Sec-5: Score the consensus to the reference
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(conBinTree != 0)
    { /*While have consensuses to compare*/
        blankSamEntry(samStruct);
        refBin = conBinTree->rightChild;

        while(refBin != 0)
        { /*While have another clusters consensus to compare*/
            if(refBin->balUChar < 0)
            { /*If this cluster has been marked to be skipped*/
                refBin = refBin->rightChild;
                continue;
            } /*If this cluster has been marked to be skipped*/

            blankSamEntry(samStruct);

            /*Prepare the minimap2 command*/
            tmpCStr = cStrCpInvsDelm(minimap2CmdCStr, minimap2CMD);
            tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
            cpParmAndArg(
                tmpCStr,
                conBin->consensusCStr,
                refBin->consensusCStr
            );

            stdinFILE = popen(minimap2CmdCStr, "r"); /*run minimap2*/

            blankSamEntry(samStruct);
            errUChar = readSamLine(samStruct, stdinFILE); /*1st line*/

            if(*samStruct->samEntryCStr != '@')
            { /*If no header*/
                pclose(stdinFILE);
                return 0;
            } /*If no header*/

            while(errUChar & 1)
            { /*While on the haeder lines*/
                blankSamEntry(samStruct);
                errUChar = readSamLine(samStruct, stdinFILE);

                if(*samStruct->samEntryCStr != '@')
                { /*If off the header entries*/
                    pclose(stdinFILE);
                    break; /*If not a header*/
                } /*If off the header entries*/
            } /*While on the haeder lines*/

            if(*samStruct->samEntryCStr == '@')
            { /*If their was only headers*/
                pclose(stdinFILE);
                continue;
            } /*If their was only headers*/

            if(!(samStruct->flagUSht & 4))
            { /*If the reads mapped to each other*/
                scoreAln(
                    minStats,
                    samStruct,
                    refStruct,   /*Reference struct to score dels with*/ 
                    &zeroUChar,  /*Mapped consensus has no Q-score*/
                    &zeroUChar   /*Mapped consensus has no Q-score*/
                ); /*Score the alignment*/

                if(checkIfKeepRead(minStats, samStruct) & 1)
                    return refBin; /*If consensus look the same*/
                    /*The file is already closed*/
            } /*If the reads mapped to each other*/

            refBin = refBin->rightChild;
        } /*While have another clusters consensus to compare*/

        conBinTree = conBinTree->leftChild;
    } /*While have consensuses to compare*/

    return 0;
} /*cmpCons*/

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
) /*Frees an base from a linked list of bases*/
{ /*freeBaseStruct*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-8 TOC: Sec-1 Sub-1: freeBaseStruct
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(*baseToFree == 0)
        return 0;
    
    if((*baseToFree)->altBase != 0)
    { /*If have an alternative base to keep as the next base*/
        (*baseToFree)->altBase->nextBase = (*baseToFree)->nextBase;

        if(lastBase != 0)
        { /*If have a previous base*/
            lastBase->nextBase = (*baseToFree)->altBase;
            free(*baseToFree);
            *baseToFree = lastBase->nextBase;
        } /*If have a previous base*/

        else
        { /*Else this was the first base*/
            lastBase = (*baseToFree)->altBase;
            free(*baseToFree);
            *baseToFree = lastBase;
        } /*Else this was the first base*/

        return *baseToFree;
    } /*If have an alternative base to keep as the next base*/

    else
    { /*Else if their is only a next base*/
        if(lastBase != 0)
        { /*If have a previous base in the list*/
            lastBase->nextBase = (*baseToFree)->nextBase;
            free(*baseToFree);
            *baseToFree = 0; /*Tell that their are no alternate bases*/
            return lastBase->nextBase;
        } /*If have a previous base in the list*/

        else
        { /*Else baseToFree was the first base*/
            lastBase = (*baseToFree)->nextBase;
            free(*baseToFree);
            *baseToFree = 0; /*Tell that their are no alternate bases*/
            return lastBase;
        } /*Else baseToFree was the first base*/
    } /*Else if their is only a next base*/
} /*freeBaseStruct*/

/*---------------------------------------------------------------------\
| Output: Modifies: majConStruct to have default settings
\---------------------------------------------------------------------*/
void initMajConStruct(
    struct majConStruct *majConSettings
    /*struct to set to default values in defaultSettings.h*/
) /*Sets input structers variables to default settings*/
{ /*initMajConStruct*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-9 TOC: Sec-1 Sub-1: initMajConStruct
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    majConSettings->useMajConBl =
       (defConMethod == defUseMajCon);
    majConSettings->minBaseQUC = majConMinBaseQ;
    majConSettings->minInsQUC = majConMinInsQ;
    majConSettings->minReadsPercBaseFlt = percBasesPerPos;
    majConSettings->minReadsPercInsFlt = percInsPerPos;

    majConSettings->lenConUL = 0; /*counter*/

    return;
} /*initMajConStruct*/

/*---------------------------------------------------------------------\
| Output: Modifies: raconStruct to have default settings
\---------------------------------------------------------------------*/
void initRaconStruct(
    struct raconStruct *raconSettings
    /*struct to set to default values in defaultSettings.h*/
) /*Sets input structers variables to default settings*/
{ /*initRaconStruct*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-10 TOC: Sec-1 Sub-1: initRaconStruct
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    raconSettings->useRaconBl =
       (defConMethod == defUseRacon);
    raconSettings->rndsRaconUC = defRoundsRacon;
    raconSettings->lenConUL = 0;

    return;
} /*initRaconStruct*/

/*---------------------------------------------------------------------\
| Output: Modifies: medakaStruct to have default settings
\---------------------------------------------------------------------*/
void initMedakaStruct(
    struct medakaStruct *medakaSettings
    /*struct to set to default values in defaultSettings.h*/
) /*Sets input structers variables to default settings*/
{ /*initMedakaStruct*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-11 TOC: Sec-1 Sub-1: initMedakaStruct
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    medakaSettings->useMedakaBl =
       (defConMethod == defUseMedaka);
    strcpy(medakaSettings->modelCStr, defMedakaModel);
    medakaSettings->condaBl = defCondaBl;
    medakaSettings->lenConUL = 0;

    return;
} /*initMedakaStruct*/

/*---------------------------------------------------------------------\
| Output: Modifies: medakaStruct to have default settings
\---------------------------------------------------------------------*/
void initConBuildStruct(
    struct conBuildStruct *conSetST
    /*struct to set to default values in defaultSettings.h*/
) /*Sets input structers variables to default settings*/
{ /*initConBuildStruct*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-12 TOC: Sec-1 Sub-1: initConBuildStruct
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    conSetST->methodAryUC[0] = defConMethod;
    conSetST->methodAryUC[1] = defNoCon;
    conSetST->lenMethodUC = 1;
    conSetST->minPercMappedReadsFlt =defMinPercMappedReads;
    
    conSetST->useStatBl = 0;
    conSetST->clustUC = 0;
    conSetST->minReadsToBuildConUL = minReadsPerBin;
    conSetST->numRndsToPolishUI = defNumPolish;
    conSetST->maxReadsToBuildConUL = defReadsPerCon;
    conSetST->minConLenUI = defMinConLen;
    conSetST->lenConUL = 0;
    conSetST->numReadsForConUL = 0;

    initMajConStruct(&conSetST->majConSet);
    initRaconStruct(&conSetST->raconSet);
    initMedakaStruct(&conSetST->medakaSet);
    initIvarSet(conSetST->ivarSetST);

    return;
} /*initConBuildStruct*/

/*--------------------------------------------------------\
| Name: ivarCon
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
    struct samEntry *samST   /*For getting sam file input*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: ivarCon
   '  - Builds a consensus using ivar from the best reads
   '    & top read
   '  o fun-13 sec-01:
   '    - Variable declerations
   '  o fun-13 sec-02:
   '    - Make the command for ivar
   '  o fun-13 sec-03:
   '    - Run ivar
   '  o fun-13 sec-04:
   '    - Clean up and trim/print consensus
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-13 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/


   ushort lenIvarCmdUS = 2048;
   char minimap2CmdStr[lenIvarCmdUS];
   char ivarCmdStr[lenIvarCmdUS];
   char *buffStr = 0;
      /*So can swap between buffers when checking the
      ` sequence. I am using the minimap2CmdStr and
      ` ivarCmdStr as buffers in the consensus check step.
      */

   int lenPrefI = 512;
   char prefixStr[lenPrefI];   /*output File name*/
   char *prefEndStr = 0;
   char *tmpStr = 0;

   uchar errUC = 0;  /*Report errors*/

   FILE *stdinFILE = 0;
   FILE *stdoutFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-13 Sec-02:
   ^  - Make the commands to call
   ^  o fun-13 sec-02 sub-01:
   ^    - Build the consensus name
   ^  o fun-13 sec-02 sub-02:
   ^    - Make the command for minimap2
   ^  o fun-13 sec-02 sub-03:
   ^    - Make the command for ivar
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-13 Sec-02 Sub-01:
   *  - Build the consensus name
   \******************************************************/

    tmpStr =
       cStrCpInvsDelm(
          conST->consensusCStr,
          conST->fqPathCStr
    ); /*Copy the consensus name*/

    /*Get to end of .fastq or .fq*/
    if(*(tmpStr - 2) == 't') tmpStr -= 6; /*.fastq case*/
    else tmpStr -= 3;                     /*.fq case*/

    tmpStr = cStrCpInvsDelm(tmpStr, "--clust-");

    /*Copy cluster number into file name*/
    tmpStr = uCharToCStr(tmpStr, *clustUC);
    tmpStr = cStrCpInvsDelm(tmpStr, "--con.fasta");

   /******************************************************\
   * Fun-13 Sec-02 Sub-02:
   *  - Make the command for minimap2
   \******************************************************/

   tmpStr = cStrCpInvsDelm(minimap2CmdStr, minimap2CMD);
   tmpStr = cpParmAndArg(tmpStr, "-t", threadsStr);

   /*Check if I am using a consensus file or best reads*/
   /*I set the best read to the last consensus when I am
   ` doing multiple iterations
   */
   if(conST->consensusCStr[0] == '\0')
      tmpStr =
         cpParmAndArg(
            tmpStr,
            conST->bestReadCStr,
            conST->consensusCStr
      ); /*Copy reads and the consensus*/

   else
      tmpStr =
         cpParmAndArg(
            tmpStr,
            conST->bestReadCStr,
            conST->fqPathCStr
      ); /*Copy the reads and the best read/consensus*/

   /******************************************************\
   * Fun-13 Sec-02 Sub-03:
   *  - Make the command for ivar
   \******************************************************/
  
    /*I will use popen(file, "w") to pass the trimmed
    ` minimap2 output to the samtools/ivar pipeline
    */
   tmpStr = cStrCpInvsDelm(ivarCmdStr, "samtools sort");
   tmpStr = cpParmAndArg(tmpStr, "-@ 1", "- |");
     /*For some odd reason the pipe has an issue with
     ` samtools sort with more than 1 thread*/
   tmpStr = cpParmAndArg(tmpStr, "samtools", "view");
   tmpStr = cpParmAndArg(tmpStr, "-@", "1");
     /*I could multi thread here, but no point*/
   tmpStr = cpParmAndArg(tmpStr, "-F", "4");
     /*No unmapped reads*/
   tmpStr = cpParmAndArg(tmpStr, "-F", "256");
     /*No secondary alignments*/
   tmpStr = cpParmAndArg(tmpStr, "-F", "512");
     /*No alignments that falied to pass filters*/
   tmpStr = cpParmAndArg(tmpStr, "-F", "2048");
     /*No supplemenatry alignments*/
   tmpStr = cpParmAndArg(tmpStr, "-b -", "|");
   tmpStr = cpParmAndArg(tmpStr, "samtools", "mpileup");
   tmpStr = cpParmAndArg(tmpStr,"-B -aa -A -d 0","-Q 0");
   tmpStr = cpParmAndArg(tmpStr, "-", "|");

   tmpStr = cpParmAndArg(tmpStr, "ivar", "consensus");

   prefEndStr =
      strCpTill(
         prefixStr,
         conST->consensusCStr,
         '.',
         &lenPrefI
   ); /*Get the name of the consensus*/
   
   tmpStr = cpParmAndArg(tmpStr, "-p", prefixStr);
   tmpStr = cpParmAndArg(tmpStr, "-i", prefixStr);
   tmpStr = cpParmAndArg(tmpStr,"-t",settings->minSupStr);

   tmpStr =
      cpParmAndArg(tmpStr, "-c", settings->minInsSupStr);
   tmpStr =
      cpParmAndArg(tmpStr, "-m", settings->minDepthStr);

   *tmpStr = ' ';
   ++tmpStr;
   *tmpStr = '-';
   ++tmpStr;
   *tmpStr = 'n';
   ++tmpStr;
   *tmpStr = ' ';
   ++tmpStr;
   *tmpStr = settings->maskC;
   ++tmpStr;
   *tmpStr = ' ';
   ++tmpStr;

   tmpStr = cpParmAndArg(tmpStr, "-q", settings->minQStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-13 Sec-03:
   ^  - Run ivar
   ^  o fun-13 sec-03 sub-01:
   ^    - Call minimap2
   ^  o fun-13 sec-03 sub-02:
   ^    - Call ivar
   ^  o fun-13 sec-03 sub-03:
   ^    - Build Consensus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-13 Sec-03 Sub-01:
   *  - Call minimap2
   \******************************************************/

   stdinFILE = popen(minimap2CmdStr, "r");
   if(stdinFILE == 0) return 32;

   blankSamEntry(samST);
   errUC = readSamLine(samST, stdinFILE);

   if(!(errUC & 1))
   { /*If: I had an error with minimap2*/
      pclose(stdinFILE);
      return 32;
   } /*If: I had an error with minimap2*/

   /******************************************************\
   * Fun-13 Sec-03 Sub-02:
   *  - Call ivar
   \******************************************************/

   stdoutFILE = popen(ivarCmdStr, "w");

   if(stdoutFILE == 0)
   { /*If: I had a pipe problem with ivar*/
      pclose(stdinFILE);
      blankSamEntry(samST);
      return 2;
   } /*If: I had a pipe problem with ivar*/

   /******************************************************\
   * Fun-13 Sec-03 Sub-03:
   *  - Build the consensus
   \******************************************************/

   while(errUC & 1)
   { /*Loop: though the entire minimap2 output*/
      if(*samST->samEntryCStr != '@')
      { /*If: This is not a header entry*/
         if(
              (samST->flagUSht & 4)    /*Unmmaped*/
           || (samST->flagUSht & 256)/*SecondaryAlignment*/
           || (samST->flagUSht & 512)  /*Failed filters*/
           || (samST->flagUSht & 2048) /*Supplemental aln*/
         ){ /*If: This is not an entry I want to keep*/
            blankSamEntry(samST);
            errUC = readSamLine(samST, stdinFILE);
            continue;
         } /*If: This is not an entry I want to keep*/

         trimSamEntry(samST);
      } /*If: This is not a header entry*/

      printSamEntry(samST, stdoutFILE);
      blankSamEntry(samST);
      errUC = readSamLine(samST, stdinFILE);
   } /*Loop: though the entire minimap2 output*/

   pclose(stdinFILE);
   stdinFILE = 0;

   pclose(stdoutFILE);
   stdoutFILE = 0;

   if(errUC & 64) return 64; /*If this was a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-13 Sec-04:
   ^  - Clean up and trim/print consensus
   ^  o fun-13 sec-04 sub-01:
   ^    - Remove the quality score file ivar outputs
   ^  o fun-13 sec-04 sub-02:
   ^    - Open consensus files & print header to out file
   ^  o fun-13 sec-04 sub-03:
   ^    - Get past anonymous bases at the start
   ^  o fun-13 sec-04 sub-04:
   ^    - Mask anonymous bases, find length, & move to end
   ^  o fun-13 sec-04 sub-05:
   ^    - Trim masked bases at end (adjust length for trim)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-13 Sec-04 Sub-01:
   *  - Remove the quality score file ivar outputs
   \******************************************************/

   /*Add the fa ending ivar output*/
   strCpTill(prefEndStr, ".qual.txt", 0, &lenPrefI);
   remove(prefixStr);

   /******************************************************\
   * Fun-13 Sec-04 Sub-02:
   *  - Open consensus files & print header to output file
   \******************************************************/

   /*Add the fa ending ivar output*/
   *prefEndStr = '.';
   *(prefEndStr + 1) = 'f';
   *(prefEndStr + 2) = 'a';
   *(prefEndStr + 3) = '\0';

   /*Reads in reference sequence to samST and removes
   ` white space. Header is ignored.
   */
   if(readInConFa(prefixStr, samST) != 1) return 4;

   stdoutFILE = fopen(conST->consensusCStr, "w");
   if(stdoutFILE == 0) return 4;

   /*Put the header to the file*/
   *prefEndStr = '\0';
   fprintf(stdoutFILE, ">%s\n", prefixStr);
   *prefEndStr = '.';

   /******************************************************\
   * Fun-13 Sec-04 Sub-03:
   *  - Get past anonymous bases at the start
   \******************************************************/

   buffStr = samST->seqCStr;

   while(*buffStr > 32) /*No white space in sequence*/
   { /*Loop: Get past masking at start*/
      switch(*buffStr & (~32))
      { /*Switch: Check if is a standard base*/
         case 'A': goto ivarConCheckSeq;
         case 'T': goto ivarConCheckSeq;
         case 'C': goto ivarConCheckSeq;
         case 'G': goto ivarConCheckSeq;
         default:
            ++buffStr;
            break;
      } /*Switch: Check if is a standard base*/
   } /*Loop: Get past masking at start*/

   /******************************************************\
   * Fun-13 Sec-04 Sub-04:
   *  - Mask anonymous bases, find length, and move to end
   \******************************************************/

   ivarConCheckSeq:

   tmpStr = buffStr;
   settings->lenConUL = 0;

   while(*tmpStr > 32) /*No white space in sequence*/
   { /*Loop: While I have a sequence to read in*/
      ++(settings->lenConUL);

      switch(*tmpStr & (~32))
      { /*Switch: Check if is a standard base*/
         case 'A': break;
         case 'T': break;
         case 'C': break;
         case 'G': break;
         default: *tmpStr = settings->maskC; break;
      } /*Switch: Check if is a standard base*/

      ++tmpStr;
   } /*Loop: While I have a sequence to read in*/

   /******************************************************\
   * Fun-13 Sec-04 Sub-05:
   *  - Trim masked bases at end (adjust length for trim)
   \******************************************************/

   --tmpStr;

   while(*tmpStr > *buffStr)
   { /*Loop: Trim off masking at end*/
      if(*tmpStr != settings->maskC) break;
      --tmpStr;
      --(settings->lenConUL);
   } /*Loop: Trim off masking at end*/

   *(tmpStr + 1) = '\0';

   fputs(buffStr, stdoutFILE);
   fputc('\n', stdoutFILE);

   fclose(stdoutFILE);
   stdoutFILE = 0;

   remove(prefixStr);

   return 1;
} /*ivarCon*/

/*--------------------------------------------------------\
| Name: initIvarSet [initIvarSet(iverSetST);]
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
