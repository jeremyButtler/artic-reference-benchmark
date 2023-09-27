/*######################################################################
# Use:
#   o buildCon is the driver file for the buildCon program.
#   o It will buid a consensus from a fastq file with or without a 
#     reference.
# Requires:
#   o cStrFun.c/h
#   o cStrToNumberFun.c/h
#   o defaultSettings.h
#   o printError.c/h
#   o samEntryStruct.c/h
#   o trimSam.c/h
#   o fqGetIdsStructs.c/h
#   o fqGetIdsAVLTree.c/h
#   o fqGetIdsHash.c/h
#   o fqGetIdsFqFun.c/h
#   o fqGetIdsSearchFq.c/h
#   buildConFun.c/h
# C libaries:
#   o string.h
#   o stdlib.h
#   o stdio.h
#   o stdint.h
######################################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' HOP: Head of program
'   header:
'     o has header files and includes
'   main:
'     o The main function
'   fun-1 getUserInput:
'     o processes usser input
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
^ Header:
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "buildConFun.h"

/*---------------------------------------------------------------------\
| Output:
|    - Modifies input variables to hold user input
|    - Returns:
|        - 0: if nothing went wrong
|        - pionter to invalid paramter
|            - if no agruments input, returns ponter to argsCStr
\---------------------------------------------------------------------*/
char * getUserInput(
    int32_t lenArgsInt,
    char *argsCStr[],
    char **fqPathCStr, /*Holds path to fastq file*/
    char **refPathCStr, /*Holds path to references*/
    char **statsPathCStr, /*Holds path to scoreReads output to use*/
    char *prefixCStr,   /*Prefix to name everything*/
    char *threadsCStr, /*Number threads for minimap2 & racon*/
    struct conBuildStruct *conSet, 
    struct minAlnStats *readToReadMinStats,/*Read pull scoring setting*/
    struct minAlnStats *minReadConStats  /*Read cluster socring set*/
); /*Reads in user input*/

int main(
    int32_t lenArgsInt, /*Number of parameters & arguments user input*/
    char *argsCStr[]  /*List of parameters & arguements user input*/
) /*Function to run everything*/
{ /*main*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Main TOC:
    '   main sec-1: Variable declerations
    '   main sec-2: Set up default settings in structures
    '   main sec-3: Get user input & set up commands
    '   main sec-4: Check user input
    '   main sec-5: Build the consensus
    '   main sec-6: Clean up, Print out version numbers, and exit
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-1: Variable declerations
    ^    main sec-1 sub-1: General variables
    ^    main sec-1 sub-2: Help message
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-1 Sub-1: General variables
    \******************************************************************/

    /*User Input or related variables (non file names)*/
    char *fqPathCStr = 0;
    char *refPathCStr = 0;
    char *statsPathCStr = 0;      /*Stats file to use in read seletion*/
    char prefixCStr[64] = defPrefix;
    char tmpPathCStr[128];               /*copy of Fastq file working on*/
    char threadsCStr[16] = defThreads; /*# of threads for system calls*/

    /*Variables for the majority consensus step*/
    struct conBuildStruct conSetting;

    /*C-strings that hold commands*/
    char tmpCmdCStr[1024];      /*Holds a quick system command*/
    char minimap2VersionCStr[256];
    char raconVersionCStr[256];
    char medakaVersionCStr[2048]; /*For gpu error messages*/
    char ivarVersionStr[256];
    char samtoolsVersionStr[256];
    char *tmpCStr = 0;         /*For string manipulation*/

    /*Miscalanious variables*/
    unsigned char errUC = 0;         /*Holds error messages*/

    /*Holds thresholds to keep alignments, errors, & matches*/
    struct samEntry samStruct; /*For reading files*/
    struct samEntry refStruct; /*For reading files*/

    struct minAlnStats readToReadMinStats; /*for read to consensus map*/
    struct minAlnStats minReadConStats; /*for read to consensus map*/

    struct readBin fastqStruct;         /*Holds my fastq file*/

    FILE *stdinFILE = 0;
    FILE *cpFILE = 0;  /*For making a copy of the input fastq file*/
        /*This avoids modifications of the input fastq file*/

    /******************************************************************\
    * Main Sec-1 Sub-2: Help message
    \******************************************************************/

    char *helpCStr = "\
            \n buildCon -fastq reads.fastq [Options ...]\
            \n Use: Builds a consensus from a fastq file.\
            \n    -fastq: [Required]\
            \n      - Fastq file with reads to search\
            \n    -ref: [None]\
            \n      - Reference start building the\
            \n        consensus with\
            \n      - This will hunt for reads if the\
            \n        reference failed.\
            \n    -prefix: [out]\
            \n      - Prefix to name output file\
            \n    -stats: [None]\
            \n      - tsv file output by score reads to\
            \n        use to select the best read by\
            \n        mapping quality.\
            \n      - Default is to use the read with the\
            \n        best median Q-score.\
            \n    -threads: [3]\
            \n      - Number of threads to use\
            \n    -min-depth: [100]\
            \n      - Min number of reads needed to build\
            \n        a consensus\
            \n    -min-perc: [0.1 or 10%]\
            \n      - Min percentage of reads needed to\
            \n        map to a read/con to keep it.\
            \n      - This is here to overide -min-depth\
            \n        when the read depth is deep.\
            \n    -max-depth: [300]\
            \n      - Max number of reads to subsample\
            \n        for buildnig a consensus.\
            \n    -extra-consensus-steps: [2]\
            \n      - Number of times to rebuild the\
            \n        consensus using a new set of best\
            \n        reads.\
            \n    -min-length: [500]\
            \n      - Discard consensuses or read mappings\
            \n        that are under the input length.\
            \n Additional Help messages:\
            \n    -h-build-consensus:\
            \n      - Help message for consensus building\
            \n    -h-read-con-map:\
            \n      - Help message for read to consensus\
            \n        mapping\
            \n    -h-read-read-map:\
            \n      - Help mesage for read to read mapping\
            \n Output:\
            \n  - File:\
            \n    o Consensus built from the fasta file\
            \n    o fastqFileName--cluster-0--consensus.fasta\
            \n  - Stdout: Program versions.\
            \n Requires:\
            \n  - Minimap2\
            \n Option specific requirments:\
            \n  - Racon\
            \n    o Required for -use-racon\
            \n    o Needs to be in the system path\
            \n  - Medaka\
            \n    o Required for -use-medaka\
            \n    o Can be installed at ~/medaka through\
            \n      the python virtual enviroment\
            \n      (git hub install) or by miniconda\
            \n  - Ivar\
            \n    o Required for -use-ivar\
            \n    o Requires: samtools in the system path\
            \n    o Needs to be in the system path\
        "; /*Help message*/

    char *conBuildHelpCStr = "\
            \n buildCon can use several different consensuses\
            \n   building methods to build a consensus. These methods\
            \n   can be run separately or can be combined together.\
            \n   The order these methods are input is what \
            \n   determines which method is used.\
            \n\
            \n    -min-con-length: [500]\
            \n       - Discard consensuses that are under\
            \n         the input length.\
            \n       - If you lower this your should also\
            \n         lower -min-read-read-map-length &\
            \n         -min-read-con-map-length\
            \n       - use -min-length to lower all three\
            \n    -max-depth:\
            \n        - Max number of reads to use in        [300]\
            \n          a consensus.\
            \n    -extra-consensus-steps:                    [2]\
            \n        - Number of times to rebuild the\
            \n          consensus using a new set of best\
            \n          reads.\
            \n  -use-ivar: [No]\
            \n    o Build a consensus using an ivar step\
            \n    -ivar-min-deth: [10]\
            \n      - Min depth for ivar not to mask.\
            \n    -ivar-min-snp-sup: [0.5] \
            \n      - Min support to not mask an snp \
            \n      - 0 to 1\
            \n    -ivar-min-ins-sup: [0.9] \
            \n      - Min support to not remove an\
            \n        insertion\
            \n      - 0 to 1\
            \n    -ivar-min-q: [10] \
            \n      - Min Q-score needed to keep a base\
            \n    -ivar-mask: [N] \
            \n      - Base to use for masking\
            \n  -use-majCon: [Yes]\
            \n   o Build a consensus using majCon\
            \n   o Disable by inputing a different method\
            \n    -maj-con-min-bases: [0.35=35%]\
            \n        - When building the majority consesus\
            \n          make a deletion in positions that\
            \n          have less than x\% of bases (35%).\
            \n    -maj-con-min-base-q: [7]\
            \n        - Minimum q-score to keep a base when\
            \n          building a majority consensus.\
            \n    -maj-con-min-ins: [0.3=30%]\
            \n        - When building the majority consesus\
            \n          ingore insertions that have support\
            \n          from less than x\% of reads (30%).\
            \n    -maj-con-min-ins-q: [5]\
            \n        - Minimum q-score to keep a insertion\
            \n          when building a majority consensus.\
            \n  -use-racon: [No]\
            \n    o Use racon to build/polish a consensus.\
            \n    o My advice is to use ivar or medaka\
            \n    -rounds-racon: [4]\
            \n      - Number of rounds to polish a\
            \n        consensus with racon\
            \n  -use-medaka: [No]\
            \n    o use medaka to build/polish a consensus\
            \n    -model: [r941_min_high_g351]\
            \n      - Model to use with medaka_consensus\
            \n\
       "; /*consensus building parameters*/

    char
        *readReadMapHelpCStr ="\
            \n The read mapping step uses a higher quality (by median\
            \n    Q-score) read in the bin to identify other high\
            \n    qaulity reads (by median Q-score) that likely came\
            \n    from the same virus. This is done by mapping reads\
            \n    to the selected read and removing all reads that are\
            \n    beneath the quality thresholds. Only the top X\
            \n    (default 300) are kept to build a consensus with.\
            \n\
            \n     -min-read-read-map-length:                     [500]\
            \n        - Minimum aligned read length needed\
            \n          to keep a read to read mapping.\
            \n        - Values less than this will not be\
            \n          extracted.\
            \n        - Change this will also require changing\
            \n          -min-read-con-map-length.\
            \n     -read-read-snps:                      [0.021 = 2.1%]\
            \n        - Minimum percentage of snps needed to\
            \n          discard a read during the read to\
            \n          read clustering step.\
            \n     -read-read-diff:                       [0.02 = 2%]\
            \n        - Minimum percent difference needed to\
            \n          discard a read during the read to\
            \n          read clustering step.\
            \n     -read-read-dels:                       [1 = 100%]\
            \n        - Minimum percentage of deletions\
            \n          needed to discard a read during\
            \n          the read to read clustering step.\
            \n     -read-read-inss:                       [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          needed to discard a read during\
            \n          the read to read clustering step.\
            \n     -read-read-indels:                     [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          and deletions needed to discard a\
            \n          read during the read to read\
            \n          clustering step.\
            \n\
            \n\
            \n The scoring settings for read mapping only counts SNPs\
            \n     and indels in reads that are at or above the user\
            \n     input thresholds. These counts are used to find the\
            \n     percentages in the comparison step.\
            \n\
            \n     -read-read-min-base-q:                 [10]\
            \n        - Minimum Q-score needed to keep an\
            \n          SNP or insertion when comparing\
            \n          reads.\
            \n     -read-read-min-mapq:                   [20]\
            \n        - Minimum mapping quality needed to\
            \n          keep a read when comparing reads.\
            \n\
            \n     -read-read-max-a-ins-homo:                [1]\
            \n        - Maximum A homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-t-ins-homo:                [1]\
            \n        - Maximum T homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-c-ins-homo:                [1]\
            \n        - Maximum C homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-g-ins-homo:                [1]\
            \n        - Maximum G homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n\
            \n     -read-read-max-a-del-homo:                [0]\
            \n        - Maximum A homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-t-del-homo:                [0]\
            \n        - Maximum T homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-c-del-homo:                [0]\
            \n        - Maximum C homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-g-del-homo:                [0]\
            \n        - Maximum G homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n\
            "; /*Read mapping help message*/


    char
        *readConMapHelpCStr = "\
            \n These paramaters are for controlling which reads are\
            \n    kept when mapped to the best read or the consensus.\
            \n\
            \n     -min-read-con-map-length:                   [500]\
            \n        - Minimum aligned read length needed to\
            \n          keep a read to consensus mapping.\
            \n        - Reads that have an alinged length less\
            \n          than this will not be extracted.\
            \n        - Change this will also require changing\
            \n          -min-read-read-map-length.\
            \n     -read-con-snps:                       [0.015 = 1.5%]\
            \n        - Minimum percentage of snps needed to\
            \n          discard a read during the read to\
            \n          consensus clustering step.\
            \n     -read-con-diff:                       [0.02 = 2%]\
            \n        - Minimum percent difference needed to\
            \n          discard a read during the read to\
            \n          consensus clustering step.\
            \n     -read-con-dels:                       [1 = 100%]\
            \n        - Minimum percentage of deletions\
            \n          needed to discard a read during\
            \n          the read to consensus clustering.\
            \n     -read-con-inss:                       [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          needed to discard a read during\
            \n          the read to consensus clustering.\
            \n     -read-con-indels:                     [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          and deletions needed to discard a\
            \n          read during the read to consensus\
            \n          clustering step.\
            \n\
            \n\
            \n The scoring settings for only counts SNPs and indels\
            \n     in reads that are at or above the input thresholds.\
            \n\
            \n     -read-con-min-base-q:                 [10]\
            \n        - Minimum Q-score needed to keep an\
            \n          SNP or insertion when clustering\
            \n          reads.\
            \n     -read-con-min-mapq:                   [20]\
            \n        - Minimum mapping quality needed to\
            \n          keep a read when clustering reads.\
            \n\
            \n     -read-con-max-a-ins-homo:             [1]\
            \n        - Maximum A homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-t-ins-homo:             [1]\
            \n        - Maximum T homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-c-ins-homo:             [1]\
            \n        - Maximum C homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-g-ins-homo:             [1]\
            \n        - Maximum G homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n\
            \n     -read-con-max-a-del-homo:             [0]\
            \n        - Maximum A homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-t-del-homo:             [0]\
            \n        - Maximum T homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-c-del-homo:             [0]\
            \n        - Maximum C homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-g-del-homo:             [0]\
            \n        - Maximum G homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n\
            "; /*Clustering help message*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-2: Set up default settings in structures
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Set default consensus step settings*/
    initConBuildStruct(&conSetting);
    initSamEntry(&samStruct);
    initSamEntry(&refStruct);

    /*Set up filter and scoring settings*/
    blankMinStats(&readToReadMinStats);
    blankMinStatsReadCon(&minReadConStats);

    /*Make sure the struct holding the fastq file does not have noise*/
    fastqStruct.refIdCStr[0] = '\0';
    fastqStruct.statPathCStr[0] = '\0';
    fastqStruct.bestReadCStr[0] = '\0';
    fastqStruct.topReadsCStr[0] = '\0';
    fastqStruct.consensusCStr[0] = '\0';
    fastqStruct.numReadsULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-3: Get user input
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr =
        getUserInput(
            lenArgsInt,
            argsCStr,
            &fqPathCStr,
            &refPathCStr,
            &statsPathCStr,
            prefixCStr,
            threadsCStr,
            &conSetting,
            &readToReadMinStats, /*to keep a read/read map*/
            &minReadConStats   /*To keep read/consensus map*/
    ); /*Get user input*/

    if(tmpCStr != 0)
    { /*If have an error*/

        if(tmpCStr == 0)
        { /*If no user input was supplied*/
            fprintf(
                stderr,
                "%s\n\nNo user input was supplied\n",
                helpCStr
            );

            exit(1);
        } /*If no user input was supplied*/

       if(
            strcmp(tmpCStr, "-v") == 0 ||
            strcmp(tmpCStr, "-V") == 0 ||
            strcmp(tmpCStr, "-version") == 0
        ) { /*If the user is requesting the version number*/
            fprintf(stdout, "buildCon built with findCoInft version:");
            fprintf(stdout, " %.8f\n", defVersion);
            exit(0);
        } /*If the user is requesting the version number*/

        if(strcmp(tmpCStr, "-h-build-consensus") == 0)
        { /*If the user wants the consensus building options*/
            fprintf(stdout, "%s", conBuildHelpCStr);
            exit(0);
        } /*If the user wants the consensus building options*/

        if(strcmp(tmpCStr, "-h-read-read-map") == 0)
        { /*If user wants to know about the read mapping parameters*/
            fprintf(stdout, "%s", readReadMapHelpCStr);
            exit(0);
        } /*If user wants to know about the read mapping parameters*/
        
        if(strcmp(tmpCStr, "-h-read-con-map") == 0)
        { /*If user wants to know about the read mapping parameters*/
            fprintf(stdout, "%s", readConMapHelpCStr);
            exit(0);
        } /*If user wants to know about the read mapping parameters*/

        if(
            strcmp(tmpCStr, "-h") == 0 ||
            strcmp(tmpCStr, "-help") == 0 ||
            strcmp(tmpCStr, "--h") == 0 ||
            strcmp(tmpCStr, "--help") == 0
        ) { /*If user wanted the help message*/
            fprintf(stdout, "%s\n", helpCStr);/*Print out help message*/
            exit(0);
        } /*If user wanted the help message*/

        if(strcmp(tmpCStr, "-ivar-min-depth") == 0)
        { /*If: Ivar min depth was invalid*/
            fputs(
               "Invalid input for -ivar-min-depth\n",
               stderr
            );/*Print out help message*/
            exit(-1);
        } /*If: Ivar min depth was invalid*/

        if(strcmp(tmpCStr, "-ivar-min-snp-sup") == 0)
        { /*If: ivar min snp support invalid*/
            fputs(
               "Invalid input for -ivar-snp-sup\n",
               stderr
            );/*Print out help message*/
            exit(-1);
        } /*If: ivar min snp support invalid*/

        if(strcmp(tmpCStr, "-ivar-min-ins-sup") == 0)
        { /*If: ivar min ins support is invalid*/
            fputs(
               "Invalid input for -ivar-ins-sup\n",
               stderr
            );/*Print out help message*/
            exit(-1);
        } /*If: ivar min ins support is invalid*/

        if(strcmp(tmpCStr, "-ivar-min-q") == 0)
        { /*If: Q-score input for ivar is invalid*/
            fputs(
               "Invalid input for -ivar-min-q\n",
               stderr
            );/*Print out help message*/
            exit(-1);
        } /*If: Q-score input for ivar is invalid*/

        fprintf(
            stderr,
            "%s\n\n%s is an invalid parameter\n",
            helpCStr,   /*Print out the help message*/
            tmpCStr     /*Print out the error*/
        ); /*Let user know about the invalid parameter*/

        exit(1);
    } /*If have an error*/

    if(
          (!conSetting.majConSet.useMajConBl)
       && (!conSetting.raconSet.useRaconBl)
       && (!conSetting.medakaSet.useMedakaBl)
       && (!conSetting.ivarSetST.useIvarBl)
    ) { /*If; No consensus method input*/
        fprintf(stderr, "No consensus method selected\n");
        fprintf(
           stderr,
           "Use: -usr-ivar, -use-medaka, or -use-racon\n"
        );
    }

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-4:
    ^  - Get minimap2, racon, & medaka versions
    ^  o main sec-4 sub-1:
    ^    - Check if fastq file exists
    ^  o main sec-4 sub-2:
    ^    - Find minimap2 version & check if exists
    ^  o main sec-4 sub-3:
    ^    - If using Racon, find version & check if exists
    ^  o main sec-4 sub-4:
    ^    - If using Medaka, find version
    ^  o main sec-4 sub-5:
    ^    - Check if ivar exists
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Main Sec-4 Sub-1:
    *  - Check if the fastq file exists
    \*****************************************************/

    stdinFILE = fopen(fqPathCStr, "r");

    if(stdinFILE == 0)
    { /*If fastq file was invalid*/
        fprintf(
            stdout,
            "Could not open provided fastq file (%s)\n",
             fqPathCStr
         ); /*Let user know about the issue*/

         exit(1);
    } /*If fastq file was invalid*/

    fclose(stdinFILE);

    if(statsPathCStr != 0)
    { /*If the user provided a stats file for read selection*/
        stdinFILE = fopen(statsPathCStr, "r");
        if(stdinFILE == 0)
        { /*If I could not open the provided file*/
            fprintf(
                stdout,
                "Could not open stats file (%s) from scoreReads\n",
                 statsPathCStr
             ); /*Let user know about the issue*/

             exit(1);
        } /*If I could not open the provided file*/

        fclose(stdinFILE);
    } /*If the user provided a stats file for read selection*/


    /*Check if the reference file exists*/
    if(refPathCStr != 0)
    { /*If the user provided a reference file*/
        stdinFILE = fopen(refPathCStr, "r");

        if(stdinFILE == 0)
        { /*If fastq file was invalid*/
            fprintf(
                stdout,
                "Could not open provided reference file (%s)\n",
                 refPathCStr
             ); /*Let user know about the issue*/

             exit(1);
        } /*If fastq file was invalid*/

        fclose(stdinFILE);
    } /*If the user provided a reference file*/

    /*****************************************************\
    * Main Sec-4 Sub-2:
    *  - Find minimap2 version & check if exists
    \*****************************************************/

    /*Set up minimap 2 check*/
    minimap2VersionCStr[0] = '\0';
    tmpCStr = cpParmAndArg(tmpCmdCStr, "minimap2", "--version");

    stdinFILE = popen(tmpCmdCStr, "r");
    fgets(minimap2VersionCStr, 256, stdinFILE);
    fclose(stdinFILE);

    if(minimap2VersionCStr[0] == '\0')
    { /*If could not find minimap2*/
        fprintf(stdout, "Minimap2 could not be found\n");
        exit(1);
    } /*If could not find minimap2*/

    /*****************************************************\
    * Main Sec-4 Sub-3:
    *  - If using Racon, find version & check if exists
    \*****************************************************/

    if(conSetting.raconSet.useRaconBl & 1)
    { /*If using racon, get the version used*/
        /*Set up racon check*/
        raconVersionCStr[0] = '\0';
        tmpCStr = cpParmAndArg(tmpCmdCStr, "racon", "--version");

        stdinFILE = popen(tmpCmdCStr, "r");
        fgets(raconVersionCStr, 256, stdinFILE);
        pclose(stdinFILE);

        if(raconVersionCStr[0] == '\0')
        { /*If racon does not exist*/
            fprintf(stderr, "Racon could not be found\n");
            exit(1);
        } /*If racon does not exist*/
    } /*If using racon, get the version used*/

    /*****************************************************\
    * Main Sec-4 Sub-4:
    *  - If using Medaka, find version & check if exists
    \*****************************************************/

    if(conSetting.medakaSet.useMedakaBl & 1)
    { /*If using medaka, check version*/
        /*Set up non-miniconda command*/
        medakaVersionCStr[0] = '\0';
        tmpCStr = cStrCpInvsDelm(tmpCmdCStr, medakaCMD);
        tmpCStr = cpParmAndArg(tmpCStr, "medaka", "--version");
        tmpCStr = cpParmAndArg(tmpCStr, medakaCMDEnd, "");

        stdinFILE = popen(tmpCmdCStr, "r");
        fgets(medakaVersionCStr, 2048, stdinFILE);
        fclose(stdinFILE);

        if(medakaVersionCStr[0] == '\0')
        { /*If python virtual enviorment medaka does not exist*/
            /*Set up miniconda medaka enviroment command*/
            tmpCStr = cStrCpInvsDelm(tmpCmdCStr, medCondCMD);
            tmpCStr = cpParmAndArg(tmpCStr, "medaka", "--version");
            tmpCStr = cStrCpInvsDelm(tmpCStr, medCondCMDEnd);

            stdinFILE = popen(tmpCmdCStr, "r");
            fgets(medakaVersionCStr, 2048, stdinFILE);
            pclose(stdinFILE);

            if(medakaVersionCStr[0] == '\0')
            { /*If medaka could not be found*/
                fprintf(stderr, "Medaka could not be found\n");
                exit(1);
            } /*If medaka could not be found*/

            /*Mark that I am using medaka from miniconda*/
            conSetting.medakaSet.condaBl = 1;
        } /*If python virtual enviorment medaka does not exist*/
    } /*If using medaka, check version*/

    /*****************************************************\
    * Main Sec-4 Sub-5:
    *  - Check if ivar exists
    \*****************************************************/

    if(conSetting.ivarSetST.useIvarBl)
    { /*If: I am using ivar*/
        ivarVersionStr[0] = '\0';

        tmpCStr = cpParmAndArg(tmpCmdCStr, "ivar", "-v");
        stdinFILE = popen(tmpCmdCStr, "r");
        fgets(ivarVersionStr, 256, stdinFILE);
        pclose(stdinFILE);

        if(ivarVersionStr[0] == '\0')
        { /*If racon does not exist*/
            fprintf(stderr, "Ivar could not be found\n");
            exit(1);
        } /*If racon does not exist*/

        samtoolsVersionStr[0] = '\0';

        tmpCStr =
           cpParmAndArg(tmpCmdCStr, "samtools", "version");

        stdinFILE = popen(tmpCmdCStr, "r");
        fgets(samtoolsVersionStr, 256, stdinFILE);
        pclose(stdinFILE);

        if(samtoolsVersionStr[0] == '\0')
        { /*If racon does not exist*/
           fprintf(stderr,"samtools could not be found\n");
           exit(1);
        } /*If racon does not exist*/
    } /*If: I am using ivar*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-5: Build the consensus
    ^   o main sec-5 sub-1: Make copy of fastq file so orignal is safe
    ^   o main sec-5 sub-2: Make copy of stats file so orignal is safe
    ^   o main sec-5 sub-3: Build the consensus
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-5 Sub-1: Make copy of fastq file so orignal is safe
    \******************************************************************/

    tmpCStr = cStrCpInvsDelm(tmpPathCStr, prefixCStr);
    tmpCStr = cStrCpInvsDelm(tmpCStr, ".fastq");
    copyFile(fqPathCStr, tmpPathCStr);
    strcpy(fastqStruct.fqPathCStr, tmpPathCStr);

    /******************************************************************\
    * Main Sec-5 Sub-2: Make copy of stats file so orignal is safe
    \******************************************************************/

    if(statsPathCStr != 0)
    { /*If I need to copy over the stats file*/
        tmpCStr = cStrCpInvsDelm(tmpPathCStr, prefixCStr);
        tmpCStr = cStrCpInvsDelm(tmpCStr, "--stats.tsv");
        stdinFILE = fopen(statsPathCStr, "r");
        cpFILE = fopen(tmpPathCStr, "w");

        while(fgets(tmpCmdCStr, 1024, stdinFILE))
            fprintf(cpFILE, "%s", tmpCmdCStr);

        fclose(stdinFILE);
        fclose(cpFILE);
        strcpy(fastqStruct.statPathCStr, tmpPathCStr);

        conSetting.useStatBl = 1;
    } /*If I need to copy over the stats file*/

    /******************************************************************\
    * Main Sec-5 Sub-2: Build the consensus
    \******************************************************************/

    errUC =
        buildCon(
            &fastqStruct,     /*Fastq file with reads*/
            refPathCStr,      /*Fasta file with reference, 0 to ignore*/
            threadsCStr,      /*Number threads to use with system calls*/
            &conSetting,
            &samStruct,
            &refStruct,
            &readToReadMinStats,  /*Min stats to keep mapped reads*/
            &minReadConStats    /*Min stats to keep mapped reads*/
    ); /*Build the consensus if possible*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-6: Clean up, Print out version numbers, and exit
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*remove(fastqStruct.bestReadCStr); *//*Make sure no extra files*/
        /*buildCon function already did this*/
    remove(fastqStruct.fqPathCStr);
    remove(fastqStruct.bestReadCStr);

    if(statsPathCStr != 0)
        remove(fastqStruct.statPathCStr); /*Make sure no extra files*/

    freeStackSamEntry(&samStruct);
    freeStackSamEntry(&refStruct);

    if(errUC & 64)
    { /*If had a memory allocation error*/
        fprintf(stdout, "Memory allocation error\n");
        exit(1);
    } /*If had a memory allocation error*/

    fprintf(stdout, "Minimap2 version: %s\n", minimap2VersionCStr);

    if(conSetting.raconSet.useRaconBl & 1)
        fprintf(stdout, "Racon version: %s\n", raconVersionCStr);
    if(conSetting.medakaSet.useMedakaBl & 1)
        fprintf(stdout, "Medaka version: %s\n", medakaVersionCStr);

    if(conSetting.ivarSetST.useIvarBl & 1)
    { /*If: I need to print out version number for ivar*/
       fprintf(stdout,"Ivar version: %s\n",ivarVersionStr);
       fprintf(
          stdout,
          "Samtools version: %s\n",
          samtoolsVersionStr
       );
    } /*If: I need to print out version number for ivar*/

    if(errUC & 16)
    { /*If had a memory allocation error*/
        fprintf(stdout, "Unable to build a consensus\n");
        exit(0);
    } /*If had a memory allocation error*/

    exit(0);
} /*main*/

/*---------------------------------------------------------------------\
| Output:
|    - Modifies input variables to hold user input
|    - Returns:
|        - 0: if nothing went wrong
|        - pionter to invalid paramter
|            - if no agruments input, returns ponter to argsCStr
\---------------------------------------------------------------------*/
char * getUserInput(
    int32_t lenArgsInt,
    char *argsCStr[],
    char **fqPathCStr, /*Holds path to fastq file*/
    char **refPathCStr, /*Holds path to references*/
    char **statsPathCStr, /*Holds path to scoreReads output to use*/
    char *prefixCStr,   /*Prefix to name everything*/
    char *threadsCStr, /*Number threads for minimap2 & racon*/
    struct conBuildStruct *conSet, 
    struct minAlnStats *readToReadMinStats,/*Read pull scoring setting*/
    struct minAlnStats *minReadConStats  /*Read cluster socring set*/
) /*Reads in user input*/
{ /*getUserInput*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: getUserInput
    '    fun-1 sec-1: variable declerations
    '    fun-1 sec-2: Get user input
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *tmpStr = 0;
   char *argStr = 0; /*Points to user input part of line*/
   char *parmStr = 0;   /*Points to argument part of parameter*/
   char firstConMethodBl = 0;
     /*Tells if user supplied a consensus method*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-2: Get user input
    ^    fun-1 sec-2 sub-1: General input
    ^    fun-1 sec-2 sub-2: Percent difference settings
    ^    fun-1 sec-2 sub-3: scoreReads read to reference settings
    ^    fun-1 sec-2 sub-4: scoreReads read to read settings
    ^    fun-1 sec-2 sub-5: scoreReads read to consensus settings
    ^    fun-1 sec-2 sub-6: scoreReads consensus to consensus settings
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-1 Sec-2 Sub-1:
    *  - General input
    \*****************************************************/

    if(lenArgsInt < 2)
        return 0;       /*Nothing input*/

    for(int32_t intArg = 1; intArg < lenArgsInt; intArg += 2)
    { /*Loop through all user arguments (for)*/
        parmStr = *(argsCStr + intArg); /*Get the parameter used*/
        argStr = *(argsCStr + intArg + 1); /*setting*/

        if(strcmp(parmStr, "-fastq") == 0)
            *fqPathCStr = argStr;  /*Fastq file to check*/

        else if(strcmp(parmStr, "-stats") == 0)
            *statsPathCStr = argStr; /*Stats file with reads to use*/

        else if(strcmp(parmStr, "-ref") == 0)
            *refPathCStr = argStr;  /*references*/

        else if(strcmp(parmStr, "-prefix") == 0)
            strcpy(prefixCStr, argStr); /*Prefix to name files with*/

        else if(strcmp(parmStr, "-threads") == 0)
            strcpy(threadsCStr, argStr); /*Number of threads*/

        else if(strcmp(parmStr, "-min-depth") == 0)
            conSet->minReadsToBuildConUL=strtoul(argStr,&tmpStr,10);

        else if(strcmp(parmStr, "-max-depth") == 0)
            conSet->maxReadsToBuildConUL=strtoul(argStr,&tmpStr,10);

        else if(strcmp(parmStr, "-min-perc") == 0)
           conSet->minPercMappedReadsFlt =
              (float) atof(argStr);
              /*Not the best way, but works*/

        else if(strcmp(parmStr, "-extra-consensus-steps") == 0)
           strToUIBase10(
              argStr,
              &conSet->numRndsToPolishUI
           ); /*Get base 10 unsigned int*/

        else if(strcmp(parmStr, "-min-length") == 0)
        { /*Else if: setting min length for everything*/
            strToUIBase10(
               argStr,
               &conSet->minConLenUI
            ); /*Get base 10 unsigned int*/

            minReadConStats->minReadLenULng =
               conSet->minConLenUI;

            readToReadMinStats->minReadLenULng =
               conSet->minConLenUI;
        } /*Else if: setting min length for everything*/

        else if(strcmp(parmStr, "-min-con-length") == 0)
            conSet->minConLenUI = strtoul(argStr, &tmpStr, 10);

        /*************************************************\
        * Fun-1 Sec-2 Sub-2:
        *  - Selecting the consensus methods
        \*************************************************/

        else if(strcmp(parmStr, "-use-racon") == 0)
        { /*Else if: using racon*/
            if(!firstConMethodBl)
            { /*If: this is the first input method*/
               conSet->lenMethodUC = 0;
               firstConMethodBl = 1;
            } /*If: this is the first input method*/

            conSet->methodAryUC[conSet->lenMethodUC] =
               defUseRacon;

            ++(conSet->lenMethodUC);

            conSet->methodAryUC[conSet->lenMethodUC] =
               defNoCon;

            --intArg;

            /*Make sure I know racon was used*/
            conSet->raconSet.useRaconBl = 1;
        } /*Else if: using racon*/

        else if(strcmp(parmStr, "-use-medaka") == 0)
        { /*Else if: using medaka*/
            if(!firstConMethodBl)
            { /*If: this is the first input method*/
               conSet->lenMethodUC = 0;
               firstConMethodBl = 1;
            } /*If: this is the first input method*/

            conSet->methodAryUC[conSet->lenMethodUC] =
               defUseMedaka;

            ++(conSet->lenMethodUC);

            conSet->methodAryUC[conSet->lenMethodUC] =
               defNoCon;

            --intArg;
            /*Make sure I know medaka was used*/
            conSet->medakaSet.useMedakaBl = 1;
        } /*Else if: using medaka*/

        else if(strcmp(parmStr, "-use-majCon") == 0)
        { /*Else if: using the majority consensus*/
            if(!firstConMethodBl)
            { /*If: this is the first input method*/
               conSet->lenMethodUC = 0;
               firstConMethodBl = 1;
            } /*If: this is the first input method*/

            conSet->methodAryUC[conSet->lenMethodUC] =
               defUseMajCon;

            ++(conSet->lenMethodUC);

            conSet->methodAryUC[conSet->lenMethodUC] =
               defNoCon;

            --intArg;

            /*So I know that majority consensus was used*/
            conSet->majConSet.useMajConBl = 1;
        } /*Else if: using the majority consensus*/

        else if(strcmp(parmStr, "-use-ivar") == 0)
        { /*Else if: using ivar*/
            if(!firstConMethodBl)
            { /*If: this is the first input method*/
               conSet->lenMethodUC = 0;
               firstConMethodBl = 1;
            } /*If: this is the first input method*/

            conSet->methodAryUC[conSet->lenMethodUC] =
               defUseIvar;

            ++(conSet->lenMethodUC);

            conSet->methodAryUC[conSet->lenMethodUC] =
               defNoCon;

            --intArg;

            /*So I know that majority consensus was used*/
            conSet->ivarSetST.useIvarBl = 1;
        } /*Else if: using ivar*/

        /*************************************************\
        * Fun-1 Sec-2 Sub-2:
        *  - Settings unique to ivar
        \*************************************************/

        else if(strcmp(parmStr, "-ivar-min-depth") == 0)
        { /*Else if: Setting the min depth to not mask*/
           if(isUInt(argStr))
             strcpy(conSet->ivarSetST.minDepthStr,argStr);
           /*else return parmStr;*/
        } /*Else if: Setting the min depth to not mask*/

        else if(strcmp(parmStr, "-ivar-min-snp-sup") == 0)
        { /*Else if: Setting the min snp support*/
           if(isUDec(argStr))
             strcpy(conSet->ivarSetST.minSupStr,argStr);
           else return parmStr;
        } /*Else if: Setting the min snp support*/

        else if(strcmp(parmStr, "-ivar-min-ins-sup") == 0)
        { /*Else if: Setting the min insterion support*/
           if(isUDec(argStr))
            strcpy(conSet->ivarSetST.minInsSupStr,argStr);
           else return parmStr;
        } /*Else if: Setting the min insterion support*/

        else if(strcmp(parmStr, "-ivar-min-q") == 0)
        { /*Else if: Setting the base Q-score*/
           if(isUInt(argStr))
             strcpy(conSet->ivarSetST.minQStr, argStr);
           else return parmStr;
        } /*Else if: Setting the base Q-score*/

        else if(strcmp(parmStr, "-ivar-mask") == 0)
           conSet->ivarSetST.maskC = *parmStr;

        /*************************************************\
        * Fun-1 Sec-2 Sub-2:
        *  - Settings unique to medaka and racon
        \*************************************************/

        else if(strcmp(parmStr, "-model") == 0)
            strcpy(conSet->medakaSet.modelCStr, argStr);

        else if(strcmp(parmStr, "-rounds-racon") == 0)
            cStrToUChar(argStr, &conSet->raconSet.rndsRaconUC);
 
        /*************************************************\
        * Fun-1 Sec-2 Sub-2:
        *  - Majority consensus settings
        \*************************************************/

        else if(strcmp(parmStr, "-maj-con-min-bases") == 0)
          sscanf(argStr,"%f",&conSet->majConSet.minReadsPercBaseFlt);

        else if(strcmp(parmStr, "-maj-con-min-ins") == 0)
           sscanf(argStr,"%f",&conSet->majConSet.minReadsPercInsFlt);

        else if(strcmp(parmStr, "-maj-con-min-base-q") == 0)
            cStrToUChar(argStr, &conSet->majConSet.minBaseQUC);

        else if(strcmp(parmStr, "-maj-con-min-ins-q") == 0)
            cStrToUChar(argStr, &conSet->majConSet.minInsQUC);
            
        /*************************************************\
        * Fun-1 Sec-2 Sub-2: Percent difference settings
        \*************************************************/

        else if(strcmp(parmStr, "-min-read-con-map-length") == 0)
        { /*Else if the user provided a minimum read length*/
            minReadConStats->minReadLenULng =
                strtoul(argStr, &tmpStr, 10);
        } /*Else if the user provided a minimum read length*/

        else if(strcmp(parmStr, "-read-con-snps") == 0)
            sscanf(argStr, "%f", &minReadConStats->minSNPsFlt);
        else if(strcmp(parmStr, "-read-con-diff") == 0)
            sscanf(argStr, "%f", &minReadConStats->minDiffFlt);
        else if(strcmp(parmStr, "-read-con-dels") == 0)
            sscanf(argStr, "%f", &minReadConStats->minDelsFlt);
        else if(strcmp(parmStr, "-read-con-inss") == 0)
            sscanf(argStr, "%f", &minReadConStats->minInssFlt);
        else if(strcmp(parmStr, "-read-con-indels") == 0)
            sscanf(argStr, "%f", &minReadConStats->minIndelsFlt);

        /**************************************************************\
        * Fun-1 Sec-2 Sub-5: scoreReads read to consensus settings
        \**************************************************************/

        else if(strcmp(parmStr, "-read-con-min-base-q") == 0)
            cStrToUChar(argStr, &minReadConStats->minQChar);

        else if(strcmp(parmStr, "-read-con-min-mapq") == 0)
            cStrToUInt(argStr, &minReadConStats->minMapqUInt);

         else if(strcmp(parmStr, "-read-con-max-a-ins-homo") == 0)
           cStrToUInt(argStr, &minReadConStats->maxHomoInsAry[0]);

         else if(strcmp(parmStr, "-read-con-max-t-ins-homo") == 0)
           cStrToUInt(argStr,&minReadConStats->maxHomoInsAry[10]);

         else if(strcmp(parmStr, "-read-con-max-c-ins-homo") == 0)
            cStrToUInt(argStr,&minReadConStats->maxHomoInsAry[1]);

         else if(strcmp(parmStr, "-read-con-max-g-ins-homo") == 0)
            cStrToUInt(argStr,&minReadConStats->maxHomoInsAry[3]);

         else if(strcmp(parmStr, "-read-con-max-a-del-homo") == 0)
           cStrToUInt(argStr, &minReadConStats->maxHomoDelAry[0]);

         else if(strcmp(parmStr, "-read-con-max-t-del-homo") == 0)
           cStrToUInt(argStr,&minReadConStats->maxHomoDelAry[10]);

         else if(strcmp(parmStr, "-read-con-max-c-del-homo") == 0)
            cStrToUInt(argStr,&minReadConStats->maxHomoDelAry[1]);

         else if(strcmp(parmStr, "-read-con-max-g-del-homo") == 0)
            cStrToUInt(argStr,&minReadConStats->maxHomoDelAry[3]);

        /**************************************************************\
        * Fun-1 Sec-2 Sub-2: Percent difference settings
        \**************************************************************/

        else if(strcmp(parmStr, "-min-read-read-map-length") == 0)
        { /*Else if the user provided a minimum read length*/
            readToReadMinStats->minReadLenULng =
                strtoul(argStr, &tmpStr, 10);
        } /*Else if the user provided a minimum read length*/

        else if(strcmp(parmStr, "-read-read-snps") == 0)
            sscanf(argStr, "%f", &readToReadMinStats->minSNPsFlt);
        else if(strcmp(parmStr, "-read-read-diff") == 0)
            sscanf(argStr, "%f", &readToReadMinStats->minDiffFlt);
        else if(strcmp(parmStr, "-read-read-dels") == 0)
            sscanf(argStr, "%f", &readToReadMinStats->minDelsFlt);
        else if(strcmp(parmStr, "-read-read-inss") == 0)
            sscanf(argStr, "%f", &readToReadMinStats->minInssFlt);
        else if(strcmp(parmStr, "-read-read-indels") == 0)
            sscanf(argStr, "%f", &readToReadMinStats->minIndelsFlt);

        /**************************************************************\
        * Fun-1 Sec-2 Sub-5: scoreReads read to readsensus settings
        \**************************************************************/

        else if(strcmp(parmStr, "-read-read-min-base-q") == 0)
            cStrToUChar(argStr, &readToReadMinStats->minQChar);

        else if(strcmp(parmStr, "-read-read-min-mapq") == 0)
            cStrToUInt(argStr, &readToReadMinStats->minMapqUInt);

         else if(strcmp(parmStr, "-read-read-max-a-ins-homo") == 0)
           cStrToUInt(argStr, &readToReadMinStats->maxHomoInsAry[0]);

         else if(strcmp(parmStr, "-read-read-max-t-ins-homo") == 0)
           cStrToUInt(argStr,&readToReadMinStats->maxHomoInsAry[10]);

         else if(strcmp(parmStr, "-read-read-max-c-ins-homo") == 0)
            cStrToUInt(argStr,&readToReadMinStats->maxHomoInsAry[1]);

         else if(strcmp(parmStr, "-read-read-max-g-ins-homo") == 0)
            cStrToUInt(argStr,&readToReadMinStats->maxHomoInsAry[3]);

         else if(strcmp(parmStr, "-read-read-max-a-del-homo") == 0)
           cStrToUInt(argStr, &readToReadMinStats->maxHomoDelAry[0]);

         else if(strcmp(parmStr, "-read-read-max-t-del-homo") == 0)
           cStrToUInt(argStr,&readToReadMinStats->maxHomoDelAry[10]);

         else if(strcmp(parmStr, "-read-read-max-c-del-homo") == 0)
            cStrToUInt(argStr,&readToReadMinStats->maxHomoDelAry[1]);

         else if(strcmp(parmStr, "-read-read-max-g-del-homo") == 0)
            cStrToUInt(argStr,&readToReadMinStats->maxHomoDelAry[3]);

        else
            return parmStr;
    } /*Loop through all user arguments (for)*/

    return 0;
} /*getUserInput*/
