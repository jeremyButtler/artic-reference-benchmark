/*#########################################################
# Name: cStrFun
# Use:
#  o Holds functions for copying or manipualting c-strings
# Libraries:
# C Standard Libraries:
#########################################################*/

#ifndef CSTRFUN_H
#define CSTRFUN_H

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' cStrFun SOH: Start Of Header
'  - Has functions that operate on c-strings
'  o Fun-01 cStrCpInvsDelim:
'    - Copy one c-string till an tab, newline, or '\0'.
'      Spaces are kept.
'  o fun-02 cpSpaceCStr:
'    - Copies two c-strings into a c-string. A space is
'      inserted before each copied c-string.
'  o fun-03 strCpTill:
'    - Copies one c-string to another c-string untill
'      either a null is reached, the specified deliminator
'      is reached, or the buffer limit is reached.
'  o macro-01 isuint:
'     - Checks if a c-string is an unsigned number
'  o macro-02 issint:
'     - Checks if a c-string is an signed number
'  o macro-03 isudec:
'     - Checks if a c-string is an unsigned decimal number
'  o macro-04 issdec:
'     - Checks if a c-string is an signed decimal number
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: cStrCpInvsDelm (Fun-01:)
| Use:
|  - Copies one c-string to another c-string untill the
|    first tab, newline, or null. The pointer to the null
|    at the end of the buffer c-string is returned.
| Input:
|  - cpToCStr:
|    o C-string to copy string to
|  - cpFromCStr:
|    o C-string with string to copy to cpToCStr
| Output:
|  - Modifies:
|    o cpToCStr to hold the copied C-string
|  - Returns:
|    o pointer to null at end of cpToCStr
\--------------------------------------------------------*/
static inline char * cStrCpInvsDelm(
    char *cpToCStr,  /*C-string to copy values to*/
    char *cpFromCStr /*C-string to copy*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: cStrCpInvsDelim
   '  - Copy one c-string till an tab, newline, or '\0'.
   '    Spaces are kept.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    while(*cpFromCStr > 31)
    { /*While have a  c-string to copy*/
        *cpToCStr = *cpFromCStr;
        ++cpToCStr;
        ++cpFromCStr;
    } /*While have a  c-string to copy*/

    *cpToCStr = '\0';
    return cpToCStr;
} /*cStrCpInvsDelim*/

/*--------------------------------------------------------\
| Name: cpParmAndArg (Fun-02:)
| Use:
|  - Copies two c-strings into a c-string. A space is
|    inserted before each copied c-string.
|  - This function was designed to copy program arguments
|    into a c-string to use for a system call.
| Input:
|  - cpToCStr:
|    o C-string to copy cpParmCStr and cpAryCStr into
|  - cpParmCStr:
|    o First item to copy into cpToCStr (the parameter)
|    o A space is inserted before and after the parameter
|  - cpArgCStr:
|    o The second item to copy into cpToCStr (the argument)
| Output:
|  - Modifies:
|    o cpToCStr to hold space, parameter, space, & argument
|  - Returns:
|    o pointer to the null at end of cpToCStr
\--------------------------------------------------------*/
static inline char * cpParmAndArg(
    char *cpToCStr, /*Holds copied parameter and argement*/
    char *cpParmCStr,/*Paramater to copy*/
    char *cpArgCStr  /*Argument to copy*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: cpSpaceCStr
   '  - Copies two c-strings into a c-string. A space is
   '    inserted before each copied c-string.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    *cpToCStr = ' ';
    ++cpToCStr;
    
    while(*cpParmCStr > 31)
    { /*While have a  c-string to copy*/
        *cpToCStr = *cpParmCStr;
        ++cpToCStr;
        ++cpParmCStr;
    } /*While have a  c-string to copy*/

    *cpToCStr = ' ';
    ++cpToCStr;

    while(*cpArgCStr > 31)
    { /*While have a  c-string to copy*/
        *cpToCStr = *cpArgCStr;
        ++cpToCStr;
        ++cpArgCStr;
    } /*While have a  c-string to copy*/

    *cpToCStr = '\0';
    return cpToCStr;
} /*cpParmAndArg*/

/*--------------------------------------------------------\
| Name: strCpTill (Fun-03:)
| Use:
|  - Copies one c-string to another c-string untill either
|    a null is reached, the specified deliminator is
|    reached, or the buffer limit is reached.
| Input:
|  - buffStr:
|    o Buffer to hold the copied c-string
|  - toCpStr:
|    o C-string to copy to buffStr
|  - tillC:
|    o Deliminator to stop at (other than '\0')
|  - limitI:
|    o Number of empty characters in buffStr
| Output:
|  - Modifies:
|    o buffStr to hold contents of toCpStr
|    o limitI to hold the empty length remaining
|  - Returns:
|    o pointer to null (end of) buffStr
\--------------------------------------------------------*/
static inline char * strCpTill(
   char *buffStr, /*Buffer to hold copy*/
   char *toCpStr, /*C-string to copy*/
   char tillC,    /*Stop copying at this character*/
   int *limitI    /*Empty space in buffer*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: strCpTill
   '  - Copy a c-string till a null ('\0') or character in
   '    tillC is found or limitI is reached
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   while(*toCpStr != tillC && *toCpStr != '\0' && limitI>0)
   { /*Loop: Copy the c-string over*/
      *buffStr = *toCpStr;
      ++buffStr;
      ++toCpStr;
      --(*limitI);
   } /*Loop: Copy the c-string over*/

   *buffStr = '\0';
   return buffStr;
} /*strCpTill*/

/*--------------------------------------------------------\
| Name: isUInt (Macro-01:)
| Use:
|  - Checks if a c-string is an unsigned number
| Input:
|  - strIn:
|    o Pointer to string to check
| Output:
|  - Returns:
|    o 1 if numeric
|    o 0 if not numeric
\--------------------------------------------------------*/
#define isUInt(strIn)({\
   char *iterStr = (strIn);\
   while(*iterStr > 47 && *iterStr < 58) ++iterStr;\
   *iterStr == '\0';\
}) /*isUInt*/

/*--------------------------------------------------------\
| Name: isSInt (Macro-02:)
| Use:
|  - Checks if a c-string is an signed number
| Input:
|  - strIn:
|    o Pointer to string to check
| Output:
|  - Returns:
|    o 1 if numeric
|    o 0 if not numeric
\--------------------------------------------------------*/
#define isSInt(strIn)({ \
   char *iterStr = (strIn) + (*(strIn) == '-'); \
   while(*iterStr > 47 && *iterStr < 58) ++iterStr;\
   *iterStr == '\0'; \
}) /*isUInt*/

/*--------------------------------------------------------\
| Name: isUDec (Macro-03:)
| Use:
|  - Checks if a c-string is an unsigned decimal number
| Input:
|  - strIn:
|    o Pointer to string to check
| Output:
|  - Returns:
|    o 1 if numeric
|    o 0 if not numeric
\--------------------------------------------------------*/
#define isUDec(strIn)({ \
   char *iterStr = (strIn); \
   while(*iterStr > 47 && *iterStr < 58) ++iterStr;\
   \
   iterStr += (*iterStr == '.'); \
   while(*iterStr > 47 && *iterStr < 58) ++iterStr;\
   \
   *iterStr == '\0'; \
})/*isUDec*/

/*--------------------------------------------------------\
| Name: isSDec (Macro-04:)
| Use:
|  - Checks if a c-string is an signed decimal number
| Input:
|  - strIn:
|    o Pointer to string to check
| Output:
|  - Returns:
|    o 1 if numeric
|    o 0 if not numeric
\--------------------------------------------------------*/
#define isSDec(strIn)({ \
   char *iterStr = (strIn) + (*(strIn) == '-'); \
   while(*iterStr > 47 && *iterStr < 58) ++iterStr;\
   \
   iterStr += (*iterStr == '.'); \
   while(*iterStr > 47 && *iterStr < 58) ++iterStr;\
   \
   *iterStr == '\0'; \
})/*isUDec*/

#endif
