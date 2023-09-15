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

#endif
