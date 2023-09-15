/*#########################################################
# Name: cStrToNumberFun
# Use:
#  - These are my numeric functions to convert numbers to
#    specific data types. 
#  - This allows me to avoid issues of overflows form the
#    base c functions, wich always convert to unsigned
#    longs or longs.
# Libraries:
#   - "dataTypeShortHand.h"
# C Standard Libraries:
#   - <stdint.h>
# Note:
#  - This uses set sizes for each data type & so is a bit
#    less universall than the base functions.
#  - Currently this is set up only for base 10.
#########################################################*/

#ifndef CSTRTONUMBERFUN_H
#define CSTRTONUMBERFUN_H

#include "dataTypeShortHand.h" /*unsigned var to uvar*/
#include <stdint.h> /*uintx_t variables*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' cStrToNumberFun SOF: Start Of Header
'  - Functions to converts c-strings to numbers
'  o fun-01 cStrToUInt:
'    - converst a base 10 c-string to an uint32_t
'  o fun-02 cStrToUSht:
'    - converst a base 10 c-string to an uint16_t
'  o fun-03 cStrToUChar:
'    - converst a base 10 c-string to an uint8_t
'  o fun-04 backwarsCStrToUInt:
'    - converst a backwards (lowest value first) base 10
'      c-string to an uint32_t
'  o fun-05 uCharToCStr:
'    - Converts an uint8_t to an c-string
'  o fun-06 cStrToInt16:
'    - converst a base 10 c-string to an int16_t
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: cStrToUInt (Fun-01:)
| Use:
|  - Converts a base 10 c-string to an 32 bit unsigned int
| Input:
|  - charUCStr:
|    o C-string with base 10 number
|  - retUInt:
|    o 32 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUInt to hold the uint32_t integer
| Note:
|  - This function will only convert the first number of
|    digits in an uint32_t unsigned integer or till
|    non-numeric character
\--------------------------------------------------------*/
static inline char * cStrToUInt(
    char *charUCStr,     /*C-string with base 10 number*/
    uint32_t *retUInt   /*Holds converted number*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: cStrToUInt
   '  - converst a base 10 c-string to an uint32_t
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*1st digit*/
    /*Convert first number*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUInt = *charUCStr - 48;
    else
    { /*Else non-numeric*/
        *retUInt = 0;
        return charUCStr;
    } /*Else non-numeric*/

    ++charUCStr;

    for(uint32_t intCnt = 0; intCnt < 7; ++intCnt)
    { /*Convert digits with no overflow concern*/
        if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
            *retUInt = *retUInt * 10 + *charUCStr - 48;
        else
            return charUCStr;
        ++charUCStr;
    } /*Convert digits with no overflow concern*/

    /*9th digit*/
    /*Convert last two digits, which could overflow*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
    { /*If have one or two more numbers*/
        if(10 * *retUInt + *charUCStr - 48 > 10000)
        { /*If can fit in one more number*/
            *retUInt = *retUInt * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If can fit in one more number*/

        else
            return charUCStr;
    } /*If have one or two more numbers*/

    else return charUCStr;

    /*10th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
    { /*If have one more number*/
        if(10 * *retUInt + *charUCStr - 48 > 10000)
        { /*If would not cause an overflow*/
            /*> 1000 will keep value as int*/
            *retUInt = *retUInt * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If would not cause an overflow*/
    } /*If have one more number*/

    return charUCStr;
} /*cStrToUInt*/

/*--------------------------------------------------------\
| Name: cStrToUSht (Fun-02:)
| Use:
|  - Converts a base 10 c-string to an 16 bit unsigned int
| Input:
|  - charUCStr:
|    o C-string with base 10 number
|  - retUSht:
|    o 16 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUInt to hold the uint16_t
| Note:
|  - This function will only convert the first number of
|    digits in an uint16_t unsigned integer or till
|    non-numeric character
\--------------------------------------------------------*/
static inline char * cStrToUSht(
    char *charUCStr,    /*C-string with base 10 number*/
    uint16_t *retUSht   /*Holds converted number*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: cStrToUSht
   '  - converst a base 10 c-string to an uint16_t
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*1st digit*/
    /*Convert first number*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUSht = *charUCStr - 48;
    else
    { /*Else non-numeric*/
        *retUSht = 0;
        return charUCStr;
    } /*Else non-numeric*/

    ++charUCStr;

    /*2nd digit*/
    /*digits were their is no risk of overlfow (2-8th)*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUSht = *retUSht * 10 + *charUCStr - 48;
    else
        return charUCStr;
    ++charUCStr;

    /*3rd digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUSht = *retUSht * 10 + *charUCStr - 48;
    else
        return charUCStr;
    ++charUCStr;

    /*Convert last two digits, which could overflow*/
    /*4th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
    { /*If have one or two more numbers*/
        if(10 * *retUSht + *charUCStr - 48 < 65536)
        { /*If can fit in one more number*/
            *retUSht = *retUSht * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If can fit in one more number*/

        else
            return charUCStr;
    } /*If have one or two more numbers*/
    else

        return charUCStr;
    /*5th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
    { /*If have one more number*/
        if(10 * *retUSht + *charUCStr - 48 > 100)
        { /*If would not cause an overflow*/
            /*> 1000 will keep value as int 16*/
            *retUSht = *retUSht * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If would not cause an overflow*/
    } /*If have one more number*/

    return charUCStr;
} /*cStrToUSht*/

/*--------------------------------------------------------\
| Name: cStrToUChar (Fun-03:)
| Use:
|  - Converts a base 10 c-string to an 8 bit unsigned int
| Input:
|  - charUCStr:
|    o C-string with base 10 number
|  - retUChar:
|    o 8 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUInt to hold the uint8_t
| Note:
|  - This function will only convert the first number of
|    digits in an uint8_t unsigned integer or till
|    non-numeric character
\--------------------------------------------------------*/
static inline char * cStrToUChar(
    char *charUCStr,       /*C-string with base 10 number*/
    unsigned char *retUChar /*Holds converted number*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: cStrToUChar
   '  - converst a base 10 c-string to an uint8_t
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(*charUCStr < 58 && *charUCStr > 47)
        *retUChar = *charUCStr - 48;
    else
    { /*Else is a non-valid character*/
        *retUChar = 0;
        return charUCStr;
    } /*Else is a non-valid character*/

    ++charUCStr;

    if(*charUCStr < 58 && *charUCStr > 47)
        *retUChar = *retUChar * 10 + *charUCStr - 48;
    else
        return charUCStr;

    ++charUCStr;

    if(*charUCStr < 58 && *charUCStr > 47)
        if(*retUChar * 10 + *charUCStr - 48 < 256)
        { /*If we can add in one more digit*/
            *retUChar = *retUChar * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If we can add in one more digit*/

    return charUCStr;
} /*cStrToUInt*/

/*--------------------------------------------------------\
| Name: backwarsCStrToUInt (Fun-04:)
| Use:
|  - Converts a backwards (lowest value first) base 10
|    c-string to an 32 bit unsigned int
| Input:
|  - charUCStr:
|    o C-string with a backwards (lowest digit first) base
|      10 number
|  - retUChar:
|    o 32 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUInt to hold the uint32_t
| Note:
|  - This function will only convert the first number of
|    digits in an uint32_t unsigned integer or till
|    non-numeric character
\--------------------------------------------------------*/
static inline char * backwarsCStrToUInt(
    char *charUCStr,    /*backwards base 10 number*/
    uint32_t *retUInt   /*Holds converted number*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: backwarsCStrToUInt
   '  - converst a backwards (lowest value first) base 10
   '    c-string to an uint32_t
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*1st digit*/
    /*Convert first number*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUInt = *charUCStr - 48;
    else
    { /*Else non-numeric*/
        *retUInt = 0;
        return charUCStr;
    } /*Else non-numeric*/

    --charUCStr;

    /*2nd digit*/
    /*digits were their is no risk of overlfow (2-8th)*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUInt = *retUInt + (*charUCStr - 48) * 10;
    else
        return charUCStr;
    --charUCStr;

    /*3rd digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUInt = *retUInt + (*charUCStr - 48) * 100;
    else
        return charUCStr;
    --charUCStr;

    /*4th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUInt = *retUInt + (*charUCStr - 48) * 1000;
    else
        return charUCStr;
    --charUCStr;

    /*5th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUInt = *retUInt + (*charUCStr - 48) * 10000;
    else
        return charUCStr;
    --charUCStr;

    /*6th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUInt = *retUInt + (*charUCStr - 48) * 100000;
    else
        return charUCStr;
    --charUCStr;

    /*7th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUInt = *retUInt + (*charUCStr - 48) * 1000000;
    else
        return charUCStr;
    --charUCStr;

    /*8th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
        *retUInt = *retUInt + (*charUCStr - 48) * 10000000;
    else
        return charUCStr;
    --charUCStr;

    /*Convert last two digits, which could overflow*/
    /*9th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
    { /*If have one or two more numbers*/
        if(*retUInt + (*charUCStr - 48) * 10000000 > 10000)
        { /*If can fit in one more number*/
            *retUInt = *retUInt + (*charUCStr - 48) * 10000000;
            --charUCStr;
        } /*If can fit in one more number*/

        else
            return charUCStr;
    } /*If have one or two more numbers*/
    else
        return charUCStr;

    /*10th digit*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
    { /*If have one more number*/
        if(*retUInt + (*charUCStr - 48) * 100000000 > 1000)
        { /*If would not cause an overflow*/
            *retUInt =
               *retUInt + (*charUCStr - 48) * 100000000;
            --charUCStr;
        } /*If would not cause an overflow*/
    } /*If have one more number*/

    return charUCStr;
} /*backwardsCStrToUInt*/

/*--------------------------------------------------------\
| Name: uCharToCStr (Fun-05:)
| Use:
|  - Converts an uint8_t to an c-string
| Input:
|  - buffCStr:
|    o Buffer to hold the base 10 number in uCharToCnvt
|    o This buffer needs at least 4 elements.
|  - uCharToCnvt:
|    o base 10 uint8_t number to convert to c-string
| Output:
|  - Returns:
|    o pointer to end of buffCStr
|  - Modifies:
|    o buffCStr to hold the base 10 number
\--------------------------------------------------------*/

static inline char * uCharToCStr(
    char *buffCStr,     /*4 element buffer to hold number*/
    unsigned char uCharToCnvt /*Char to convert*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: uCharToCStr
   '  - Converts an uint8_t to an c-string
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char onesUC = 0;
    char tensUC = 0;

    onesUC = (uCharToCnvt % 10) + 48;

    /*Move to next digit in character*/
    uCharToCnvt = uCharToCnvt / 10;

    if(uCharToCnvt > 0)
        tensUC = (uCharToCnvt % 10) + 48;
    else
    { /*Else their is only a ones position*/
       *buffCStr = onesUC;
       ++buffCStr;

       *buffCStr = '\0';
       return buffCStr;
    } /*Else their is only a ones position*/

    /*Move to next digit in character*/
    uCharToCnvt = uCharToCnvt / 10;

    if(uCharToCnvt > 0)
    { /*If looking at the 100s position*/
        *buffCStr = (uCharToCnvt % 10) + 48;
        ++buffCStr;
    } /*If looking at the 100s position*/

    *buffCStr = tensUC;
    ++buffCStr;

    *buffCStr = onesUC;
    ++buffCStr;

    *buffCStr = '\0';
    return buffCStr;
} /*uCharToCStr*/

/*--------------------------------------------------------\
| Name: cStrToInt16 (Fun-06:)
| Use:
|  - Converts a base 10 c-string to an 16 bit signed int
| Input:
|  - charUCStr:
|    o C-string with base 10 number
|  - retUSht:
|    o 16 bit signed integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUInt to hold the int16_t
| Note:
|  - This function will only convert the first number of
|    digits in an int16_t integer or till a non-numeric
|    character
\--------------------------------------------------------*/
static inline char * cStrToInt16(
    char *cStrToCnvt, /*CString with base 10 number*/
    int16_t *numS     /*Holds the converted signed short*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: cStrToInt16
   '  - converst a base 10 c-string to an int16_t
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    *numS = 0;

    /*Check if it is a negative number*/
    if(*cStrToCnvt == '-')
    { /*If have a negative symbol*/
        ++cStrToCnvt; /*Move of the negative symbol*/

        if(*(cStrToCnvt)-58 < 0 && *(cStrToCnvt)-47 > 0)
            *numS = -1 * (*cStrToCnvt - 48);

        /*- alone is not a numeric value*/
        else return (cStrToCnvt - 1);
    } // If have a negative symbol

    else if(*cStrToCnvt - 58 < 0 && *cStrToCnvt - 47 > 0)
        *numS = *cStrToCnvt - 48;

    else return cStrToCnvt;

    ++cStrToCnvt;    /*Move to the next character*/

    /*convert the next two digits (no overflow risk*/
    for(uint32_t intCnt = 0; intCnt < 2; ++intCnt)
    { /*Convert digits with no overflow concern*/
        if(*cStrToCnvt - 58 < 0 && *cStrToCnvt - 47 > 0)
            *numS = *numS * 10 + *cStrToCnvt - 48;

        else return cStrToCnvt;

        ++cStrToCnvt;
    } /*Convert digits with no overflow concern*/

    /*4th digit*/
    /*Convert last two digits, which could overflow*/
    if(*cStrToCnvt - 58 < 0 && *cStrToCnvt - 47 > 0)
    { /*If have one or two more numbers*/
        if(10*(*numS & (~(1 << 15)))+*cStrToCnvt-48 > 1000)
        { /*If can fit in one more number*/
            *numS = *numS * 10 + *cStrToCnvt - 48;
            ++cStrToCnvt;
        } /*If can fit in one more number*/

        else
            return cStrToCnvt;
    } /*If have one or two more numbers*/
    else
        return cStrToCnvt;

    /*5th digit*/
    if(*cStrToCnvt - 58 < 0 && *cStrToCnvt - 47 > 0)
    { /*If have one more number*/
        if(10*(*numS & (~(1 << 15)))+*cStrToCnvt-48 > 1000)
        { /*If would not cause an overflow*/
            *numS = *numS * 10 + *cStrToCnvt - 48;
            ++cStrToCnvt;
        } /*If would not cause an overflow*/
    } /*If have one more number*/

    return cStrToCnvt;
} /*cStrToInt16*/

/*--------------------------------------------------------\
| Name: strToULBase10 (Fun-07:)
| Use:
|  - Converts a base 10 c-string to an unsigned long
| Input:
|  - numericStr:
|    o C-string with base 10 number to convert
|  - retUL:
|    o unsigned long to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUL to hold the converted number
| Note:
|  - This function will only convert the first number of
|    digits in an unsigned long or till a non-numeric
|    character
\--------------------------------------------------------*/
static inline char * strToULBase10(
    char *numericStr,    /*C-string with base 10 number*/
    ulong *retUL         /*Holds converted number*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: cStrToUInt
   '  - converst a base 10 c-string to an uint32_t
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    ulong tmpUL = 0;
    *retUL = 0;

    for(
       uchar ucDig = 0;
       ucDig < ((sizeof(ulong) << 3) / 3.333333333) - 1;
      ++ucDig
    ) { /*Convert digits with no overflow concern*/
       /* Logic behind: ((sizeof(ulong)<<3)/3.333333333)+1
       `  o sizeof(ulong) << 3:
       `    - This gives the number of bits (powers of two)
       `  o / 3.3333333...:
       `    - This converts powers of 2 to powers of 10
       `  o -1:
       `    - This undercounts by two digits.
       `    - The ones digit is missed without -1.
       `    - The second digit is lost with - 1.
       */

        if(*numericStr - 58 < 0 && *numericStr - 47 > 0)
            *retUL = *retUL * 10 + *numericStr - 48;
        else return numericStr;

        ++numericStr;
    } /*Convert digits with no overflow concern*/

    /*Convert last two digits, which could overflow*/
    /*2nd to last digit*/
    if(*numericStr - 58 < 0 && *numericStr - 47 > 0)
    { /*If have one or two more numbers*/
        tmpUL = (10 * *retUL + *numericStr - 48);

        if(tmpUL < *retUL) return numericStr;

        else
        { /*Else: I can fit in one more number*/
            *retUL = tmpUL;
            ++numericStr;
        } /*Else: I can fit in one more number*/
    } /*If have one or two more numbers*/

    else return numericStr;

    /*10th digit*/
    if(*numericStr - 58 < 0 && *numericStr - 47 > 0)
    { /*If have one more number*/
        tmpUL = (10 * *retUL + *numericStr - 48);

        if(tmpUL < *retUL) return numericStr;

        else
        { /*Else: I can fit in one more number*/
            *retUL = tmpUL;
            ++numericStr;
        } /*Else: I can fit in one more number*/
    } /*If have one more number*/

    return numericStr;
} /*cStrToUInt*/

/*This is just an branchless attempt. I need to benchmark
` to see if it is faster. I suspect it will be slower. It
` has to convert the entire number, not a portion
*/
#define bitStrToULBase10(base10Str, retStr, retUL)\
   char keepBl = 0;\
   ulong tmpUL = 0;\
   (retStr) = (base10Str);\
   \
   for(\
       uchar ucDig = 0;\
       ucDig < ((sizeof(ulong) << 3) / 3.333333333) - 1;\
      ++ucDig\
   ){ /*Convert safe digits in string*/\
      /* Logic behind: ((sizeof(ulong)<<3)/3.333333333)+1
      `  o sizeof(ulong) << 3:
      `    - This gives the number of bits (powers of two)
      `  o / 3.3333333...:
      `    - This converts powers of 2 to powers of 10
      `  o -1:
      `    - This undercounts by two digits.
      `    - The ones digit is missed without -1.
      `    - The second digit is lost with - 1.
      */\
      keepBl &= (*(retStr) > 47) & (*(retStr) < 58);\
      (retUL)=((-keepBl) & (10 * (retUL) + *(retStr)-47));\
      (retStr) += keepBl;\
   } /*Convert safe digits in string*/\
   \
   /*At this point I need to worry about overflows*/\
   /*See if I have a 2nd to last digit to convert*/\
   keepBl &= (*(retStr) > 47) & (*(retStr) < 58);\
   tmpUL = ((-keepBl) & (10 * (retUL) + *(retStr) - 47));\
   keepBl &= ((retUL) < tmpUL);\
   (retUL) = ((retUL) & (keepBl-1)) | (tmpUL & (-keepBl));\
   (retStr) += keepBl;\
   \
   /*See if I have a last digit to convert*/\
   keepBl &= (*(retStr) > 47) & (*(retStr) < 58);\
   tmpUL = ((-keepBl) & (10 * (retUL) + *(retStr) - 47));\
   keepBl &= ((retUL) < tmpUL);\
   (retUL) = ((retUL) & (keepBl-1)) | (tmpUL & (-keepBl));\
   (retStr) += keepBl;\
}
  
#endif
