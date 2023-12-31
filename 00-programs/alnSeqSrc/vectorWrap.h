/*#########################################################
# Name: vectorWrap
# Use:
#  - Is a wraper for SSE2 equivlent instruction sets
#    for SSE2, AVX2, and NEON SMID instructions.
#  - This is not designed to replace better options, such
#    as SIMDe, but it is designed to provide single
#    functions that can be complied for all supported
#    instruction sets.
#  - All vector functions wrapped by macros to ensure all
#    code is inlined
#  - WARNING: Not all comparisons returen the same result.
#    o This can be avoided for bit counts by calling
#      fixIXCmpMaskCnt, which corrects the bit count.
#  - This code does not compile to scalar, but the
#    vectorBits and vectorBytes do do compile in scalar.
#    This is here to not cause errors when compiling.
# Compiling for gcc:
#  - sse: -DSSE
#  - avx2: -DAVX2 -mavx2
#  - neon: -mfpu=neon-vfpv4 or -mfpu=neon
# Supporting functions:
#  - Are here to to provide support for the vectors or do
#    things that would be difficult for vecotors. You can 
#    find these in sec-04.
#    o These are mostly for int8 vectors.
#  - These functions may be slower than doing scalar
#    operations instead.
# Datatype:
#  - vectIX is an alias for the vector data type. The X in
#    vectIX is the size of the datatype working on with the
#    vecotor. For example a vectors of int8_t's would be
#    vectI8.
#    o Currnetly datatypes for vectI8, vectI16, and
#      vectI32.
#    o 64 bit numbers are not supported, due to intel
#      having no comaparisons for SSE2.
# Definitions:
#   - vectorBits:
#     o Is the number of bits in the vector.
#     o This goes to 8 when compiled to scalar
#   - vectorBytes:
#     o Is the number of bytes a vector holds
#     o This goes to 1 when compiled to scalar
# Wrapper functions:
#  - For each function there are 8 bit, 16 bit, and 23 bit
#    variants. This is represented by X in the function
#    namnes listed here.
#  - Functions either take ro return a minVectIX or vactIX.
#  - Input/Output
#    o mmLoadIX
#      - Loads a an array into the vector
#      - Wrapper for _mmx_load_siX((__mxi *) array)
#    o mmLoadUIX
#      - mmLoad for unaligned data
#    o mmSetZeroIX
#      - Sets all values in a vectors with zeros
#      - Wrapper for _mmx_setzero_siX()
#      - Note: this is equivlent to mmSet1IX for NEON.
#        So, do not depend on it always being a latency
#        of 1, like for SSE/AVX.
#    o mmSet1IX
#      - Sets all values in a vectors to input number
#      - Wrapper for _mmx_set1_epiX()
#    o mmStoreIX
#      - Stores a vector into an array
#      - Use vectorBytes to determine the array size in
#        bytes
#      - Wrapper for mmx_store_siX((__mxi *) array, vector)
#      - This function casts input array, so it is only
#        datatype specific for the vector type, not the
#        array.
#  - Compaisons (returns mmaskx)
#    o fixIXCmpMaskCnt:
#      - These functions make sure that bit counts are
#        corrected for the differences in return
#        comparisons between comparisons.
#    o mmStoreIXCmpMask:
#      - This stores a minVectIX into an unsigned long
#      - Note, because of differences between SSE and NEON,
#        this function will have slighty different output.
#        o NEON64 returns the acutal vector full of ones
#          or zeros.
#        o NEON returns a long with multiple bits per
#          datatype (number bits in datatype / 2).
#        o SSE and AVX2 returns an long (techincally int)
#          with 1 bit for every 8 bits
#          (number of bits in datatype / 8).
#    o mmCmpGtIX
#      - Do a greater then comparison on vectors
#    o mmCmpLtIX
#      - Do a lesser then comparison on vectors
#    o mmCmpEqIX
#      - Do an is equal comparison on vectors
#  - Logical
#    o mmAndNotIX
#      - Not the first vector and then do an & (and) on the
#        first and second vector
#      - Wrapper for mmX_andnot_siX(vectOne, vectTwo)
#    o mmAndIX
#      - Perform an & (and) opteration on both vectors
#      - Wrapper for mmX_and_siX(vectOne, vectTwo)
#    o mmOrIX
#      - Perform an | (or) opteration on both vectors
#      - Wrapper for mmX_or_siX(vectOne, vectTwo)
#    o mmXOrIX
#      - Perform an ^ (xor) opteration on both vectors
#      - Wrapper for mmX_or_siX(vectOne, vectTwo)
#  - Bit mainpulation
#    o mmShiftRightIX
#      - shift a vector right by input bytes
#    o mmShiftLeftIX
#      - shift vector left by input bytes
#  - Math
#    o mmAddIX
#      - adds two vectors together
#    o mmAddSatIX
#      - Saturation add two vectors
#    o mmSubIX
#      - subtact two vectors
#    o mmSubSatIX
#      - Satuturation Subtact two vectors
#    o cnvtIXToIY
#      - Converts one vectIX type to a vectIY type.
#      - EX: convert vectI8 to vectI16
# Libraries:
# C Standard Libraries:
#   - <limits.h>
#   - <stdint.h>
#   - <immintrin.h>
#   - <arm_neon.h>
#########################################################*/

/*
 gcc SMID flags
   - SSE is default gcc (no flag needed)
   - AVX is -mavx
   - AVX2 is -mavx2
*/

#ifndef VECTORWRAP_H
#define VECTORWRAP_H

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o sec-01:
'    - Included libraries
'  o sec-02:
'    - variable defintions 
'  o sec-03:
'    - Structures (Currently not used)
'  o sec-04:
'    - Functions
'  o sec-05:
'    - Macros defing vector functions
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-01:
^  - Included libraries
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

/*This is to avoid makefile errors when compiling with
` vectorWrap, but not using vectors. This would be a
` case were the user is donig a scalar fallback
*/

#if defined AVX2 || SSE || NEON || NEON64

#ifdef DEBUG
<immintrin_dbg.h> /*For debugging*/
#endif

#ifdef AVX2                  /*256 bit registers*/
  #include <immintrin.h>
  #include "vectWrapAvx512Fun.h"

  #define vectorBits 256
  #define vectorBytes 32

  typedef __m256i vectI8;  /*vector of bytes (8bits)*/
  typedef __m256i vectI16; /*vector of shorts (16 bits)*/
  typedef __m256i vectI32; /*vector of ints (32 bits)*/

  typedef __m256i vectPrefU8;  /*Prefer unsigned bytes*/
  typedef __m256i vectPrefU16; /*Prefer unsigned shorts*/
  typedef __m256i vectPrefU32; /*Prefer unsinged ints*/

#elif SSE || SSE4           /*128 bit registers*/
  #ifdef SSE4
     #define SSE
  #endif 

  #include <immintrin.h>

  #define vectorBits 128
  #define vectorBytes 16

  typedef __m128i vectI8;  /*vector of bytes (8bits)*/
  typedef __m128i vectI16; /*vector of shorts (16 bits)*/
  typedef __m128i vectI32; /*vector of ints (32 bits)*/

  typedef __m128i vectPrefU8;  /*Prefer unsigned bytes*/
  typedef __m128i vectPrefU16; /*Prefer unsigned shorts*/
  typedef __m128i vectPrefU32; /*Prefer unsinged ints*/

#elif NEON
  #include <arm_neon.h>  /*Arm neon support*/

  #define vectorBits 128
  #define vectorBytes 16

  typedef int8x16_t vectI8;  /*vector of bytes (8bits)*/
  typedef int16x8_t vectI16; /*vector of bytes (8bits)*/
  typedef int32x4_t vectI32; /*vector of bytes (8bits)*/

  typedef uint8x16_t vectPrefU8; /*Prefer unsigned bytes*/
  typedef uint16x8_t vectPrefU16;/*Prefer unsigned shorts*/
  typedef uint32x4_t vectPrefU32;/*Prefer unsigned ints*/

#elif NEON64
  #include <arm_neon.h>  /*Arm neon support*/

  #define vectorBits 64
  #define vectorBytes 8

  typedef int8x8_t vectI8;   /*vector of bytes (8bits)*/
  typedef int16x4_t vectI16; /*vector of bytes (8bits)*/
  typedef int32x2_t vectI32; /*vector of bytes (8bits)*/

  typedef uint8x8_t vectPrefU8;  /*Prefer unsigned bytes*/
  typedef uint16x4_t vectPrefU16;/*Prefer unsigned shorts*/
  typedef uint32x2_t vectPrefU32;/*Prefer unsigned ints*/

#else
  #error No vector support (use AVX2, SSE4, or SSE)
#endif

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-02:
^  - variable defintions 
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*Probably in stdint, but must makking sure I have the
`maximum value for each uintx_t value
*/
#define defMaxUI64 0xFFFFFFFFFFFFFFFF
#define defMaxUI32 0xFFFFFFFF
#define defMaxUI16 0xFFFF
#define defMaxUI8  0xFF

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-03:
^  - Structures (Not used for this program)
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-04:
^  - Functions (Not used for this program)
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-05:
^  - Macros defing vector functions
^  o sec-05 sub-02:
^    - Macros for AVX support (256 bit vectors)
^  o sec-05 sub-03:
^    - Macros for SEE support (128 bit vectors)
*  o sec-05 sub-04:
*    - Macros for arm neon 128bit support
*  o sec-05 sub-05:
*    - Macros for arm 64bit support
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*********************************************************\
* Sec-05 Sub-02:
*  - Macros for AVX support (256 bit vectors)
*  o sec-05 sub-02 cat-01:
*    - Input/output functions (load, set, store)
*  o sec-05 Sub-02 Cat-02:
*    - Comparison functions (if)
*  o sec-05 Sub-02 Cat-03:
*    - Logical functions (and, or, xor)
*  o sec-05 Sub-02 Cat-04:
*    - Bit manipulation functions (shifts)
*  o sec-05 sub-02 cat-05:
*    - Math functions (add and sub)
*  o sec-05 sub-02 cat-06:
*    - Conversion (casting)
\*********************************************************/

#elif AVX2

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-02 Cat-01:
  +  - Input/output functions (load, set, store)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  // Aligned loading
  #define mmLoadI8(inVectI8, array) { \
    (inVectI8) = _mm256_load_si256((__m256i *) (array)); \
  }

  #define mmLoadI16(inVectI16, array) { \
    (inVectI16) = _mm256_load_si256((__m256i *) (array)); \
  }

  #define mmLoadI32(inVectI32, array) { \
    (inVectI32) = _mm256_load_si256((__m256i *) (array)); \
  }

  // Unaligned loading
  #define mmLoadUI8(inVectI8, array) { \
    (inVectI8) = _mm256_loadu_si256((__m256i *) (array)); \
  }

  #define mmLoadUI16(inVectI16, array) { \
    (inVectI16) = _mm256_loadu_si256((__m256i *) (array));\
  }

  #define mmLoadUI32(inVectI32, array) { \
    (inVectI32) = _mm256_loadu_si256((__m256i *) (array));\
  }

  // Making a zero vector
  #define mmSetZeroI8(inVectI8) {\
    (inVectI8) = _mm256_setzero_si256(); \
  }

  #define mmSetZeroI16(inVectI16) {\
    (inVectI16) = _mm256_setzero_si256(); \
  }

  #define mmSetZeroI32(inVectI32) {\
    (inVectI32) = _mm256_setzero_si256(); \
  }

  // Set a single element for all vector positions
  #define mmSet1I8(inVectI8, valC) { \
    (inVectI8) = _mm256_set1_epi8((valC)); \
  }

  #define mmSet1I16(inVectI16, valC) { \
    (inVectI16) = _mm256_set1_epi16((valC)); \
  }

  #define mmSet1I32(inVectI32, valC) { \
    (inVectI32) = _mm256_set1_epi32((valC)); \
  }

  // Store vector output into an array
  #define mmStoreI8(array, inVectI8) { \
    _mm256_store_si256((__m256i *) (array), (inVectI8)); \
  }

  #define mmStoreI16(array, inVectI16) { \
    _mm256_store_si256((__m256i *) (array), (inVectI16)); \
  }

  #define mmStoreI32(array, inVectI32) { \
    _mm256_store_si256((__m256i *) (array), (inVectI32)); \
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-02 Cat-02:
  +  - Comparison functions (if)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  // Functions for converting mask to longs
  // retUL is an unsigned long
  #define storeI8CmpMask(retUL, inVectI8) { \
    (retUL) = _mm256_movemask_epi8((inVectI8)); \
  } // Mask is from comparing 8 bit integers

  // retUL is an unsinged long
  #define storeI16CmpMask(retUL, inVectI16) { \
    (retUL) = _mm256_movemask_epi8((inVectI16)); \
  } // Mask is from comparing 16 bit integers

  // retUL is an unsinged long
  #define storeI32CmpMask(retUL, inVectI32) { \
    (retUL) = _mm256_movemask_epi8((inVectI32)); \
  } // Mask is from comparing 32 bit integers

  // Correct the bitcounts in the comparison masks
  // Intel AVX2 returns a 256 bit vector that is converted
  // to characters
  #define fixI8CmpMaskCnt(retUL,inUL) ((retUL) = (inUL))
  #define fixI16CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>1)
  #define fixI32CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>2)

  // Comparison functions
  #define mmCmpGtI8(outMaskI8, vectOneI8, vectTwoI8) { \
    (outMaskI8) = \
      _mm256_cmpgt_epi8((vectOneI8), (vectTwoI8)); \
  }

  #define mmCmpGtI16(outMaskI16, vectOneI16, vectTwoI16){ \
    (outMaskI16) = \
      _mm256_cmpgt_epi16((vectOneI16), (vectTwoI16)); \
  }

  #define mmCmpGtI32(outMaskI32, vectOneI32, vectTwoI32){ \
    (outMaskI32) = \
      _mm256_cmpgt_epi32((vectOneI32), (vectTwoI32)); \
  }

  #define mmCmpLtI8(outMaskI8, vectOneI8, vectTwoI8) { \
    (outMaskI8) = \
      _mm256_cmplt_epi8((vectOneI8), (vectTwoI8)); \
  }

  #define mmCmpLtI16(outMaskI16, vectOneI16, vectTwoI16){ \
    (outMaskI16) = \
      _mm256_cmplt_epi16((vectOneI16), (vectTwoI16)); \
  }

  #define mmCmpLtI32(outMaskI32, vectOneI32, vectTwoI32){ \
    (outMaskI32) = \
      _mm256_cmplt_epi32((vectOneI32), (vectTwoI32)); \
  }

  #define mmCmpEqI8(outMaskI8, vectOneI8, vectTwoI8) { \
    (outMaskI8) = \
       _mm256_cmpeq_epi8((vectOneI8), (vectTwoI8)); \
  }

  #define mmCmpEqI16(outMaskI16, vectOneI16, vectTwoI16){ \
    (outMaskI16) =\
      _mm256_cmpeq_epi16((vectOneI16), (vectTwoI16)); \
  }

  #define mmCmpEqI32(outMaskI32, vectOneI32, vectTwoI32){ \
    (outMaskI32) = \
      _mm256_cmpeq_epi32((vectOneI32), (vectTwoI32)); \
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-02 Cat-03:
  +  - Logical functions (and, or, xor)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  // and not function. The extra functions here are to
  // support the NEON functions
  #define mmAndNotI8(outVectI8, notVectI8, inVectI8) { \
    (outVectI8) = \
      _mm256_andnot_si256((notVectI8), (inVectI8)); \
  }

  #define mmAndNotI16(outVectI16, notVectI16, inVectI16){ \
    (outVectI16) = \
      _mm256_andnot_si256((notVectI16), (inVectI16)); \
  }

  #define mmAndNotI32(outVectI32, notVectI32, inVectI32){ \
    (outVectI32) = \
      _mm256_andnot_si256((notVectI32), (inVectI32)); \
  }

  // And functions (extra functions for NEON support)
  #define mmAndI8(outVectI8, vectOneI8, vectTwoI8) { \
    (outVectI8) = \
      _mm256_and_si256((vectOneI8), (vectTwoI8)); \
  }

  #define mmAndI16(outVectI16, vectOneI16, vectTwoI16) { \
    (outVectI16) = \
      _mm256_and_si256((vectOneI16), (vectTwoI16)); \
  }

  #define mmAndI32(outVectI32, vectOneI32, vectTwoI32) { \
    (outVectI32) = \
      _mm256_and_si256((vectOneI32), (vectTwoI32)); \
  }

  // Or functions
  #define mmOrI8(outVectI8, vectOneI8, vectTwoI8) { \
    (outVectI8) = \
      _mm256_or_si256((vectOneI8),(vectTwoI8)); \
  }

  #define mmOrI16(outVectI16, vectOneI16, vectTwoI16) { \
    (outVectI16) = \
      _mm256_or_si256((vectOneI16), (vectTwoI16)); \
  }

  #define mmOrI32(outVectI32, vectOneI32, vectTwoI32) { \
    (outVectI32) = \
      _mm256_or_si256((vectOneI32), (vectTwoI32)); \
  }

  // Or functions
  #define mmXOrI8(outVectI8, vectOneI8, vectTwoI8) { \
    (outVectI8) = \
      _mm256_xor_si256((vectOneI8), (vectTwoI8)); \
  }

  #define mmXOrI16(outVectI16, vectOneI16, vectTwoI16) { \
    (outVectI16) = \
      _mm256_xor_si256((vectOneI16), (vectTwoI16));\
  }

  #define mmXOrI32(outVectI32, vectOneI32, vectTwoI32) { \
    (outVectI32) = \
      _mm256_xor_si256((vectOneI32), (vectTwoI32)); \
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-02 Cat-04:
  +  - Bit manipulation functions (shifts)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  #define mmShiftRightI8(outVectI8, inVectI8, numBytesI){\
    (outVectI8) = \
       _mm256_srli_si256((inVectI8), (numBytesI)); \
  }

  // the << 1 is to account for two bytes per short
  // This is to ensure that this 
  #define mmShiftRightI16(outVectI16,inVectI16,numShorts){\
    (outVectI16) = \
      _mm256_srli_si256((inVectI16), (numShorts) <<1);\
  }

  // The << 2 is to account for four bytes per int
  #define mmShiftRightI32(outVectI32, inVectI32, numInts){\
    (outVectI) = \
      _mm256_srli_si256((inVectI32), (numInts) << 2);\
  }

  // Shift left
  #define mmShiftLeftI8(outVectI8, inVectI8, numBytesI){\
    (outVectI8) = \
      _mm256_slli_si256((inVectI8), (numBytesI)); \
  }

  // the << 1 is to account for two bytes per short
  // This is to ensure that this 
  #define mmShiftLeftI16(outVectI16,inVectI16,numShorts){\
    (outVectI16) =\
      _mm256_slli_si256((inVectI16), (numShorts) <<1);\
  }

  // The << 2 is to account for four bytes per int
  #define mmShiftLeftI32(outVectI32, inVectI32, numInts){\
    (outVectI) = \
      _mm256_slli_si256((inVectI32), (numInts) << 2);\
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-02 Cat-05:
  +  - Math functions (add and sub)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   #define mmAddI8(outVectI8, vectOneI8, vectTwoI8) {\
     (outVectI8)=_mm256_add_epi8((vectOneI8),(vectTwoI8));\
   }

   #define mmAddI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = \
       _mm256_add_epi16((vectOneI16), (vectTwoI16)); \
   }

   #define mmAddI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = \
       _mm256_add_epi32((vectOneI32), (vectTwoI32)); \
   }

   // Saturation addition
   #define mmAddSatI8(outVectI8, vectOneI8, vectTwoI8) {\
     (outVectI8) = \
       _mm256_adds_epi8((vectOneI8), (vectTwoI8)); \
   }

   #define mmAddSatI16(outVectI16,vectOneI16,vectTwoI16){ \
     (outVectI16) = \
       _mm256_adds_epi16((vectOneI16), (vectTwoI16)); \
   }

   #define mmAddSatI32(outVectI32,vectOneI32,vectTwoI32){ \
     (outVectI32) = \
       _mm256_adds_epi32((vectOneI32), (vectTwoI32)); \
   }

   // Subtraction

   #define mmSubI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = \
        _mm256_sub_epi8((vectOneI8), (vectTwoI8)); \
   }

   #define mmSubI16(outVectI16, vectOneI16, vectTwoI16){ \
     (outVectI16) = \
       _mm256_sub_epi16((vectOneI16), (vectTwoI16)); \
   }

   #define mmSubI32(outVectI32, vectOneI32, vectTwoI32){ \
     (outVectI32) = \
       _mm256_sub_epi32((vectOneI32), (vectTwoI32)); \
   }

   // Satuturation Subtraction
   #define mmSubSatI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = \
       _mm256_subs_epi8((vectOneI8), (vectTwoI8)); \
   }

   #define mmSubSatI16(outVectI16,vectOneI16,vectTwoI16){ \
     (outVectI16) = \
       _mm256_subs_epi16((vectOneI16), (vectTwoI16)); \
   }

   #define mmSubSatI32(outVectI32,vectOneI32,vectTwoI32){ \
     (outVectI32) = \
       _mm256_subs_epi32((vectOneI32), (vectTwoI32)); \
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-02 Cat-06:
   +  - Conversion (casting)
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    // These are here to support NEON

   #define cnvtI8ToI16(outVectI16, inVectI8) { \
     (outVectI16) = (inVectI8); \
   }

   #define cnvtI8ToI32(outVectI32, inVectI8) { \
     (outVectI32) = (inVectI8); \
   }

   #define cnvtI16ToI8(outVectI8, inVectI16) { \
     (outVectI8) = (inVectI16); \
   }

   #define cnvtI16ToI32(outVectI32, inVectI16) { \
     (outVectI32) = (inVectI16); \
   }

   #define cnvtI16ToI32(outVectI32, inVectI16) { \
     (outVectI32) = (inVectI16); \
   }

   #define cnvtI32ToI8(outVectI8, inVectI32) { \
     (outVectI8) = (inVectI32); \
   }

   #define cnvtI32ToI16(outVectI16, inVectI32) { \
     (outVectI16) = (inVectI32); \
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-04 Sub-02 Cat-07:
   +  - Max
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /* Byte maximize
   `  SSE2 does not support I8 max
   */
   #define mmMaxPrefI8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = \
           _mm256_max_epi8((firstVectI8), (secVectI8)); \
   }

   #define mmMaxU8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = \
        _mm256_max_epu8((firstVectI8), (secVectI8)); \
   }

   /*16 bit maximizes*/
   #define mmMaxI16(outVectI16, firstVectI16, secVectI16){\
      (outVectI16) = \
         _mm256_max_epi16((firstVectI16), (secVectI16)); \
   }

   #define mmMaxPrefU16( \
         outVectU16, \
         firstVectU16, \
         secVectU16 \
      ){\
        (outVectU16) = \
           _mm256_max_epu16((firstVectU16), (secVectU16)); \
   }

   /*32 unsigned bit maximizes (if possible)*/
   #define mmMaxPrefU32( \
      outVectU32, \
      firstVectU32, \
      secVectU32 \
   ){\
      (outVectU32) = \
         _mm256_max_epu32((firstVectU32), (secVectU32)); \
   }

   /*32 signed bit maximizes (if possible)*/
   #define mmMaxPrefI32( \
      outVectI32, \
      firstVectI32, \
      secVectI32 \
   ){\
      (outVectI32) = \
         _mm256_max_epi32((firstVectI32), (secVectI32)); \
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-02 Cat-08:
   +  - Min
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /* Byte minimize
   `  SSE2 does not support I8 min
   */
   #define mmMinPrefI8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = \
         _mm256_min_epi8((firstVectI8), (secVectI8)); \
   }

   #define mmMinU8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = \
        _mm256_min_epu8((firstVectI8), (secVectI8)); \
   }

   /*16 bit minimize*/
   #define mmMinI16(outVectI16, firstVectI16, secVectI16){\
      (outVectI16) = \
         _mm256_min_epi16((firstVectI16), (secVectI16)); \
   }

   #define mmMinPrefU16( \
      outVectU16, \
      firstVectU16, \
      secVectU16 \
   ){\
      (outVectU16) = \
         _mm256_min_epu16((firstVectU16), (secVectU16)); \
   }

   /*32 unsigned bit minimize (if possible)*/
   #define mmMinPrefU32( \
      outVectU32, \
      firstVectU32, \
      secVectU32 \
   ){\
     (outVectU32) = \
        _mm256_min_epu32((firstVectU32), (secVectU32)); \
   }

   /*32 signed bit minimize (if possible)*/
   #define mmMinPrefI32( \
      outVectI32, \
      firstVectI32, \
      secVectI32 \
   ){\
     (outVectI32) = \
        _mm256_min_epi32((firstVectI32), (secVectI32)); \
   }

/*********************************************************\
* Sec-05 Sub-03:
*  - Macros for SSE support (128 bit)
*  o sec-05 sub-03 cat-01:
*    - Input/output functions (load, set, store)
*  o sec-05 sub-03 cat-02:
*    - Comparison functions (if)
*  o sec-05 sub-03 cat-03:
*    - Logical functions (and, or, xor)
*  o sec-05 sub-03 cat-04:
*    - Bit manipulation functions (shifts)
*  o sec-05 sub-03 cat-05:
*    - Math functions (add and sub)
*  o sec-05 sub-03 cat-06:
*    - Conversion (casting)
*  o sec-05 sub-03 cat-07:
*    - Max and min
\*********************************************************/

#elif SSE

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-03 Cat-01:
  +  - Input/output functions (load, set, store)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  /*Aligned loading*/
  #define mmLoadI8(inVectI8, array) { \
    (inVectI8) = _mm_load_si128((__m128i *) (array)); \
  }

  #define mmLoadI16(inVectI16, array) { \
    (inVectI16) = _mm_load_si128((__m128i *) (array)); \
  }

  #define mmLoadI32(inVectI32, array) { \
    (inVectI32) = _mm_load_si128((__m128i *) (array)); \
  }

  /*Unaligned loading*/
  #define mmLoadUI8(inVectI8, array) { \
    (inVectI8) = _mm_loadu_si128((__m128i *) (array)); \
  }

  #define mmLoadUI16(inVectI16, array) { \
    (inVectI16) = _mm_loadu_si128((__m128i *) (array)); \
  }

  #define mmLoadUI32(inVectI32, array) { \
    (inVectI32) = _mm_loadu_si128((__m128i *) (array)); \
  }

  /*Making a zero vector*/
  #define mmSetZeroI8(inVectI8) {\
    (inVectI8) = _mm_setzero_si128(); \
  }

  #define mmSetZeroI16(inVectI16) {\
    (inVectI16) = _mm_setzero_si128(); \
  }

  #define mmSetZeroI32(inVectI32) {\
    (inVectI32) = _mm_setzero_si128(); \
  }

  // Set a single element for all vector positions
  #define mmSet1I8(inVectI8, valC) { \
    (inVectI8) = _mm_set1_epi8((valC)); \
  }

  #define mmSet1I16(inVectI16, valC) { \
    (inVectI16) = _mm_set1_epi16((valC)); \
  }

  #define mmSet1I32(inVectI32, valC) { \
    (inVectI32) = _mm_set1_epi32((valC)); \
  }

  // Store vector output into an array
  #define mmStoreI8(array, inVectI8) { \
    _mm_store_si128((__m128i *) (array), (inVectI8)); \
  }

  #define mmStoreI16(array, inVectI16) { \
    _mm_store_si128((__m128i *) (array), (inVectI16)); \
  }

  #define mmStoreI32(array, inVectI32) { \
    _mm_store_si128((__m128i *) (array), (inVectI32)); \
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-03 Cat-02:
  +  - Comparison functions (if)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  // Functions for converting mask to longs
  // retUL is an unsigned long
  #define storeI8CmpMask(retUL, inVectI8) { \
    (retUL) = _mm_movemask_epi8((inVectI8)); \
  } // Mask is from comparing 8 bit integers

  // retUL is an unsinged long
  #define storeI16CmpMask(retUL, inVectI16) { \
    (retUL) = _mm_movemask_epi8((inVectI16)); \
  } // Mask is from comparing 16 bit integers

  // retUL is an unsinged long
  #define storeI32CmpMask(retUL, inVectI32) { \
    (retUL) = _mm_movemask_epi8((inVectI32)); \
  } // Mask is from comparing 32 bit integers

  // Correct the bitcounts in the comparison masks
  // Intel SSE returns a 128 bit vector that is converted
  // to characters
  #define fixI8CmpMaskCnt(retUL,inUL) ((retUL) = (inUL))
  #define fixI16CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>1)
  #define fixI32CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>2)

  // Comparison functions
  #define mmCmpGtI8(outMaskI8, vectOneI8, vectTwoI8) { \
    (outMaskI8) = _mm_cmpgt_epi8((vectOneI8),(vectTwoI8));\
  }

  #define mmCmpGtI16(outMaskI16, vectOneI16, vectTwoI16){ \
    (outMaskI16) = \
      _mm_cmpgt_epi16((vectOneI16),(vectTwoI16));\
  }

  #define mmCmpGtI32(outMaskI32, vectOneI32, vectTwoI32){ \
    (outMaskI32) = \
      _mm_cmpgt_epi32((vectOneI32),(vectTwoI32));\
  }

  #define mmCmpLtI8(outMaskI8, vectOneI8, vectTwoI8) { \
    (outMaskI8) = _mm_cmplt_epi8((vectOneI8),(vectTwoI8));\
  }

  #define mmCmpLtI16(outMaskI16, vectOneI16, vectTwoI16){ \
    (outMaskI16) = \
      _mm_cmplt_epi16((vectOneI16),(vectTwoI16));\
  }

  #define mmCmpLtI32(outMaskI32, vectOneI32, vectTwoI32){ \
    (outMaskI32) = \
      _mm_cmplt_epi32((vectOneI32),(vectTwoI32));\
  }

  #define mmCmpEqI8(outMaskI8, vectOneI8, vectTwoI8) { \
    (outMaskI8) = _mm_cmpeq_epi8((vectOneI8),(vectTwoI8));\
  }

  #define mmCmpEqI16(outMaskI16, vectOneI16, vectTwoI16){ \
    (outMaskI16) = \
      _mm_cmpeq_epi16((vectOneI16),(vectTwoI16));\
  }

  #define mmCmpEqI32(outMaskI32, vectOneI32, vectTwoI32){ \
    (outMaskI32) = \
      _mm_cmpeq_epi32((vectOneI32),(vectTwoI32));\
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-03 Cat-03:
  +  - Logical functions (and, or, xor)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  // and not function. The extra functions here are to
  // support the NEON functions
  #define mmAndNotI8(outVectI8, notVectI8, inVectI8) { \
    (outVectI8) =_mm_andnot_si128((notVectI8),(inVectI8));\
  }

  #define mmAndNotI16(outVectI16, notVectI16, inVectI16){ \
    (outVectI16) = \
      _mm_andnot_si128((notVectI16), (inVectI16));\
  }

  #define mmAndNotI32(outVectI32, notVectI32, inVectI32){ \
    (outVectI32) = \
      _mm_andnot_si128((notVectI32), (inVectI32));\
  }

  // And functions (extra functions for NEON support)
  #define mmAndI8(outVectI8, vectOneI8, vectTwoI8) { \
    (outVectI8) = _mm_and_si128((vectOneI8), (vectTwoI8));\
  }

  #define mmAndI16(outVectI16, vectOneI16, vectTwoI16) { \
    (outVectI16) = \
      _mm_and_si128((vectOneI16), (vectTwoI16)); \
  }

  #define mmAndI32(outVectI32, vectOneI32, vectTwoI32) { \
    (outVectI32) = \
      _mm_and_si128((vectOneI32), (vectTwoI32)); \
  }

  // Or functions
  #define mmOrI8(outVectI8, vectOneI8, vectTwoI8) { \
    (outVectI8) = _mm_or_si128((vectOneI8), (vectTwoI8)); \
  }

  #define mmOrI16(outVectI16, vectOneI16, vectTwoI16) { \
    (outVectI16) =_mm_or_si128((vectOneI16),(vectTwoI16));\
  }

  #define mmOrI32(outVectI32, vectOneI32, vectTwoI32) { \
    (outVectI32)= _mm_or_si128((vectOneI32),(vectTwoI32));\
  }

  // Or functions
  #define mmXOrI8(outVectI8, vectOneI8, vectTwoI8) { \
    (outVectI8)= _mm_xor_si128((vectOneI8), (vectTwoI8));\
  }

  #define mmXOrI16(outVectI16, vectOneI16, vectTwoI16) { \
    (outVectI16)=_mm_xor_si128((vectOneI16),(vectTwoI16));\
  }

  #define mmXOrI32(outVectI32, vectOneI32, vectTwoI32) { \
    (outVectI32)=_mm_xor_si128((vectOneI32),(vectTwoI32));\
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-03 Cat-04:
  +  - Bit manipulation functions (shifts)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  #define mmShiftRightI8(outVectI8, inVectI8, numBytesI){\
    (outVectI8)=_mm_srli_si128((inVectI8), (numBytesI));\
  }

  // the << 1 is to account for two bytes per short
  // This is to ensure that this 
  #define mmShiftRightI16(outVectI16,inVectI16,numShorts){\
    (outVectI16) =\
      _mm_srli_si128((inVectI16), (numShorts) <<1);\
  }

  // The << 2 is to account for four bytes per int
  #define mmShiftRightI32(outVectI32, inVectI32, numInts){\
    (outVectI32) = \
      _mm_srli_si128((inVectI32), (numInts) << 2);\
  }

  // Shift left
  #define mmShiftLeftI8(outVectI8, inVectI8, numBytesI){\
    (outVectI8) = _mm_slli_si128((inVectI8), (numBytesI));\
  }

  // the << 1 is to account for two bytes per short
  // This is to ensure that this 
  #define mmShiftLeftI16(outVectI16,inVectI16,numShorts){\
    (outVectI16) =\
      _mm_slli_si128((inVectI16), (numShorts) <<1);\
  }

  // The << 2 is to account for four bytes per int
  #define mmShiftLeftI32(outVectI32, inVectI32, numInts){\
    (outVectI) = \
      _mm_slli_si128((inVectI32), (numInts) << 2);\
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-03 Cat-05:
  +  - Math functions (add and sub)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   #define mmAddI8(outVectI8, vectOneI8, vectTwoI8) {\
     (outVectI8) = _mm_add_epi8((vectOneI8), (vectTwoI8));\
   }

   #define mmAddI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = \
       _mm_add_epi16((vectOneI16), (vectTwoI16)); \
   }

   #define mmAddI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = \
       _mm_add_epi32((vectOneI32), (vectTwoI32)); \
   }

   // Saturation addition
   #define mmAddSatI8(outVectI8, vectOneI8, vectTwoI8) {\
     (outVectI8) = _mm_adds_epi8((vectOneI8),(vectTwoI8));\
   }

   #define mmAddSatI16(outVectI16,vectOneI16,vectTwoI16){ \
     (outVectI16) = \
       _mm_adds_epi16((vectOneI16), (vectTwoI16)); \
   }

   #define mmAddSatI32(outVectI32,vectOneI32,vectTwoI32){ \
     (outVectI32) = \
       _mm_adds_epi32((vectOneI32), (vectTwoI32)); \
   }

   // Subtraction

   #define mmSubI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = _mm_sub_epi8((vectOneI8), (vectTwoI8));\
   }

   #define mmSubI16(outVectI16, vectOneI16, vectTwoI16){ \
     (outVectI16) = \
       _mm_sub_epi16((vectOneI16), (vectTwoI16)); \
   }

   #define mmSubI32(outVectI32, vectOneI32, vectTwoI32){ \
     (outVectI32) = \
       _mm_sub_epi32((vectOneI32), (vectTwoI32)); \
   }

   // Satuturation Subtraction
   #define mmSubSatI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = _mm_subs_epi8((vectOneI8),(vectTwoI8));\
   }

   #define mmSubSatI16(outVectI16,vectOneI16,vectTwoI16){ \
     (outVectI16) = \
       _mm_subs_epi16((vectOneI16), (vectTwoI16)); \
   }

   #define mmSubSatI32(outVectI32,vectOneI32,vectTwoI32){ \
     (outVectI32) = \
       _mm_subs_epi32((vectOneI32), (vectTwoI32)); \
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-03 Cat-06:
   +  - Conversion (casting)
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    /* These are here to support NEON*/

   #define cnvtI8ToI16(outVectI16, inVectI8) { \
     (outVectI16) = (inVectI8); \
   }

   #define cnvtI8ToI32(outVectI32, inVectI8) { \
     (outVectI32) = (inVectI8); \
   }

   #define cnvtI16ToI8(outVectI8, inVectI16) { \
     (outVectI8) = (inVectI16); \
   }

   #define cnvtI16ToI32(outVectI32, inVectI16) { \
     (outVectI32) = (inVectI16); \
   }

   #define cnvtI16ToI32(outVectI32, inVectI16) { \
     (outVectI32) = (inVectI16); \
   }

   #define cnvtI32ToI8(outVectI8, inVectI32) { \
     (outVectI8) = (inVectI32); \
   }

   #define cnvtI32ToI16(outVectI16, inVectI32) { \
     (outVectI16) = (inVectI32); \
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-03 Cat-07:
   +  - Max
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /* Byte maximize
   `  SSE2 does not support I8 max
   */
   #ifdef SSE4
     #define mmMaxPrefI8(outVectI8,firstVectI8,secVectI8){\
        (outVectI8) = \
           _mm_max_epu8((firstVectI8), (secVectI8)); \
     }
   #else
     #define mmMaxPrefI8(outVectI8,firstVectI8,secVectI8){\
        (outVectI8) = \
           _mm_max_epu8((firstVectI8), (secVectI8)); \
     }
   #endif

   #define mmMaxU8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = \
        _mm_max_epu8((firstVectI8), (secVectI8)); \
   }

   /*16 bit maximizes*/
   #define mmMaxI16(outVectI16, firstVectI16, secVectI16){\
      (outVectI16) = \
         _mm_max_epi16((firstVectI16), (secVectI16)); \
   }

   /*16 bit maximizes*/
   #ifdef SSE4
      #define mmMaxPrefU16( \
         outVectU16, \
         firstVectU16, \
         secVectU16 \
      ){\
        (outVectU16) = \
           _mm_max_epu16((firstVectU16), (secVectU16)); \
      }
   #else
      #define mmMaxPrefU16( \
         outVectU16, \
         firstVectU16, \
         secVectU16 \
      ){\
        (outVectU16) = \
           _mm_max_epi16((firstVectU16), (secVectU16)); \
      }
   #endif

   /*32 unsigned bit maximizes (if possible)*/
   #ifdef SSE4
      #define mmMaxPrefU32( \
         outVectU32, \
         firstVectU32, \
         secVectU32 \
      ){\
        (outVectU32) = \
           _mm_max_epu32((firstVectU32), (secVectU32)); \
      }
   #else
      #define mmMaxPrefU32( \
         outVectU32, \
         firstVectU32, \
         secVectU32 \
      ){\
        (outVectU32) = \
           _mm_max_epi16((firstVectU32), (secVectU32)); \
      }
   #endif

   /*32 signed bit maximizes (if possible)*/
   #ifdef SSE4
      #define mmMaxPrefI32( \
         outVectI32, \
         firstVectI32, \
         secVectI32 \
      ){\
        (outVectI32) = \
           _mm_max_epi32((firstVectI32), (secVectI32)); \
      }
   #else
      #define mmMaxPrefI32( \
         outVectI32, \
         firstVectI32, \
         secVectI32 \
      ){\
        (outVectI32) = \
           _mm_max_epi16((firstVectI32), (secVectI32)); \
      }
   #endif

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-03 Cat-08:
   +  - Min
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /* Byte minimize
   `  SSE2 does not support I8 min
   */
   #ifdef SSE4
     #define mmMinPrefI8(outVectI8,firstVectI8,secVectI8){\
        (outVectI8) = \
           _mm_min_epu8((firstVectI8), (secVectI8)); \
     }
   #else
     #define mmMinPrefI8(outVectI8,firstVectI8,secVectI8){\
        (outVectI8) = \
           _mm_min_epu8((firstVectI8), (secVectI8)); \
     }
   #endif

   #define mmMinU8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = \
        _mm_min_epu8((firstVectI8), (secVectI8)); \
   }

   /*16 bit minimize*/
   #define mmMinI16(outVectI16, firstVectI16, secVectI16){\
      (outVectI16) = \
         _mm_min_epi16((firstVectI16), (secVectI16)); \
   }

   /*16 bit minimize*/
   #ifdef SSE4
      #define mmMinPrefU16( \
         outVectU16, \
         firstVectU16, \
         secVectU16 \
      ){\
        (outVectU16) = \
           _mm_min_epu16((firstVectU16), (secVectU16)); \
      }
   #else
      #define mmMinPrefU16( \
         outVectU16, \
         firstVectU16, \
         secVectU16 \
      ){\
        (outVectU16) = \
           _mm_min_epi16((firstVectU16), (secVectU16)); \
      }
   #endif

   /*32 unsigned bit minimize (if possible)*/
   #ifdef SSE4
      #define mmMinPrefU32( \
         outVectU32, \
         firstVectU32, \
         secVectU32 \
      ){\
        (outVectU32) = \
           _mm_min_epu32((firstVectU32), (secVectU32)); \
      }
   #else
      #define mmMinPrefU32( \
         outVectU32, \
         firstVectU32, \
         secVectU32 \
      ){\
        (outVectU32) = \
           _mm_min_epi16((firstVectU32), (secVectU32)); \
      }
   #endif

   /*32 signed bit minimize (if possible)*/
   #ifdef SSE4
      #define mmMinPrefI32( \
         outVectI32, \
         firstVectI32, \
         secVectI32 \
      ){\
        (outVectI32) = \
           _mm_min_epi32((firstVectI32), (secVectI32)); \
      }
   #else
      #define mmMinPrefI32( \
         outVectI32, \
         firstVectI32, \
         secVectI32 \
      ){\
        (outVectI32) = \
           _mm_min_epi16((firstVectI32), (secVectI32)); \
      }
   #endif

/*********************************************************\
* Sec-05 Sub-04:
*  - Macros for arm neon support (NEED TO FILL)
*  o sec-04 sub-04 cat-01:
*    - Input/output functions (load, set, store)
*  o sec-04 sub-04 cat-02:
*    - Comparison functions (if)
*  o sec-04 sub-04 cat-03:
*    - Logical functions (and, or, xor)
*  o sec-04 sub-04 cat-04:
*    - Bit manipulation functions (shifts)
*  o sec-04 sub-04 cat-05:
*    - Addition and subtraction
*  o sec-04 sub-04 cat-06:
*    - Conversion (casting)
\*********************************************************/

/*Insert
   inVectI16 = vsetq_lane_s16(valI16, inVectI16, posI)
*/

#elif NEON
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-04 Cat-01:
  +  - Input/output functions (load, set, store)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   #define mmLoadI8(outVectI8, arrayI8) { \
    (outVectI8) = vld1q_s8((int8_t *) (arrayI8));\
   }

   #define mmLoadI16(outVectI16, arrayI16) { \
    (outVectI8) = vld1q_s16((int16_t *) (arrayI16));\
   }

   #define mmLoadI32(outVectI32, arrayI32) { \
    (outVectI8) = vld1q_s32((int32_t *) (arrayI32));\
   }

   // Unaliged load options.
   // This is here to support intel, which does have this
   // function
   #define mmLoadUI8(outVectI8, arrayI8) { \
    (outVectI8) = vld1q_s8((int8_t *) (arrayI8));\
   }

   #define mmLoadUI16(outVectI16, arrayI16) { \
    (outVectI16) = vld1q_s16((int16_t *) (arrayI16));\
   }

   #define mmLoadUI32(outVectI32, arrayI32) { \
    (outVectI32) = vld1q_s32((int32_t *) (arrayI32));\
   }

   // single value loads
   #define mmSet1I8(inVectI8, valC) { \
     (inVectI8) = vdupq_n_s8((valC)); \
   }

   #define mmSet1I16(inVectI16, valC) { \
     (inVectI16) = vdupq_n_s16((valC)); \
   }

   #define mmSet1I32(inVectI32, valC) { \
     (inVectI32) = vdupq_n_s32((valC)); \
   }

   // Set zero values (here to support intel)
   #define mmSeqZeroI8(inVectI8) { \
     (inVectI8) = vsubq_s8(0); \
   }

   #define mmSeqZeroI16(inVectI16) { \
     (inVectI16) = vsubq_s16(0); \
   }

   #define mmSeqZeroI32(inVectI32) { \
     (inVectI32) = vsubq_s32(0); \
   }

  // Output functions; storing vectors in arrays
  #define mmStoreI8(array, inVectI8) { \
    vstlq_s8((int8_t *) (array), (inVectI8)); \
  }

  #define mmStoreI16(array, inVectI16) { \
    vstlq_s16((int16_t *) (array), (inVectI16)); \
  }

  #define mmStoreI32(array, inVectI32) { \
    vstlq_s32((int32_t *) (array), (inVectI32)); \
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-04 Cat-04:
  +  - Comparison functions (if)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// Store the mask functions (unsigned long returns)
   // The logic for these store functions came from
   // https://community.arm.com/arm-community-blogs/b/infrastructure-solutions-blog/posts/porting-x86-vector-bitmask-optimizations-to-arm-neon
   #define storeI8CmpMask(retUL, inVectI8) { \
     uint16x8_t res16x8 = vreinterpretq_u16_s8((inVectI8)); \
     uint8x8_t res8x8 = vshrn_n_u16(res16x8, 4); \
     uint64x1_t res64x1 = vreinterpret_u64_u8(res8x8); \
     (retUL) = vget_lane_u64(res64x1, 0); \
   }

   #define storeI16CmpMask(retUL, inVectI16) { \
     uint32x4_t res32x4 =vreinterpretq_u32_s16((inVectI16));\
     uint16x4_t res16x4 = vshrn_n_u32(res32x4, 4); \
     uint64x1_t res64x1 = vreinterpret_u64_u16(res16x4); \
     (retUL) = vget_lane_u64(res64x1, 0); \
   }

   #define storeI32CmpMask(retUL, inVectI32) { \
     uint64x2_t res64x2 =vreinterpretq_u64_s32((inVectI32));\
     uint32x2_t res32x2 = vshrn_n_u64(tmp64x2, 4); \
     uint64x1_t res64x1 = vreinterpret_u64_u32(tmp32x2); \
     (retUL) = vget_lane_u64(res64x1, 0); \
   }

   // Correct the bitcounts in the comparison masks
   #define fixI8CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>2)
   #define fixI16CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>3)
   #define fixI32CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>4)

   // greater then comparisions
   #define mmCmpGtI8(outMaskI8, vectOneI8, vectTwoI8) { \
     (outMaskI8) = vcgtq_s8((vectOneI8), (vectTwoI8)); \
   };

   #define mmCmpGtI16(outMaskI16,vectOneI16,vectTwoI16) { \
     (outMaskI16) = vcgtq_s16((vectOneI16), (vectTwoI16));\
   };

   #define mmCmpGtI32(outMaskI32,vectOneI32,vectTwoI32) { \
     (outMaskI32) = vcgtq_s32((vectOneI32), (vectTwoI32));\
   };

   // Less then comparisions
   #define mmCmpLtI8(outMaskI8, vectOneI8, vectTwoI8) { \
     (outMaskI8) = vcltq_s8((vectOneI8), (vectTwoI8)); \
   };

   #define mmCmpLtI16(outMaskI16,vectOneI16,vectTwoI16) { \
     (outMaskI16) = vcltq_s16((vectOneI16), (vectTwoI16));\
   };

   #define mmCmpLtI32(outMaskI32,vectOneI32,vectTwoI32) { \
     (outMaskI32) = vcltq_s32((vectOneI32), (vectTwoI32));\
   };

   // Less then comparisions
   #define mmCmpEqI8(outMaskI8, vectOneI8, vectTwoI8) { \
     (outMaskI8) = vceqq_s8((vectOneI8), (vectTwoI8)); \
   };

   #define mmCmpEqI16(outMaskI16,vectOneI16,vectTwoI16) { \
     (outMaskI16) = vceqq_s16((vectOneI16), (vectTwoI16));\
   };

   #define mmCmpEqI32(outMaskI32,vectOneI32,vectTwoI32) { \
     (outMaskI32) = vceqq_s32((vectOneI32), (vectTwoI32));\
   };

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-04 Cat-03:
  +  - Logical functions (and, or, xor)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   #define mmAndNotI8(outVectI8, notVectI8, inVect8) { \
     (outVectI8) = \
       vandq_s8(vmvnq_s8((notVectI8)), (inVect8)); \
   }

   #define mmAndNotI16(outVectI16, notVectI16, inVect16){ \
     (outVectI16) = \
       vandq_s16(vmvnq_s16((notVectI16)), (inVect16)); \
   }

   #define mmAndNotI32(outVectI32, notVectI32, inVect32){ \
     (outVectI32) = \
       vandq_s32(vmvnq_s32((notVectI32)), (inVect32)); \
   }

   #define mmAndI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vandq_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmAndI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = vandq_s16((vectOneI16), (vectTwoI16));\
   }

   #define mmAndI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = vandq_s32((vectOneI32), (vectTwoI32));\
   }

   #define mmOrI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vorrq_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmOrI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI160 = vorrq_s16((vectOneI16), (vectTwoI16));\
   }

   #define mmOrI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = vorrq_s32((vectOneI32), (vectTwoI32));\
   }

   #define mmXorI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = veorq_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmXorI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = veorq_s16((vectOneI16), (vectTwoI16));\
   }

   #define mmXorI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = veorq_s32((vectOneI32), (vectTwoI32));\
   }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-04 Cat-04:
  +  - Bit manipulation functions (shifts)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  #define mmShiftRightI8(outVectI8, inVectI8, numBytesI){\
    (outVectI8) = vshrq_s8((inVectI8), (numBytesI)); \
  }

  #define mmShiftRightI16(outVectI16,inVectI16,numBytesI){\
    (outVectI16) = vshrq_s16((inVectI16), (numBytesI)); \
  }

  #define mmShiftRightI32(outVectI32,inVectI32,numBytesI){\
    (outVectI32) = vshrq_s32((inVectI32), (numBytesI)); \
  }

  #define mmShiftLeftI8(outVectI8, inVectI8, numBytesI){\
    (outVectI8) = vshlq_s8((inVectI8), (numBytesI)); \
  }

  #define mmShiftLeftI16(outVectI16,inVectI16,numBytesI){\
    (outVectI16) = vshlq_s16((inVectI16), (numBytesI)); \
  }

  #define mmShiftLeftI32(outVectI32,inVectI32,numBytesI){\
    (outVectI32) = vshlq_s32((inVectI32), (numBytesI)); \
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-04 Cat-05:
  +  - Math functions (add and sub)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   // Addition (non-saturating)
   #define mmAddI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vaddq_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmAddI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = vaddq_s16((vectOneI16), (vectTwoI16));\
   }

   #define mmAddI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = vaddq_s32((vectOneI32),(vectTwoI32)); \
   }

   // Addition (saturating)
   #define mmAddSatI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vqaddq_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmAddSatI16(outVectI16,vectOneI16,vectTwoI16){ \
     (outVectI16) = vqaddq_s16((vectOneI16),(vectTwoI16));\
   }

   #define mmAddSatI32(outVectI32,vectOneI32,vectTwoI32){ \
     (outVectI32) = vqaddq_s32((vectOneI32),(vectTwoI32));\
   }

   // subtraction (non-saturating)
   #define mmsubI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vsubq_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmsubI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = vsubq_s16((vectOneI16),(vectTwoI16)); \
   }

   #define mmsubI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = vsubq_s32((vectOneI32),(vectTwoI32)); \
   }

   // subtraction (saturating)
   #define mmSubSatI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vqsubq_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmSubSatI16(outVectI16,vectOneI16,vectTwoI16){ \
     (outVectI16) = vqsubq_s16((vectOneI16),(vectTwoI16));\
   }

   #define mmSubSatI32(outVectI32,vectOneI32,vectTwoI32){ \
     (outVectI32) = vqsubq_s32((vectOneI32),(vectTwoI32));\
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-04 Cat-06:
   +  - Conversion (casting)
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   #define cnvtI8ToI16(outVectI16, inVectI8) { \
     (outVectI16) = vreinterpretq_s16_s8((inVectI8)); \
   }

   #define cnvtI8ToI32(outVectI32, inVectI8) { \
     (outVectI32) = vreinterpretq_s32_s8((inVectI8)); \
   }

   #define cnvtI16ToI8(outVectI8, inVectI16) { \
     (outVectI8) = vreinterpretq_s8_s16((inVectI16)); \
   }

   #define cnvtI16ToI32(outVectI32, inVectI16) { \
     (outVectI32) = vreinterpretq_s32_s16((inVectI16)); \
   }

   #define cnvtI16ToI32(outVectI32, inVectI16) { \
     (outVectI32) = vreinterpretq_s32_s16((inVectI16)); \
   }

   #define cnvtI32ToI8(outVectI8, inVectI32) { \
     (outVectI8) = vreinterpretq_s8_s32((inVectI32)); \
   }

   #define cnvtI32ToI16(outVectI16, inVectI32) { \
     (outVectI16) = vreinterpretq_s16_s32((inVectI32)); \
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-04 Cat-07:
   +  - Max
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /* Byte maximize
   `  SSE2 does not support I8 max
   */
   #define mmMaxPrefI8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = vmaxq_s8((firstVectI8), (secVectI8)); \
   }

   #define mmMaxU8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = vmaxq_u8((firstVectI8), (secVectI8)); \
   }

   /*16 bit maximizes*/
   #define mmMaxI16(outVectI16, firstVectI16, secVectI16){\
      (outVectI16)=vmaxq_s16((firstVectI16),(secVectI16));\
   }

   #define mmMaxPrefU16( \
         outVectU16, \
         firstVectU16, \
         secVectU16 \
   ){\
      (outVectU16)=vmaxq_u16((firstVectU16),(secVectU16));\
   }

   /*32 unsigned bit maximizes (if possible)*/
   #define mmMaxPrefU32( \
      outVectU32, \
      firstVectU32, \
      secVectU32 \
   ){\
      (outVectU32)=vmaxq_u32((firstVectU32),(secVectU32));\
   }

   /*32 signed bit maximizes (if possible)*/
   #define mmMaxPrefI32( \
      outVectI32, \
      firstVectI32, \
      secVectI32 \
   ){\
      (outVectI32)=vmaxq_s32((firstVectI32),(secVectI32));\
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-04 Cat-08:
   +  - Min
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /* Byte minimize
   `  SSE2 does not support I8 min
   */
   #define mmMinPrefI8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = vminq_s8((firstVectI8), (secVectI8)); \
   }

   #define mmMinU8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = vminq_u8((firstVectI8), (secVectI8)); \
   }

   /*16 bit minimize*/
   #define mmMinI16(outVectI16, firstVectI16, secVectI16){\
      (outVectI16)=vminq_s16((firstVectI16),(secVectI16));\
   }

   #define mmMinPrefU16( \
      outVectU16, \
      firstVectU16, \
      secVectU16 \
   ){\
      (outVectU16)=vminq_u16((firstVectU16),(secVectU16));\
   }

   /*32 unsigned bit minimize (if possible)*/
   #define mmMinPrefU32( \
      outVectU32, \
      firstVectU32, \
      secVectU32 \
   ){\
     (outVectU32)=vminq_u32((firstVectU32),(secVectU32));\
   }

   /*32 signed bit minimize (if possible)*/
   #define mmMinPrefI32( \
      outVectI32, \
      firstVectI32, \
      secVectI32 \
   ){\
     (outVectI32)=vminq_s32((firstVectI32),(secVectI32));\
   }

/*********************************************************\
* Sec-05 Sub-05:
*  - Macros for arm (64bit) support (NEED TO FILL)
*  o sec-05 sub-05 cat-01:
*    - Input/output functions (load, set, store)
*  o sec-05 sub-05 cat-02:
*    - Comparison functions (if)
*  o sec-05 sub-05 cat-03:
*    - Logical functions (and, or, xor)
*  o sec-05 sub-05 cat-04:
*    - Bit manipulation functions (shifts)
*  o sec-04 sub-05 cat-06:
*    - Conversion (casting)
\*********************************************************/

#elif NEON64
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-05 Cat-01:
  +  - Input/output functions (load, set, store)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   // Regular loads
   #define mmLoadI8(outVectI8, arrayI8) { \
    (outVectI8) = vld1_s8((int8_t *) (arrayI8));\
   }

   #define mmLoadI16(outVectI16, arrayI16) { \
    (outVectI8) = vld1_s16((int16_t *) (arrayI16));\
   }

   #define mmLoadI32(outVectI32, arrayI32) { \
    (outVectI8) = vld1_s32((int32_t *) (arrayI32));\
   }

   // Unligned load options (here to support intel)
   #define mmLoadUI8(outVectI8, arrayI8) { \
    (outVectI8) = vld1_s8((int8_t *) (arrayI8));\
   }

   #define mmLoadUI16(outVectI16, arrayI16) { \
    (outVectI16) = vld1_s16((int16_t *) (arrayI16));\
   }

   #define mmLoadUI32(outVectI32, arrayI32) { \
    (outVectI32) = vld1_s32((int32_t *) (arrayI32));\
   }

   // single value loads
   #define mmSet1I8(inVectI8, valC) { \
     (inVectI8) = vdup_n_s8((valC)); \
   }

   #define mmSet1I16(inVectI16, valC) { \
     (inVectI16) = vdup_n_s16((valC)); \
   }

   #define mmSet1I32(inVectI32, valC) { \
     (inVectI32) = vdup_n_s32((valC)); \
   }

   // Set zero values (here to support intel)
   #define mmSeqZeroI8(inVectI8) { \
     (inVectI8) = vsub_s8(0); \
   }

   #define mmSeqZeroI16(inVectI16) { \
     (inVectI16) = vsub_s16(0); \
   }

   #define mmSeqZeroI32(inVectI32) { \
     (inVectI32) = vsub_s32(0); \
   }

  // Output functions; storing vectors in arrays
  #define mmStoreI8(array, inVectI8) { \
    vstl_s8((int8_t *) (array), (inVectI8)); \
  }

  #define mmStoreI16(array, inVectI16) { \
    vstl_s16((int16_t *) (array), (inVectI16)); \
  }

  #define mmStoreI32(array, inVectI32) { \
    vstl_s32((int32_t *) (array), (inVectI32)); \
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-05 Cat-05:
  +  - Comparison functions (if)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /*Store the mask functions (unsigned long returns)
   ` The logic for these store functions came from
   ` https://community.arm.com/arm-community-blogs/b/infrastructure-solutions-blog/posts/porting-x86-vector-bitmask-optimizations-to-arm-neon
   */
   #define storeI8CmpMask(retUL, inVectI8) { \
     (retUL) = \
       vget_lane_u64(vreinterpret_u64_s8((inVectI8)), 0);\
   }

   #define storeI16CmpMask(retUL, inVectI16) { \
     (retUL) = \
      vget_lane_u64(vreinterpret_u64_s16((inVectI16)), 0);\
   }

   #define storeI32CmpMask(retUL, inVectI32) { \
     (retUL) = \
      vget_lane_u64(vreinterpret_u64_s32((inVectI32)), 0);\
   }

   // Correct the bitcounts in the comparison masks
   #define fixI8CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>3)
   #define fixI16CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>4)
   #define fixI32CmpMaskCnt(retUL,inUL) ((retUL)=(inUL)>>5)

   // greater then comparisions
   #define mmCmpGtI8(outMaskI8, vectOneI8, vectTwoI8) { \
     (outMaskI8) = vcgt_s8((vectOneI8), (vectTwoI8)); \
   };

   #define mmCmpGtI16(outMaskI16,vectOneI16,vectTwoI16) { \
     (outMaskI16) = vcgt_s16((vectOneI16), (vectTwoI16));\
   };

   #define mmCmpGtI32(outMaskI32,vectOneI32,vectTwoI32) { \
     (outMaskI32) = vcgt_s32((vectOneI32), (vectTwoI32)); \
   };

   // Less then comparisions
   #define mmCmpLtI8(outMaskI8, vectOneI8, vectTwoI8) { \
     (outMaskI8) = vclt_s8((vectOneI8), (vectTwoI8)); \
   };

   #define mmCmpLtI16(outMaskI16,vectOneI16,vectTwoI16) { \
     (outMaskI16) = vclt_s16((vectOneI16), (vectTwoI16)); \
   };

   #define mmCmpLtI32(outMaskI32,vectOneI32,vectTwoI32) { \
     (outMaskI32) = vclt_s32((vectOneI32), (vectTwoI32)); \
   };

   // Less then comparisions
   #define mmCmpEqI8(outMaskI8, vectOneI8, vectTwoI8) { \
     (outMaskI8) = vceq_s8((vectOneI8), (vectTwoI8)); \
   };

   #define mmCmpEqI16(outMaskI16,vectOneI16,vectTwoI16) { \
     (outMaskI16) = vceq_s16((vectOneI16), (vectTwoI16)); \
   };

   #define mmCmpEqI32(outMaskI32,vectOneI32,vectTwoI32) { \
     (outMaskI32) = vceq_s32((vectOneI32), (vectTwoI32)); \
   };

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-05 Cat-03:
  +  - Logical functions (and, or, xor)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   #define mmAndNotI8(outVectI8, notVectI8, inVect8) { \
     (outVectI8) = \
       vand_s8(vmvn_s8((notVectI8)), (inVect8)); \
   }

   #define mmAndNotI16(outVectI16, notVectI16, inVect16){ \
     (outVectI16) = \
       vand_s16(vmvn_s16((notVectI16)), (inVect16)); \
   }

   #define mmAndNotI32(outVectI32, notVectI32, inVect32){ \
     (outVectI32) = \
       vand_s32(vmvn_s32((notVectI32)), (inVect32)); \
   }

   #define mmAndI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vand_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmAndI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = vand_s16((vectOneI16), (vectTwoI16)); \
   }

   #define mmAndI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = vand_s32((vectOneI32), (vectTwoI32)); \
   }

   #define mmOrI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vorr_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmOrI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = vorr_s16((vectOneI16), (vectTwoI16)); \
   }

   #define mmOrI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = vorr_s32((vectOneI32), (vectTwoI32)); \
   }

   #define mmXorI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = veor_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmXorI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = veor_s16((vectOneI16), (vectTwoI16)); \
   }

   #define mmXorI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = veor_s32((vectOneI32), (vectTwoI32)); \
   }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-05 Cat-04:
  +  - Bit manipulation functions (shifts)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  #define mmShiftRightI8(outVectI8, inVectI8, numBytesI){\
    (outVectI8) = vshr_s8((inVectI8), (numBytesI)); \
  }

  #define mmShiftRightI16(outVectI16,inVectI16,numBytesI){\
    (outVectI16) = vshr_s16((inVectI16), (numBytesI)); \
  }

  #define mmShiftRightI32(outVectI32,inVectI32,numBytesI){\
    (outVectI32) = vshr_s32((inVectI32), (numBytesI)); \
  }

  #define mmShiftLeftI8(outVectI8, inVectI8, numBytesI){\
    (outVectI8) = vshl_s8((inVectI8), (numBytesI)); \
  }

  #define mmShiftLeftI16(outVectI16, inVectI16, numBytesI){\
    (outVectI16) = vshl_s16((inVectI16), (numBytesI)); \
  }

  #define mmShiftLeftI32(outVectI32, inVectI32, numBytesI){\
    (outVectI32) = vshl_s32((inVectI32), (numBytesI)); \
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++\
  + Sec-05 Sub-05 Cat-05:
  +  - Math functions (add and sub)
  \*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   // Addition (non-saturating)
   #define mmAddI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vadd_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmAddI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = vadd_s16((vectOneI16), (vectTwoI16)); \
   }

   #define mmAddI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = vadd_s32((vectOneI32), (vectTwoI32)); \
   }

   // Addition (saturating)
   #define mmAddSatI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vqadd_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmAddSatI16(outVectI16,vectOneI16,vectTwoI16){ \
     (outVectI16) = vqadd_s16((vectOneI16), (vectTwoI16));\
   }

   #define mmAddSatI32(outVectI32,vectOneI32,vectTwoI32){ \
     (outVectI32) = vqadd_s32((vectOneI32), (vectTwoI32));\
   }

   // subtraction (non-saturating)
   #define mmsubI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vsub_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmsubI16(outVectI16, vectOneI16, vectTwoI16) { \
     (outVectI16) = vsub_s16((vectOneI16), (vectTwoI16)); \
   }

   #define mmsubI32(outVectI32, vectOneI32, vectTwoI32) { \
     (outVectI32) = vsub_s32((vectOneI32), (vectTwoI32)); \
   }

   // subtraction (saturating)
   #define mmSubSatI8(outVectI8, vectOneI8, vectTwoI8) { \
     (outVectI8) = vqsub_s8((vectOneI8), (vectTwoI8)); \
   }

   #define mmSubSatI16(outVectI16,vectOneI16,vectTwoI16){ \
     (outVectI16) = vqsub_s16((vectOneI16), (vectTwoI16));\
   }

   #define mmSubSatI32(outVectI32,vectOneI32,vectTwoI32){ \
     (outVectI32) = vqsub_s32((vectOneI32), (vectTwoI32));\
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-05 Cat-06:
   +  - Conversion (casting)
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   #define cnvtI8ToI16(outVectI16, inVectI8) { \
     (outVectI16) = vreinterpret_s16_s8((inVectI8)); \
   }

   #define cnvtI8ToI32(outVectI32, inVectI8) { \
     (outVectI32) = vreinterpret_s32_s8((inVectI8)); \
   }

   #define cnvtI16ToI8(outVectI8, inVectI16) { \
     (outVectI8) = vreinterpret_s8_s16((inVectI16)); \
   }


   #define cnvtI16ToI32(outVectI32, inVectI16) { \
     (outVectI32) = vreinterpret_s32_s16((inVectI16)); \
   }

   #define cnvtI16ToI32(outVectI32, inVectI16) { \
     (outVectI32) = vreinterpret_s32_s16((inVectI16)); \
   }

   #define cnvtI32ToI8(outVectI8, inVectI32) { \
     (outVectI8) = vreinterpret_s8_s32((inVectI32)); \
   }

   #define cnvtI32ToI16(outVectI16, inVectI32) { \
     (outVectI16) = vreinterpret_s16_s32((inVectI32)); \
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-05 Cat-07:
   +  - Max
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /* Byte maximize
   `  SSE2 does not support I8 max
   */
   #define mmMaxPrefI8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = vmax_s8((firstVectI8), (secVectI8)); \
   }

   #define mmMaxU8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = vmax_u8((firstVectI8), (secVectI8)); \
   }

   /*16 bit maximizes*/
   #define mmMaxI16(outVectI16, firstVectI16, secVectI16){\
      (outVectI16)=vmax_s16((firstVectI16),(secVectI16));\
   }

   #define mmMaxPrefU16( \
         outVectU16, \
         firstVectU16, \
         secVectU16 \
   ){\
      (outVectU16)=vmax_u16((firstVectU16),(secVectU16));\
   }

   /*32 unsigned bit maximizes (if possible)*/
   #define mmMaxPrefU32( \
      outVectU32, \
      firstVectU32, \
      secVectU32 \
   ){\
      (outVectU32)=vmax_u32((firstVectU32),(secVectU32));\
   }

   /*32 signed bit maximizes (if possible)*/
   #define mmMaxPrefI32( \
      outVectI32, \
      firstVectI32, \
      secVectI32 \
   ){\
      (outVectI32)=vmax_s32((firstVectI32),(secVectI32));\
   }

   /*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Sec-05 Sub-05 Cat-08:
   +  - Min
   \+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /* Byte minimize
   `  SSE2 does not support I8 min
   */
   #define mmMinPrefI8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = vmin_s8((firstVectI8), (secVectI8)); \
   }

   #define mmMinU8(outVectI8,firstVectI8,secVectI8){\
      (outVectI8) = vmin_u8((firstVectI8), (secVectI8)); \
   }

   /*16 bit minimize*/
   #define mmMinI16(outVectI16, firstVectI16, secVectI16){\
      (outVectI16)=vmin_s16((firstVectI16),(secVectI16));\
   }

   #define mmMinPrefU16( \
      outVectU16, \
      firstVectU16, \
      secVectU16 \
   ){\
      (outVectU16)=vmin_u16((firstVectU16),(secVectU16));\
   }

   /*32 unsigned bit minimize (if possible)*/
   #define mmMinPrefU32( \
      outVectU32, \
      firstVectU32, \
      secVectU32 \
   ){\
     (outVectU32)=vmin_u32((firstVectU32),(secVectU32));\
   }

   /*32 signed bit minimize (if possible)*/
   #define mmMinPrefI32( \
      outVectI32, \
      firstVectI32, \
      secVectI32 \
   ){\
     (outVectI32)=vmin_s32((firstVectI32),(secVectI32));\
   }
#endif /*Commands for AVX2, SSE2, SSE4, and NEON*/

#else
  /*This only happens for scaler code. This is here so
  ` the user can use these values to modify their code,
  ` but also have then got to one byte for scaler
  */
  #define vectorBits 8
  #define vectorBytes 1

#endif /*User is using scalar side of their pogram, this
       ` is to disable compiler errrors when vectorWarp.c
       ` is included in build for vector versions of code
       */
#endif /*if this header has not already been defined*/
