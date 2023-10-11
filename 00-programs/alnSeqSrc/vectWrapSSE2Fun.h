#ifndef VECTWRAPSSE2_H
#define VECTWRAPSSE2_h 

#include <immintrin.h>

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-01:
^  - Definitions and variable declerations
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define vectorBits 128
#define vectorBytes 16

typedef __m128i vectI8;  /*vector of bytes (8bits)*/
typedef __m128i vectI16; /*vector of shorts (16 bits)*/
typedef __m128i vectI32; /*vector of ints (32 bits)*/

typedef __m128i vectU8;  /*Prefer unsigned bytes*/
typedef __m128i vectU16; /*Prefer unsigned shorts*/
typedef __m128i vectU32; /*Prefer unsinged ints*/

typedef __m128i mask8;  /*mask of 8 bit values*/
typedef __m128i mask16; /*Mask of 16 bit values*/
typedef __m128i mask32; /*mask of 32 bit values*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-02:
^  - Macros
^  o sec-02 sub-01:
^    - Input/output functions (load, set, store)
^  o sec-02 sub-02:
^    - Output functions
^  o sec-02 sub-03:
^    - Comparison functions (if)
^  o sec-02 sub-04:
^    - Logical functions (andNot, and, or, xor)
^  o sec-02 sub-05:
^    - Bit manipulation functions (shifts)
^  o sec-02 sub-06:
^    - Math functions
^  o sec-02 sub-07:
^    - Casts [for NEON & masks] and Conversions
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*********************************************************\
* Sec-02 Sub-01:
*  - Input functions
*  o sec-02 sub-01 cat-01:
*    - aligned loading
*  o sec-02 sub-01 cat-02:
*    - Unaligned loading
*  o sec-02 sub-01 cat-03:
*    - make zero vectors
*  o sec-02 sub-01 cat-04:
*    - Make vectors of one element
*  o sec-02 sub-01 cat-05:
*    - Insert an element into an vector
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-01:
+  - Aligned loading
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define load_I8Ary_retVectI8(arrayI8) { \
  _mm_load_si128((__m128i *) (arrayI8)); \
}

#define load_I16Ary_retVectI16(arrayI16) { \
  _mm_load_si128((__m128i *) (arrayI16));\
}

#define load_I32Ary_retVectI32(arrayI32) { \
  _mm_load_si128((__m128i *) (arrayI32));\
}


#define load_U8Ary_retVectU8(arrayU8) { \
  _mm_load_si128((__m128i *) (arrayU8)); \
}

#define load_U16Ary_retVectU16(arrayU16) { \
  _mm_load_si128((__m128i *) (arrayU16));\
}

#define load_U32Ary_retVectU32(arrayU32) { \
  _mm_load_si128((__m128i *) (arrayU32));\
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-01:
+  - Unaligned loading
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#define loadu_I8Ary_retVectI8(arrayI8) { \
  _mm_loadu_si128((__m128i *) (arrayI8)); \
}

#define loadu_I16Ary_retVectI16(arrayI16) { \
  _mm_loadu_si128((__m128i *) (arrayI16));\
}

#define loadu_I32Ary_retVectI32(arrayI32) { \
  _mm_loadu_si128((__m128i *) (arrayI32));\
}


#define loadu_U8Ary_retVectU8(arrayU8) { \
  _mm_loadu_si128((__m128i *) (arrayU8)); \
}

#define loadu_U16Ary_retVectU16(arrayU16) { \
  _mm_loadu_si128((__m128i *) (arrayU16));\
}

#define loadu_U32Ary_retVectU32(arrayU32) { \
  _mm_loadu_si128((__m128i *) (arrayU32));\
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-03:
+  - make zero vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define zero_retVectI8()  (_mm_setzero_si128();)
#define zero_retVectI16() (_mm_setzero_si128();)
#define zero_retVectI32() (_mm_setzero_si128();)

#define zero_retVectU8()  (_mm_setzero_si128();)
#define zero_retVectU16() (_mm_setzero_si128();)
#define zero_retVectU32() (_mm_setzero_si128();)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-04:
+  - Make vectors of one element
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define set1_I8_retVectI8(valC) (mm_set1_epi8((valC));)

#define set1_I16_retVectI16(valI16) \
   (mm_set1_epi16((valI16));)

#define set1_I32_retVectI32(valI32) \
   (_mm128_set1_epi32((valI32));)


#define set1_U8_retVectU8(valUC) (mm_set1_epi8((valUC));)

#define set1_U16_retVectU16(valC) \
   (mm_set1_epi16((valU16));)

#define set1_U32_retVectU32(valU32) \
   (_mm_set1_epi32((valU32));)


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-05:
+  - Insert an element into an vector
+  - https://stackoverflow.com/questions/58303958/how-to-implement-16-and-32-bit-integer-insert-and-extract-operations-with-avx-51
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*TODO: find the sse insert functions*/
#define insert_I8_retVectI8(insI8, posI){\
   _mm128_mask_set1_epi8(insVectI8,1UL<<posI,insI8);\
}

#define insert_I16_retVectI16(insI16, posI){\
   _mm128_mask_set1_epi16(insVectI16,1UL<<posI,insI16);\
}         

#ifdef SSE4
   #define insert_I32_retVectI32(insI32, posI){\
      _mm128_mask_set1_epi32(insVectI32,1UL<<posI,insI32);\
   }
#elif
  /*SSE2 replacement (call insert twice)*/
   #define insert_I32_retVectI32(insI32, posI){\
   }
#endif


#define insert_U8_retVectU8(insU8, posU){\
   _mm128_mask_set1_epi8(insVectU8,1UL<<posU,insU8);\
}

#define insert_U16_retVectU16(insU16, posU){\
   _mm128_mask_set1_epi16(insVectU16,1UL<<posU,insU16);\
}         

#ifdef SSE4
   #define insert_U32_retVectU32(insU32, posU){\
      _mm128_mask_set1_epi32(insVectU32,1UL<<posU,insU32);\
   }
#elif
  /*SSE2 replacement (call insert twice)*/
   #define insert_U32_retVectU32(insU32, posU){\
   }
#endif

/*********************************************************\
* Sec-02 Sub-02:
*  - Output functions
*  - sec-02 sub-02 cat-01:
*    - Store vector ouput into an array
*  - sec-02 sub-02 cat-02:
*    - Store masks into longs
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-01:
+  - Store vector ouput into an array
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define store_vectI8_retAryI8(array, inVectI8) { \
  _mm_store_si128((__m128i *) (array), (inVectI8)); \
}

#define store_vectI16_retAryI16(array, inVectI16) { \
  _mm_store_si128((__m128i *) (array), (inVectI16)); \
}

#define store_vectI32_retAryI32(array, inVectI32) { \
  _mm_store_si128((__m128i *) (array), (inVectI32)); \
}


#define store_vectU8_retAryU8(array, inVectU8) { \
  _mm_store_si128((__m128i *) (array), (inVectU8)); \
}

#define store_vectU16_retAryU16(array, inVectU16) { \
  _mm_store_si128((__m128i *) (array), (inVectU16)); \
}

#define store_vectU32_retAryU32(array, inVectU32) { \
  _mm_store_si128((__m128i *) (array), (inVectU32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-02:
+  - Store masks into longs
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define store_mask8_retUL(retUL, mask8) { \
   (retUL) = _mm_movemask_epi8((mask8)) \
} /*Mask is 8 bit integers*/

#define store_mask16_retUL(retUL, mask16) { \
   (retUL) = _mm_movemask_epi8((mask16)) \
} /*Mask is 16 bit integers*/

#define store_mask32_retUL(retUL, mask32) { \
   (retUL) = _mm_movemask_epi8((mask32)) \
} /*Mask is 32 bit integers*/

/*********************************************************\
* Sec-02 Sub-03:
*  - Comparison functions (if)
*  o sec-02 sub-03 cat-01:
*    - Comparisons (return masks)
*  o sec-02 sub-03 cat-02:
*    - Fix differences in population counts (total 1's)
*    - This is here because I am also supporting SSE2 and
*      NEON.
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-01:
+  - Comparisons (return masks)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*8-bit comparisions; returns a mask8*/
#define cmpeq_vectI8_retMask8(vectOneI8, vectTwoI8) { \
   _mm_cmpeq_epi8((vectOneI8), (vectTwoI8)); \
}

#define cmpeq_vectU8_retMask8(vectOneI8, vectTwoI8) { \
   _mm_cmpeq_epi8((vectOneI8), (vectTwoI8)); \
}

#define cmpgt_vectI8_retMask8(vectOneI8, vectTwoI8) {\
   _mm_cmpgt_epi8((vectOneI8), (vectTwoI8)); \
}

#define cmplt_vectI8_retMask8(vectOneI8, vectTwoI8) { \
   _mm_cmplt_epi8((vectOneI8), (vectTwoI8)); \
}

/*16-bit comparisions; returns a mask16*/

#define cmpeq_vectI16_retMask16(vectOneI16, vectTwoI16) { \
   _mm_cmpeq_epi16((vectOneI16), (vectTwoI16)); \
}

#define cmpeq_vectU16_retMask16(vectOneI16, vectTwoI16) { \
   _mm_cmpeq_epi16((vectOneI16), (vectTwoI16)); \
}

#define cmpgt_vectI16_retMask16(vectOneI16, vectTwoI16) {\
   _mm_cmpgt_epi16((vectOneI16), (vectTwoI16)); \
}

#define cmplt_vectI16_retMask16(vectOneI16, vectTwoI16) { \
   _mm_cmplt_epi16((vectOneI16), (vectTwoI16)); \
}

/*32-bit comparisions; returns a mask32*/

#define cmpeq_vectI32_retMask32(vectOneI32, vectTwoI32) { \
   _mm_cmpeq_epi32((vectOneI32), (vectTwoI32)); \
}

#define cmpeq_vectU32_retMask32(vectOneI32, vectTwoI32) { \
   _mm_cmpeq_epi32((vectOneI32), (vectTwoI32)); \
}

#define cmpgt_vectI32_retMask32(vectOneI32, vectTwoI32) {\
   _mm_cmpgt_epi32((vectOneI32), (vectTwoI32)); \
}

#define cmpLt_vectI32_retMask32(vectOneI32, vectTwoI32) { \
   _mm_cmplt_epi32((vectOneI32), (vectTwoI32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-02:
+  - Fix differences in population counts (total 1's)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* This is here to support SSE2, AVX2, and NEON. 
`  The returned mask in AVX128 has only one bit per
`  data type in comparison
`  SSE (all) and AVX2 returns a 128 bit vector that is
`   converted to characters
*/
#define fix_mask8_popcount(inUL)  ((inUL))
#define fix_mask16_popcount(inUL) ((inUL) >> 1)
#define fix_mask32_popcount(inUL) ((inUL) >> 2)

/*********************************************************\
* Sec-02 Sub-04:
*  - Logical functions (andNot, and, or, xor)
*  o sec-02 sub-02 cat-01:
*    - and not functions
*  o sec-02 sub-02 cat-02:
*    - and functions
*  o sec-02 sub-02 cat-03:
*    - or functions
*  o sec-02 sub-02 cat-04:
*    - xor functions
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-01:
+  - and not functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define andNot_vectI8_retVectI8(notVectI8, inVectI8) { \
   _mm_andnot_si128((notVectI8), (inVectI8)); \
}

#define andNot_vectU8_retVectU8(notVectU8, inVectU8) { \
   _mm_andnot_si128((notVectU8), (inVectU8)); \
}

#define andNot_vectI16_retVectI16(notVectI16, inVectI16) {\
   _mm_andnot_si128((notVectI16), (inVectI16)); \
}

#define andNot_vectU16_retVectU16(notVectU16, inVectU16) {\
   _mm_andnot_si128((notVectU16), (inVectU16)); \
}

#define andNot_vectI32_retVectI32(notVectI32, inVectI32) {\
   _mm_andnot_si128((notVectI32), (inVectI32)); \
}

#define andNot_vectU32_retVectU32(notVectU32, inVectU32) {\
   _mm_andnot_si128((notVectU32), (inVectU32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-02:
+  - and functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define and_vectI8_retVectI8(vectOneI8, vectTwoI8) { \
   _mm_and_si128((vectOneI8), (vectTwoI8)); \
}

#define and_vectU8_retVectU8(vectOneU8, vectTwoU8) { \
   _mm_and_si128((vectOneU8), (vectTwoU8)); \
}

#define and_vectI16_retVectI16(vectOneI16, vectTwoI16) { \
   _mm_and_si128((vectOneI16), (vectTwoI16)); \
}

#define and_vectU16_retVectU16(vectOneU16, vectTwoU16) { \
   _mm_and_si128((vectOneU16), (vectTwoU16)); \
}

#define and_vectI32_retVectI32(vectOneI32, vectTwoI32) { \
   _mm_and_si128((vectOneI32), (vectTwoI32)); \
}

#define and_vectU32_retVectU32(vectOneU32, vectTwoU32) { \
   _mm_and_si128((vectOneU32), (vectTwoU32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-03:
+  - or functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define or_vectI8_retVectI8(vectOneI8, vectTwoI8) { \
   _mm_or_si128((vectOneI8), (vectTwoI8)); \
}

#define or_vectU8_retVectU8(vectOneU8, vectTwoU8) { \
   _mm_or_si128((vectOneU8), (vectTwoU8)); \
}

#define or_vectI16_retVectI16(vectOneI16, vectTwoI16) { \
   _mm_or_si128((vectOneI16), (vectTwoI16)); \
}

#define or_vectU16_retVectU16(vectOneU16, vectTwoU16) { \
   _mm_or_si128((vectOneU16), (vectTwoU16)); \
}

#define or_vectI32_retVectI32(vectOneI32, vectTwoI32) { \
   _mm_or_si128((vectOneI32), (vectTwoI32)); \
}

#define or_vectU32_retVectU32(vectOneU32, vectTwoU32) { \
   _mm_or_si128((vectOneU32), (vectTwoU32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-03:
+  - or functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define xor_vectI8_retVectI8(vectOneI8, vectTwoI8) { \
   _mm_xor_si128((vectOneI8), (vectTwoI8)); \
}

#define xor_vectU8_retVectU8(vectOneU8, vectTwoU8) { \
   _mm_xor_si128((vectOneU8), (vectTwoU8)); \
}

#define xor_vectI16_retVectI16(vectOneI16, vectTwoI16) { \
   _mm_xor_si128((vectOneI16), (vectTwoI16)); \
}

#define xor_vectU16_retVectU16(vectOneU16, vectTwoU16) { \
   _mm_xor_si128((vectOneU16), (vectTwoU16)); \
}

#define xor_vectI32_retVectI32(vectOneI32, vectTwoI32) { \
   _mm_xor_si128((vectOneI32), (vectTwoI32)); \
}

#define xor_vectU32_retVectU32(vectOneU32, vectTwoU32) { \
   _mm_xor_si128((vectOneU32), (vectTwoU32)); \
}

/*********************************************************\
* Sec-02 Sub-05:
*  - Bit manipulation functions (shifts)
*  o sec-02 sub-05 cat-01:
*    - Shift numbers in vectors right by x bits
*  o sec-02 sub-05 cat-03:
*    - Shift a vector right by x bytes
*  o sec-02 sub-05 cat-02:
*    - Shift numbers in vectors left by x bits
*  o sec-02 sub-05 cat-04:
*    - Shift a vector left by x bytes
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-01:
+  - Shift each element right by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sr_vectI8_retVectI8(inVectI8, numBitsI){\
    _mm_srli_epi8((inVectI8), (numBitsI)); \
}

#define sr_vectU8_retVectU8(inVectU8, numBitsI){\
    _mm_srli_epi8((inVectU8), (numBitsI)); \
}

#define sr_vectI16_retVectI16(inVectI16, numBitsI){\
    _mm_srli_epi16((inVectI16), (numBitsI)); \
}

#define sr_vectU16_retVectU16(inVectU16, numBitsI){\
    _mm_srli_epi16((inVectU16), (numBitsI)); \
}

#define sr_vectI32_retVectI32(inVectI32, numBitsI){\
    _mm_srli_epi32((inVectI32), (numBitsI)); \
}

#define sr_vectU32_retVectU32(inVectU32, numBitsI){\
    _mm_srli_epi32((inVectU32), (numBitsI)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-02:
+  - Shift a vector right by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*arm idea from
https://stackoverflow.com/questions/11854819/attempt-to-convert-sse2-fast-corner-score-code-to-arm-neon
   _mm_srli_si128(q0, 4) = vextq_s16(q0, zeroVect, 2);
      copies first two 0's from 0 vect to q0
   _mm_srli_si128(q0, 4) = vextq_s16(zeroVect, q0, 2);
      copies last two numbers from 0 vect to q0
      This is my own logic (not tested)
   _mm_srli_si128(q0, 2) = vextq_s16(q0, zeroVect, 1);
      copies first 0 from 0 vect to q0
*/

#define srvect_vectI8_retVectI8(inVectI8, numBytesI){\
    _mm_srli_si128((inVectI8), (numBytesI)); \
}

#define srvect_vectU8_retVectU8(inVectU8, numBytesI){\
    _mm_srli_si128((inVectU8), (numBytesI)); \
}

#define srvect_vectI16_retVectI16(inVectI16, numBytesI){\
    _mm_srli_si128((inVectI16)); \
}

#define srvect_vectU16_retVectU16(inVectU16, numBytesI){\
    _mm_srli_si128((inVectU16), (numBytesI)); \
}

#define srvect_vectI32_retVectI32(inVectI32, numBytesI){\
    _mm_srli_si128((inVectI32), (numBytesI)); \
}

#define srvect_vectU32_retVectU32(inVectU32, numBytesI){\
    _mm_srli_si128((inVectU32), (numBytesI)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-03:
+  - Shift numbers in vectors left by x bits
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sl_vectI8_retVectI8(inVectI8, numBitsI){\
    _mm_slli_epi8((inVectI8), (numBitsI)); \
}

#define sl_vectU8_retVectU8(inVectU8, numBitsI){\
    _mm_slli_epi8((inVectU8), (numBitsI)); \
}

#define sl_vectI16_retVectI16(inVectI16, numBitsI){\
    _mm_slli_epi16((inVectI16), (numBitsI)); \
}

#define sl_vectU16_retVectU16(inVectU16, numBitsI){\
    _mm_slli_epi16((inVectU16), (numBitsI)); \
}

#define sl_vectI32_retVectI32(inVectI32, numBitsI){\
    _mm_slli_epi32((inVectI32), (numBitsI)); \
}

#define sl_vectU32_retVectU32(inVectU32, numBitsI){\
    _mm_slli_epi32((inVectU32), (numBitsI)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-04:
+  - Shift a vector left by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sllvect_vectI8_retVectI8(inVectI8, numBytesI){\
   _mm_slli_si128((inVectI8), (numBytesI)); \
}

#define sllvect_vectU8_retVectU8(inVectU8, numBytesI){\
   _mm_slli_si128((inVectU8), (numBytesI)); \
}

#define sllvect_vectI16_retVectI16(inVectI16, numBytesI){\
   _mm_slli_si128((inVectI16), (numBytesI)); \
}

#define sllvect_vectU16_retVectU16(inVectU16, numBytesI){\
   _mm_slli_si128((inVectU16), (numBytesI)); \
}

#define sllvect_vectI32_retVectI32(inVectI32, numBytesI){\
   _mm_slli_si128((inVectI32), (numBytesI)); \
}

#define sllvect_vectU32_retVectU32(inVectU32, numBytesI){\
   _mm_slli_si128((inVectU32), (numBytesI)); \
}

/*********************************************************\
* Sec-02 Sub-06:
*  - Math functions
*  o sec-02 sub-06 cat-01:
*    - addition max 64 bit; [sse2-avx128]; no epu
*  o sec-02 sub-06 cat-02:
*    - staturation addition; max 16 bit [sse2-avx128]; epu
*  o sec-02 sub-06 cat-03:
*    - subtraction; max 64bit [sse2-avx128]; no epu
*  o sec-02 sub-06 cat-04:
*    - staturation subtraction; max 16 bit; [sse2-avx128];u
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-01:
+  - addition max 64 bit; [sse2-avx128]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define add_vectI8_retVectI8(vectOneI8, vectTwoI8) {\
   _mm_add_epi8((vectOneI8), (vectTwoI8)); \
}

#define add_vectI16_retVectI16(vectOneI16, vectTwoI16) {\
   _mm_add_epi16((vectOneI16), (vectTwoI16)); \
}

#define add_vectI32_retVectI32(vectOneI32, vectTwoI32) {\
   _mm_add_epi32((vectOneI32), (vectTwoI32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-02:
+  - staturation addition; max 16 bit [sse2 to avx128]; epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*TODO: Check if hav esigned saturation adding*/

#define addSat_vectU8_retVectU8(vectOneU8, vectTwoU8) {\
   _mm_adds_epu8((vectOneU8), (vectTwoU8)); \
}

#define addSat_vectU16_retVectU16(vectOneU16, vectTwoU16){\
   _mm_adds_epu16((vectOneU16), (vectTwoU16)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-03:
+  - subtraction; max 64bit [sse2-avx128]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sub_vectI8_retVectI8(vectOneI8, vectTwoI8) { \
    _mm_sub_epi8((vectOneI8), (vectTwoI8)); \
}

#define sub_vectI16_retVectI16(vectOneI16, vectTwoI16) { \
    _mm_sub_epi16((vectOneI16), (vectTwoI16)); \
}

#define sub_vectI32_retVectI32(vectOneI32, vectTwoI32) { \
    _mm_sub_epi32((vectOneI32), (vectTwoI32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-04:
+  - staturation subtraction; max 16 bit; [sse2-avx128]); u
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define subSat_vectU8_retVectU8(vectOneU8, vectTwoU8) {\
   _mm_subs_epu8((vectOneU8), (vectTwoU8)); \
}

#define subSat_vectU16_retVectU16(vectOneU16, vectTwoU16) {\
   _mm_subs_epu16((vectOneU16), (vectTwoU16)); \
}

/*********************************************************\
* Sec-02 Sub-07:
*  - Conversions [Here for NEON support]
*  o sec02 sub-07 cat-01:
*    - Convert signed vectors to other signed vectors
*  o sec02 sub-07 cat-02:
*    - Convert unsigned vectors to unsigned vectors
*  o sec02 sub-07 cat-02:
*    - Convert unsigned vectors to unsigned vectors
*  o sec02 sub-07 cat-06:
*    - Convert masks to vectors
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-01:
+  - Convert signed vectors to other signed vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectI8_to_VectI16(inVectI8) { \
  (vectI16) (inVectI8); \
}

#define cast_vectI8_to_VectI32(inVectI8) { \
  (vectI32) (inVectI8); \
}

#define cast_vectI16_to_VectI8(inVectI16) { \
  (vectI8) (inVectI16); \
}

#define cast_vectI16_to_VectI32(inVectI16) { \
  (vectI32) (inVectI16); \
}

#define cast_vectI32_to_VectI8(inVectI32) { \
  (vectI8) (inVectI32); \
}

#define cast_vectI32_to_VectI16(inVectI32) { \
  (vectI16) (inVectI32); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-02:
+  - Convert unsigned vectors to unsigned vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectU8_to_VectU16(inVectU8) { \
  (vectU16) (inVectU8); \
}

#define cast_vectU8_to_VectU32(inVectU8) { \
  (vectU32) (inVectU8); \
}

#define cast_vectU16_to_VectU8(inVectU16) { \
  (vectU8) (inVectU16); \
}

#define cast_vectU16_to_VectU32(inVectU16) { \
  (vectU32) (inVectU16); \
}

#define cast_vectU32_to_VectU8(inVectU32) { \
  (vectU8) (inVectU32); \
}

#define cast_vectU32_to_VectU16(inVectU32) { \
  (vectU16) (inVectU32); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-02:
+  - Convert masks to vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define mask8_retVectI8(inMask8){ \
   _mm_movemask_epi8((inMask8)); \
}

#define mask16_retVectI16(inMask16){ \
   _mm_movemask_epi8((inMask8)); \
}

#define mask32_retVectI32(inMask32){ \
   _mm_movemask_epi8((inMask8)); \
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-07:
+  - Max
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#ifdef SSE4
   #define max_vectI8_retVectI8(firstVectI8,secVectI8){\
      _mm_max_epi8((firstVectI8), (secVectI8)); \
   }
#else
   /*SSE2 does not support I8 max, So I need to do a
   ` branchless max
   */
   #define max_vectI8_retVectI8(firstVectI8,secVectI8){\
      cmpVectI8 = _mm_cmplt_epi8(firstVectI8,secVectI8); \
      xorVectI8 = _mm_xor_epi8(firstVectI8, secVectI8); \
      cmpVectI8 = _mm_and_epi8(cmpVectI8, xorVectI8); \
      /*At this point the xor value is only kept if
      ` the frist vector is less than the second vector
      */\
      cmpVectI8 = _mm_xor_epi8(cmpVectI8, firstVectI8); \
      /*The first vector is only kept if it is >= the
      ` second vector, else the second vector is kept.
      ` For min, change firstVectI8 to secVectI8.
      */
#endif

#define max_vectU8_retVectU8(firstVectU8,secVectU8){\
   _mm128_max_epu8((firstVectU8), (secVectU8)); \
}

#define max_vectI16_retVectI16(firstVectI16,secVectI16){\
   _mm128_max_epi16((firstVectI16), (secVectI16)); \
}

#ifdef SSE4
  #define max_vectU16_retVectU16(firstVectU16,secVectU16){\
     _mm128_max_epu16((firstVectU16), (secVectU16)); \
  }
#else
  #define max_vectU16_retVectU16(firstVectU16,secVectU16){\
     vect1No1 = _mm_srli_epi16(firstVectU16, 1); \
     vect2No1 = _mm_srli_epi16(secVectU16, 1); \
     vec2No1 = _mm_max_epi16(vect1No1, vect2No1); \
       /*The only problem case is when |vect1 - vect2| = 1
       ` This next part is to handle this problem.
       */ \
     \
     oneVect = _mm_set1_epi16(1); \
     vect1No1 = _mm_xor_si128(firstVectU16, secVectU16); \
       /*This only keeps the differnt bits*/ \
     vect1No1 = _mm_cmpeq_epi16(vect1No1, oneVect); \
       /*This tells me if the difference was only by 1*/\
     vect1No1 = _mm_and_si128(vect1No1, oneVect); \
       /*vect1No1 is either -1 (all ones) or 0 for each
       ` short. I am changing to have either be 0 or 1
       */\
    vect2No1 | vect1No1; \
       /*Add in the ones*/\
  }
#endif

#ifdef SSE4
  #define max_vectI32_retVectI32(firstVectI32,secVectI32){\
     _mm128_max_epi32((firstVectI32), (secVectI32)); \
  } /*SSE2 will use branchless max function*/

#else
  /*SS#2 does not support 32 bit max*/
  #define min_vectI32_retVectI32(firstVectI32,secVectI32){\
    cmpVectI32 = _mm_cmplt_epi32(firstVectI32,secVectI32);\
    xorVectI32 = _mm_xor_epi32(firstVectI32, secVectI32);\
    cmpVectI32 = _mm_and_epi32(cmpVectI32, xorVectI32); \
      /*At this point the xor value is only kept if
      ` the frist vector is less than the second vector
      */\
    cmpVectI32 = _mm_xor_epi32(cmpVectI32, firstVectI32);\
      /*The first vector is only kept if it is >= the
      ` second vector, else the second vector is kept.
      ` For min, change firstVectI8 to secVectI8.
      */\
  } /*SSE2 will use branchless min function*/
#endif

#ifdef SSE4
  #define max_vectU32_retVectU32(firstVectU32,secVectU32){\
     _mm_max_epu32((firstVectU32), (secVectU32)); \
  }

#else
   /*SSE2 does not support 32 bit axes*/

  #define max_vectU32_retVectU32(firstVectU32,secVectU32){\
     /*Remove the 1/0 bit and shit sign one back*/ \
     vect1No1 = _mm_srli_epi16(firstVectU16, 1); \
     vect2No1 = _mm_srli_epi16(secVectU16, 1); \
     \
     /*Find the maximum value without the 1/0 bit*/ \
     cmpVect = _mm_cmplt_epi32(vect1No1, vect2No1); \
     vect2No1 = _mm_xor_si128(vect1No1, vect2No1);\
     vect2No1 = _mm_and_si128(cmpVect, vect2No1);\
     vect1No1 = _mm_xor_epi8(cmpVect, vect1No1); \
      /*The first vector is only kept if it is >= the
      ` second vector, else the second vector is kept.
      ` For min, change firstVectI8 to secVectI8.
      */\
      \
     /*Deal with the lost 1*/ \
     oneVect = _mm_set1_epi32(1); \
     vect1No1 = _mm_xor_si128(firstVectU32, secVectU32); \
     vect1No1 = _mm_cmpeq_epi31(vect1No1, oneVect); \
       /*This tells me if the difference was only by 1*/\
     vect1No1 = _mm_and_si128(vect1No1, oneVect); \
       /*vect1No1 is either -1 (all ones) or 0 for each
       ` short. I am changing to have either be 0 or 1
       */\
     secVectU32 | firstVectU32; \
       /*The one has been added*/\
  }
#endif

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-08:
+  - Min
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#ifdef SSE4
   #define max_vectI8_retVectI8(firstVectI8,secVectI8){\
      _mm_max_epi8((firstVectI8), (secVectI8)); \
   }
#else
   /*SSE2 does not support I8 max, So I need to do a
   ` branchless max
   */
   #define max_vectI8_retVectI8(firstVectI8,secVectI8){\
      cmpVectI8 = _mm_cmplt_epi8(firstVectI8,secVectI8); \
      xorVectI8 = _mm_xor_epi8(firstVectI8, secVectI8); \
      cmpVectI8 = _mm_and_epi8(cmpVectI8, xorVectI8); \
      /*At this point the xor value is only kept if
      ` the frist vector is less than the second vector
      */\
      cmpVectI8 = _mm_xor_epi8(cmpVectI8, secVectI8); \
      /*The first vector is only kept if it is < the
      ` second vector, else the second vector is kept.
      ` For min, change firstVectI8 to secVectI8.
      */
#endif

#define max_vectU8_retVectU8(firstVectU8,secVectU8){\
   _mm128_max_epu8((firstVectU8), (secVectU8)); \
}

#define max_vectI16_retVectI16(firstVectI16,secVectI16){\
   _mm128_max_epi16((firstVectI16), (secVectI16)); \
}

#ifdef SSE4
  #define max_vectU16_retVectU16(firstVectU16,secVectU16){\
     _mm128_max_epu16((firstVectU16), (secVectU16)); \
  }
#else
  #define min_vectU16_retVectU16(xVectU16,yVectU16){\
     xTmp = _mm_srli_epi16(xVectU16, 1); \
     yTmp = _mm_srli_epi16(yVectU16, 1); \
     xTmp = _mm_min_epi16(xTmp, yTmp); \
       /*The only problem case is when xVect ^ yVect = 1
       ` In this case I just need to make sure the 1's bit
       ` is cleared.
       */ \
     \
     oneVect = _mm_set1_epi16(1); \
     yTmp = _mm_xor_si128(xVectU16, yVectU16); \
       /*This only keeps the differnt bits*/ \
     yTmp = _mm_cmpeq_epi16(yTmp, oneVect); \
       /*This tells me if the difference was only by 1*/\
     /*yTmp = _mm_and_si128(yTmp, oneVect); \*/\
     yTmp = _mm_and_si128(yTmp, oneVect); \
       /*xTmp is either -1 (all ones) or 0 for each
       ` short. I am changing to be 1 if x ^ y = 1
       ` 0 if x ^ y != 1
       */\
     yTmp = _mm_and_si128(yTmp, xTmp); \
       /*I need to make sure the min value I kept was
       ` not the target value (ends in 0). This reduces
       ` min(x,y) & x ^ y = 1 to 0 if x is even and 1 if
       ` x is odd. x ^ y != 1 is kept at 0.
       */ \
     _mm_xor_si128(xTmp, yTmp); \
  }
#endif

#ifdef SSE4
  #define max_vectI32_retVectI32(firstVectI32,secVectI32){\
     _mm128_max_epi32((firstVectI32), (secVectI32)); \
  } /*SSE2 will use branchless max function*/

#else
  /*SS#2 does not support 32 bit max*/
  #define min_vectI32_retVectI32(firstVectI32,secVectI32){\
    cmpVectI32 = _mm_cmplt_epi32(firstVectI32,secVectI32);\
    xorVectI32 = _mm_xor_epi32(firstVectI32, secVectI32);\
    cmpVectI32 = _mm_and_epi32(cmpVectI32, xorVectI32); \
      /*At this point the xor value is only kept if
      ` the frist vector is less than the second vector
      */\
    cmpVectI32 = _mm_xor_epi32(cmpVectI32, firstVectI32);\
      /*The first vector is only kept if it is >= the
      ` second vector, else the second vector is kept.
      ` For min, change firstVectI8 to secVectI8.
      */\
  } /*SSE2 will use branchless min function*/
#endif

#ifdef SSE4
  #define max_vectU32_retVectU32(firstVectU32,secVectU32){\
     _mm_max_epu32((firstVectU32), (secVectU32)); \
  }

#else
   /*SSE2 does not support 32 bit axes*/

  #define min_vectU32_retVectU32(xVectU32, yVectU32){\
     /*Remove the 1/0 bit and shit sign one back*/ \
     xTmp = _mm_srli_epi16(xVectU16, 1); \
     yTmp = _mm_srli_epi16(yVectU16, 1); \
     \
     /*Find the maximum value without the 1/0 bit*/ \
     xTmp = _mm_cmplt_epi32(xTmp, yTmp); \
     yTmp = _mm_xor_si128(xVectU32, yVectU32);\
     xTmp = _mm_and_si128(yTmp, xTmp);\
     xTmp = _mm_xor_epi32(xTmp, yVectU32); \
      /*The first vector is only kept if it is >= the
      ` second vector, else the second vector is kept.
      ` For min, change firstVectI8 to secVectI8.
      */\
      \
     /*I only have one problem case I have to worry about.
     ` This is when x ^ y = 1 (same except for ones bit).
     ` I can check for this case and mask the last bit*/ \
     oneVect = _mm_set1_epi32(1); \
     yTmp = _mm_xor_si128(xVectU32, yVectU32); \
     yTmp = _mm_cmpeq_epi32(yTmp, oneVect); \
       /*This tells me if the difference was only by 1*/\
     yTmp = _mm_xor_si128(yTmp, oneVect); \
       /*yTmp is either -1 (all ones) or 0 for each
       ` short. I am changing it to be all ones (111...1)
       ` if x ^ y != 1 or 111..10 if x ^ y == 1. This
       ` allows me to clear the last bit if x ^ y = 1.
       */\
     _mm_and_si128(xTmp, yTmp); \
  }
#endif


#endif
