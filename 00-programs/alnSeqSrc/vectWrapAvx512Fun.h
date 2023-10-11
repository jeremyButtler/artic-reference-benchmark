#ifndef VECTWRAPAVX512_H
#define VECTWRAPAVX512_h 

#include <immintrin.h>

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-01:
^  - Definitions and variable declerations
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define vectorBits 512
#define vectorBytes 64

typedef __m512i vectI8;  /*vector of bytes (8bits)*/
typedef __m512i vectI16; /*vector of shorts (16 bits)*/
typedef __m512i vectI32; /*vector of ints (32 bits)*/

typedef __m512i vectU8;  /*Prefer unsigned bytes*/
typedef __m512i vectU16; /*Prefer unsigned shorts*/
typedef __m512i vectU32; /*Prefer unsinged ints*/

typedef __mask64 mask8;  /*mask of 8 bit values*/
typedef __mask32 mask16; /*Mask of 16 bit values*/
typedef __mask16 mask32; /*mask of 32 bit values*/

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
  _mm512_load_si512((__m512i *) (arrayI8)); \
}

#define load_I16Ary_retVectI16(arrayI16) { \
  _mm512_load_si512((__m512i *) (arrayI16));\
}

#define load_I32Ary_retVectI32(arrayI32) { \
  _mm512_load_si512((__m512i *) (arrayI32));\
}



#define load_U8Ary_retVectU8(arrayU8) { \
  _mm512_load_si512((__m512i *) (arrayU8)); \
}

#define load_U16Ary_retVectU16(arrayU16) { \
  _mm512_load_si512((__m512i *) (arrayU16));\
}

#define load_U32Ary_retVectU32(arrayU32) { \
  _mm512_load_si512((__m512i *) (arrayU32));\
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-01:
+  - Unaligned loading
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define loadu_I8Ary_retVectI8(arrayI8) { \
  _mm512_loadu_si512((__m512i *) (arrayI8)); \
}

#define loadu_I16Ary_retVectI16(arrayI16) { \
  _mm512_loadu_si512((__m512i *) (arrayI16));\
}

#define loadu_I32Ary_retVectI32(arrayI32) { \
  _mm512_loadu_si512((__m512i *) (arrayI32));\
}



#define loadu_U8Ary_retVectU8(arrayU8) { \
  _mm512_loadu_si512((__m512i *) (arrayU8)); \
}

#define loadu_U16Ary_retVectU16(arrayU16) { \
  _mm512_loadu_si512((__m512i *) (arrayU16));\
}

#define loadu_U32Ary_retVectU32(arrayU32) { \
  _mm512_loadu_si512((__m512i *) (arrayU32));\
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-03:
+  - make zero vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define zero_retVectI8() (_mm512_setzero_si512();)
#define zero_retVectI16()(_mm512_setzero_si512();)
#define zero_retVectI32()(_mm512_setzero_si512();)

#define zero_retVectU8() (_mm512_setzero_si512();)
#define zero_retVectU16()(_mm512_setzero_si512();)
#define zero_retVectU32()(_mm512_setzero_si512();)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-04:
+  - Make vectors of one element
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define set1_I8_retVectI8(valC) (mm512_set1_epi8((valC));)

#define set1_I16_retVectI16(valI16) \
   (mm512_set1_epi16((valI16));)

#define set1_I32_retVectI32(valI32) \
   (_mm512_set1_epi32((valI32));)


#define set1_U8_retVectU8(valUC)(mm512_set1_epi8((valUC));)

#define set1_U16_retVectU16(valC) \
   (mm512_set1_epi16((valU16));)

#define set1_U32_retVectU32(valU32) \
   (_mm512_set1_epi32((valU32));)


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-05:
+  - Insert an element into an vector
+  - https://stackoverflow.com/questions/58303958/how-to-implement-16-and-32-bit-integer-insert-and-extract-operations-with-avx-51
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define insert_I8_retVectI8(insI8, posI){\
   _mm512_mask_set1_epi8(insVectI8,1UL<<posI,insI8);\
}

#define insert_I16_retVectI16(insI16, posI){\
   _mm512_mask_set1_epi16(insVectI16,1UL<<posI,insI16);\
}         

#define insert_I32_retVectI32(insI32, posI){\
   _mm512_mask_set1_epi32(insVectI32,1UL<<posI,insI32);\
}


#define insert_U8_retVectU8(insU8, posU){\
   _mm512_mask_set1_epi8(insVectU8,1UL<<posU,insU8);\
}

#define insert_U16_retVectU16(insU16, posU){\
   _mm512_mask_set1_epi16(insVectU16,1UL<<posU,insU16);\
}         

#define insert_U32_retVectU32(insU32, posU){\
   _mm512_mask_set1_epi32(insVectU32,1UL<<posU,insU32);\
}

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
  _mm512_store_si512((__m512i *) (array), (inVectI8)); \
}

#define store_vectI16_retAryI16(array, inVectI16) { \
  _mm512_store_si512((__m512i *) (array), (inVectI16)); \
}

#define store_vectI32_retAryI32(array, inVectI32) { \
  _mm512_store_si512((__m512i *) (array), (inVectI32)); \
}


#define store_vectU8_retAryU8(array, inVectU8) { \
  _mm512_store_si512((__m512i *) (array), (inVectU8)); \
}

#define store_vectU16_retAryU16(array, inVectU16) { \
  _mm512_store_si512((__m512i *) (array), (inVectU16)); \
}

#define store_vectU32_retAryU32(array, inVectU32) { \
  _mm512_store_si512((__m512i *) (array), (inVectU32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-02:
+  - Store masks into longs
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define store_mask8_retUL(retUL, mask8) { \
  _store_mask64((__mask64 *) (retUL), (mask8)); \
} /*Mask is 8 bit integers*/

#define store_mask16_retUL(retUL, mask16) { \
  _store_mask64((__mask32 *) (retUL), (mask16)); \
} /*Mask is 16 bit integers*/

#define store_mask32_retUL(retUL, mask32) { \
  _store_mask16((__mask16 *) (retUL), (mask32)); \
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
   _mm512_cmpeq_epi8_mask((vectOneI8), (vectTwoI8)); \
}

#define cmpeq_vectU8_retMask8(vectOneI8, vectTwoI8) { \
   _mm512_cmpeq_epi8_mask((vectOneI8), (vectTwoI8)); \
}

#define cmpgt_vectI8_retMask8(vectOneI8, vectTwoI8) {\
   _mm512_cmpgt_epi8_mask((vectOneI8), (vectTwoI8)); \
}

#define cmplt_vectI8_retMask8(vectOneI8, vectTwoI8) { \
   _mm512_cmplt_epi8_mask((vectOneI8), (vectTwoI8)); \
}

/*16-bit comparisions; returns a mask16*/

#define cmpeq_vectI16_retMask16(vectOneI16, vectTwoI16) { \
   _mm512_cmpeq_epi16_mask((vectOneI16), (vectTwoI16)); \
}

#define cmpeq_vectU16_retMask16(vectOneI16, vectTwoI16) { \
   _mm512_cmpeq_epi16_mask((vectOneI16), (vectTwoI16)); \
}

#define cmpgt_vectI16_retMask16(vectOneI16, vectTwoI16) {\
   _mm512_cmpgt_epi16_mask((vectOneI16), (vectTwoI16)); \
}

#define cmplt_vectI16_retMask16(vectOneI16, vectTwoI16) { \
   _mm512_cmplt_epi16_mask((vectOneI16), (vectTwoI16)); \
}

/*32-bit comparisions; returns a mask32*/

#define cmpeq_vectI32_retMask32(vectOneI32, vectTwoI32) { \
   _mm512_cmpeq_epi32_mask((vectOneI32), (vectTwoI32)); \
}

#define cmpeq_vectU32_retMask32(vectOneI32, vectTwoI32) { \
   _mm512_cmpeq_epi32_mask((vectOneI32), (vectTwoI32)); \
}

#define cmpgt_vectI32_retMask32(vectOneI32, vectTwoI32) {\
   _mm512_cmpgt_epi32_mask((vectOneI32), (vectTwoI32)); \
}

#define cmpLt_vectI32_retMask32(vectOneI32, vectTwoI32) { \
   _mm512_cmplt_epi32_mask((vectOneI32), (vectTwoI32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-02:
+  - Fix differences in population counts (total 1's)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* This is here to support SSE2, AVX2, and NEON. 
`  The returned mask in AVX512 has only one bit per
`  data type in comparison
*/
#define fix_mask8_popcount(inUL) ((inUL))
#define fix_mask16_popcount(inUL) ((inUL))
#define fix_mask32_popcount(inUL) ((inUL))

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
   _mm512_andnot_si512((notVectI8), (inVectI8)); \
}

#define andNot_vectU8_retVectU8(notVectU8, inVectU8) { \
   _mm512_andnot_si512((notVectU8), (inVectU8)); \
}

#define andNot_vectI16_retVectI16(notVectI16, inVectI16) {\
   _mm512_andnot_si512((notVectI16), (inVectI16)); \
}

#define andNot_vectU16_retVectU16(notVectU16, inVectU16) {\
   _mm512_andnot_si512((notVectU16), (inVectU16)); \
}

#define andNot_vectI32_retVectI32(notVectI32, inVectI32) {\
   _mm512_andnot_si512((notVectI32), (inVectI32)); \
}

#define andNot_vectU32_retVectU32(notVectU32, inVectU32) {\
   _mm512_andnot_si512((notVectU32), (inVectU32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-02:
+  - and functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define and_vectI8_retVectI8(vectOneI8, vectTwoI8) { \
   _mm512_and_si512((vectOneI8), (vectTwoI8)); \
}

#define and_vectU8_retVectU8(vectOneU8, vectTwoU8) { \
   _mm512_and_si512((vectOneU8), (vectTwoU8)); \
}

#define and_vectI16_retVectI16(vectOneI16, vectTwoI16) { \
   _mm512_and_si512((vectOneI16), (vectTwoI16)); \
}

#define and_vectU16_retVectU16(vectOneU16, vectTwoU16) { \
   _mm512_and_si512((vectOneU16), (vectTwoU16)); \
}

#define and_vectI32_retVectI32(vectOneI32, vectTwoI32) { \
   _mm512_and_si512((vectOneI32), (vectTwoI32)); \
}

#define and_vectU32_retVectU32(vectOneU32, vectTwoU32) { \
   _mm512_and_si512((vectOneU32), (vectTwoU32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-03:
+  - or functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define or_vectI8_retVectI8(vectOneI8, vectTwoI8) { \
   _mm512_or_si512((vectOneI8), (vectTwoI8)); \
}

#define or_vectU8_retVectU8(vectOneU8, vectTwoU8) { \
   _mm512_or_si512((vectOneU8), (vectTwoU8)); \
}

#define or_vectI16_retVectI16(vectOneI16, vectTwoI16) { \
   _mm512_or_si512((vectOneI16), (vectTwoI16)); \
}

#define or_vectU16_retVectU16(vectOneU16, vectTwoU16) { \
   _mm512_or_si512((vectOneU16), (vectTwoU16)); \
}

#define or_vectI32_retVectI32(vectOneI32, vectTwoI32) { \
   _mm512_or_si512((vectOneI32), (vectTwoI32)); \
}

#define or_vectU32_retVectU32(vectOneU32, vectTwoU32) { \
   _mm512_or_si512((vectOneU32), (vectTwoU32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-03:
+  - or functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define xor_vectI8_retVectI8(vectOneI8, vectTwoI8) { \
   _mm512_xor_si512((vectOneI8), (vectTwoI8)); \
}

#define xor_vectU8_retVectU8(vectOneU8, vectTwoU8) { \
   _mm512_xor_si512((vectOneU8), (vectTwoU8)); \
}

#define xor_vectI16_retVectI16(vectOneI16, vectTwoI16) { \
   _mm512_xor_si512((vectOneI16), (vectTwoI16)); \
}

#define xor_vectU16_retVectU16(vectOneU16, vectTwoU16) { \
   _mm512_xor_si512((vectOneU16), (vectTwoU16)); \
}

#define xor_vectI32_retVectI32(vectOneI32, vectTwoI32) { \
   _mm512_xor_si512((vectOneI32), (vectTwoI32)); \
}

#define xor_vectU32_retVectU32(vectOneU32, vectTwoU32) { \
   _mm512_xor_si512((vectOneU32), (vectTwoU32)); \
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
    _mm512_srli_epi8((inVectI8), (numBitsI)); \
}

#define sr_vectU8_retVectU8(inVectU8, numBitsI){\
    _mm512_srli_epi8((inVectU8), (numBitsI)); \
}

#define sr_vectI16_retVectI16(inVectI16, numBitsI){\
    _mm512_srli_epi16((inVectI16), (numBitsI)); \
}

#define sr_vectU16_retVectU16(inVectU16, numBitsI){\
    _mm512_srli_epi16((inVectU16), (numBitsI)); \
}

#define sr_vectI32_retVectI32(inVectI32, numBitsI){\
    _mm512_srli_epi32((inVectI32), (numBitsI)); \
}

#define sr_vectU32_retVectU32(inVectU32, numBitsI){\
    _mm512_srli_epi32((inVectU32), (numBitsI)); \
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
    _mm512_srli_si512((inVectI8), (numBytesI)); \
}

#define srvect_vectU8_retVectU8(inVectU8, numBytesI){\
    _mm512_srli_si512((inVectU8), (numBytesI)); \
}

#define srvect_vectI16_retVectI16(inVectI16, numBytesI){\
    _mm512_srli_si512((inVectI16)); \
}

#define srvect_vectU16_retVectU16(inVectU16, numBytesI){\
    _mm512_srli_si512((inVectU16), (numBytesI)); \
}

#define srvect_vectI32_retVectI32(inVectI32, numBytesI){\
    _mm512_srli_si512((inVectI32), (numBytesI)); \
}

#define srvect_vectU32_retVectU32(inVectU32, numBytesI){\
    _mm512_srli_si512((inVectU32), (numBytesI)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-03:
+  - Shift numbers in vectors left by x bits
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sl_vectI8_retVectI8(inVectI8, numBitsI){\
    _mm512_slli_epi8((inVectI8), (numBitsI)); \
}

#define sl_vectU8_retVectU8(inVectU8, numBitsI){\
    _mm512_slli_epi8((inVectU8), (numBitsI)); \
}

#define sl_vectI16_retVectI16(inVectI16, numBitsI){\
    _mm512_slli_epi16((inVectI16), (numBitsI)); \
}

#define sl_vectU16_retVectU16(inVectU16, numBitsI){\
    _mm512_slli_epi16((inVectU16), (numBitsI)); \
}

#define sl_vectI32_retVectI32(inVectI32, numBitsI){\
    _mm512_slli_epi32((inVectI32), (numBitsI)); \
}

#define sl_vectU32_retVectU32(inVectU32, numBitsI){\
    _mm512_slli_epi32((inVectU32), (numBitsI)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-04:
+  - Shift a vector left by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sllvect_vectI8_retVectI8(inVectI8, numBytesI){\
   _mm512_slli_si512((inVectI8), (numBytesI)); \
}

#define sllvect_vectU8_retVectU8(inVectU8, numBytesI){\
   _mm512_slli_si512((inVectU8), (numBytesI)); \
}

#define sllvect_vectI16_retVectI16(inVectI16, numBytesI){\
   _mm512_slli_si512((inVectI16), (numBytesI)); \
}

#define sllvect_vectU16_retVectU16(inVectU16, numBytesI){\
   _mm512_slli_si512((inVectU16), (numBytesI)); \
}

#define sllvect_vectI32_retVectI32(inVectI32, numBytesI){\
   _mm512_slli_si512((inVectI32), (numBytesI)); \
}

#define sllvect_vectU32_retVectU32(inVectU32, numBytesI){\
   _mm512_slli_si512((inVectU32), (numBytesI)); \
}

/*********************************************************\
* Sec-02 Sub-06:
*  - Math functions
*  o sec-02 sub-06 cat-01:
*    - addition max 64 bit; [sse2-avx512]; no epu
*  o sec-02 sub-06 cat-02:
*    - staturation addition; max 16 bit [sse2-avx512]; epu
*  o sec-02 sub-06 cat-03:
*    - subtraction; max 64bit [sse2-avx512]; no epu
*  o sec-02 sub-06 cat-04:
*    - staturation subtraction; max 16 bit; [sse2-avx512];u
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-01:
+  - addition max 64 bit; [sse2-avx512]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define add_vectI8_retVectI8(vectOneI8, vectTwoI8) {\
   _mm512_add_epi8((vectOneI8), (vectTwoI8)); \
}

#define add_vectI16_retVectI16(vectOneI16, vectTwoI16) {\
   _mm512_add_epi16((vectOneI16), (vectTwoI16)); \
}

#define add_vectI32_retVectI32(vectOneI32, vectTwoI32) {\
   _mm512_add_epi32((vectOneI32), (vectTwoI32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-02:
+  - staturation addition; max 16 bit [sse2 to avx512]; epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define addSat_vectU8_retVectU8(vectOneU8, vectTwoU8) {\
   _mm512_adds_epu8((vectOneU8), (vectTwoU8)); \
}

#define addSat_vectU16_retVectU16(vectOneU16, vectTwoU16){\
   _mm512_adds_epu16((vectOneU16), (vectTwoU16)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-03:
+  - subtraction; max 64bit [sse2-avx512]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sub_vectI8_retVectI8(vectOneI8, vectTwoI8) { \
    _mm512_sub_epi8((vectOneI8), (vectTwoI8)); \
}

#define sub_vectI16_retVectI16(vectOneI16, vectTwoI16) { \
    _mm512_sub_epi16((vectOneI16), (vectTwoI16)); \
}

#define sub_vectI32_retVectI32(vectOneI32, vectTwoI32) { \
    _mm512_sub_epi32((vectOneI32), (vectTwoI32)); \
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-04:
+  - staturation subtraction; max 16 bit; [sse2-avx512]); u
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define subSat_vectI8_retVectI8(vectOneI8, vectTwoI8) {\
   _mm512_subs_epi8((vectOneI8), (vectTwoI8)); \
}

#define subSat_vectU8_retVectU8(vectOneU8, vectTwoU8) {\
   _mm512_subs_epu8((vectOneU8), (vectTwoU8)); \
}

#define subSat_vectI16_retVectI16(vectOneI16, vectTwoI16) {\
   _mm512_subs_epi16((vectOneI16), (vectTwoI16)); \
}

#define subSat_vectU16_retVectU16(vectOneU16, vectTwoU16) {\
   _mm512_subs_epu16((vectOneU16), (vectTwoU16)); \
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
   _mm512_movm_epi8((inMask8)); \
}

#define mask16_retVectI16(inMask16){ \
   _mm512_movm_epi16((inMask16)); \
}

#define mask32_retVectI32(inMask32){ \
   _mm512_movm_epi32((inMask32)); \
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-07:
+  - Max
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define max_vectI8_retVectI8(firstVectI8,secVectI8){\
   _mm512_max_epi8((firstVectI8), (secVectI8)); \
} /*SSE2 will use branchless max function*/

#define max_vectU8_retVectU8(firstVectU8,secVectU8){\
   _mm512_max_epu8((firstVectU8), (secVectU8)); \
} /*SSE2 will use vector max function*/

#define max_vectI16_retVectI16(firstVectI16,secVectI16){\
   _mm512_max_epi16((firstVectI16), (secVectI16)); \
} /*SSE2 will use vector max function*/

#define max_vectU16_retVectU16(firstVectU16,secVectU16){\
   _mm512_max_epu16((firstVectU16), (secVectU16)); \
} /*SSE2 will use banchless max function*/
/*  See min for how this might work for SSE2*/

#define max_vectI32_retVectI32(firstVectI32,secVectI32){\
   _mm512_max_epi32((firstVectI32), (secVectI32)); \
} /*SSE2 will use branchless max function*/

#define min_vectU32_retVectU32(firstVectU32,secVectU32){\
   _mm512_min_epu32((firstVectU32), (secVectU32)); \
} /*SSE2 will use branchless min function*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-08:
+  - Min
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define min_vectI8_retVectI8(firstVectI8,secVectI8){\
   _mm512_min_epi8((firstVectI8), (secVectI8)); \
} /*SSE2 will use branchless min function*/

#define min_vectU8_retVectU8(firstVectU8,secVectU8){\
   _mm512_min_epu8((firstVectU8), (secVectU8)); \
} /*SSE2 will use vector min function*/

#define min_vectI16_retVectI16(firstVectI16,secVectI16){\
   _mm512_min_epi16((firstVectI16), (secVectI16)); \
} /*SSE2 will use vector min function*/

#define min_vectU16_retVectU16(firstVectU16,secVectU16){\
   _mm512_min_epu16((firstVectU16), (secVectU16)); \
} /*SSE2 will use banchless min function*/
/*//This is how an SSE2 min16U might look
  // 9 OP + three extra vector registers
  // A branchless non-vect min would be 4 OP + comparison
  // which would take a least 3 operations

  // 16 bit (8 OP on 8 items)
  oneVect = set1(1);
  vect1No1 = vect1 >> 1;
  vect2No1 = vect2 >> 1;
  vec2No1 = min(vect1No1, vect2No1);
    // The only problem case is when |vect1 - vect2| = 1.
  vect1No1 = vect1 ^ vect2;
  vect1No1 = ifeq(vect1No1, oneVect);
  vect1No1 = vect1No1 ^ oneVect;
  vect2No1 & vec1No1;

  // 32 bit (unsiged) (11 OP on 4 items) + 1 load
  // Without vectors (branchless), this would be 5 OP per
  // item, but have little load cost (20OP total).

  // Remove the 1/0 bit and shit sign one back
  vect1No1 = vect1 >> 1;
  vect2No1 = vect2 >> 1;

  // Find the minimum value without the 1/0 bit
  vect2No1 = iflt(vect1No1, vect2No1);
  vect1No1 = vect1 ^ vect2;
  vect1No1 = vect1No1 & vect2No1;
  vect1No1 = vect1 ^ vect2No1; // Make sure this is min

  // Account for the 1 bit
  oneVect = set1(1);
  vect2No1 = vect1 ^ vect2;
  vect2No1 = ifeq(vect1No1, oneVect);
  vect2No1 = vect2No1 ^ oneVect;
    // This will allow me to remove the 1 from the values
    // that are one differnt. This also avoids the error
    // of subtracing a negative value
  minVect & vec2No1;
*/


#define min_vectI32_retVectI32(firstVectI32,secVectI32){\
   _mm512_min_epi32((firstVectI32), (secVectI32)); \
} /*SSE2 will use branchless min function*/

#define min_vectU32_retVectU32(firstVectU32,secVectU32){\
   _mm512_min_epu32((firstVectU32), (secVectU32)); \
} /*SSE2 will use branchless min function*/

#endif
