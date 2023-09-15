/*#########################################################
# Name: twoBitArrays
# Use:
#   o Holds functions to handle two bit arrays
# Includes:
# C Standard Includes:
#   - <stdint.h>
#########################################################*/

#ifndef TWOBITARRAYS_H
#define TWOBITARRAYS_H

#include <stdint.h>
#include <stdlib.h>

/*Get the negative flag for a character*/
/*--------------------------------------------------------\
| Use:
|  - This is just to avoid (negChar < 0), which can have
|    branching issues on some cpus (likely very old
|    machines). It takes the same amount of time. The
|    (sizeof(unsigned char) << 3) - 1 statement will be
|    converted to a constant during compile time.
| Output:
|  - Returns:
|    o 1: if negChar was negative
|    o 0: if negChar was positive
\--------------------------------------------------------*/
static inline char getCNegBit(
   char negChar /*Character to get the negative flag from*/
){
   return 
         (((unsigned char) negChar)
      >> (((sizeof(unsigned char) << 3) - 1)));
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-01 getTwoBitElm:
'     o Get an element from a two bit array
'  - fun-02 twoBitAryShiftBytsForNewElm:
'     o Make room in unit8_t for two more two-bit elements
'  - fun-03 twoBitMvToNextElm:
'     o Moves to the next element in a two-bit array
'  - fun-04 twoBitMvForXElm:
'     o Moves two bit array pointer by x to next element
'  - fun-05 twoBitMvBackOneElm:
'     o Moves back one element in a 2-bit array
'  - fun-06 twoBitMvBackXElm:
'     o Moves back x elements in a 2-bit array
'  - fun-08 changeTwoBitElm:
'     o Changes a single two bit value in a two bit array.
'  - fun-09 blankTwoBitLimb:
'     o Sets a uint8_t (a limb) in a two bit array to 0.
'       Each limb holds four two bit elments.
'  - fun-10 twoBitMvToNextLimb:
'     o Moves to start of the next limb (uint8_t) in two
'       bit array
'  - fun-11 twoBitMvToLastLimb:
'     o Moves to start of the previous limb (uint8_t) in
'       two bit array
'  - fun-12 makeTwoBit:
'     o Make a two bit array struct
'  - fun-13 cpTwoBitPos:
'     o copy pointers of cpTwoBitST to dupTwoBitST
'  - fun-14 freeTwoBit:
'     o Frees a two bit array
'  - fun-15 twoBitMvXElmFromStart:
'     o Change the two bit array postion to x elements
'       from the starting position
'  - fun-16 twoBitGetLen:
'     o Get the length of a two-bit array
'  - fun-17 twoBitGetIndex:
'    - Returns the index of the current two bit element
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Struct-01: twoBitArry
|  - Holds the two bit arrays and pointers to manage it
\--------------------------------------------------------*/
typedef struct twoBitAry
{ // twoBitArry structure
  uint8_t *firstLimbUCPtr; // First limb in two bit array
  uint8_t *limbOnUCPtr;    // Limb currently working on
  int8_t elmOnC;         // element on in the limb
  unsigned long lenAryUL;  // Number of limbs in the array
}twoBitAry;

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Two bits of interest from the two bit array
| Note:
|  - Each limb has four two bit elements
\--------------------------------------------------------*/
static inline uint8_t getTwoBitElm(
  struct twoBitAry *twoBitST // Array to get element from
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: getElmFromToBitUCAry
   '  - Get an element from a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*This was a branchless attempt. It had little affect
      on the speed. For this to work, the case must
      be reversed
    */
    char shiftC = twoBitST->elmOnC << 1;
    return
        (*twoBitST->limbOnUCPtr & (3 << shiftC)) >> shiftC;
} // getTwoBitElm

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitST to point to next element in two-bit array
\--------------------------------------------------------*/
static inline void twoBitMvToNextElm(
  struct twoBitAry *twoBitST
    // Two bit array to change index
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-1 Sub-1: twoBitMvToNextElm
   '  - Moves to the next element in a two-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   ++twoBitST->elmOnC;
   twoBitST->limbOnUCPtr += (twoBitST->elmOnC >> 2);
   twoBitST->elmOnC &= 3;
   return;
} // twoBitMvToNextElm

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements ahead
\--------------------------------------------------------*/
static inline void twoBitMvForXElm(
  struct twoBitAry *twoBitST,// Array to change index for
  unsigned long shiftByUL    // Number elements to shift by
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-1 Sub-1: twoBitMvForXElm
   '  - Moves forward in two bit array by x elements
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Get whole shifts to perform
   shiftByUL += twoBitST->elmOnC;
   twoBitST->limbOnUCPtr += (shiftByUL >> 2);
   twoBitST->elmOnC = (shiftByUL & 3);
   return;
} // twoBitMvForXElm

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitST to point to the previous element
\--------------------------------------------------------*/
static inline void twoBitMvBackOneElm(
  struct twoBitAry *twoBitST // array to move back in
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-1 Sub-1: twoBitMvBackOneElm
   '  - Moves back one element in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*Does not alter speed from a switch*/
   --twoBitST->elmOnC;
   twoBitST->limbOnUCPtr -= getCNegBit(twoBitST->elmOnC);
   /*Adds an extra 1 if I am moving back extra elements*/
   twoBitST->elmOnC = 3 & twoBitST->elmOnC;
    /*-1 goes to 3, and all other values are 3 or less and
    ` so remain the same
    */
   return; 
} // twoBitMvToNextElm

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBotST to point to X elements back    
\--------------------------------------------------------*/
static inline void twoBitMvBackXElm(
  struct twoBitAry *twoBitST, // To bit array to move back
  long shiftL  // number elements to shift back
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-1 Sub-1: twoBitMvBackXElm
   '  - Moves back X elements back in a 2-bit array
   '  - THIS HAS NO BEEN TESTED
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*Figure out how many elements I am moving back*/
   twoBitST->elmOnC = twoBitST->elmOnC - (shiftL & 3);

   twoBitST->limbOnUCPtr -=
      ((shiftL >> 2) + getCNegBit(twoBitST->elmOnC));
   /*Adds an extra 1 if I am moving back extra elements*/

   twoBitST->elmOnC =
        twoBitST->elmOnC
      ^ (  -(getCNegBit(twoBitST->elmOnC))
         & (twoBitST->elmOnC ^ (twoBitST->elmOnC + 4))
        );
     /*This is a minimize function, but instead of
     ` Returning the minimum if negative, it returns
     ` the minimum + 4.
     ` -(twoBitST->elmOnC < 0)
     `   Is -1 (111....) if elmOnC is negative
     `   Is 0 other if elmOnC is not negative.
     ` elmOnC + 4
     `   This is the correct ending index if elmOnC is
     `   negative.
     ` elmOnC ^ (-1 & (elmOnC ^ (elmOnC + 4)))
     `   This returns elmOnC + 4
     ` elmOnC ^ (0 & (elmOnC ^ (elmOnC + 4)))
     `   This returns elmOnC
     */
   return;
} // twoBitMvBackXElm

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o Changes the taget element in twoBitUCArray to
|      input value
\--------------------------------------------------------*/
static inline void changeTwoBitElm(
  struct twoBitAry *twoBitST, //Points to element to change
  uint8_t newValueUC          // New value; only the first
                              // two bits can be set to 1
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: Sec-1 Sub-1: changeTwoBitElm
   '  - Changes a single two bit value in a two bit array.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /* This was a branchless attempt. It had little affect
    `  on the acutal speed and may have had a slight slow
    `  down.*/
    char valC = twoBitST->elmOnC << 1;
    char clearC = 3 << valC;

    valC = newValueUC << valC;
    /*This shifts the new value to its correct position
    ` in the two bit array
    */

    *twoBitST->limbOnUCPtr =
      (*twoBitST->limbOnUCPtr & (~clearC)) | valC;
    /*This clears and then sets the value to its correct
    `  postion:
    ` & (~valC):
    `   ~valC removes
    ` and then
    ` sets the postion to the input value (| valC).
    */

    return;
} // changeTwoBitElm

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    - Sets the current limb in a two-bit array to 0
\--------------------------------------------------------*/
static inline void blankTwoBitLimb(
  struct twoBitAry *twoBitST // two bit array to restart
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: Sec-1 Sub-1: blankTwoBitLimb
   '  - Sets a limb in a two bit array to 0.
   '    Each limb holds four two bit elments.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   twoBitST->limbOnUCPtr = 0;
   twoBitST->elmOnC = 0;
} // blankTwoBitLimb

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o Moves to the next limb in a two bit array
|  - Returns:
|    o 0: For a succesfull move
|    o 1: For memory out of bounds error
\--------------------------------------------------------*/
static inline void twoBitMvToNextLimb(
  struct twoBitAry *twoBitST // two-bit array to move back
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: Sec-1 Sub-1: twoBitMvToNextLimb
   '  - Moves back to start of limb in a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   ++twoBitST->limbOnUCPtr;
   twoBitST->elmOnC = 0;

   return;
} // twoBitMvToNextLimb

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o lastLimbUCPtr to point to the previous limb
\--------------------------------------------------------*/
static inline void twoBitMvToLastLimb(
    struct twoBitAry *twoBitST // array to move back in
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: Sec-1 Sub-1: twoBitMvToLastLimb
   '  - Moves to the start of the previous limb (uint8_t)
   '    in the two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   --twoBitST->limbOnUCPtr;
   twoBitST->elmOnC = 0;

   return;
} // twoBitMvToLastLimb

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    - twoBitAry structure with an unit8_t array of limbs
|    - blankAryBl = 1: returns a blank twoBitAry structer
\--------------------------------------------------------*/
static inline struct twoBitAry * makeTwoBit(
  unsigned long lenArryUL, // Length of new array
  char blankAryBl          // 1: Make a blank strucuter
                           // 0: Make structure with array
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: makeTwoBit
   '  - Make a two bit array struct
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   struct twoBitAry *twoBitST =
     malloc(sizeof(struct twoBitAry));

   if(twoBitST == 0) return 0; // Memory error
   twoBitST->elmOnC = 0;

   switch(blankAryBl)
   { // Switch check if making an array
     case 1:
       twoBitST->firstLimbUCPtr = 0;
       twoBitST->limbOnUCPtr = 0;
       twoBitST->lenAryUL = 0;
       return twoBitST;
   } // Switch check if making an array

   twoBitST->firstLimbUCPtr =
     calloc((lenArryUL >> 2) + 1, sizeof(uint8_t));
     // Each limb (uint8_t) has four elements, so I need
     // to divide the number of elements by 4 (>> 2)

   if(twoBitST->firstLimbUCPtr == 0) return 0;

   twoBitST->lenAryUL = ((lenArryUL >> 2) << 2) + 4;
   twoBitST->limbOnUCPtr = twoBitST->firstLimbUCPtr;

   return twoBitST;
} // makeTwoBitArray

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o dupTwoBitST to hold same pointers/values as
|      cpTwoBitST
\--------------------------------------------------------*/
static inline void cpTwoBitPos(
  struct twoBitAry *cpTwoBitST,  // Structer to copy
  struct twoBitAry *dupTwoBitST  // Struct to copy to
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: Sec-1 Sub-1: cpTwoBitPos
   '  - copy pointers of cpTwoBitST to dupTwoBitST
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   dupTwoBitST->firstLimbUCPtr =cpTwoBitST->firstLimbUCPtr;
   dupTwoBitST->limbOnUCPtr = cpTwoBitST->limbOnUCPtr;
   dupTwoBitST->elmOnC = cpTwoBitST->elmOnC;
   dupTwoBitST->lenAryUL = cpTwoBitST->lenAryUL;

   return;
} // cpTwoBitPos

/*--------------------------------------------------------\
| Output:
|  Frees: the input towBitAry strucuter
\--------------------------------------------------------*/
static inline void freeTwoBit(
  struct twoBitAry *stToFree, // Two bit array to free
  char twoBitOnStackBl,       // 0: free stToFree
                              // 1: stToFree on stack
  char doNotFreeAryBl         // 0: Free everthing
      // 1: Do not free the structer in the array. Only
      //    use this if this stuct is a copy
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-14 TOC: Sec-1 Sub-1: freeTwoBit
   '  - Frees a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(doNotFreeAryBl)
   { // switch, check if freeing the array
     case 0:
       if(stToFree->firstLimbUCPtr != 0)
         free(stToFree->firstLimbUCPtr);
   } // switch, check if freeing the array
   
   switch(twoBitOnStackBl)
   { // switch; check if freeing the twobit structer
     case 1: free(stToFree);
   } // switch; check if freeing the twobit structer

   return;
} // freeTwoBit

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements after the start of
|      the array
\--------------------------------------------------------*/
static inline void twoBitMvXElmFromStart(
  struct twoBitAry *twoBitST,// Array to change index for
  unsigned long shiftByUL     // Number elements to shift by
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-15 TOC: Sec-1 Sub-1: twoBitMvXElmFromStart
   '  - Change the two bit array postion to x elements
   '    from the starting position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Get whole shifts to perform
   twoBitST->limbOnUCPtr =
     twoBitST->firstLimbUCPtr + (shiftByUL >> 2);

   twoBitST->elmOnC = (shiftByUL & 3);
   return;
} // twoBitMvXElmFromStart

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Length of the two bit array
\--------------------------------------------------------*/
static inline unsigned long twoBitGetLen(
  struct twoBitAry *twoBitST
    // two-bit array to get length for
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-16 TOC: Sec-1 Sub-1: twoBitGetLen
   '  - Get the length of a two-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   return twoBitST->lenAryUL << 2;
} // twoBitGetLen

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o The index of the elemen on in the twoBitST struct
\--------------------------------------------------------*/
static inline unsigned long twoBitGetIndex(
  struct twoBitAry *twoBitST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-17 TOC: Sec-01: twoBitGetIndex
   '  - Returns the index of the current two bit element
   '    on
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  return (
     ((twoBitST->limbOnUCPtr-twoBitST->firstLimbUCPtr) <<2)
    + twoBitST->elmOnC
  ); // Return the index
} // twoBitGetIndex

#endif
