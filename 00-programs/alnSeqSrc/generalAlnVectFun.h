#include "vectorWrap.h"

/*Deltions are the last best scores*/
/*snps are the previous rounds ins vector (not updated)*/
/*ins are the new vector*/
/*gapVect is 0 if no gap selected & > 0 if gap selected*/
/*By gap I mean that I am not recording the direction*/
/*After this, only the ins vector needs to be reloaded
` for the next round*/
#define vectI16MaxScoreAndGap( \
   insScVectI16, \
   snpScVectI16, \
   delScVectI16, /*has best scores*/\
   gapVectI16, \ /*This is also used as a temp vector*/\
   gapOpenVectI16, \ /*Gap opening penalty*/\
   gapExtendVectI16 \ /*Gap extension penalty*/\
) { \
   mmMaxI16((gapVectI16),(insScVectI16),(snpScVectI16));\
   mmMaxI16((delScVectI16),(gapScVectI16),(delScVectI16));\
   mmCmpEqVectI16(\
      (gapVectI16),\
      (delScVectI16),\
      (snpScVectI16) \
   ); /*See if the kept score was an SNP or gap*/\
   \
   /*Find the gap penalties for the next round*/\
   mmAndI16(snpScVectI16, gapVectI16, gapExtendVectI16);\
   mmXorI16(gapVectI16, gapVectI16, gapOpenVectI16);\
   mmAddI16(gapVectI16, gapVectI16, snpScVectI16);\
   \
   /*Get the next snp scores*/\
   snpScVectI16 = insScVectI16; /*For next round*/\
}

#define vectI16GetIndelScore( \
   indelScVectI16, \
   gapVectI16 \
){\
   (indelScVectI16)=mmAddI16((indelScVectI16, gapVectI16);\
}
