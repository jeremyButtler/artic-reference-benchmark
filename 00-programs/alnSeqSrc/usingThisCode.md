# Use:

This is here to give you an idea of the functions and data
  structures you will need to use to use this code in your
  own code.

This is not the best guide and does not cover everything,
  but I hope it will help.

## Reading in sequences

Sequences are stored in the seqStruct, which is struct-01
  of seqStruct.h. seqStruct stores the sequence id,
  sequence, and q-score entries. It also stores the length
  of the id, sequence, and q-score and the size of each
  buffer storing the id, sequence, and q-score.

An seqStruct should be initialized with the initSeqST
  function (fun-08 in seqStruct.h) only for the first use.
  After initialization use the blankSeqST() function 
  (fun-10 seqStruct.h) to reset the structure to no
  sequence.

Use the freeSeqST function (fun-07 seqStruct.h) to clean
  up. freeSeqST takes a pointer to an seqStruct structure
  and a variable to mark if the structure is on the heap
  or stack. If heapBl is set to 1, then the input structure
  will be freed, but if not, then only the buffers will be
  freed. freeSeqST(seqStruct, 1) does not set the pointer
  to the seqStruct to 0, you must do this.

A sequence can be read in from a fastq file using the
  readFqSeq function (fun-03 in seqStruct.h) or fasta file
  using the readFaSeq function (fun-04 in seqStruct.h).
  Both functions take a pointer to the FILE with the
  sequence and a pointer to an seqStruct structure. The
  FILE pointer will be set to the next entry in the
  fasta/fastq file.

If you compiled alnSetStruct.c without -DNOSEQCNVT, then
  you will also need to run seqToLookupIndex()
  (Fun-06 in alnSetStruct.c/h) to convert the sequence to
  each base to the look up index's for base comparisons.
  You can reconvert the sequence back to nucleotides with
  lookupIndexToSeq() (Fun-07 in alnSetStruct.c/h). However,
  all output bases will be in uppercase. These functions
  do nothing when -DNOSEQCNVT is used.

You can reverse complement a sequence using the
  reverseComplementSeq() function (fun-01 in seqStruct.h).

```
/*Here is an example.*/

FILE faFILE = fopen("sequence.fasta", "r");
struct seqStruct *sequence = 0;

if(faFILE == 0) return 0;

sequence=malloc(sizeof(struct seqStruct));

if(sequence == 0)
{
  fclose(faFILE);
  return 0;
}

initSeqST(&sequence);

switch(readFaSeq(faFILE, &sequence))
{ /*Switch: Check what kind of error*/
  case 0:
    fclose(faFILE);
    freeSeqST(sequence, 1);
    sequence = 0;
    return 0;        /*EOF*/

  case 1:
    fclose(faFILE);
    return sequence;

  case 2:
    fclose(faFILE);
    freeSeqST(sequence, 1);
    sequence = 0;
    return 0;        /*invalid file*/

  case 64:
    fclose(faFILE);
    freeSeqST(sequence, 1);
    sequence = 0;
    return 0;       /*memory error*/
} /*Switch: Check what kind of error*/

/*Only if aligning sequences*/
seqToLookupIndex(&sequence);

/*Do something with the index's here*/

lookupIndexToSeq(&sequence)
```

## Alignment settings

Settings for the alignment settings are stored in the
  alnSet structure (struct-01 alnSetStruct.h). This
  structure holds the output line wrap size, which
  alignment to run, if you are doing a reference/query
  scan, the scoring matrix, the gap extension penalty, the
  gap starting penalty, the min score to keep an
  alternative alignment, and the order to prefer a single
  direction. The default values for these settings can all
  be modified in alnSeqDefaults.h.
  
You can initialize an alnSet structure using the
  initAlnSet() function (fun-01 in alnSetStruct.h). This 
  will set all values in the alnSet structure to the
  default values in alnSeqDefaults.h.

For changing the scoring matrix you can use the
  setBasePairScore() function (fun-02 in alnSetStruct.h) to
  change the value for individual base pairs in a matrix or
  use the readInScoreFile() function
  (fun-04 alnSetStruct.h) to read in a scoring matrix. This
  scoring matrix is different than the normal format, so
  see scoring-matrix.txt for an example of my format.

Changing the non-scoring values for an alnSet structure
  requires you to manually change them. If the value ends
  in Bl, then it should only ever be a 1 or 0.

Finally you can free the alnSet structure using the 
  freeAlnSet() function (fun-03 in alnSetStruct.h). This
  takes in a pointer to an alnSet structure and a 1 or 0.
  If you use a 1, then the alnSet structure will not be
  freed, but if you use a 0 it will be freed. Just remember
  to set the pointer to null calling freeAlnSet(alnSet, 0);

```
An example

struct alnSet *settings = 0;
FILE *matrixFILE = fopen("scoring-matrix.txt", "r");

if(matrixFILE == 0) return ERROR;

settings = malloc(sizeof(struct alnSet));

if(settings == 0)
{
  fclose(matrixFILE);
  return ERROR;
}

initAlnSet(settings);
readInScoreFile(alnSetST, matrixFILE);

/*Do something*/

fclose(matrixFILE);
freeAlnSet(settings, 0);

return 0;
```
  
## Doing an alignment

alnSeq can do either an Needleman Wunsch, Hirschberg, or an
  Waterman Smith alignment. The Needleman Wunsch alignment
  is done by the NeedlemanAln function
  (fun-01 in needle.c/h), the Hirschberg is done by the 
  Hirschberg function (fun-01 in hirschberg.c/h), the
  Waterman Smith alignment is done by the WatermanAln
  function (fun-01 waterman.c/h), and the reference and
  query scan is done with WatermanAltAln
  (Fun-02 waterman.h/c).

All alignment functions take an pointer to an seqStruct
  structure with the reference sequence, an pointer to an
  seqStruct structure with the query sequence, and an
  pointer to a alignment settings structure.

The Needleman Wunsch, Waterman Smith, and query/reference
  scan waterman alignments will return an alnMatrixStruct
  structure (struct-01 alnMatrixStruct.h). The
  alnMatrixStruct has the directional matrix,
  best score (Waterman)/ending score (Needleman), and the
  extra scores kept for the reference/query scan
  (if was called for).

The Hirschberg returns an alnStruct datatype, which has the
  alignment. The alignment is stored as a series of gaps,
  SNPs, and matches for the reference and query separately.
  To get a printable alignment you would have to call
  alnSTToSeq (Fun-01 alnStruct.c/h) (I have not tested this
  function yet).

Once you are finished with the matrix you can free it with
  the freeAlnMatrixST() function (fun-02 in
  alnMatrixStruct.h). Just remember to set it to null after
  freeing.

```
struct seqStruct querySeq;
struct seqStruct refSeq;
struct alnSet settings;

struct alnMatrixStruct *matrix = 0;

FILE *faFILE = 0;

initAlnSet(&settings);
initSeqST(&querySeq);
initSeqST(&refSeq);

faFILE = fopen("query.fasta", "r");
if(faFILE == 0) return 0;
fclose(faFILE);

switch(readFaSeq(faFILE, &querySeq))
{
/*Check to see if read in sequence & handle errors here*/
}

faFILE = fopen("ref.fasta", "r");

if(faFILE == 0)
{
/*Handle file error here*/
}

fclose(faFILE);

switch(readFaSeq(faFILE, &refSeq))
{
/*Check to see if read in sequence & handle errors here*/
}

/* This is to avoid errors. Originally I wanted to be able
`  to align in the middle of sequences, but I never set
`  that option up fully. It would crash the current code.
`  Note the Hirschberg is set up for changing these values.
*/
queryST.endAlnUI = queryST.lenSeqUI;
queryST.offsetUI = 0;

refST.endAlnUI = refST.lenSeqUI;
refST.offsetUI = 0;

matrix = WatermanAln(querySeq, refSeq, settings, "");
/* or matrix = NeedlemanAln(querySeq, refSeq, settings);
`  or matrix = Hirschberg(querySeq, refSeq, settings);
`    For Hirschberg do nothing else except freeing the
`    other variables. I will change the printing methods
`    for this later.
*/

/*Do something with matrix here*/

freeSeqST(&refSeq, 0);
freeSeqST(&qrySeq, 0);
freeAlnMatrixST(alignmentMatrix);
freeAlnSet(settings, 0);
```
  
## Printing alignments

### Printing alternative alignments

Alternative alignments can be printed with the
  printAltWaterAlns function (fun-04 waterman.c/h). This
  functions takes the alnMatrixStruct returned by
  WatermanAltAln (not WatermanAln), a min score to keep an
  alternative alignment, and an output file.

### Printing primary alignments

To print out the alignment you first have to create an
  alnStruct. This is already done for the Hirschberg 
  alignment, but not for the Needleman or Waterman
  alignments. 

Converting the directional matrix to an alignment array 
  is done with the cnvtDirMatrixToAlnAry() function
  (fun-02 in alnStruct.h). This function takes in a pointer
  to an seqStruct structure with the reference sequence,
  an pointer to an seqStruct structure with the query
  sequence, a scoreStruct structure
  (alnMatrixStruct->bestScoreST), and the direction matrix
  (alnMatrixStruct->dirMatrixST). The returned value is
  an alnStruct with the alignment.

After building the alignment array you can then free the
  alnMatrixStruct structure. Just remember to save the best
  score before freeing the alnMatrixStruct.

You can print out the alignment stored in the alnStruct
  using the printAln function (fun-03 in alnStruct.c/h).

Just make sure to free the alignment array structure
  (alnStruct) with freeAlnSt(alnStruct, 1);. Were 1 tells
  the function to free the structure and variables in the
  structure. A 0 would just free the variable inside the
  structure.

```
char *alignedQuery = 0;
char *alignedRef = 0;
long bestScoreL = 0;
struct alnStruct *alignment = 0;

/*Other variables, such as the output file*/

seqToLookupIndex(&sequence);

matrix = WatermanAln(refSeqST, querySeqST, settings);

/*Get the alignment array*/
alignment =
  dirMatrixToAlnST(
    refSeq,
    querySeq,
    matrix->dirMatrixST,
    matrix->besScoreST,
    1
); /*Build the alignment array*/

bestScoreL->matrix->bestScoreST.scoreL;

freeAlnMatrixST(matrix);
matrix=0;

if(alignment == 0)
{
  /*Handle memory errors*/
}

/*Printing out the structure*/

lookupIndexToSeq(&sequence);

printAln(
  outFILE,              /*File to output to*/
  outFileStr,           /*File name*/
  refSeqST,
  querySeqST,
  alignment 
  bestScoreL,
  settings,             /*unsigned short from alnSetStruct*/
  alignment             /*alnStruct*/
  scoreMatrixFileNamea  /*Name of user input score matrix*/
    /*Use defMatrixNameStr if using the default matrix*/
);

fclose(outFILE);

freeSeqST(refSeqST);
freeSeqST(querySeqST);

freeAlnST(alignment, 1);
freeAlnSet(settings, 0);

return 0;
```

## Two bit arrays

All twoBit functions are static inlined (only a
  twoBitArrays.h file).

The two bit array is an array of uint8_t's that each hold
  four two bit elements. This array is stored in a
  twoBitAry structure (struct-01 twoBitArrays.h). The
  twoBitAry structure has a pointer to the first unit8_t
  (limb) in the array (firstLimbUCPtr), a pointer to the 
  limb currently one (limbOnUCPtr), a counter telling which
  element we are on in the limb, and the length of the
  array.

You can create a twoBitAry structure using the makeTwoBit
  function (fun-12 in twoBitArrys.h). This function will
  return a two bit array structure allocated on the head
  with a two bit array of the requested size. The inputs
  are the number of elements in the two bit array and 0.
  However, if you wish to not make an array you can use
  twoBitStruct = makeTowBit(0, 1).

When finished the two bit array can be freed with
  the freeTwoBit function (fun-14 twoBitArrays.h). This
  takes a pointer to the twoBitAry structure to free,
  a 1 or 0 (free structure) to tell if to not free the
  structure, and a 1 or 0 (free array) to tell if to free
  the array stored in the structure. If you free the
  structure, make sure you set the pointer to null.

Their are several functions that are used to move around
  the two bit array, most of which only take a pointer to
  a two bit array structure or a pointer to the structure
  and a additional value.

To get the current two bit element on in the array use
  the getTwoBitElm() function
  (fun-01 twoBitArrays.h). To move forward one element use

You can get the index of the element you are on in the two
  bit array using the twoBitIndex() function (fun-17
  twoBitArrays.h).

You can get the length of the two bit array by using the
  twoBitGetLen() function (fun-16 twoBitArrays.c).
 
To change the current element on in a two bit array use the
  changeTwoBitElm() function (fun-08 in twoBitArrays.h)

You can do a shallow copy of the two bit array by using
  the cpTwoBitPos() function (fun-13 twoBitArrays.h). This
  will allow you to work on the same two bit array with two
  separate pointers.
  
You can move forward in a two bit array using the
  twoBitMvToNextElm() (fun-03 twoBitArrays.h) and
  twoBitMvForXElm functions (fun-04 twoBitArrays.h).
  You can also use the twoBitMvXElmFromStart() function
  (fun-15 twoBitArrays.h) to move from the starting index
  instead of the current position.

You can move backwards in a two bit array using the
  twoBitMvBackOnElm() (fun-05 twoBitArrays.h) and
  twoBitMvBackXElm functions (fun-06 twoBitArrays.h).

Finally if you need to blank a limb to 0 you can use the
  blankTwoBitLimb() function (fun-09 twoBitArrays.h).

```
/*An example of using two bit arrays*/

struct twoBitAry *twoBitST = makeTwoBit(10, 0);

/*Add some values into the array*/
changeTwoBitElm(twoBitST, 1);
twoBitMvToNextElm(twoBitSt);

changeTwoBitElm(twoBitST, 0);
twoBitMvToNextElm(twoBitSt);

changeTwoBitElm(twoBitST, 3);
twoBitMvToNextElm(twoBitSt);

/*Move back to the start of the array*/
twoBitMvXElmFromStart(twoBitST, 0);

/*Print out the first element*/
printf("%u\n", getETwoBitElm(twoBitST);

freeTwoBit(twoBitST, 0, 0);
twoBitST = 0;
return 0;
```
