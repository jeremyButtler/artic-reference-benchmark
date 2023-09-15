#define defVersion 20230914 /*Version number*/

#define defPrefix "scaffold"
#define defStichOverwrite 0 /*Overwrite output file*/
#define defThreads 3    /*Number of threads to use*/
#define defMaskBase 'N'/*Mask to use on low support bases*/
#define defMinStichDepth 3 /*Min depth to do voting at*/
   /*Under this depth I require 100% support*/
#define defMinStichSup 66
  /*Min % support needed to select one disagreement over
  ` another
  */
#define defMinimap2Bl 1/*1:replace alnSeq with minimap2*/
  /*I need to look at the mem waterman/Hirschberg
  ` combination in alnSeq. For some odd reason the
  ` Hirschberg is adding a few extra bases to the ref
  */

/*Command to run mininmap2*/
#define stichMinimap2CMD "minimap2 --eqx --secondary=no -a -x map-ont"
