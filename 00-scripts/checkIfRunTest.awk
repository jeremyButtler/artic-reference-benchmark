##########################################################
# Name: checkIfRunTest.awk
# Use:
#   - Checks to see if I have run a benchAll.sh test
# Input:
#   - See the BEGIN block (everything is needed)
# Output:
#   - Returns
#     o 1 if this test has already been run
#     o 0 if this test has not been run
##########################################################

BEGIN{
   matchBl = 0;
   if(progamStr == "") exit; # Program used
   if(depthI == "") exit;    # Subsampled read depth
   if(percMutI == "") exit;  # percent muated reference by
   if(seedI == "") exit;     # Seed used in test
   if(modelStr == "") exit;  # Model of medaka used
   if(usedMedakaBl == "") exit; # Did you use medaka
   if(usedMajConBl == "") exit; # Did you use majCon
   if(minDepthI == "") exit;    # Min read depth (buildCon)
   if(minLenI == "") exit;      # min read/consensus length
   if(fqDirStr == "") exit;     # dir with fastq files
   if(schemeDirStr == "") exit; # scheme directory
   if(schemeStr == "") exit;    # scheme name
   if(verStr == "") exit;       # scheme version
   if(prefixStr == "") exit;      # prefix
} # BEGIN

{ # MAIN
  if($1 != progamStr) next;
  if($3 != depthI) next;
  if($5 != percMutI) next;
  if($24 != seedI) next;
  if($25 != modelStr) next;
  if($26 != usedMedakaBl) next;
  if($27 != usedMajConBl) next;
  if($28 != minLenI) next;
  if($29 != minDepthI) next;
  if($30 != fqDirStr) next;
  if($31 != schemeDirStr) next;
  if($32 != schemeStr) next;
  if($33 != verStr) next;
  if($34 != prefixStr) next;

  matchBl = 1;
  exit;
} # MAIN

END{
   if(matchBl == 1) print "1";
   else print "0";
} # END
