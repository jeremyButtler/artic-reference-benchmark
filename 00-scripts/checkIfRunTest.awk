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
} # BEGIN

{ # MAIN
  # I am using && = "", because some paratemters are
  # unique to specific programs
  # General settings
  if($1 != programStr && programStr != "") next;
  if($3 != depthI && depthI != "") next;
  if($5 != percMutI && percMutI != "") next;
  if($22 != seedI && seedI != "") next;
  if($23 != modelStr && modelStr != "") next;
  if($36 != prefixStr && prefixStr != "") next;

  # Specific to buildcon
  if($24 != usedMedakaBl && usedMedakaBl != "") next;
  if($25 != usedMajconBl && usedMajconBl != "") next;
  if($26 != usedIvarBl && usedIvarBl != "") next;
  if($27 != usedRaconBl && usedRaconBl != "") next;
  if($28 != medakaPolishBl && medakaPolishBl != "") next;
  if($29 != ivarPolishBl && ivarPolishBl != "") next;
  if($30 != minLenI && minLenI != "") next;
  if($31 != minDepthI && minDepthI != "") next;

  # File paths and schemes
  if($32 != fqDirStr && fqDirStr != "") next;
  if($33 != schemeDirStr && schemeDirStr != "") next;
  if($34 != schemeStr && schemeStr != "") next;
  if($35 != verStr && verStr != "") next;
     # Scheme version (not program)

  matchBl = 1;
  exit;
} # MAIN

END{
   if(matchBl == 1) print "1";
   else print "0";
} # END
