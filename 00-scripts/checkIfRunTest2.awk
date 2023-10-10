#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# checkIfRunTest.awk TOC: Table Of contents
#  - Check if ran a test with benchAll.sh already
#  o sec-01:
#    - Help message
#  o sec-02:
#    - Find the positions of the columns to check
#  o sec-03:
#    - Check for duplicates
#  o sec-04:
#    - Print out results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

BEGIN{
   if(helpStr != "")
   { # If the user wanted the help message
      printf "Name: checkIfRunTest.awk\n";
      printf "Use: \n";
      printf "  - checks if two entires in from bechAll.sh"
      printf " are the same\n";
      printf "  - All non filled inputs are ignored\n"

      printf "Input:\n";
      printf "  file.tsv: [Required]\n";
      printf "    - file with entries to check for";
      printf " duplicates.\n";

      printf "  -v programStr=\"$program\":\n";
      printf "     - Program used to build the consensus";
      printf "\n";

      printf "  -v StatsProgram=\"$statsProgram\":\n";
      printf "     - Program used to get stats for the";
      printf " consensus\n";

      printf "  -v depthI=\"$depthI\":\n";
      printf "     - Depth amplicons were subsampled to\n";

      printf "  -v percMutI=\"$refMutI\":\n";
      printf "     - Percentage of mutations made in";
      printf " reference\n";

      printf "  -v seedI=\"$seedI\":\n";
      printf "     - Seed used in benchmarking\n";

      printf "  -v modelStr=\"$modelStr\":\n";
      printf "     - Model used with Medaka\n";

      printf "  -v prefixStr=\"$prefixStr\":\n";
      printf "     - Prefix used to name the consensus\n";

      printf "  -v usedMedakaBl=\"$medUsedBl\":\n";
      printf "     - TRUE = used medaka; FALSE = no\n";

      printf "  -v usedMajConBl=\"$majUsedBl\":\n";
      printf "     - TRUE = used majcon; FALSE = no\n";

      printf "  -v usedIvarBl=\"$ivarUsedBl\":\n";
      printf "     - TRUE = used ivar; FALSE = no\n";

      printf "  -v usedRaconBl=\"$raconUsedBl\":\n";
      printf "     - TRUE = used racon; FALSE = no\n";

      printf "  -v medakaPolishbl=\"$medPolbl\":\n";
      printf "     - TRUE = did final polish with Medaka;";
      printf "\n     - FALSE = no polish.\n";

      printf "  -v ivarPolishbl=\"$ivarPolbl\":\n";
      printf "     - TRUE = did final polish with ivar;\n";
      printf "     - FALSE = no polish.\n";

      printf "  -v minLenI=\"$minLenI\":\n";
      printf "     - Min consensus length for buildcon\n";

      printf "  -v minDepthI=\"$minDepthI\":\n";
      printf "     - Min read depth for buildcon\n";

      printf "  -v fqDirStr=\"$fastqDirectoryStr\":\n";
      printf "     - Location of directory with amplicon";
      printf " fastq files\n";

      printf "  -v schemeDirStr=\"$schemeDirStr\":\n";
      printf "     - Directroy with scheme used\n";

      printf "  -v schemeStr=\"$schemeStr\":\n";
      printf "     - Name of scheme used\n";

      printf "  -v verStr=\"$versionStr\":\n";
      printf "     - Version of the scheme used\n";

      printf "  -v helpStr=\"help\":\n";
      printf "     - Print this help message\n";
      exit;
   } # If the user wanted the help message

   getline; # Get the header

   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Sec-02:
   #  - Find the positions of the columns to check
   #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   # Make sure not program parameter input. I sometimes do
   # this to make life easier
   sub(/^-/, "", statsProgamStr);

   # Find all the columns I want to check
   for(iCol = 1; iCol <= NF; ++iCol)
   { # Loop: Find the column of each header entry
     # General variables
     if($iCol == "program") progColI = iCol;
     if($iCol == "StatsProgram") statsProgColI = iCol;
     if($iCol == "depth") depthColI = iCol;
     if($iCol == "percMutate") mutColI = iCol;
     if($iCol == "usedSeed") seedColI = iCol;
     if($iCol == "medakaModel") modelColI = iCol;
     if($iCol == "prefix") prefColI = iCol;

     # Consensus building methods (most buildcon only)
     if($iCol == "usedMedaka") usedMedColI = iCol;
     if($iCol == "usedMajcon") usedMajColI = iCol;
     if($iCol == "usedIvar") usedIvarColI = iCol;
     if($iCol == "usedRacon") usedRacColI = iCol;
     if($iCol == "medakaPolish") medPolColI = iCol;
     if($iCol == "ivarPolish") ivarPolColI = iCol;

     # buildcon specific
     if($iCol == "minLen") minLenColI = iCol;
     if($iCol == "minDepth") minDepthColI = iCol;

     # File paths and schemes
     if($iCol == "fastq") fqColI = iCol;
     if($iCol == "schemeDir") schemeDirColI = iCol;
     if($iCol == "scheme") schemeColI = iCol;
     if($iCol == "schemeVer") verColI = iCol;
   } # Loop: Find the column of each header entry

   matchBl = 0;
}; # BEGIN


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Check for duplicates
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

{ # MAIN
  # I am using && = "", because some paratemters are
  # unique to specific programs
  # General settings
  if($progColI != programStr && programStr != "") next;
  if($depthColI != depthI && depthI != "") next;
  if($mutColI != percMutI && percMutI != "") next;
  if($seedColI != seedI && seedI != "") next;
  if($modelColI != modelStr && modelStr != "") next;
  if($prefColI != prefixStr && prefixStr != "") next;

  if($statsProgColI!=statsProgamStr && statsProgramStr!="")
     next;

  # Consensus building methods (most buildcon only)
  if($usedMedColI!=usedMedakaBl && usedMedakaBl!="") next;
  if($usedMajColI!=usedMajconBl && usedMajconBl!="") next;
  if($usedIvarColI != usedIvarBl && usedIvarBl != "") next;
  if($usedRacColI!= usedRaconBl && usedRaconBl != "") next;

  if($medPolColI != medakaPolishBl && medakaPolishBl != "")
     next;

  if($ivarPolColI != ivarPolishBl && ivarPolishBl != "")
     next;

   # buildcon specific
  if($minLenColI != minLenI && minLenI != "") next;
  if($minDepthColI != minDepthI && minDepthI != "") next;

  # File paths and schemes
  if($fqColI != fqDirStr && fqDirStr != "") next;

  if($schemeDirColI != schemeDirStr && schemeDirStr != "")
     next;

  if($schemeColI != schemeStr && schemeStr != "") next;
  if($verColI != verStr && verStr != "") next;
     # Scheme version (not program)

  matchBl = 1;
  exit;
}; # MAIN

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-04:
#  - Print out results
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

END{print matchBl;};
