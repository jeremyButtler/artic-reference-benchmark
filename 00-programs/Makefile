CC=cc
CXX=c++

# Find the directory of this make file
# from: https://stackoverflow.com/questions/18136918/how-to-get-current-relative-directory-of-your-makefile
# For more recent version of GMAKE
curDirStr=$(shell pwd)
#makeDirStr := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
condaStr=$(shell conda info | grep -i 'base environment' | sed 's/base.*: //; s/  *.read only.//; s/ //g')
condaCallStr=$(condaStr)/etc/profile.d/conda.sh

#  MAKEFILE_LIST changes as include files come in and out
#    of scope. The last item is the current file. 
#  lastword plucks the last item; Makefile name, relative
#    to pwd
#  realpath is built-in to make, and resolves to a
#    canonical path from filesystem root
#  dir trims off the filename, leaving just the directory.

all: buildAlnSeq buildFqGetIds buildScoreReads buildTrimPrimers buildTrimSamFile buildBuildCon buildStich makeLilo

buildAlnSeq:
	make -C alnSeqSrc CC=$(CC)|| printf "Failed to make alnSeq\n" && exit;
	mv alnSeqSrc/alnSeq ./ || printf "Failed to move alnSeq\n" && exit;
	make -C alnSeqSrc mid CC=$(CC)|| printf "Failed to make alnSeqMid\n" && exit;
	mv alnSeqSrc/alnSeqMid ./ || printf "Failed to move alnSeqMid\n" && exit;

buildFqGetIds:
	# fqGetIds works better with vector support.
	make -C fqGetIdsSrc sse CC=$(CC) || make -C fqGetIdsSrc neon CC=$(CC) || make -C fqGetIdsSrc CC=$(CC) || printf "Failed to complile fqGetIds\n" && exit;
	mv fqGetIdsSrc/fqGetIds ./ || printf "Failed to move fqGetIds\n" && exit;

buildScoreReads:
	make -C findCoInfctSrc scoreReads CC=$(CC) || printf "Failed to complile scoreReads\n" && exit;
	mv findCoInfctSrc/scoreReads ./ || printf "Failed to move scoreReads\n" && exit;

buildTrimPrimers:
	make -C findCoInfctSrc trimPrimers CC=$(CC) || printf "Failed to complile trimPrimers\n" && exit;
	mv findCoInfctSrc/trimPrimers ./ || printf "Failed to move trimPrimers\n" && exit;

buildTrimSamFile:
	make -C findCoInfctSrc trimSam CC=$(CC) || printf "Failed to complile trimSamFile\n" && exit;
	mv findCoInfctSrc/trimSamFile ./ || printf "Failed to move trimSamFile\n" && exit;

buildBuildCon:
	make -C findCoInfctSrc buildCon CC=$(CC) || printf "Failed to complile buildCon\n" && exit;
	mv findCoInfctSrc/buildCon ./ || printf "Failed to move buildCon\n" && exit;

buildStich:
	make -C stichSrc CC=$(CC) || printf "Failed to complile stich\n" && exit;
	mv stichSrc/stich ./ || printf "Failed to move stich\n" && exit;

makeLilo:
	# Install Porechop for Lilo (installed in 00-programs/bin)
	# Porchop needs to be in the path for Lilo. Have the
	# benchmarking script add it to the path.
	#   --home=$(curDirStr) allows me to do a local install
	source $(condaStr)/etc/profile.d/conda.sh || printf "Is conda activated?\n";
	ls Porechop-1 || git clone https://github.com/sclamons/Porechop-1;
	#cd Porechop-1 && python3 setup.py install --home=$(curDirStr);
	# Can not do local install because messes up LILO
	cd Porechop-1 && sudo python3 setup.py install;
	# install Lilo
	ls Lilo || git clone https://github.com/amandawarr/Lilo;
	cd Lilo && conda env create --file LILO.yaml;
	cd Lilo && conda env create --file scaffold_builder.yaml ;

removeprograms:
	source $(condaCallStr) || printf "Is conda activated?\n";
	conda init bash || printf "Failed to activate bash\n";
	conda remove -n LILO --all -y LILO || printf "No LILO enviroment in conda\n";
	conda remove --all -n scaffold_builder -y scaffold_builder || printf "No scaffold_builder conda envrioment\n";
	rm -f -r bin || printf "Porchop not compiled in 00-programs\n";
	rm -f -r lib || printf "Porchop not compiled in 00-programs\n";
	rm alnSeq || printf "alnSeq not compiled\n";
	rm alnSeqMid || printf "alnSeqMid not compiled\n";
	rm fqGetIds || printf "fqGetIds not compiled\n";
	rm buildCon || printf "buildCon not compiled\n";
	rm trimPrimers || printf "trimPrimers not compiled\n";
	rm trimSamFile || printf "trimSamFile not compiled\n";
	rm scoreReads || printf "scoreReads not compiled\n";

# Cleans everything not downloaded from github
# (everything that is not my code)
removeall:
	rm -r -f Lilo || printf "Lilo source already cleaned\n";
	rm -r -f Porechop-1 || printf "Porchop source already cleaned\n";
	source $(condaCallStr) || printf "Is conda activated?\n";
	conda init bash || printf "Failed to activate bash\n";
	conda remove -n LILO --all -y LILO || printf "No LILO enviroment in conda\n";
	conda remove --all -n scaffold_builder -y scaffold_builder || printf "No scaffold_builder conda envrioment\n";
	rm -f -r bin || printf "Porchop not compiled in 00-programs\n";
	rm -f -r lib || printf "Porchop not compiled in 00-programs\n";
	rm alnSeq || printf "alnSeq not compiled\n";
	rm alnSeqMid || printf "alnSeqMid not compiled\n";
	rm fqGetIds || printf "fqGetIds not compiled\n";
	rm buildCon || printf "buildCon not compiled\n";
	rm trimPrimers || printf "trimPrimers not compiled\n";
	rm trimSamFile || printf "trimSamFile not compiled\n";
	rm scoreReads || printf "scoreReads not compiled\n";

# This is to make sure pwd works across board.
test:
	pwd;
	printf "%s\n" $(curDirStr);
