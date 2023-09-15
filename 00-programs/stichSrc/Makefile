PREFIX=/usr/local/bin

CC=cc

# These are flags I do not want the user to over write,
# acidentally.
COREFLAGS=\
  -Wall\
  --std=gnu99 \
  -static\
  -O3
  #-DINSSNPDEL
# I need --std=gnu99 for popen (popen is not in --std=c99)
# This is likely do to popen being POSIX. So, this may not
# work well on windows

# These are here for the user to overwrite
CFLAGS=-DBLANK

DEBUGFLAGS=\
   DNOSEQCNVT
   #-DBYTEMATRIX
   #-DNOGAPOPEN
# Current options: -DNOSEQCNVT, -DHIRSCHTWOBIT

SOURCE=\
  samFunSrc/trimSam.c \
  samFunSrc/seqStruct.c\
  stichMinimapFun.c\
  stichInputAndHelp.c\
  stich.c
 
# Build findCoInfct
all:
	$(CC) $(COREFLAGS) $(CFLAGS) $(SOURCE) -o stich

debug:
	$(CC) -Wall  -static --std=gnu99 -O0 -ggdb $(DEBUGFLAGS) $(SOURCE) -o stichDebug
	# Used to use -g, but -ggdb provides more info for gdb
	bash debug.sh
	# edit debugCMDs.txt to change the gdb commands

egcc:
	#egcc $(COREFLAGS) $(CFLAGS) $(SOURCE) -o alnSeq
	egcc $(COREFLAGS) $(CFLAGS) $(SOURCE) -o stich
gcc:
	gcc  $(COREFLAGS) $(CFLAGS) $(SOURCE) -o stich
cc:
	cc   $(COREFLAGS) $(CFLAGS) $(SOURCE) -o stich

clean:
	rm stichDebug || printf ""; # Only thing to clean up

install:
	mv stich $(PREFIX) || printf "Unable to install stich at %s\n Change this with make PREFIX=/path/to/install install\n" $(PREFIX) && exit;
	chmod a+x $(PREFIX)/stich;
