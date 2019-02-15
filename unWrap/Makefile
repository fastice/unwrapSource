CC =		gcc
CFLAGS  =	 $(FLAGS) -c -I$(INCLUDEPATH)

OBJS= unwrap.o unwrapPhase.o phaseImage.o residueImage.o branchCuts.o \
      addPhaseRamp.o bcut.o labelRegions.o interpolatePhase.o

.KEEP_STATE:
all:	$(OBJS)





	mv *.o $(MACHTYPE)-$(OSTYPE)
