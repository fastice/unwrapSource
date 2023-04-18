C =		gcc
ROOTDIR =	/Users/ian
PROGDIR =       $(ROOTDIR)/progs/GIT
INCLUDEPATH =	$(ROOTDIR)/progs/GIT
IHOME =		~
BINDIR =	$(IHOME)/bin/$(MACHTYPE)
#
CFLAGS =	'-O3 -m32 -I$(INCLUDEPATH) $(COMPILEFLAGS)'
CCFLAGS =  '-O3 -m32 -D$(MACHTYPE) $(COMPILEFLAGS) '
#-Wunused-variable'

CCFLAGS1= '-O3'
# uncomment to debug
#CFLAGS =	'-g -m32 -I$(INCLUDEPATH) $(COMPILEFLAGS)'
#CCFLAGS =  '-g -m32 -D$(MACHTYPE) $(COMPILEFLAGS)'
#CCFLAGS1= '-g'

COMMON=	$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/llToImageNew.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/julianDay.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/parseInputFile.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/polintVec.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/vectorFunc.o


ERS1CODE =	ers1Code/$(MACHTYPE)-$(OSTYPE)/initComplexImage.o ers1Code/$(MACHTYPE)-$(OSTYPE)/rowIncrement.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/rowPtr.o ers1Code/$(MACHTYPE)-$(OSTYPE)/initPowerImage.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/fillBuffer.o ers1Code/$(MACHTYPE)-$(OSTYPE)/openImage.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/writeBuffer.o ers1Code/$(MACHTYPE)-$(OSTYPE)/pixSize.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/getRadiometricParams.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/initOffsets.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/fractShiftBuffer.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/integerShift.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/initShifts.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/readShifts.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/freeImage.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/computePhase.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/initPhaseImage.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/getAzVaryingShifts.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/interpolateShifts.o \
                ers1Code/$(MACHTYPE)-$(OSTYPE)/earthRadius.o

STANDARD =	$(PROGDIR)/clib/$(MACHTYPE)-$(OSTYPE)/standard.o

RECIPES  =	$(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/polint.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/nrutil.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/ratint.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/four1.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svdfit.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svdcmp.o \
		$(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svdvar.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svbksb.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/pythag.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/hunt.o

UNWRAPCODE =	unWrap/$(MACHTYPE)-$(OSTYPE)/phaseImage.o unWrap/$(MACHTYPE)-$(OSTYPE)/residueImage.o \
                unWrap/$(MACHTYPE)-$(OSTYPE)/branchCuts.o unWrap/$(MACHTYPE)-$(OSTYPE)/unwrapPhase.o \
                unWrap/$(MACHTYPE)-$(OSTYPE)/addPhaseRamp.o unWrap/$(MACHTYPE)-$(OSTYPE)/labelRegions.o \
                unWrap/$(MACHTYPE)-$(OSTYPE)/interpolatePhase.o


TARGETS = unwrap

all: $(TARGETS)

UNWRAPDIRS =	unWrap ers1Code

unwrap:	
	@for i in ${UNWRAPDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=1; \
			cd $(PROGDIR); \
		); done
		gcc -m32 $(CCFLAGS1)   \
                unWrap/$(MACHTYPE)-$(OSTYPE)/unwrap.o $(UNWRAPCODE) $(ERS1CODE) $(COMMON)  $(STANDARD) $(RECIPES) \
                -lm -o $(BINDIR)/unwrap

