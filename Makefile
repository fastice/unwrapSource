#
# Modify this section to point to where stuff is.
# Current names are for a specific file system.
# ROOTDIR: root directory for code (e.g. /home/username directory).
# PROGDIR: location for top source code directory (default ROOTDIR/progs/GIT64)
# BINHOME: root directory for binaries
# BINNAME: archetecture dependent basename for bin dir
# BINDIR: directory for binaries (default BINHOME/bin/BINNAME) (will create if doesn't exist)
# INCLUDEPATH: include path (default PROGDIR anything else could cause a problem)
# Various directors can be overridden with environment variable or from make command
# make BINHOME=/a/path/to/binhome
#
# Base directory for code
USER =	$(shell id -u -n)
$(info  "USER is $(USER)")
MACHTYPE = $(shell uname -m)
OSTYPE = $(shell uname -s)
$(info "MACHTYPE/OSTYPE $(MACHTYPE) $(OSTYPE)")
#
# Default rootdir
ifneq ($(ROOTDIR)),)
	ROOTDIR =	$(dir $(CURDIR))
endif
$(info ROOTDIR="$(ROOTDIR)")
#
# Default root for source code
ifneq ($(PROGDIR)),)
	PROGDIR =       $(dir $(CURDIR))
endif
$(info PROGDIR ="$(PROGDIR)")
#
# Default location root for compiled programs
ifneq ($(BINHOME)),)
	BINHOME =		~$(USER)
endif
# For historical reasons, can compile with 32-bit memory model using MEM=-m32
# In almost all cases, should be compiled as 64bit.
ifneq ($(MEM), -m32)
	BINNAME=	$(MACHTYPE)
	FFTDIR = $(MACHTYPE)-$(OSTYPE)
else
	BINNAME =	i386
	FFTDIR = i386-$(OSTYPE)
endif
#
# Default binary directory
ifneq ($(BINDIR)),)
	BINDIR =	$(BINHOME)/bin/$(BINNAME)
endif
$(info "BINDIR = $(BINDIR)")
#
# Create bin dir if it doesn't exist
$(shell mkdir -p $(BINDIR))
#
# Default include path
ifneq ($(INCLUDEPATH)),)
	INCLUDEPATH =	$(PROGDIR)
endif
$(info INCLUDEPATH ="$(INCLUDEPATH)")
#
# Compiler stuff
#
C =		gcc
#
CFLAGS =	'-O3 $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS)'
CCFLAGS =  '-O3 $(MEM) $(COMPILEFLAGS) '
GDAL = -lgdal -lcurl  -lsqlite3 -llzma -lpoppler -lopenjp2 -lssh2 -llcms2
#
CCFLAGS1= -O3 
#-no-pie
# uncomment to debug
#CFLAGS =	'-g $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS)'
#CCFLAGS =  '-g $(MEM) -D$(MACHTYPE) $(COMPILEFLAGS)'
#CCFLAGS1= '-g'
#
ifneq ($(OSTYPE),darwin)
	NOPIE =	-no-pie
endif
$(info NOPIE ="$(NOPIE)")
#
# ******** SHOULD NOT NEED TO MODIFY BELOW HERE *********
#
COMMON=	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/geojsonCode.o \
                $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/llToImageNew.o \
		$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/julianDay.o \
                $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/lltoxy1.o \
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

UNWRAPDIRS =	unWrap ers1Code $(PROGDIR)/clib $(PROGDIR)/cRecipes $(PROGDIR)/mosaicSource/common

unwrap:	
	@for i in ${UNWRAPDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=1; \
			cd $(PROGDIR); \
		); done
		g++ $(MEM) $(CCFLAGS1)   \
                unWrap/$(MACHTYPE)-$(OSTYPE)/unwrap.o $(UNWRAPCODE) $(ERS1CODE) $(COMMON)  $(STANDARD) $(RECIPES) \
                -lm $(GDAL) -o $(BINDIR)/unwrap

