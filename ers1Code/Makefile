CC =		gcc
CFLAGS =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS=  initComplexImage.o rowIncrement.o rowPtr.o initPowerImage.o \
        fillBuffer.o openImage.o writeBuffer.o pixSize.o \
        getRadiometricParams.o initOffsets.o \
        fractShiftBuffer.o integerShift.o initShifts.o readShifts.o \
        freeImage.o computePhase.o initPhaseImage.o getAzVaryingShifts.o \
        interpolateShifts.o earthRadius.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shelf uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)
	mv *.o $(MACHTYPE)-$(OSTYPE)
