/*
    Include file for ers1 images.
*/
#include "clib/standard.h"
#include "cRecipes/cRecipes.h"

/*
           Constants
*/
#define RANGEONLY 0
#define RANGEANDAZIMUTH 1
#define INIT 2
#define UPDATE 3

#define RETURN 0
#define CONTINUE 1
#define OK 0
#define MAXLAYERS 10000
#define COMPLEX 0
#define POWER 1
#define AMPLITUDE 2
#define CORRELATION 3
#define PHASE 4
#define FLOAT 5
#define UNWRAPPED 6
#define EMINOR 6356.7523142
#define EMAJOR 6378.137
#define F 1.0 / 298.2572236
#if PAF == 1
#define RANGESIZE 4992
#define AZIMUTHSIZE 26368
#define BUFFERSIZE 10223616
#define AMPLITUDELINESIZE 4992
#elif PAF == 2
#define RANGESIZE 1248
#define AZIMUTHSIZE 1318
#define BUFFERSIZE 6579456
#define AMPLITUDELINESIZE 1248
#elif PAF == 3
#define RANGESIZE 4992
#define AZIMUTHSIZE 5376
#define BUFFERSIZE 10223616
#define AMPLITUDELINESIZE 4992
#else
#define RANGESIZE 2048
#define AZIMUTHSIZE 12800
#define BUFFERSIZE 8388608
#define AMPLITUDELINESIZE 2048
#endif
#define OVERLAP 80
#define READIMAGE 0
#define WRITEIMAGE 1
#define DSHIFT 0.05

#define KERNELSIZE 10
#define NKERNEL 20
#define AZIMUTHPIXELSIZE 3.8999999
#define RANGEPIXELSIZE 7.8999

#define RANGERESOLUTION 7.9
#define LAMBDAERS1 0.05656
#define LARGEINT 2000000000

#define DESCENDING 0
#define ASCENDING 1

/*
   Global variables declarations
*/
extern int32_t RangeSize;       /* Range size of complex image */
extern int32_t AzimuthSize;     /* Azimuth size of complex image */
extern int32_t BufferSize;      /* Size of nonoverlap region of the buffer */
extern int32_t BufferLines;     /* # of lines of nonoverlap in buffer */
extern double RangePixelSize;   /* Size in m of range pixel */
extern double AzimuthPixelSize; /* Size in m of azimuth pixel */

/*
     Typedefs.
*/

typedef struct ers1OffsetType
{
    int32_t startRange;
    int32_t endRange;
    int32_t startAzimuth;
    int32_t endAzimuth;
    float *a1;
    float *a2;
    float *a3;
    float *a4;
    int32_t *columnLeft;
    int32_t *columnRight;
    int32_t *rowUp;
    int32_t *rowDown;
    float *deltaR;
    float *deltaA;
} ers1Offset;

typedef struct ers1ShiftType
{
    int32_t startRange;
    int32_t endRange;
    int32_t startAzimuth;
    int32_t endAzimuth;
    int32_t *columnShift;
    int32_t *rowShift;
    float *deltaR;
    float *deltaA;
    int32_t *aKernel;
    int32_t *rKernel;
    /*
        Stuff for azimuth varying shifts
    */
    int32_t shiftFlag;
    int32_t gRo;
    int32_t gAo;
    int32_t dRg;
    int32_t dAg;
    int32_t gNr;
    int32_t gNa;
    float **rangeShifts;
    float **azimuthShifts;
} ers1Shift;

typedef int16_t ers1ComplexElement;

typedef float (*ers1Kernel)[KERNELSIZE];

typedef struct ers1ComplexType
{
    ers1ComplexElement r;
    ers1ComplexElement i;
} ers1Complex;

typedef struct ers1FloatComplexType
{
    float r;
    float i;
} ers1FloatComplex;

typedef struct ers1TimeType
{
    int16_t hour;
    int16_t minute;
    int16_t second;
} ers1Time;

typedef struct ers1DateType
{
    int16_t year;
    int16_t day;
    int16_t month;
} ers1Date;

typedef struct ers1BaselineType
{
    float minute;
    float bParallel;
    float bPerpendicular;
    float omegaA;
    float omegaR;
    float deltaBn;
    int16_t firstOrbit;
    int16_t secondOrbit;
} ers1Baseline;

typedef char ers1Amplitude;

typedef float ers1Power;

typedef float ers1Phase;

typedef struct headerType
{
    FILE *infp;
    FILE *outfp;
    int32_t rs; /* Image size */
    int32_t as;
    int32_t currentRow;
    int32_t firstBufferRow;
    int32_t lastBufferRow;
    int32_t bufferSize;
    int32_t imageType;
    int32_t imageMode;
} headerInfo;

typedef struct ers1GenericImageType
{
    headerInfo header;
    char *image;
} ers1GenericImage;

typedef struct ers1ComplexImageType
{
    headerInfo header;
    ers1Complex *image;
} ers1ComplexImage;

typedef struct ers1AmplitudeImageType
{
    headerInfo header;
    ers1Amplitude *image;
} ers1AmplitudeImage;

typedef struct ers1PowerImageType
{
    headerInfo header;
    ers1Power *image;
} ers1PowerImage;

typedef struct ers1PhaseImageType
{
    headerInfo header;
    ers1Phase *image;
} ers1PhaseImage;

typedef ers1Complex (*ers1Complex2D)[RANGESIZE];

/*
    Create and init full size (standard size) complex image.
*/
ers1ComplexImage *initFullSizeComplexImage();

/*
    Create and init complex image.
*/
ers1ComplexImage *initComplexImage(int32_t rs, int32_t as, int32_t currentLine,
                                   int32_t bufferSize);
/*
    Increment row pointer.
*/
void rowIncrement(ers1GenericImage *image);
/*
   Return pointer to current row.
*/
void *rowPtr(ers1GenericImage *image);

/*
    Create and init full size (standard size) power image.
*/
ers1PowerImage *initFullSizePowerImage();

/*
    Create and init power image.
*/
ers1PowerImage *initPowerImage(int32_t rs, int32_t as, int32_t currentLine,
                               int32_t bufferSize);
/*
    Fill current buffer starting at current row-overlap.
*/
void fillBuffer(ers1GenericImage *image);
/*
    Open file for image
*/
void openImage(ers1GenericImage *image, char *filename, int32_t rwFlag);
/*
   Write current buffer to output file.
*/
void writeBuffer(ers1GenericImage *image);
/*
    Compute pixel size for a given image type
*/
int32_t pixSize(int32_t imageType);
/*
    Get radiometric params from ERS1 .RAD_1 file.
*/
void getRadiometricParams(char *radiometricFile, float *a1,
                          float *a2, float *a3);

/*
    Compute offsets as function of range using either deltaR and deltaA
    or range dependent offsets in offsetFile.
*/
void initOffsets(ers1Offset *offsets, float *omega, char *offsetFile,
                 float deltaR, float deltaA, ers1ComplexImage *imageA);

/*
   Interpolate complex Image

    void interpolateComplexImage(ers1ComplexImage *imageA,
                               ers1Complex2D rowB,ers1Complex *interpolatedRow,
                               ers1Offset offsets, int32_t rowIndex);
*/
/*
   Shift entire buffer by fraction of a pixels.
*/
void fractShiftBuffer(ers1GenericImage *image, ers1Shift shifts,
                      ers1Kernel kernel);

/*
   Shift complex Image by integer valued shift.
*/
void integerShift(ers1ComplexImage *imageA,
                  ers1Complex *rowB,
                  ers1Complex *interpolatedRow,
                  ers1Shift shifts, int32_t rowIndex);
/*
    Compute shifts as function of range using either deltaR and deltaA
    or range dependent shifts in shiftFile.
*/
void initShifts(ers1Shift *shifts, char *shiftFile,
                float deltaR, float deltaA,
                ers1ComplexImage *imageA, ers1Kernel *kernel,
                float *dR1, float *dA1, int32_t modeFlag);

/*
   Read shift file. Returns dR, dA, omega as function of range.
*/
void readShifts(char *shiftFile, float *dR, float *dA, float *omega);
/*
    Free complex image.
*/
void freeImage(ers1GenericImage *image);
/*
    Compute phase for flattening image.
*/
double *computePhase(double omegaR, double bn, double bp, double rc,
                     double h, double lat, double *phase);
/*
    Create and init phase image.
*/
ers1PhaseImage *initPhaseImage(int32_t rs, int32_t as, int32_t currentRow,
                               int32_t bufferSize);
/*
    Get array of shifts from shiftfile to allow shifts to be updated along
    track.
*/
void getAzVaryingShifts(ers1Shift *shifts, char *shiftFile);
/*
   Interpolate registration data to get shifts for given azimuth (dr,da).
*/

void interpolateShifts(int32_t az, ers1Shift *shifts, float *dR1, float *dA1);
/*
   Earth Radius
*/
double earthRadius(double lat, double rp, double re);
