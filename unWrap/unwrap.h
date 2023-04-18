#define VALIDPIXEL 0x1
#define PRESIDUE 0x2
#define NRESIDUE 0x4
#define RESIDUE 0x6
#define UNCHARGED 0x8
#define BRANCHCUT 0x10
#define LOWDENSITY 0x20
#define REGION 0x40
#define MASK 0x80
/* NOTE THIS IS THE SAME FLAG AS LOWDENSITY BUT
                       LOWDENSITY IS NOT NEEDED BY THE TIME BORDER IS */

#define BORDER 0x20
#define NOTBORDER 0xDF
#define CHARGED 0xF7
#define NORESIDUE 0xF9
#define INVALIDPIXEL 0xFE
#define HIGHDENSITY 0xDF
#define MAXCOMPLEX 16000

#define OPPOSITE 0
#define SAME 1
#define NEUTRAL 2
#define IGNORECHARGE 3
#define validPixel(A) ((A)&VALIDPIXEL ? TRUE : FALSE)
#define maskPixel(A) ((A)&MASK ? TRUE : FALSE)

#define charged(A) ((A)&UNCHARGED ? FALSE : TRUE)
#define notBranch(A) ((A)&BRANCHCUT ? FALSE : TRUE)
#define noFPNeighbors(A) ((A).left < LARGEINT || (A).upper < LARGEINT ? FALSE : TRUE)

#define noSPNeighbors(A) ((A).right < LARGEINT || (A).lower < LARGEINT ? FALSE : TRUE)

typedef struct rListType
{
    int32_t r;
    int32_t a;
} rList;

typedef struct distType
{
    int32_t r1;
    int32_t r2;
    int32_t a1;
    int32_t a2;
    int32_t d;
} dist;

typedef struct regionElement
{
    int32_t label;
    int32_t nTotal;
    int32_t nLeft;
    int32_t charge;
    int32_t mask;
    rList *pts;
    dist *distances;
    struct regionElement *next;
} regionList;

typedef struct ptType
{
    int32_t i;
    int32_t j;
} pt;

typedef struct nodeListType
{
    int32_t node;
    float distance;
} nodeListElement;

typedef struct fpNeighborsType
{
    int32_t left;
    int32_t upper;
} fpNeighbors;

typedef struct spNeighborsType
{
    int32_t right;
    int32_t lower;
} spNeighbors;

typedef struct columnType
{
    int32_t start;
    int32_t end;
    float phaseOffset;
} columnStruct;

typedef struct eqTableType
{
    int32_t a;
    int32_t b;
    float dp;
} eqTableElement;

typedef struct unwrapImageType
{
    int32_t rangeSize;
    int32_t azimuthSize;
    int32_t nPlus;
    int32_t nMinus;
    int32_t nLeft;
    int32_t left;
    int32_t right;
    int32_t up;
    int32_t down;
    int32_t blockSize;
    int32_t residueDensity;
    int32_t thresh;
    char *mask;
    char *cutFile;
    char *labelFile;
    float **phase;
    float **dpdr;
    float **dpda;
    char **bCuts;
    int32_t **labels;
    int32_t nLabels;
    int32_t *labelCount;
    FILE *fpReport;
} unwrapImageStructure;

void phaseImage(unwrapImageStructure *unwrapImage,
                ers1ComplexImage *complexImage, float omegaR, float omegaA);

void residueImage(unwrapImageStructure *unwrapImage);

void branchCuts(unwrapImageStructure *unwrapImage);

void unwrapPhase(unwrapImageStructure *unwrapImage);

void interpolatePhase(unwrapImageStructure *unwrapImage);

void addPhaseRamp(unwrapImageStructure *unwrapImage,
                  float omegaR, float omegaA, int32_t noCenter);

void labelRegions(unwrapImageStructure *unwrapImage);
