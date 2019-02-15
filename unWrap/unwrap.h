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
#define validPixel(A) ((A) & VALIDPIXEL ? TRUE : FALSE)
#define maskPixel(A) ((A) & MASK ? TRUE : FALSE)

#define charged(A) ((A) & UNCHARGED ? FALSE : TRUE)
#define notBranch(A) ((A) & BRANCHCUT ? FALSE : TRUE)
#define noFPNeighbors(A) ((A).left < LARGEINT || (A).upper < LARGEINT ? FALSE : TRUE)

#define noSPNeighbors(A) ((A).right < LARGEINT || (A).lower < LARGEINT ? FALSE : TRUE)

typedef struct rListType {
    int r;
    int a;
} rList;

typedef struct distType {
   int r1;
   int r2;
   int a1;
   int a2;
   int d;
} dist;

typedef struct regionElement {
    int label;
    int nTotal;
    int nLeft;
    int charge;
    int mask;
    rList *pts;
    dist *distances;
    struct regionElement *next;    
} regionList;

typedef struct ptType {
    int i;
    int j;
} pt;

typedef struct nodeListType {
    int node;
    float distance;
} nodeListElement;

typedef struct fpNeighborsType {
    int left;
    int upper;
} fpNeighbors;

typedef struct spNeighborsType {
    int right;
    int lower;
} spNeighbors;

typedef struct columnType {
     int start;
     int end;
     float phaseOffset;
} columnStruct;

typedef struct eqTableType {
     int a;
     int b;
     float dp;
} eqTableElement;


typedef struct unwrapImageType {
    int rangeSize;
    int azimuthSize;
    int nPlus;
    int nMinus;
    int nLeft;
    int left;
    int right;
    int up;
    int down;   
    int blockSize; 
    int residueDensity; 
    long int thresh;
    char *mask; 
    char *cutFile; 
    char *labelFile;
    float **phase;
    float **dpdr;
    float **dpda;
    char **bCuts;
    int **labels;
    int nLabels;
    int *labelCount;
    FILE  *fpReport;
 } unwrapImageStructure;

void phaseImage(unwrapImageStructure *unwrapImage,
                ers1ComplexImage *complexImage, float omegaR, float omegaA);

void residueImage( unwrapImageStructure *unwrapImage );

void branchCuts( unwrapImageStructure *unwrapImage );

void unwrapPhase( unwrapImageStructure *unwrapImage );

void interpolatePhase( unwrapImageStructure *unwrapImage );

void addPhaseRamp( unwrapImageStructure *unwrapImage,
                   float omegaR, float omegaA, int noCenter);


void labelRegions( unwrapImageStructure *unwrapImage );
