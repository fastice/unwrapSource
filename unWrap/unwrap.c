#include "stdio.h"
#include "string.h"
#include "unwrapInclude.h"
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
/*#include "ers1/getLocC_p/parfile.h"*/
#include "mosaicSource/common/geocode.h"
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
int32_t llConserveMem = 999;

typedef struct OffsetsType
{
    int32_t nr;
    int32_t na;
    int32_t rO;
    int32_t aO;
    float deltaA;
    float deltaR;
    /* Azimuth params */
    double c1;
    double dbcds;
    double dbhds;
    double doffdx;
    /* Range params */
    double bn;
    double bp;
    double dBn;
    double dBp;
    double rConst;
    /* Pixels spacing */
    float **da;
    float **dr;
    char *file;
    char *azParamsFile;
    char *rFile;
    char *rParamsFile;
} Offsets;

int32_t HemiSphere = NORTH;
int32_t DemType = LATLON;
double Rotation = 45.;
double SLat = -91.0;

/*
    Phase unwrapping program.
*/

static void readArgs(int32_t argc, char *argv[],
                     int32_t *rangeSize, int32_t *azimuthSize,
                     float *omegaA, float *omegaR, int32_t *blockSize,
                     int32_t *residueDensity, int32_t *noInterpolate, char **mask,
                     int32_t *noCenter, int32_t *thresh, char **labelFile,
                     char **offsetsFile, char **geodat, char **baselineFile,
                     double *ambThresh, char **reportFile);

void resolveAmb(unwrapImageStructure *unwrapImage, char *offsetsFile,
                char *geodat, char *baselineFile, double ambThresh);
static void readBaseParams(char *file, double *bn, double *bp, double *dbn, double *dbp);
static void readRangeOffsets(Offsets *offsets);
static void usage();

static void doStats(int32_t *nDiff, double *diff, double *diffSq, double *mean, double *errorThresh,
                    inputImageStructure inputImage, unwrapImageStructure *unwrapImage,
                    double bnc, double dbn, double bpc, double dbp,
                    double thetaC, double ReH, double Re, Offsets offsets, double thresh, int32_t *refLabel, double ambThresh);

/*
   Global variables definitions
*/
int32_t RangeSize = RANGESIZE;              /* Range size of complex image */
int32_t AzimuthSize = AZIMUTHSIZE;          /* Azimuth size of complex image */
int32_t BufferSize = BUFFERSIZE;            /* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;                  /* # of lines of nonoverlap in buffer */
double RangePixelSize = RANGEPIXELSIZE;     /* Range PixelSize */
double AzimuthPixelSize = AZIMUTHPIXELSIZE; /* Azimuth PixelSize */

int main(int32_t argc, char *argv[])
{
    FILE *fp;
    ers1ComplexImage *complexImage;   /* Input image */
    unwrapImageStructure unwrapImage; /* Unwrapped image */
    int32_t rangeSize, azimuthSize;   /* Image size */
    int32_t blockSize, residueDensity;
    int32_t noInterpolate;
    char *geodat, *baselineFile, *offsetsFile, *reportFile;
    int32_t thresh;
    int32_t noCenter;
    char *mask;
    char *labelFile;
    double ambThresh;
    float omegaR, omegaA; /* Phase ramp parameters */
    int32_t i, j;         /* LCV */
                          /*
                              Read command line args and compute filenames
                          */
    readArgs(argc, argv, &rangeSize, &azimuthSize, &omegaA, &omegaR,
             &blockSize, &residueDensity, &noInterpolate, &mask,
             &noCenter, &thresh, &labelFile,
             &offsetsFile, &geodat, &baselineFile, &ambThresh, &reportFile);

    if (offsetsFile != NULL)
        fprintf(stderr, "Using Offsets\n%s\n%s\n%s\n",
                offsetsFile, geodat, baselineFile);
    /*
        Init and Open images
    */
    complexImage = initComplexImage(rangeSize, azimuthSize, 0, BUFFERSIZE);
    openImage((ers1GenericImage *)complexImage, NULL, READIMAGE);
    unwrapImage.rangeSize = rangeSize;
    unwrapImage.azimuthSize = azimuthSize;
    unwrapImage.blockSize = blockSize;
    unwrapImage.labelFile = labelFile;
    unwrapImage.thresh = thresh;
    unwrapImage.mask = mask;
    unwrapImage.cutFile = NULL; /* Set this to a file name for debugging */

    unwrapImage.residueDensity = residueDensity;
    unwrapImage.phase = (float **)malloc(azimuthSize * sizeof(float *));
    unwrapImage.bCuts = (char **)malloc(azimuthSize * sizeof(char *));
    unwrapImage.labels = (int32_t **)malloc(azimuthSize * sizeof(int32_t *));

    for (i = 0; i < azimuthSize; i++)
    {
        unwrapImage.phase[i] = (float *)malloc(rangeSize * sizeof(float));
        unwrapImage.bCuts[i] = (char *)malloc(rangeSize * sizeof(char));
        unwrapImage.labels[i] = (int32_t *)malloc(rangeSize * sizeof(int));
    }
    if (reportFile != NULL)
    {
        fprintf(stderr, "Opening report file %s\n", reportFile);
        unwrapImage.fpReport = fopen(reportFile, "w");
    }
    else
        unwrapImage.fpReport = NULL;

    if (unwrapImage.mask != NULL)
    {
        fprintf(stderr, "USING MASK IMAGE %s\n", mask);
        if (unwrapImage.fpReport != NULL)
            fprintf(unwrapImage.fpReport, ";  Masked with %s\n", mask);
    }
    /*
        Compute phaseimage and free complex image.
    */
    fprintf(stderr, "phaseImage\n");
    phaseImage(&unwrapImage, complexImage, omegaR, omegaA);

    freeImage((ers1GenericImage *)complexImage);
    /*
        Compute residues
    */
    fprintf(stderr, "residueImage\n");
    residueImage(&unwrapImage);
    /*
        Draw branch cuts
    */
    fprintf(stderr, "branchCuts %f %f\n", (float)LARGEINT, (double)LARGEINT);
    /*
       skip branch cuts if no residues and no mask
    */
    if (unwrapImage.mask != NULL ||
        (unwrapImage.nPlus > 0) ||
        (unwrapImage.nPlus > 0))
        branchCuts(&unwrapImage);
    /*
       Unwrap phase
    */
    fprintf(stderr, "unwrapPhase\n");
    unwrapPhase(&unwrapImage);

    if (offsetsFile != NULL)
        resolveAmb(&unwrapImage, offsetsFile, geodat, baselineFile, ambThresh);
    /*
       Interpolate phase.
    */
    if (noInterpolate == FALSE)
    {
        fprintf(stderr, "InterpolatePhase\n");
        interpolatePhase(&unwrapImage);
    }
    /*
       Add phase ramp
    */
    fprintf(stderr, "addPhaseRamp\n");
    addPhaseRamp(&unwrapImage, omegaR, omegaA, noCenter);

    /*
        Output total phase image
    */
    for (i = 0; i < azimuthSize; i++)
        fwriteBS(unwrapImage.phase[i], rangeSize * sizeof(float), 1, stdout, FLOAT32FLAG);
    /*
        Output label file
    */
    if (labelFile != NULL)
    {
        fprintf(stderr, "Writing label file:  %s %i %i\n", labelFile, azimuthSize, rangeSize);
        fp = fopen(labelFile, "w");
        for (i = 0; i < azimuthSize; i++)
            fwriteBS(unwrapImage.labels[i], rangeSize * sizeof(int32_t), 1, fp, INT32FLAG);
    }

    /*return; */
}

void resolveAmb(unwrapImageStructure *unwrapImage, char *offsetsFile,
                char *geodat, char *baselineFile, double ambThresh)
{
    Offsets offsets;
    double bnc, bpc, dbn, dbp;
    double bn, bp, bsq;
    inputImageStructure inputImage;
    int32_t azSlcSize;
    int32_t r, a;
    int32_t r1, a1;
    int32_t i, j;
    int32_t intCorrection;
    double correction;
    double range, x;
    double delta, deltaOffset;
    double *diff, *diffSq;
    double refDif;
    double *errorThresh;
    int32_t refLabel, maxCount;
    int32_t *nDiff;
    double Re, ReH, theta, thetaC, thetaD;
    double phase;
    double scalePhase;
    double thresh;
    double *mean, sq, sigma;
    int32_t currentLabel;
    int32_t labelCount;
    fprintf(stderr, "Read geodat file \n");
    inputImage.stateFlag = TRUE;
    inputImage.passType = DESCENDING;

    parseInputFile(geodat, &inputImage);
    initllToImageNew(&inputImage);

    offsets.rFile = offsetsFile;
    fprintf(stderr, "Read range offsets \n");
    readRangeOffsets(&offsets);
    fprintf(stderr, "Read baseline params\n");
    readBaseParams(baselineFile, &bnc, &bpc, &dbn, &dbp);

    ReH = inputImage.cpAll.Re + inputImage.par.H;
    Re = inputImage.cpAll.Re;
    thetaC = acos((pow(inputImage.cpAll.RCenter, 2.0) + ReH * ReH - Re * Re) /
                  (2.0 * inputImage.cpAll.RCenter * ReH));
    scalePhase = (1.0 / (2.0 * PI)) * (0.5) * inputImage.par.lambda;
    fprintf(stderr, "Re = %lf\n", Re);
    fprintf(stderr, "ReH = %lf\n", ReH);
    fprintf(stderr, "ThetaD = %lf\n", thetaC * RTOD);
    fprintf(stderr, "Scale Phase %lf\n", scalePhase);

    /*
       Malloc and initialize
    */
    nDiff = (int32_t *)malloc(sizeof(int) * (unwrapImage->nLabels + 1));
    diff = (double *)malloc(sizeof(double) * (unwrapImage->nLabels + 1));
    diffSq = (double *)malloc(sizeof(double) * (unwrapImage->nLabels + 1));
    errorThresh = (double *)malloc(sizeof(double) * (unwrapImage->nLabels + 1));
    mean = (double *)malloc(sizeof(double) * (unwrapImage->nLabels + 1));
    for (i = 0; i <= unwrapImage->nLabels; i++)
    {
        nDiff[i] = 0;
        diff[i] = 0.0;
        diffSq[i] = 0.;
        errorThresh[i] = 1e30;
        mean[i] = 0;
    }
    fprintf(stderr, "Max Label %i\n", unwrapImage->nLabels);
    /*
        Compute differences
     */
    thresh = -1.0;
    fprintf(stderr, "Starting do stats\n");
    doStats(nDiff, diff, diffSq, mean, errorThresh, inputImage, unwrapImage,
            bnc, dbn, bpc, dbp, thetaC, ReH, Re, offsets, thresh, &refLabel, ambThresh);
    fprintf(stderr, "Finished do stats\n");
    /*
       Fix phase
    */
    refDif = mean[refLabel] / scalePhase;
    fprintf(stderr, "%f %f\n", mean[refLabel], refDif);
    for (i = 0; i < inputImage.azimuthSize; i++)
    {
        for (j = 0; j < inputImage.rangeSize; j++)
        {
            /* Only change labels that are valid */
            if (unwrapImage->labels[i][j] >= 0 &&
                unwrapImage->labels[i][j] <= unwrapImage->nLabels)
            {
                /* Change or invalidate */
                if (mean[unwrapImage->labels[i][j]] > (-LARGEINT + 10000))
                {
                    /* Correction, everything references to refDif */
                    correction = (refDif - mean[unwrapImage->labels[i][j]] / scalePhase) / (2.0 * PI);
                    if (correction > 0)
                        intCorrection = (int)(correction + 0.5);
                    else
                        intCorrection = (int)(correction - 0.5);
                    correction = intCorrection * 2.0 * PI; /* to nearest ambiguity */
                    unwrapImage->phase[i][j] += correction;
                }
                else
                    unwrapImage->phase[i][j] = -(float)LARGEINT;
            }
        }
    }
}

/* ******************************************************************************************/
static void doStats(int32_t *nDiff, double *diff, double *diffSq, double *mean, double *errorThresh,
                    inputImageStructure inputImage, unwrapImageStructure *unwrapImage,
                    double bnc, double dbn, double bpc, double dbp,
                    double thetaC, double ReH, double Re, Offsets offsets, double thresh, int32_t *refLabel, double ambThresh)
{
    double sq, sigma;
    double x;
    int32_t a, r;
    int32_t i, j;
    double bn, bp, bsq;
    int32_t labelCount, currentLabel, maxCount;
    double theta, thetaD;
    double delta, phase, deltaOffset, scalePhase;
    double range;
    float tmp, tmp1, tmp2;
    int32_t totalGood, maxLabel, maxGood;
    FILE *fp, *fp1, *fp2;
    float zero;
    zero = 0.0;

    scalePhase = (1.0 / (2.0 * PI)) * (0.5) * inputImage.par.lambda;
    for (i = 0; i <= unwrapImage->nLabels; i++)
    {
        nDiff[i] = 0;
        diff[i] = 0.0;
        diffSq[i] = 0.;
    }
    for (i = 0; i < offsets.na; i++)
    {
        a = (int)(offsets.aO + i * offsets.deltaA) / inputImage.nAzimuthLooks;
        x = (a - (inputImage.azimuthSize / 2.)) / (double)inputImage.azimuthSize;
        /* Baseline */
        bn = bnc + x * dbn;
        bp = bpc + x * dbp;
        bsq = bn * bn + bp * bp;
        /* Loop on range */
        for (j = 0; j < offsets.nr; j++)
        {
            r = (int)(offsets.rO + j * offsets.deltaR) / inputImage.nRangeLooks;

            /* Get label count */
            if (unwrapImage->labels[a][r] >= 0 &&
                unwrapImage->labels[a][r] <= unwrapImage->nLabels)
                labelCount = unwrapImage->labelCount[unwrapImage->labels[a][r]];
            else
                labelCount = -1;

            /* Only proceed if enough hits */
            if (labelCount > 25000)
            {
                /* Compute flat earth phase */
                range = inputImage.cpAll.RNear + r * inputImage.rangePixelSize;
                theta = acos((range * range + ReH * ReH - Re * Re) / (2.0 * (ReH)*range));
                thetaD = theta - thetaC;
                delta = -bn * sin(thetaD) - bp * cos(thetaD) + (0.5 * bsq / range);
                if (a > 0 && a < inputImage.azimuthSize &&
                    r > 0 && r < inputImage.rangeSize)
                {
                    phase = unwrapImage->phase[a][r];
                    currentLabel = unwrapImage->labels[a][r];
                }
                else
                    phase = -(float)LARGEINT;
                /* Convert  valid phase to meters */
                if (phase > (-LARGEINT + 10000))
                    delta = delta + scalePhase * phase;
                else
                    delta = -(double)LARGEINT;
                /* Convert valid offset to meters */
                if (offsets.dr[i][j] > (-LARGEINT + 10000))
                {
                    deltaOffset = offsets.dr[i][j] *
                                  inputImage.rangePixelSize / inputImage.nRangeLooks;
                }
                else
                    deltaOffset = -(double)LARGEINT;
                /* Accumulate difference for pixels with both valid phase and offsets */
                if (delta > (-LARGEINT + 10000) && deltaOffset > (-LARGEINT + 10000))
                {
                    if ((fabs(delta - deltaOffset - mean[currentLabel]) < thresh || thresh < 0))
                    {
                        diff[currentLabel] += delta - deltaOffset;
                        diffSq[currentLabel] += (delta - deltaOffset) * (delta - deltaOffset);
                        nDiff[currentLabel]++;
                    }
                }
                tmp = (float)delta;
                tmp1 = (float)deltaOffset;
                tmp2 = (float)range;
            }
        }
    }
    /*
      Go through results, compute mean, kill of ones with bad sigma
    */
    maxCount = 0;
    totalGood = 0;
    maxGood = 0;
    maxLabel = 0;
    if (unwrapImage->fpReport != NULL)
        fprintf(unwrapImage->fpReport, "\n; lambda = %f\n;Component      size    relative Size    nDiff    mean         sq       sigma     thresh\n", inputImage.par.lambda);
    for (i = 0; i <= unwrapImage->nLabels; i++)
    {
        if (nDiff[i] > 1000)
        {
            mean[i] = diff[i] / (double)nDiff[i];
            sq = diffSq[i] / (double)nDiff[i];
            sigma = sqrt(sq - mean[i] * mean[i]);
            errorThresh[i] = sigma;
            /*  Don't use if sigma not beter than .2  */
            if (sigma > ambThresh)
            {
                nDiff[i] = 0;
                mean[i] = -(double)LARGEINT;
            }
            else
            {
                totalGood += unwrapImage->labelCount[i];
                maxGood = max(maxGood, unwrapImage->labelCount[i]);
            }
            maxLabel = max(maxLabel, unwrapImage->labelCount[i]);
            if (nDiff[i] > maxCount)
            {
                maxCount = nDiff[i];
                *refLabel = i;
            }
            if (unwrapImage->fpReport != NULL)
                fprintf(unwrapImage->fpReport, "%10i %10i %10.2f %% %10i %10.5f %10.5f %10.5f %10.5f\n", i, unwrapImage->labelCount[i],
                        (float)unwrapImage->labelCount[i] / (float)(unwrapImage->azimuthSize * unwrapImage->rangeSize) * 100.0, nDiff[i], mean[i], sq, sigma, ambThresh);

            fprintf(stderr, "%i %i %i %lf %lf %lf %lf\n", i, unwrapImage->labelCount[i], nDiff[i], mean[i], sq, sigma, ambThresh);
        }
        else
        {
            mean[i] = -(double)LARGEINT;
        }
    }
    if (unwrapImage->fpReport != NULL)
    {
        fprintf(unwrapImage->fpReport, "Unwrapped %%       %10.2f %%\n", (float)totalGood / (float)(unwrapImage->azimuthSize * unwrapImage->rangeSize) * 100.0);
        fprintf(unwrapImage->fpReport, "Largest Good %%    %10.2f %%\n", (float)maxGood / (float)(unwrapImage->azimuthSize * unwrapImage->rangeSize) * 100.0);
        fprintf(unwrapImage->fpReport, "Largest Overall %% %10.2f %%\n", (float)maxLabel / (float)(unwrapImage->azimuthSize * unwrapImage->rangeSize) * 100.0);
    }
}

/*
   Read baseline information
*/

static void readBaseParams(char *file, double *bn, double *bp, double *dbn, double *dbp)
{
    int32_t lineCount = 0, eod;
    char line[256];
    double dum;
    FILE *fp;

    fprintf(stderr, "\nREADING BASELINE INFORMATION %s\n", file);

    if (file != NULL)
    {
        fp = openInputFile(file);
        lineCount = getDataString(fp, lineCount, line, &eod);
        /*
            Read baseline params
        */
        if (sscanf(line, "%lf%lf%lf%lf%lf", bn, bp, dbn, &dum, dbp) != 5)
            error("%s  %i", "parseBase -- Missing baseline at line:", lineCount);
    }
    else if (file != NULL)
    {
        error("parseBase: baseFile not implemented yet");
    }
    else
        error("parseBase: missing baseline estimate/parameter file");

    fprintf(stderr, "Bn,Bp,dBn,dBp %f %f %f %f\n", *bn, *bp, *dbn, *dbp);
}

static void readRangeOffsets(Offsets *offsets)
{
    char *datFile, buf[1024];
    char line[1024];
    int32_t lineCount, eod;
    FILE *fp;
    int32_t rO, aO, nr, na;
    float deltaA, deltaR;
    double c1;
    float *fBufR;
    char *file;
    int32_t i;
    /*
       Read inputfile
    */
    file = offsets->rFile;
    datFile = &(buf[0]);
    buf[0] = '\0';

    datFile = strcat(datFile, file);
    datFile = strcat(datFile, ".dat");

    fp = openInputFile(datFile);
    lineCount = 0;
    lineCount = getDataString(fp, lineCount, line, &eod);

    if (sscanf(line, "%i%i%i%i%f%f", &rO, &aO, &nr, &na, &deltaR, &deltaA) != 6)
        error("%s  %i of %s",
              "readOffsets -- Missing image parameters at line:",
              lineCount, datFile);
    fclose(fp);

    offsets->nr = nr;
    offsets->na = na;
    offsets->rO = rO;
    offsets->aO = aO;
    offsets->deltaA = deltaA;
    offsets->deltaR = deltaR;
    offsets->dr = (float **)malloc(sizeof(float *) * na);
    /*
       Malloc buffers
     */
    fBufR = (float *)malloc(sizeof(float) * na * nr);
    for (i = 0; i < na; i++)
        offsets->dr[i] = &(fBufR[i * nr]);
    /*
        Read files
     */
    fprintf(stderr, "--- Reading offset file %s", offsets->rFile);

    fp = openInputFile(offsets->rFile);
    for (i = 0; i < na; i++)
        freadBS(offsets->dr[i], sizeof(float), nr, fp, FLOAT32FLAG);
    fclose(fp);
}

static void readArgs(int32_t argc, char *argv[], int32_t *rangeSize, int32_t *azimuthSize,
                     float *omegaA, float *omegaR, int32_t *blockSize, int32_t *residueDensity, int32_t *noInterpolate, char **mask,
                     int32_t *noCenter, int32_t *thresh, char **labelFile, char **offsets, char **geodat, char **baselineFile,
                     double *ambThresh, char **reportFile)
{
    int32_t filenameArg;
    char *argString;
    int32_t i;
    if (argc == 1)
        fprintf(stderr, "\nFor usage: unwrap -help\n");
    if (argc < 1 || argc > 25)
        usage();         /* Check number of args */
    *azimuthSize = 1280; /* Default values */
    *rangeSize = 1024;
    *omegaA = 0.0;
    *omegaR = 0.0;
    *labelFile = NULL;
    *blockSize = 5;
    *noInterpolate = FALSE;
    *residueDensity = 0;
    *noCenter = FALSE;
    *mask = NULL;
    *thresh = -1;
    *ambThresh = 0.2;
    *offsets = NULL;
    *geodat = NULL;
    *baselineFile = NULL;
    *reportFile = NULL;
    if (argc == 2)
        argc++;
    for (i = 1; i < argc - 1; i += 2)
    {
        argString = strchr(argv[i], '-');
        if (argString != NULL)
        {
            if (strstr(argString, "rs") != NULL)
                sscanf(argv[i + 1], "%i", rangeSize);
            else if (strstr(argString, "as") != NULL)
                sscanf(argv[i + 1], "%i", azimuthSize);
            else if (strstr(argString, "omegaA") != NULL)
                sscanf(argv[i + 1], "%f", omegaA);
            else if (strstr(argString, "omegaR") != NULL)
                sscanf(argv[i + 1], "%f", omegaR);
            else if (strstr(argString, "blockSize") != NULL)
                sscanf(argv[i + 1], "%i", blockSize);
            else if (strstr(argString, "thresh") != NULL)
                sscanf(argv[i + 1], "%i", thresh);
            else if (strstr(argString, "ambThresh") != NULL)
                sscanf(argv[i + 1], "%lf", ambThresh);
            else if (strstr(argString, "reportFile") != NULL)
            {
                *reportFile = argv[i + 1];
                fprintf(stderr, "Report File %s\n", *reportFile);
            }
            else if (strstr(argString, "labelFile") != NULL)
            {
                *labelFile = argv[i + 1];
                fprintf(stderr, "%s\n", *labelFile);
            }
            else if (strstr(argString, "residueDensity") != NULL)
                sscanf(argv[i + 1], "%i", residueDensity);
            else if (strstr(argString, "msk") != NULL)
            {
                *mask = argv[i + 1];
            }
            else if (strstr(argString, "noInterpolate") != NULL)
            {
                *noInterpolate = TRUE;
                i--;
            }
            else if (strstr(argString, "noCenter") != NULL)
            {
                *noCenter = TRUE;
                i--;
            }
            else if (strstr(argString, "useOffsets") != NULL)
            {
                *offsets = argv[i + 1];
                *geodat = argv[i + 2];
                *baselineFile = argv[i + 3];
                i += 2;
            }
            else if (strstr(argString, "help") != NULL)
                usage();
            else
                usage();
        }
        else
            usage();
    }
}

static void usage()
{
    error("\n%s\n%s\n\n%s\n%s\n%s\n%s\n%s\n\n%s\n%s\n",
          "Computes unwrapped phase image from complex interferogram",
          "Complex image input from stdio, unwrapped phase image output to stdio",
          "Usage:  ",
          "unwrap  -msk maskfile -rs rangeSize -as azimuthSize\n"
          "useOffsets offsets geodat baselineFile (params)\n"
          "-omegaA  omegaA -omegaR omegaR -noCenter -ambThresh ambThresh",
          "-labelFile labelFile -reportFile reportFile -thresh thresh",
          "        -blockSize blockSize -residueDensity residueDensity",
          "        -noInterpolate",
          "Defaults: -as 1280 -rs 1024 -omegaA 0.0 -omegaR 0.0 -blockSize 5",
          "          -residueDensity 0");
}
