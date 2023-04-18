#include "stdio.h"
#include "string.h"
#include "unwrapInclude.h"
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>

/*
    Phase unwrapping program.
*/

static void readArgs(int32_t argc, char *argv[], int32_t *azimuthSize,
                     int32_t *rangeSize, float *omegaA, float *omegaR,
                     int32_t *blockSize);

static void usage();

int main(int argc, char *argv[])
{
    ers1ComplexImage *complexImage;   /* Input image */
    unwrapImageStructure unwrapImage; /* Unwrapped image */
    int32_t rangeSize, azimuthSize;   /* Image size */
    char b[10000];
    float tmp[10000];
    int32_t blockSize;
    float omegaR, omegaA; /* Phase ramp parameters */
    int32_t i, j;         /* LCV */
                          /*
                              Read command line args and compute filenames
                          */
    readArgs(argc, argv, &rangeSize, &azimuthSize, &omegaA, &omegaR, &blockSize);
    /*
        Init and Open images
    */
    complexImage = initComplexImage(rangeSize, azimuthSize, 0, BUFFERSIZE);
    openImage((ers1GenericImage *)complexImage, NULL, READIMAGE);
    unwrapImage.rangeSize = rangeSize;
    unwrapImage.azimuthSize = azimuthSize;
    unwrapImage.blockSize = blockSize;
    unwrapImage.phase = (float **)malloc(azimuthSize * sizeof(float *));
    unwrapImage.bCuts = (char **)malloc(azimuthSize * sizeof(char *));
    unwrapImage.labels = (int32_t **)malloc(azimuthSize * sizeof(int32_t *));

    for (i = 0; i < azimuthSize; i++)
    {
        unwrapImage.phase[i] = (float *)malloc(rangeSize * sizeof(float));
        unwrapImage.bCuts[i] = (char *)malloc(rangeSize * sizeof(char));
        unwrapImage.labels[i] = (int32_t *)malloc(rangeSize * sizeof(int));
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
    fprintf(stderr, "branchCuts\n");

    branchCuts(&unwrapImage);
    for (i = 0; i < unwrapImage.azimuthSize; i++)
    {
        for (j = 0; j < unwrapImage.rangeSize - 1; j++)
        {
            b[j] = 0;
            tmp[j] = 0;
            if (unwrapImage.bCuts[i][j] & BRANCHCUT)
                b[j] = 50;
            if (unwrapImage.bCuts[i][j] & PRESIDUE)
                b[j] = 100;
            if (unwrapImage.bCuts[i][j] & NRESIDUE)
                b[j] = (char)150;
            if (unwrapImage.labels[i][j] != LARGEINT)
                tmp[j] = (float)unwrapImage.labels[i][j];
        }
        fwriteBS(b, rangeSize * sizeof(char), 1, stdout, BYTEFLAG);
    }
}

static void readArgs(int32_t argc, char *argv[],
                     int32_t *rangeSize, int32_t *azimuthSize,
                     float *omegaA, float *omegaR, int32_t *blockSize)
{
    int32_t filenameArg;
    char *argString;
    int32_t i;
    if (argc == 1)
        fprintf(stderr, "\nFor usage: unwrap -help\n");
    if (argc < 1 || argc > 13)
        usage();         /* Check number of args */
    *azimuthSize = 1280; /* Default values */
    *rangeSize = 1024;
    *omegaA = 0.0;
    *omegaR = 0.0;
    *blockSize = 5;
    if (argc == 2)
        argc++;
    for (i = 1; i < argc - 1; i += 2)
    {
        argString = strchr(argv[i], '-');
        if (argString != NULL)
        {
            if (strstr(argString, "sr") != NULL)
                sscanf(argv[i + 1], "%i", rangeSize);
            else if (strstr(argString, "sa") != NULL)
                sscanf(argv[i + 1], "%i", azimuthSize);
            else if (strstr(argString, "omegaA") != NULL)
                sscanf(argv[i + 1], "%f", omegaA);
            else if (strstr(argString, "omegaR") != NULL)
                sscanf(argv[i + 1], "%f", omegaR);
            else if (strstr(argString, "blockSize") != NULL)
                sscanf(argv[i + 1], "%i", blockSize);
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
    error("\n%s\n%s\n\n%s\n%s\n\n%s\n",
          "Computes unwrapped phase image from complex interferogram",
          "Complex image input from stdio, unwrapped phase image output to stdio",
          "Usage:  ",
          "unwrap  -sr rangeSize -omegaA -sa azimuthSize omegaA -omegaR omegaR",
          "Defaults: -sa 1280 -sr 1024 -omegaA 0.0 -omegaR 0.0 -blockSize blockSize");
}
