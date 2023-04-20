#include "ers1.h"
#include "cRecipes/cRecipes.h"
#include <math.h>
#include <stdlib.h>
/*
     Initialize kernel.
*/
static void initKernel(ers1Kernel *kernel);

/*
    Compute shifts as function of range using either deltaR and deltaA
    or range dependent shifts in shiftFile.
*/

void initShifts(ers1Shift *shifts, char *shiftFile, float deltaR, float deltaA,
                ers1ComplexImage *imageA, ers1Kernel *kernel,
                float *dR1, float *dA1, int32_t modeFlag)
{
    float rFract, aFract; /* Float parts of shift */
    int32_t rInt, aInt;   /* Integer parts of shifts */
    float minAzimuth, maxAzimuth;
    float *dR, *dA;
    float firstDeltaR;
    int32_t i; /* LCV */
               /*
                   Allocate arrays
               */
    dR = (float *)malloc(RangeSize * sizeof(float));
    dA = (float *)malloc(RangeSize * sizeof(float));

    if (modeFlag == INIT)
    {
        shifts->columnShift = (int32_t *)malloc(RangeSize * sizeof(int));
        shifts->rowShift = (int32_t *)malloc(RangeSize * sizeof(int));
        shifts->aKernel = (int32_t *)malloc(RangeSize * sizeof(int));
        shifts->rKernel = (int32_t *)malloc(RangeSize * sizeof(int));
        /*
            Initialize kernel
        */
        initKernel(kernel);
        if (shiftFile != NULL && shifts->shiftFlag == RANGEONLY)
            readShifts(shiftFile, dR, dA, NULL);
    } /* End if init */
      /*
          Init shifts for each range.
      */
    minAzimuth = 9999;
    maxAzimuth = -9999;
    shifts->endRange = -1;
    shifts->startRange = 0;
    /*
        Loop to compute shifts
    */
    for (i = 0; i < RangeSize; i++)
    {
        if (shiftFile != NULL && shifts->shiftFlag == RANGEONLY)
        {
            deltaR = dR[i];
            deltaA = dA[i];
        }
        else if (shifts->shiftFlag == RANGEANDAZIMUTH)
        {
            deltaR = dR1[i];
            deltaA = dA1[i];
        }
        /*
            Compute range shifts
        */
        rInt = (int32_t)deltaR;
        rFract = deltaR - (float)rInt;
        if ((float)fabs((double)rFract) < 0.51 * DSHIFT)
            rFract = 0.0;
        else if (rFract < 0)
        { /* Always use neg fract shift */
            rInt--;
            rFract = -(1.0 + rFract);
        }
        else
            rFract = -rFract;
        shifts->columnShift[i] = -rInt; /* Column, range */
        if (i == 0)
            firstDeltaR = deltaR;
        shifts->rKernel[i] = (int)((-rFract + 0.5 * DSHIFT) * NKERNEL);
        if (shifts->rKernel[i] >= NKERNEL)
            shifts->rKernel[i] = NKERNEL - 1;
        /*
            Compute Azimuth shifts
        */
        aInt = (int32_t)deltaA;
        aFract = deltaA - (float)aInt;
        if (fabs((double)aFract) < 0.51 * DSHIFT)
            aFract = 0.0;
        else if (aFract < 0)
        { /* Always use neg fract shift */
            aInt--;
            aFract = -(1.0 + aFract);
        }
        else
            aFract = -aFract;
        shifts->rowShift[i] = -aInt; /* Row, azimuth */
        shifts->aKernel[i] = (int)((-aFract + 0.5 * DSHIFT) * NKERNEL);
        if (shifts->aKernel[i] >= NKERNEL)
            shifts->aKernel[i] = NKERNEL - 1;
        /*
           Get min and max azimuth shifts.
        */
        minAzimuth = min(minAzimuth, deltaA);
        maxAzimuth = max(maxAzimuth, deltaA);
        /*
           Precompute whether row column coordinates for interpolation
        */
    }
    /*
        Compute start and Stop range
    */
    if (firstDeltaR > 0)
        shifts->startRange = (int)firstDeltaR + 1;
    else
        shifts->startRange = 0;
    if (firstDeltaR < 0)
        shifts->endRange = (int)(imageA->header.rs + firstDeltaR);
    else
        shifts->endRange = imageA->header.rs;

    /*
       Compute start and end azimuth
    */
    if (maxAzimuth > 0)
        shifts->startAzimuth = (int)maxAzimuth + 1;
    else
        shifts->startAzimuth = 0;

    if (minAzimuth < 0)
        shifts->endAzimuth =
            (int)(imageA->header.as + minAzimuth);
    else
        shifts->endAzimuth = imageA->header.as;

    fprintf(stderr, "Column, row = %i %i\n", shifts->columnShift[0], shifts->rowShift[0]);
    fprintf(stderr, "rFract,aFract %f %f \n", rFract, aFract);
    fprintf(stderr, "start,endRange = %i  %i\n", shifts->startRange, shifts->endRange);
    fprintf(stderr, "start,endAzimuth = %i  %i\n\n", shifts->startAzimuth, shifts->endAzimuth);

    /*
        Free memory and return
    */
    free(dR);
    free(dA);
    return;
}

/*
     Initialize kernel.
*/
static void initKernel(ers1Kernel *kernel)
{
    float shift;
    float *a;
    struct fComplex
    {
        float r;
        float i;
    } data[32];
    int32_t i, j, k, count;
    double exponent, over32;

    over32 = 1.0 / 32.0;
    count = 0;
    *kernel = (ers1Kernel)malloc(KERNELSIZE * NKERNEL * sizeof(float));
    *kernel = (ers1Kernel) & ((*kernel)[0][KERNELSIZE / 2]);
    for (shift = 0; shift < 1.0; shift += DSHIFT)
    {
        for (i = 0; i <= 16; i++)
        { /* Load freq values of shift */
            exponent = (double)i * 6.28319 * (double)shift * over32;
            data[i].r = (float)cos(exponent);
            data[i].i = (float)sin(exponent);
            if (i > 0 && i < 16)
            {
                j = 32 - i;
                data[j].r = (float)cos(exponent);
                data[j].i = (float)sin(-exponent);
            }
        }
        /*
           Set up point for numerical recipes start with one convention
        */
        a = (float *)data;
        a--;
        four1(a, 32, 1);
        /*
           Initialize kernel;
        */
        for (k = 0; k < KERNELSIZE / 2; k++)
        {
            (*kernel)[count][k] = data[k].r * over32;
            (*kernel)[count][-k - 1] = data[31 - k].r * over32;
        }
        count++;
    }
    return;
}
