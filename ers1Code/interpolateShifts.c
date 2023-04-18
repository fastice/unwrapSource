#include "ers1.h"
#include <math.h>
#include <stdlib.h>
/*
   Interpolate registration data to get shifts for given azimuth (dr,da).
*/

void interpolateShifts(int32_t az, ers1Shift *shifts, float *dR1, float *dA1)
{
    int32_t minAz, maxAz, minRg, maxRg;
    float *tmpA;
    float *tmpR;
    float a1, a2;
    int32_t i1, i;
    /*
        Malloc arrays
    */
    tmpA = (float *)malloc(RangeSize * sizeof(float));
    tmpR = (float *)malloc(RangeSize * sizeof(float));
    /*
       Compute min and max range and azimuth values
    */
    minAz = shifts->gAo;
    maxAz = minAz + shifts->dAg * (shifts->gNa - 1);
    minRg = shifts->gRo;
    maxRg = minRg + shifts->dRg * (shifts->gNr - 1);
    /*
       Use edge values if of grid
    */
    if (az <= minAz)
    {
        for (i = 0; i < shifts->gNa; i++)
        {
            tmpA[i] = shifts->azimuthShifts[0][i];
            tmpR[i] = shifts->rangeShifts[0][i];
        }
    }
    else if (az >= maxAz)
    {
        for (i = 0; i < shifts->gNa; i++)
        {
            tmpA[i] = shifts->azimuthShifts[shifts->gNa - 1][i];
            tmpR[i] = shifts->rangeShifts[shifts->gNa - 1][i];
        }
    }
    else
    { /* Interpolate between grid rows */
        a1 = (az - minAz) / (float)shifts->dAg;
        i1 = (int)a1;
        a2 = a1 - i1;
        a1 = 1.0 - a2;
        /* fprintf(stderr,"%f  %f    %i %i\n",a1,a2,i1,az); */
        if (i1 < 0 || i1 > (shifts->gNa - 2))
            error("*** getAzVaryingShifts ***");
        for (i = 0; i < shifts->gNa; i++)
        {
            tmpA[i] = shifts->azimuthShifts[i1][i] * a1 +
                      shifts->azimuthShifts[i1 + 1][i] * a2;
            tmpR[i] = shifts->rangeShifts[i1][i] * a1 +
                      shifts->rangeShifts[i1 + 1][i] * a2;
        }
    }
    /*
       Off left of grid
    */
    for (i = 0; i < minRg; i++)
    {
        dR1[i] = tmpR[0];
        dA1[i] = tmpA[0];
    }
    /*
      On grid, interpolate between grid points
    */
    for (i = minRg; i < maxRg; i++)
    {
        a1 = (i - minRg) / shifts->dRg;
        i1 = (int)a1;
        a2 = a1 - i1;
        a1 = 1.0 - a2;
        if (i1 < 0 || i1 > (shifts->dRg - 2))
            error("*** getAzVaryingShifts 1***");
        dR1[i] = tmpR[i1] * a1 + tmpR[i1 + 1] * a2;
        dA1[i] = tmpA[i1] * a1 + tmpA[i1 + 1] * a2;
    }
    /*
        Off right side  of grid
    */
    for (i = maxRg; i < RangeSize; i++)
    {
        dR1[i] = tmpR[shifts->gNr - 1];
        dA1[i] = tmpA[shifts->gNr - 1];
    }
    /*
        Free memory and return
    */
    free(tmpR);
    free(tmpA);
    return;
}
