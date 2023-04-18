#include "unwrapInclude.h"
#include <math.h>
#include <stdlib.h>
/*
     Compute phase image for unwrapping program.
*/

void phaseImage(unwrapImageStructure *unwrapImage,
                ers1ComplexImage *complexImage,
                float omegaR, float omegaA)
{
    ers1Complex *rowA;
    int32_t *rCount, *aCount;
    char **mask;
    struct
    {
        double r;
        double i;
    } tmp;
    double pr, pi;
    double phaseRamp;
    FILE *fp;
    int32_t i, j;
    int32_t left, right, up, down;
    int32_t minUpi, maxDowni, maxLeftj, minRightj;
    left = 0.05 * complexImage->header.rs;
    right = 0.95 * complexImage->header.rs;
    down = 0.05 * complexImage->header.as;
    up = 0.95 * complexImage->header.as;
    rCount = (int32_t *)malloc(complexImage->header.rs * sizeof(int));
    for (i = 0; i < complexImage->header.rs; i++)
        rCount[i] = 0;
    aCount = (int32_t *)malloc(complexImage->header.as * sizeof(int));
    for (i = 0; i < complexImage->header.as; i++)
        aCount[i] = 0;

    /*
        Input mask image if there is one
    */
    if (unwrapImage->mask != NULL)
    {
        /*
           Malloc space for mask image
        */
        mask = (char **)malloc(complexImage->header.as * sizeof(char *));
        for (i = 0; i < complexImage->header.as; i++)
            mask[i] = (char *)malloc(complexImage->header.rs * sizeof(char));
        /*
           Read file
        */
        fp = openInputFile(unwrapImage->mask);
        for (i = 0; i < complexImage->header.as; i++)
            freadBS(mask[i], complexImage->header.rs * sizeof(char), 1, fp, BYTEFLAG);
        fprintf(stderr, "Input Mask\n");
    }

    for (i = 0; i < complexImage->header.as; i++)
    {
        /*
            Get row Pointers.
        */
        rowA = (ers1Complex *)rowPtr((ers1GenericImage *)complexImage);
        for (j = 0; j < complexImage->header.rs; j++)
        {
            unwrapImage->bCuts[i][j] = 0;
            /*
               Compute phase for valid pixels
            */
            if ((abs(rowA[j].r) < 1) && (abs(rowA[j].i) < 1))
            {

                unwrapImage->bCuts[i][j] |= VALIDPIXEL; /* Mark valid pixel */
                pr = cos(phaseRamp);
                pi = sin(phaseRamp);
                unwrapImage->phase[i][j] = (float)atan2(pi, pr);
                rCount[j]++;
                aCount[i]++; /* Counters for finding edge */
            }
            else
            { /* Compute phase */

                unwrapImage->bCuts[i][j] |= VALIDPIXEL; /* Mark valid pixel */
                phaseRamp = (double)(omegaR * (float)j + omegaA * (float)i);
                pr = cos(phaseRamp);
                pi = sin(phaseRamp);
                tmp.r = (double)rowA[j].r * pr - (double)rowA[j].i * pi;
                tmp.i = (double)rowA[j].r * pi + (double)rowA[j].i * pr;
                unwrapImage->phase[i][j] = (float)atan2(tmp.i, tmp.r);
                /*
                   Apply mask
                */
                if ((unwrapImage->mask != NULL) && (mask[i][j] == 0))
                    unwrapImage->bCuts[i][j] |= MASK;
            }
        }
        rowIncrement((ers1GenericImage *)complexImage);
    }
    /*
       Find borders
    */
    maxLeftj = -1;
    minRightj = complexImage->header.rs;
    maxDowni = -1;
    minUpi = complexImage->header.as;
    i = 0;
    while (rCount[i] >= 0.9 * complexImage->header.as && i <
                                                             complexImage->header.rs)
    {
        maxLeftj++;
        i++;
    }
    i = 0;
    while (aCount[i] >= 0.9 * complexImage->header.rs && i <
                                                             complexImage->header.as)
    {
        maxDowni++;
        i++;
    }
    i = complexImage->header.rs - 1;
    while (rCount[i] >= 0.9 * complexImage->header.as && i > 0)
    {
        minRightj--;
        i--;
    }
    i = complexImage->header.as - 1;
    while (aCount[i] >= 0.9 * complexImage->header.rs && i > 0)
    {
        minUpi--;
        i++;
    }
    maxLeftj++;
    minRightj--;
    maxDowni++;
    minUpi--;
    unwrapImage->left = maxLeftj;
    unwrapImage->right = minRightj;
    unwrapImage->up = minUpi;
    unwrapImage->down = maxDowni;
    free(rCount);
    free(aCount);
    fprintf(stderr, "l,r,d,u  %i  %i  %i  %i\n", maxLeftj, minRightj, maxDowni, minUpi);
    /*
         Throw out border pixels.
    */

    for (i = 0; i < complexImage->header.as; i++)
    {
        for (j = 0; j <= maxLeftj; j++)
        {
            unwrapImage->bCuts[i][j] &= INVALIDPIXEL; /* Mark invalid pix */
            unwrapImage->bCuts[i][j] &= ~MASK;
            unwrapImage->phase[i][j] = 0.0;
        }
        for (j = minRightj; j < complexImage->header.rs; j++)
        {
            unwrapImage->bCuts[i][j] &= INVALIDPIXEL; /* Mark invalid pix */
            unwrapImage->bCuts[i][j] &= ~MASK;
            unwrapImage->phase[i][j] = 0.0;
        }
    }

    for (j = 0; j < complexImage->header.rs; j++)
    {
        for (i = 0; i <= maxDowni; i++)
        {
            unwrapImage->bCuts[i][j] &= INVALIDPIXEL; /* Mark invalid pix */
            unwrapImage->bCuts[i][j] &= ~MASK;
            unwrapImage->phase[i][j] = 0.0;
        }
    }

    for (j = 0; j < complexImage->header.rs; j++)
    {
        for (i = minUpi; i < complexImage->header.as; i++)
        {
            unwrapImage->bCuts[i][j] &= INVALIDPIXEL; /* Mark invalid pix */
            unwrapImage->bCuts[i][j] &= ~MASK;
            unwrapImage->phase[i][j] = 0.0;
        }
    }

    if (unwrapImage->mask != NULL)
    {
        for (i = 0; i < complexImage->header.as; i++)
            free(mask[i]);
        free(mask);
    }
    return;
}
