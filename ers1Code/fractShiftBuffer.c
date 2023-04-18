#include "ers1.h"
#include <stdlib.h>
/*typedef ers1Complex (*tmpBuffer)[RangeSize + KERNELSIZE];
typedef ers1Complex (*iBuffer)[RangeSize];
*/
/*
   Shift entire buffer by fraction of a pixels.
*/
void fractShiftBuffer(ers1GenericImage *image, ers1Shift shifts,
                      ers1Kernel kernel)
{
    headerInfo *header;       /* Image header */
    ers1Complex *imageBuffer; /* Image Buffer */
    ers1Complex *tImage, *tSave;
    float resultR, resultI;
    int32_t index2D;
    int32_t tImageSize;
    int32_t pixelSize;
    int32_t nRows, nCols; /* Number of rows and columns in buffer */
    int32_t kI;
    int32_t colOffset;
    ers1Complex zero;
    int32_t i, j, k;

    header = &(image->header);

    if (header->currentRow < header->lastBufferRow)
        return; /* Return if buffer not full */
                /*
                    Do only if output image
                */
    if (header->imageMode == READIMAGE)
        error("*** fractShiftBuffer: image READ type ***");
    zero.r = 0;
    zero.i = 0;
    pixelSize = sizeof(ers1Complex);
    imageBuffer = (ers1Complex *)image->image;
    nCols = RangeSize;
    nRows = header->lastBufferRow - header->firstBufferRow + 1;
    /*
        Start of Buffer
    */
    index2D = OVERLAP * RangeSize;
    imageBuffer = (ers1Complex *)&(imageBuffer[index2D]);

    tImageSize = nRows * (RangeSize + KERNELSIZE) * pixelSize;
    tImage = (ers1Complex *)malloc(tImageSize);
    tSave = tImage;
    index2D = 0 * (RangeSize + KERNELSIZE) + KERNELSIZE / 2;
    tImage = (ers1Complex *)&(tImage[index2D]);
    colOffset = RangeSize + KERNELSIZE / 2;
    /*
            Apply interpolation filter kernel in azimuth direction.
    */

    for (j = 0; j < nCols; j++)
    { /* Do Column */
        kI = shifts.aKernel[j];
        for (i = 0; i < nRows; i++)
        {
            resultR = 0.0;
            resultI = 0.0;
            for (k = -KERNELSIZE / 2; k < KERNELSIZE / 2; k++)
            {
                index2D = (i + k) * RangeSize + j;
                /*
                if( j==512 && i == 0) fprintf(stderr,"%i %i %i\n",i,
                    imageBuffer[index2D].r,imageBuffer[index2D].i);
                if( j==512 && i == 511) fprintf(stderr,"%i %i %i\n",i,
                    imageBuffer[index2D].r,imageBuffer[index2D].i);
                */
                resultR += (float)imageBuffer[index2D].r * kernel[kI][k];
                resultI += (float)imageBuffer[index2D].i * kernel[kI][k];
            }
            index2D = i * (RangeSize + KERNELSIZE) + j;
            tImage[index2D].r = (ers1ComplexElement)resultR;
            tImage[index2D].i = (ers1ComplexElement)resultI;
        }
    }
    /*
            Zero pad extra columns
    */
    for (i = 0; i < nRows; i++) /* Do Rows */
        for (j = -KERNELSIZE / 2; j < 0; j++)
        {
            index2D = i * (RangeSize + KERNELSIZE) + j;
            tImage[index2D] = zero;
            tImage[index2D + colOffset] = zero;
        }
    /*
        Apply interpolation filter in the range direction
    */
    for (j = 0; j < nCols; j++)
    { /* Do Column */
        kI = shifts.rKernel[j];
        for (i = 0; i < nRows; i++)
        {
            resultR = 0.0;
            resultI = 0.0;
            for (k = -KERNELSIZE / 2; k < KERNELSIZE / 2; k++)
            {
                index2D = i * (RangeSize + KERNELSIZE) + j + k;
                resultR += (float)tImage[index2D].r * kernel[kI][k];
                resultI += (float)tImage[index2D].i * kernel[kI][k];
            }
            index2D = i * RangeSize + j;
            imageBuffer[index2D].r = (ers1ComplexElement)resultR;
            imageBuffer[index2D].i = (ers1ComplexElement)resultI;
        }
    }

    /*
        Free buffer and return
    */
    free(tSave);
    return;
}
