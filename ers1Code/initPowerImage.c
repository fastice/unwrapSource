#include "ers1.h"
#include <stdlib.h>

/*
    Create and init full size (standard size) power image.
*/
ers1PowerImage *initFullSizePowerImage()
{
    return initPowerImage(RangeSize, AzimuthSize, 0, BufferSize);
}

/*
    Create and init power image.
*/
ers1PowerImage *initPowerImage(int32_t rs, int32_t as, int32_t currentRow,
                               int32_t bufferSize)
{
    ers1PowerImage *image;
    int32_t size, imageSize;
    /*
        Malloc space for image structure
    */
    image = (ers1PowerImage *)malloc(sizeof(ers1PowerImage));
    imageSize = rs * as * sizeof(ers1Power);
    if (bufferSize > imageSize)
        bufferSize = imageSize;

    image->header.rs = rs; /* Init variables */
    image->header.as = as;
    image->header.currentRow = currentRow;
    image->header.lastBufferRow = -1;
    image->header.firstBufferRow = -1;
    image->header.bufferSize = bufferSize;
    image->header.imageType = POWER;
    /*
       Malloc image buffer.
    */
    size = bufferSize + OVERLAP * 2 * rs * sizeof(ers1Power);
    image->image = (ers1Power *)malloc(size);
    return image;
}
