#include "ers1.h"
#include <stdlib.h>
/*
    Create and init full size (standard size) complex image.
*/
ers1ComplexImage *initFullSizeComplexImage()
{
    return initComplexImage(RangeSize, AzimuthSize, 0, BufferSize);
}

/*
    Create and init complex image.
*/
ers1ComplexImage *initComplexImage(int32_t rs, int32_t as, int32_t currentRow,
                                   int32_t bufferSize)
{
    ers1ComplexImage *image;
    int32_t size;
    int32_t imageSize;
    /*
        Malloc space for image structure
    */
    image = (ers1ComplexImage *)malloc(sizeof(ers1ComplexImage));
    imageSize = rs * as * sizeof(ers1Complex);
    if (bufferSize > imageSize)
        bufferSize = imageSize;
    image->header.rs = rs; /* Init variables */
    image->header.as = as;
    image->header.currentRow = currentRow;
    image->header.lastBufferRow = -1;
    image->header.firstBufferRow = -1;
    image->header.bufferSize = bufferSize;
    image->header.imageType = COMPLEX;
    /*
       Malloc image buffer.
    */
    size = bufferSize + OVERLAP * 2 * rs * sizeof(ers1Complex);
    fprintf(stderr, "bs,rs,as %i %i %i\n", size, rs, as);
    image->image = (ers1Complex *)malloc(size);

    return image;
}
