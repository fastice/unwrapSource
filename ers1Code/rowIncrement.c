#include "ers1.h"
/*
    increment row ptr
*/

void rowIncrement(ers1GenericImage *image)
{
    headerInfo *header;
    int32_t pixelSize;
    int32_t nBufferRows, nRowsLeft;
    header = &(image->header);
    header->currentRow++; /* Increment row pointer */

    if (header->currentRow > header->lastBufferRow &&
        header->imageMode == WRITEIMAGE)
    {
        pixelSize = pixSize(header->imageType);
        writeBuffer(image);
        /*
                Update buffer if not endof image
        */
        if (header->currentRow >= header->as)
            return;
        nBufferRows = header->bufferSize / (header->rs * pixelSize);
        nRowsLeft = header->as - header->currentRow;
        header->firstBufferRow = header->currentRow;
        header->lastBufferRow = header->currentRow +
                                min(nRowsLeft, nBufferRows) - 1;
    }
}
