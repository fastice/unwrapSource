#include "ers1.h"

/*
    Compute ptr to current row
*/

    void *rowPtr( ers1GenericImage *image ) 
{
    headerInfo *header;
    char *imageBuffer;
    int pixelSize;
    int imageType;
    int rowInBuffer;    /* Row in buffer of current row */
    int nBytesPerRow;     /* Number of bytes per row */
  

    header = &(image->header);
    imageType = header->imageType;
    if( header->currentRow >= header->as || header->currentRow < 0 ) 
        error("*** rowPtr: invalid current row  %i***\n",header->currentRow);
    pixelSize = pixSize(imageType);
  
    if(header->currentRow > header->lastBufferRow &&
        header->imageMode == READIMAGE)  fillBuffer(image);
/*
    Compute buffer address for current row.
*/
    nBytesPerRow = pixelSize * header->rs;
    rowInBuffer = header->currentRow - header->firstBufferRow;
    if( rowInBuffer < 0 ) 
        error("*** rowPtr: invalid current row or start of Buffer ***\n");
    imageBuffer = (char *) image->image;
    imageBuffer += (rowInBuffer + OVERLAP) * nBytesPerRow;
    return imageBuffer;
}
