#include "ers1.h"

/*
    Open file for image
*/

void openImage(ers1GenericImage *image, char *filename, int32_t rwFlag)
{
    FILE *fp;
    int32_t nBufferRows; /* Number of rows in a buffer */
    int32_t pixelSize;

    if (rwFlag == READIMAGE)
    { /* Open image for read */
        if (filename == NULL)
            fp = stdin;
        else
            fp = fopen(filename, "r");
        if (fp == NULL)
            error("*** openImage: Error opening %s ***\n", filename);
        image->header.infp = fp;
        image->header.imageMode = READIMAGE;
    }
    else if (rwFlag == WRITEIMAGE)
    {
        if (filename == NULL)
            fp = stdout;
        else
            fp = fopen(filename, "w");
        if (fp == NULL)
            error("*** openImage: Error opening %s ***\n", filename);
        pixelSize = pixSize(image->header.imageType);
        image->header.outfp = fp;
        image->header.firstBufferRow = 0; /* Init buffer for write operations */
        nBufferRows = image->header.bufferSize / (image->header.rs * pixelSize);
        image->header.lastBufferRow = min(nBufferRows, image->header.as) - 1;
        image->header.imageMode = WRITEIMAGE;
    }
    else
        error("*** openImage: Invalid mode flag ***\n");
}
