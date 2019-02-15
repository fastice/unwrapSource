#include "ers1.h"

/*
   Write current buffer to output file.
*/

    void writeBuffer(ers1GenericImage *image) 
{
    headerInfo *header;  /* Image header */
    char *imageBuffer; /* Image Buffer */ 
    FILE *fp;
    int imageType;
    int bstype;
    int nBytesPerRow;
    int pixelSize;
    int nRows;        /* Number of rows to read from file */

    header = &(image->header);
    imageType = header->imageType;
    imageBuffer = (char *) image->image; 
    pixelSize = pixSize(imageType);  
    nBytesPerRow = pixelSize * header->rs;
    
    bstype=FLOAT32FLAG;
    if(imageType==COMPLEX) bstype=INT16FLAG; 
    else if(imageType==AMPLITUDE) bstype=BYTEFLAG;

    imageBuffer += nBytesPerRow * OVERLAP;  /* Start of data */
   
    nRows = header->lastBufferRow - header->firstBufferRow + 1;
    fp = header->outfp;
    fwriteBS(imageBuffer, nBytesPerRow, nRows,fp,bstype);
    fprintf(stderr,"Finished write Buffer, nRows = %i\n",nRows);
    return;
}
