#include "ers1.h"

/*
    Fill current buffer starting at current row-overlap.
*/

    void fillBuffer(ers1GenericImage *image) 
{
    headerInfo *header;  /* Image header */
    char *imageBuffer, *tmpBuffer; /* Image Buffer */
    FILE *fp;
    int firstRow;     /* First row to read in file, from current - overlap */
    int nRows;        /* Number of rows to read from file */
    int nBufferRows;  /* Number of rows not including overlap regions */
    int nOver;        /* If end of image, number of overlap rows to 0 pad */
    int beginSecondOverlap;  /* Byte offset begining of second overlap region */
    int nRowsLeft;    /* Number of rows from current to end of image */
    int pixelSize;    /* Pixel size for image type */
    int imageType;    /* Image type lable */
    int nBytes;       /* Number of bytes to operate on */
    int nPixels;      /* Number of pixels to read into buffer */
    int n;
    int i;            /* LCV */
    int bstype;
    int filePosition;
    int nBytesPerRow;
    header = &(image->header);
    fp = header->infp;
    imageType = header->imageType;
    imageBuffer = (char *) image->image;   
    pixelSize = pixSize(imageType);
    bstype=FLOAT32FLAG;
    if(imageType==COMPLEX) bstype=INT16FLAG; else if(imageType==AMPLITUDE) bstype=BYTEFLAG;
/*
    Number of rows that can fit in buffer.
*/
    nBufferRows = header->bufferSize/(header->rs * pixelSize);
/*
    Number of row that have not been read yet.
*/
    nRowsLeft = header->as - header->currentRow;
/*
    IF input stdin read buffer with no overlap region
*/

    if(fp == stdin) {
          nBytesPerRow = pixelSize * header->rs;
          imageBuffer += nBytesPerRow * OVERLAP;  /* Start of data */  
          nRows = min(nBufferRows,nRowsLeft);
          header->lastBufferRow = header->currentRow + 
                                  min(nBufferRows,nRowsLeft) -1;
          header->firstBufferRow = header->currentRow;
          nPixels = nRows * header->rs;
          n = (int) freadBS(imageBuffer,pixelSize,nPixels,fp,bstype);
          fprintf(stderr,"Completing stdin read %i %i\n",
               header->firstBufferRow,header->lastBufferRow);
          return;
    }
/*
    Initial Row to read (including overlap region).
*/
    firstRow = header->currentRow - OVERLAP;
/*
    If near beginning zero fill overlap region.
*/    
    if (firstRow < 0) {
        nBytes = pixelSize * header->rs * (-firstRow);
        for(i=0; i < nBytes; i++) {       /* Zero fill initial overlap */
                *imageBuffer = 0;
                imageBuffer++;
        }
fprintf(stderr,"Begin overlap \n");
        nRows = firstRow + OVERLAP;
        firstRow = 0;
    } else nRows = OVERLAP;
/*
    Compute nRows
*/
    nRows += min(nBufferRows,nRowsLeft);    
/*
    Take care of second overlap region
*/
    nRows += OVERLAP;                        /* Add rows for overlap region */
    if( (firstRow + nRows) > header->as) {   /* If overlap past endof img */
        nOver = (firstRow + nRows) - header->as; /* Then zero fill overlap */
        nRows = header->as - firstRow;          /* And set nRows to end of img*/
        nBytes = pixelSize * header->rs * nOver;
        beginSecondOverlap = 
            header->bufferSize + pixelSize * header->rs * OVERLAP;
fprintf(stderr,"End overlap \n");
        if(firstRow ==0) tmpBuffer = (char *) image->image;
        else tmpBuffer = imageBuffer;
        for(i = beginSecondOverlap; i < beginSecondOverlap+nBytes; i++ ) 
            *(tmpBuffer+i) = 0;
    }
/* 
    Position file pointer 
*/
    filePosition = firstRow * header->rs * pixelSize; 
    fp = header->infp;
    fseek(fp, filePosition ,SEEK_SET);                  
/*
    fill Buffer
*/  
    nPixels = nRows * header->rs;
    freadBS(imageBuffer,pixelSize,nPixels,fp,bstype);
    header->lastBufferRow = header->currentRow + min(nBufferRows,nRowsLeft) -1;
    header->firstBufferRow = header->currentRow;

    fprintf(stderr,"Completing file read %i %i %i %i %i\n",
        header->firstBufferRow,header->lastBufferRow,filePosition,firstRow,nPixels);
    return;
}
