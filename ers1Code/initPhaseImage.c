#include "ers1.h"
#include <stdlib.h>

/*
    Create and init phase image.
*/
    ers1PhaseImage *initPhaseImage(int rs, int as, int currentRow,
    int bufferSize )
{
    ers1PhaseImage *image;
    int size,imageSize;
/*
    Malloc space for image structure 
*/
   image = (ers1PhaseImage *) malloc(sizeof(ers1ComplexImage));
   imageSize = rs *as * sizeof(ers1Phase);
   if(bufferSize > imageSize) bufferSize = imageSize;
 
   image->header.rs = rs ; /* Init variables */
   image->header.as = as;
   image->header.currentRow = currentRow;
   image->header.lastBufferRow = -1;
   image->header.firstBufferRow = -1;
   image->header.bufferSize = bufferSize;
   image->header.imageType = POWER;
/*
   Malloc image buffer.
*/
   size = bufferSize + OVERLAP * 2 * rs * sizeof(ers1Phase);
   image->image = (ers1Phase *) malloc(size);
   return image;
}
