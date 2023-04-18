#include "ers1.h"

/*
    Free complex image.
*/
void freeImage(ers1GenericImage *image)
{
   /*
      Free image buffer.
   */
   free(image->image);
   /*
      Free space for image structure
   */
   free(image);
   return;
}
