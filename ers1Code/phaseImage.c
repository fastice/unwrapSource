#include "ers1.h"

/*
     Compute phase image for unwrapping program.
*/

    void phaseImage(unwrapImageStructure *unwrapImage,
                           ers1ComplexImage *complexImage,
                           float omegaR, float omegaA)
{
    ers1Complex *rowA;
    ers1Complex tmp;
    float pr, pi, phaseRamp;
    int i,j;
    int left,right, up, down;
    int minUpi,maxDowni,maxLeftj,minRightj;
    left = 0.25 * complexImage->header.rs;
    right = 0.75 * complexImage->header.rs;
    down = 0.25 * complexImage->header.as;
    up = 0.75 * complexImage->header.as;
    maxLeftj = -1;
    minRightj = complexImage->header.rs;
    maxDowni = -1;
    minUpi = complexImage->header.as;

    for(i = 0; i < complexImage->header.as; i++ ) {
/*
    Get row Pointers.
*/ 
      rowA = (ers1Complex *) rowPtr((ers1GenericImage *) complexImage);      
      for(j=0; j < complexImage->header.rs; j++) {
          unwrapImage->bCuts[i][j] = 0;
          if( (abs(rowA[j].r) < 1) && (abs(rowA[j].i) < 1) ) { 
               unwrapImage->bCuts[i][j] &= INVALIDPIXEL; /* Mark invalid pix */
               unwrapImage->phase[i][j] = 0.0;
               if( j < left && i > down && i < up) 
                   maxLeftj = max(maxLeftj,j);
               else if( j > right  && i > down && i < up)
                   minRightj = min(minRightj,j);
               if( i < down && j > left && j < right)
                    maxDowni = max(maxDowni,i);
               else if( i > up  && j > left && j < right) 
                    minUpi = min(minUpi,i);
          } else { /* Compute phase */
              unwrapImage->bCuts[i][j] |= VALIDPIXEL; /* Mark valid pixel */
              phaseRamp = (double)(omegaR * (float)j + omegaA * (float)i);
              pr = (float)cos(phaseRamp); pi = (float)sin(phaseRamp);
              tmp.r = rowA[j].r * pr - rowA[j].i * pi;
              tmp.i = rowA[j].r * pi + rowA[j].i * pr;
              unwrapImage->phase[i][j] = (float) atan2((double)tmp.r,
                                                      (double)tmp.i);
          }
      } 
      rowIncrement((ers1GenericImage *) complexImage);
    }
fprintf(stderr,"l,r,d,u  %i  %i  %i  %i\n",maxLeftj,minRightj,maxDowni,minUpi);
/*
     Throw out border pixels.
*/
    for(i = 0; i < complexImage->header.as; i++ ) {
        for(j=0; j < maxLeftj+2; j++) {
            unwrapImage->bCuts[i][j] &= INVALIDPIXEL; /* Mark invalid pix */
            unwrapImage->phase[i][j] = 0.0;
        }
        for(j=minRightj-1; j < complexImage->header.rs; j++) {
            unwrapImage->bCuts[i][j] &= INVALIDPIXEL; /* Mark invalid pix */
            unwrapImage->phase[i][j] = 0.0;
        }
    }
    for(j = 0; j < complexImage->header.rs; j++ ) {
        for(i=0; i < maxDowni+2; i++) {
            unwrapImage->bCuts[i][j] &= INVALIDPIXEL; /* Mark invalid pix */
            unwrapImage->phase[i][j] = 0.0;
        }
        for(i=minUpi-1; i < complexImage->header.as; i++) {
            unwrapImage->bCuts[i][j] &= INVALIDPIXEL; /* Mark invalid pix */
            unwrapImage->phase[i][j] = 0.0;
        }
    }

    return;
}  

