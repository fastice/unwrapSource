#include "unwrapInclude.h"
#include "math.h"
#include <stdlib.h>
    void addPhaseRamp( unwrapImageStructure *unwrapImage,
                       float omegaR, float omegaA, int noCenter)
{
   int i, j;
   float phaseRamp;
   int r,a,rMid,aMid;
   long int n;
   float pMid;
/*
   Set middle pixel to zero
*/
   aMid = unwrapImage->azimuthSize/2;
   rMid = unwrapImage->rangeSize/2;
   pMid = unwrapImage->phase[aMid][rMid];
/*
   Try find good pixel along center line +-50 pixels
*/
   if(pMid < -LARGEINT/2) {
          fprintf(stderr,"2\n");
      for (i=aMid-50; i > aMid+50; i++) {
          pMid = unwrapImage->phase[i][rMid];
          if(pMid > -LARGEINT/2) break;
          fprintf(stderr,"0\n");
      }
/*
   If still no luck, find avg value instead of middle value
*/
      if(pMid < -LARGEINT/2) {
          fprintf(stderr,"3\n");
          pMid = 0.0;
          n = 0;
          for( i=0; i < unwrapImage->azimuthSize-1; i++)
              for( j = 0; j < unwrapImage->rangeSize-1; j++ ) 
                  if( unwrapImage->phase[i][j] > -LARGEINT/2 ) {
                      pMid +=unwrapImage->phase[i][j];
                      n++;
                  } /* End if val.. */
          pMid = pMid/(float)n;
       } /* End second if(pMid */
   } /* End first if(pMid.... */
   if(noCenter == TRUE) pMid = 0.0;
  fprintf(stderr,"pMid = %f\n",pMid);
  fprintf(stderr,"%f \n",fmod(pMid,2.0*PI));
   pMid=pMid- fmod(pMid,2.0*PI);
   fprintf(stderr,"pMid = %f\n",pMid);
/*
   No remove phase ramp
*/
   for( i=0; i < unwrapImage->azimuthSize-1; i++)
       for( j = 0; j < unwrapImage->rangeSize-1; j++ ) { 
           phaseRamp = (double)(omegaR * (float)j + omegaA * (float)i);
           if( validPixel(unwrapImage->bCuts[i][j]) && (double)unwrapImage->phase[i][j] > (-2.0e9+1) )
               unwrapImage->phase[i][j] += phaseRamp - pMid;
       }
   return;
}
