#include "unwrapInclude.h"
#include <stdlib.h>

    void residueImage( unwrapImageStructure *unwrapImage )
{
   float sumdp, dp;
   float twoPi;
   int nPlus,nMinus;
   int i, j;
   twoPi = 2.0 * PI;
   unwrapImage->nPlus = 0;
   unwrapImage->nMinus = 0;
   for( i=0; i < unwrapImage->azimuthSize-1; i++)
       for( j = 0; j < unwrapImage->rangeSize-1; j++ ) { 
            if( validPixel( unwrapImage->bCuts[i][j] ) &&
                validPixel( unwrapImage->bCuts[i][j+1] ) &&
                validPixel( unwrapImage->bCuts[i+1][j] ) &&
                validPixel( unwrapImage->bCuts[i+1][j+1] )  ){               
		sumdp = 0.0;            
		dp = unwrapImage->phase[i][j+1]   - unwrapImage->phase[i][j]; 
		sumdp =  (dp <= -PI) ? dp + twoPi : (dp > PI) ? dp - twoPi : dp;
		dp = unwrapImage->phase[i+1][j+1] - unwrapImage->phase[i][j+1];
		sumdp += (dp <= -PI) ? dp + twoPi : (dp > PI) ? dp - twoPi : dp;
		dp = unwrapImage->phase[i+1][j] - unwrapImage->phase[i+1][j+1];
		sumdp += (dp <= -PI) ? dp + twoPi : (dp > PI) ? dp - twoPi : dp;
		dp = unwrapImage->phase[i][j]     - unwrapImage->phase[i+1][j];
		sumdp += (dp <= -PI) ? dp + twoPi : (dp > PI) ? dp - twoPi : dp;
        	if(sumdp > 6.) {
	            unwrapImage->bCuts[i][j] |= PRESIDUE;
		    unwrapImage->nPlus++;
		} else if(sumdp < -6.) {
		    unwrapImage->bCuts[i][j] |= NRESIDUE;
		    unwrapImage->nMinus++;
		}
            }
       }
       unwrapImage->nLeft = unwrapImage->nPlus + unwrapImage->nMinus;
     fprintf(stderr,"\nnPlus = %i  nMinus = %i\n",
      unwrapImage->nPlus,unwrapImage->nMinus);
   return;
}
