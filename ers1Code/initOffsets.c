#include "ers1.h"
#include <stdlib.h>
/*
    Compute offsets as function of range using either deltaR and deltaA
    or range dependent offsets in offsetFile.
*/
    void initOffsets(ers1Offset *offsets, float *omega, char *offsetFile,
                            float deltaR,float deltaA,
                            ers1ComplexImage *imageA)
{
    float rFract, aFract;    /* Float parts of offset */
    int rInt, aInt;          /* Integer parts of offsets */
    float minAzimuth, maxAzimuth;
    float *dR,*dA;
    float firstDeltaR;
    float t,u;
    int i;                   /* LCV */
/*
    Allocate arrays 
*/
    dR = (float *)malloc(RangeSize * sizeof(float));
    dA = (float *)malloc(RangeSize * sizeof(float));

    offsets->columnLeft =  (int *) malloc(RangeSize*sizeof(int));         
    offsets->columnRight = (int *) malloc(RangeSize*sizeof(int));
    offsets->rowUp =       (int *) malloc(RangeSize*sizeof(int));
    offsets->rowDown =     (int *) malloc(RangeSize*sizeof(int));
    offsets->a1 = (float *) malloc(RangeSize*sizeof(float));         
    offsets->a2 = (float *) malloc(RangeSize*sizeof(float));
    offsets->a3 = (float *) malloc(RangeSize*sizeof(float));
    offsets->a4 = (float *) malloc(RangeSize*sizeof(float));
    offsets->deltaR = (float *) malloc(RangeSize*sizeof(float));
    offsets->deltaA = (float *) malloc(RangeSize*sizeof(float));

/*
    Init offsets for each range.
*/ 
    minAzimuth = 9999;
    maxAzimuth = -9999;
    offsets->endRange =   -1;
    offsets->startRange = 0;
    if(offsetFile != NULL) readShifts(offsetFile,dR,dA,omega); 
    for(i = 0; i < RangeSize; i++ ) {              
        if(offsetFile != NULL) {  /* Range independent offsets */
            deltaR = dR[i];
            deltaA = dA[i];    
        } else omega[i] = 0.0;
        rInt = (int) deltaR;
        rFract = deltaR - (float)rInt;
        aInt = (int) deltaA;
        aFract = deltaA - (float)aInt;
        if(i == 0) firstDeltaR = deltaR;
/*
   Get min and max azimuth shifts.
*/
         minAzimuth = min(minAzimuth,deltaA);
         maxAzimuth = max(maxAzimuth,deltaA);
/*
    Compute interpolation constants
*/
         if( deltaR > 0 ) t = 1.0 - rFract;  /* Interpolation constants */
         else t = -rFract;
         if( deltaA > 0 ) u = 1.0 - aFract; 
         else u = -aFract;

         offsets->a1[i] = (1.0 - t) * (1.0 - u);
         offsets->a2[i] = t * (1.0 - u);
         offsets->a3[i] = t * u;
         offsets->a4[i] = (1.0 - t) * u;
         offsets->deltaR[i] = deltaR;
         offsets->deltaA[i] = deltaA;
/*
   Precompute whether row column coordinates for interpolation
*/
         if( deltaR > 0 ) {
             offsets->columnLeft[i] = -rInt -1;     /* Column, range */
             offsets->columnRight[i] = -rInt;              
         } else {
             offsets->columnLeft[i] = -rInt;
             offsets->columnRight[i] = -rInt + 1;
         }
         if( deltaA > 0 ) {
             offsets->rowDown[i] = -aInt - 1;    /* Row, azimuth */
             offsets->rowUp[i] = -aInt;              
         } else {
             offsets->rowDown[i] = -aInt;
             offsets->rowUp[i] = -aInt + 1;
         }

    }
/*
    Compute start and Stop range 
*/
     if(firstDeltaR > 0) offsets->startRange=(int)firstDeltaR + 1;
     else offsets->startRange = 0;
     if(firstDeltaR < 0) 
         offsets->endRange = (int)(imageA->header.rs + firstDeltaR);
     else offsets->endRange = imageA->header.rs; 

/*
   Compute start and end azimuth
*/
     if(maxAzimuth > 0) offsets->startAzimuth = (int) maxAzimuth + 1;
     else offsets->startAzimuth = 0;
           
     if(minAzimuth < 0) offsets->endAzimuth = 
          (int)(imageA->header.as + minAzimuth);
     else offsets->endAzimuth = imageA->header.as;
fprintf(stderr,"rowUp, rowDown = %i %i\n",
offsets->rowUp[0], offsets->rowDown[0]);
fprintf(stderr,"columnLeft,columRight = %i %i\n",
offsets->columnLeft[0],offsets->columnRight[0]);
fprintf(stderr,"rFract,aFract %f %f \n",rFract,aFract);
fprintf(stderr,"start,endRange = %i  %i\n",offsets->startRange,offsets->endRange);
fprintf(stderr,"start,endAzimuth = %i  %i\n\n",
offsets->startAzimuth,offsets->endAzimuth);
fprintf(stderr,"t u %f  %f\n",t,u);
fprintf(stderr,"a1,a2,a3,a4 \n %f \n %f \n %f \n %f \n",
offsets->a1[0],offsets->a2[0],offsets->a3[0],offsets->a4[0]);
/*
    Free memory and return
*/
    free(dR);
    free(dA);
    return;
}

