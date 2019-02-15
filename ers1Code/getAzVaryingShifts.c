#include "ers1.h"
#include <math.h>
#include <stdlib.h>
/*
    Get array of shifts from shiftfile to allow shifts to be updated along
    track. File format ro,ao,dr,da,nr,na, rshifts,ashifts
*/
    void getAzVaryingShifts(ers1Shift *shifts, char *shiftFile) 
{
    FILE *fp;
    float tmp[6];
    float **rShifts,**aShifts;
    int i;

    fp = fopen(shiftFile,"r");
/*
   Read header info
*/
    freadBS(tmp,sizeof(float)*6,1,fp,FLOAT32FLAG);
    shifts->gRo = (int)tmp[0];
    shifts->gAo = (int)tmp[1];
    shifts->dRg = (int)tmp[2];
    shifts->dAg = (int)tmp[3];
    shifts->gNr = (int)tmp[4];
    shifts->gNa = (int)tmp[5]; 
/*
    Malloc memory for shift matrices    
*/
    rShifts = (float **) malloc( shifts->gNa * sizeof(float) );
    aShifts = (float **) malloc( shifts->gNa * sizeof(float) ); 
    for(i=0; i < shifts->gNa; i++) {
        rShifts[i] = (float *) malloc( shifts->gNr * sizeof(float) );
        aShifts[i] = (float *) malloc( shifts->gNr * sizeof(float) ); 
    }
    shifts->rangeShifts = rShifts;
    shifts->azimuthShifts = aShifts;
/*
   Input data range then azimuth shifts
*/
   for(i=0; i < shifts->gNa; i++)
      freadBS(rShifts[i],sizeof(float),shifts->gNr,fp,FLOAT32FLAG);
   for(i=0; i < shifts->gNa; i++)
      freadBS(aShifts[i],sizeof(float),shifts->gNr,fp,FLOAT32FLAG);
}
