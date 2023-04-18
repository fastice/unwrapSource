#include "ers1.h"
#include <stdlib.h>
/*
   Read shift file and interpolate shifts
*/
void readShifts(char *shiftFile, float *dR, float *dA, float *omega)
{
     FILE *fp;
     float *r;
     float freq;
     float range, deltaR, deltaA;
     int32_t lineCount, ncolumns, eod, count, i;
     char line[LINEMAX];
     /*
          Malloc arrays
     */
     r = (float *)malloc(RangeSize * sizeof(float));
     /*
         If no shiftFile initialize omega to zero
     */
     if (shiftFile == NULL && omega != NULL)
     {
          for (i = 0; i < RangeSize; i++)
               omega[i] = 0.0;
          return;
     }
     /*
        Open shift file and read data.
     */
     fp = openInputFile(shiftFile);
     lineCount = 0;
     count = 0;
     lineCount = goToBeginningOfData(fp, lineCount, &ncolumns);

     /*
        Loop until end of data.
     */
     while (TRUE)
     {
          lineCount = getDataString(fp, lineCount, line, &eod);
          if (eod == TRUE)
               break;
          sscanf(line, "%f%f%f%f", &range, &deltaR, &deltaA, &freq);
          r[count] = range;
          if (dR != NULL)
               dR[count] = deltaR;
          if (dA != NULL)
               dA[count] = deltaA;
          if (omega != NULL)
               omega[count] = freq;
          count++;
          if (count > RangeSize)
               error("readShifts: Invalid number of data points %i/%i",
                     count, RangeSize);
     }
     if (count != RangeSize)
          error("readShifts: Invalid number of data points %i/%i",
                count, RangeSize);
     /*
         Free memory and return
     */
     free(r);
     return;
}
