#include "ers1.h"


/*
   Shift complex Image by integer valued shift.
*/
    void integerShift(ers1ComplexImage *imageA,
                            ers1Complex *rowA, 
                            ers1Complex *interpolatedRow,
                            ers1Shift shifts, int rowIndex)
{
     int col, row;               /* Column and row shift */
     int j;                       /* LCV */
     int index2D;
     ers1Complex zero;
     zero.r = 0;
     zero.i = 0;
/*
    Zero rows that can be interpolated.
*/
    if(rowIndex < shifts.startAzimuth || rowIndex >= shifts.endAzimuth ) {
        for(j=0; j < imageA->header.rs; j++) interpolatedRow[j] = zero;
        return;
    }
/*
    Interpolated rows that can be.
*/
    for(j = 0; j < shifts.startRange; j++ ) interpolatedRow[j] = zero;
    for(j = shifts.startRange; j < shifts.endRange; j++ ) {
            row = shifts.rowShift[j];
            col =  j + shifts.columnShift[j];
            index2D = row*RangeSize+col;
            interpolatedRow[j] = rowA[index2D];
    } 
    for(j = shifts.endRange; j < imageA->header.rs; j++ ) 
         interpolatedRow[j] = zero;
}

