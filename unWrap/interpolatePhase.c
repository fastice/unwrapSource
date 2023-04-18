#include "unwrapInclude.h"
#include <math.h>
#include <stdlib.h>
/*
    Do simple interpolation on phase unwrapping output.
*/
void interpolatePhase(unwrapImageStructure *unwrapImage)
{
    int32_t i, j;
    float dpa, dpr;

    for (i = 0; i < unwrapImage->azimuthSize - 1; i++)
        for (j = 0; j < unwrapImage->rangeSize - 1; j++)
        {
            if (unwrapImage->bCuts[i][j] & BRANCHCUT)
            {
                dpa = 1000.0;
                dpr = 1000.0;
                if (i > 0 && i < unwrapImage->azimuthSize - 1)
                    if (!(unwrapImage->bCuts[i - 1][j] & BRANCHCUT) &&
                        !(unwrapImage->bCuts[i + 1][j] & BRANCHCUT))
                        dpa = (float)fabs((double)(unwrapImage->phase[i + 1][j] -
                                                   unwrapImage->phase[i - 1][j]));
                if (j > 0 && j < unwrapImage->rangeSize - 1)
                    if (!(unwrapImage->bCuts[i][j - 1] & BRANCHCUT) &&
                        !(unwrapImage->bCuts[i][j + 1] & BRANCHCUT))
                        dpr = (float)fabs((double)(unwrapImage->phase[i][j + 1] -
                                                   unwrapImage->phase[i][j - 1]));
                if (dpa < PI)
                {
                    unwrapImage->phase[i][j] = 0.5 * (unwrapImage->phase[i + 1][j] +
                                                      unwrapImage->phase[i - 1][j]);
                }
                else if (dpr < PI)
                {
                    unwrapImage->phase[i][j] = 0.5 * (unwrapImage->phase[i][j + 1] +
                                                      unwrapImage->phase[i][j - 1]);
                }
            }
        }
    return;
}
