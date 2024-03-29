#include "ers1.h"
/*
    Compute pixel size for a given image type
*/

int32_t pixSize(int32_t imageType)
{
    int32_t pixelSize;

    switch (imageType)
    {
    case COMPLEX:
        pixelSize = sizeof(ers1Complex);
        break;
    case POWER:
        pixelSize = sizeof(ers1Power);
        break;
    case AMPLITUDE:
        pixelSize = sizeof(ers1Amplitude);
        break;
    default:
        error("\n%s\n", "pixSize: Invalid imageType");
    }
    return pixelSize;
}
