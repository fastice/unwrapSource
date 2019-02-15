#include "ers1.h" 
#include "math.h"
/*
    Compute radial distance from earth center to elipsoid surface

    lat should be in radians
*/
    double earthRadius(double lat, double rp, double re)
{    
    double earthRad;
    double n,x,z;
/*    fprintf(stderr,"rp,re,lat %f %f %f",rp,re,lat); */


    n = (re*re)/sqrt( pow(re*cos(lat),2.0) +  pow(rp*sin(lat),2.0) );

    x = n*cos(lat);
    z = pow((rp/re),2.0)* n * sin(lat);
    earthRad = sqrt(x*x + z*z);
    
    return earthRad;
}




/*
    Wrong, fixed 02/07/96 
    earthRad = re*re*rp*rp /(rp*rp*cos(lat)*cos(lat) + re*re*sin(lat)*sin(lat));
    earthRad = sqrt(earthRad);
*/
