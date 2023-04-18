#include "ers1.h"
#include "math.h"
#include <stdlib.h>
/*
    Free complex image.
    Modified 07/17/95 to use earthRadius.
*/
double *computePhase(double omegaR, double bn, double bp, double rc,
                     double h, double lat, double *phase)
{
    double theta, thetaD, ptod;
    double r0, delta, bsq;
    double c1, c2, Re, thetaC;
    int32_t i, j;
    /*
        If first time, malloc mem
    */
    j = 0;
    if (phase == NULL)
    {
        phase = (double *)malloc(RangeSize * sizeof(double));
        j = -1;
    }
    /*
        Use linear range ramp
    */
    if (bn < -9998.)
    {
        for (i = 0; i < RangeSize; i++)
            phase[i] = (double)(omegaR * (float)i);
    }
    else
    {
        /*
           Else use curved earth
        */
        bsq = bn * bn + bp * bp;
        Re = earthRadius((lat * DTOR), EMINOR, EMAJOR) * 1000.0;
        /*
                  Re = sqrt(pow(EMAJOR*cos(lat*DTOR),2.0) +
                            pow(EMINOR*sin(lat*DTOR),2.0) ) * 1000.0;
        */
        c1 = 2.0 * h * Re + pow(h, 2.0);
        c2 = 2.0 * (Re + h);
        thetaC = acos((pow(rc, 2.0) + c1) / (c2 * rc));

        if (j < 0)
            fprintf(stderr, "ThetaC, Re,h,rc, c1,c2 = %f %f  %f  %f  %f %f\n", thetaC * RTOD, Re, h, rc, c1, c2);

        ptod = -(4.0 * PI / LAMBDAERS1);

        for (i = -RangeSize / 2; i < RangeSize / 2; i++)
        { /* Compute Phase */
            r0 = rc + RangePixelSize * (double)i;
            theta = acos((pow(r0, 2.0) + c1) / (c2 * r0));
            thetaD = (theta - thetaC);

            /*
               Add 0.5 factor on last term 10/14/94
            */
            delta = -bn * sin(thetaD) - bp * cos(thetaD) + (0.5 * bsq / r0);

            /*if(j < 0) fprintf(stdout,"%f %f %f %f %f %f  %f  %f\n",r0,theta,thetaD,delta,ptod*delta,bn,bp);
             */
            phase[i + RangeSize / 2] = ptod * delta;
            /*
            if(j < 0) fprintf(stderr,"%i td,%f r0,%f bn,bp %f %f d,p %f %f   \n",i,thetaD,r0,bn,bp,delta,ptod*delta);
            */
        }
    }
    return phase;
}
