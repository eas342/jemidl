#include <stdio.h>
#include "idl_export.h"

int drizzle1d_natural(double *influx,
                      double *inivar,
                      double *inbins,
                      IDL_LONG inbins_n,
                      double *outbins,
                      IDL_LONG outbins_n,
                      double *preflux,
                      double *preivar,
                      double *outflux,
                      double *outivar)
{
    long inbin, outbin;
    for(inbin=0; inbin < inbins_n-1; inbin++) {
        for(outbin=0; outbin < outbins_n-1; outbin++) {
            if (outbins[outbin] >= inbins[inbin+1]) break; //not overlapping, try next inbin
            if (inbins[inbin] >= outbins[outbin+1]) continue; //not overlapping, try next outbin
            //case 1: input bin overhangs the output bin on the left
            double fraction = 1.0;
            if (inbins[inbin] <= outbins[outbin] && inbins[inbin+1] < outbins[outbin+1]) {
                fraction = (inbins[inbin+1]-outbins[outbin])/
                    (inbins[inbin+1]-inbins[inbin]);
            }
            //case 2: input bin overhangs the output bin on the right
            else if (inbins[inbin+1] >= outbins[outbin+1] && inbins[inbin] > outbins[outbin]) {
                fraction = (outbins[outbin+1]-inbins[inbin])/
                    (inbins[inbin+1] - inbins[inbin]);
            }
            //case 3: input bin is totally encompassed by the output bin
            else if (inbins[inbin] > outbins[outbin] && inbins[inbin+1] < outbins[outbin+1]) {
                fraction = 1.0;
            }
            //case 4: the output bin is encompassed by the input bin...
            else if (inbins[inbin] <= outbins[outbin] && inbins[inbin+1] >= outbins[outbin+1]) {
                fraction = (outbins[outbin+1]-outbins[outbin])/
                    (inbins[inbin+1] - inbins[inbin]);
            }
            double newbinivar = fraction * inivar[inbin] + outivar[outbin];
            if (newbinivar == 0.) continue;
            outflux[outbin] = (influx[inbin]*fraction*inivar[inbin]
                               + outflux[outbin]*outivar[outbin]) / newbinivar;
            outivar[outbin] = newbinivar;
        }
    }
    return 1;
}

int drizzle1d(int argc, void* argv[])
{
    if (argc != 10) return 0.0;

    return drizzle1d_natural((double *) argv[0],
                             (double *) argv[1],
                             (double *) argv[2],
                             (IDL_LONG) argv[3],
                             (double *) argv[4],
                             (IDL_LONG) argv[5],
                             (double *) argv[6],
                             (double *) argv[7],
                             (double *) argv[8],
                             (double *) argv[9]);
}
