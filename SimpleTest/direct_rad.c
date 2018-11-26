/*
 * Program to do the direct radiation calculation.
 * Usage: direct_rad dRadSoft < simple.tbin > direct.radFlux
 * where "dRadSoft" is the softening of the radiation field around the
 * stars.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <assert.h>

inline
void SPLINE_RAD(double r2, double twoh2, double *b)
{
    double u, h2inv;
    h2inv = (4.0/(twoh2));
    if ((r2) < (twoh2)) {
        u = sqrt((r2)*h2inv);
        if (u < 1.0) {
            *b = u*h2inv*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
        }
        else {
            *b = -1.0/15.0/r2 + u*h2inv*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u
                                        - 1.0/6.0*u*u*u);
        }
    }
    else {
        *b = 1.0/r2;
    }
}

int
main(argc, argv)
     int argc;
     char **argv;
{
  XDR xdrs;
  struct dump h;
  int i, j;

  assert(argc == 2);
  double dRadSoft = atof(argv[1]);

  xdrstdio_create(&xdrs, stdin, XDR_DECODE);

  xdr_header(&xdrs, &h);

  int nSph = h.nsph;
  struct gas_particle *pGas;
  pGas = malloc(nSph*sizeof(*pGas));
  
  for(i = 0; i < nSph; i++) {
      xdr_gas(&xdrs, &pGas[i]);
  }

  int nStar = h.nstar;
  struct star_particle *pStar;
  pStar = malloc(nStar*sizeof(*pStar));
  
  for(i = 0; i < nStar; i++) {
      xdr_star(&xdrs, &pStar[i]);
  }
  
  printf("%d\n", nSph);
  double twoh2 = 4*dRadSoft*dRadSoft;
  for(i = 0; i < nSph; i++) {
      double dFlux = 0.0;
      for(j = 0; j < nStar; j++) {
  
          double dx = pGas[i].pos[0] - pStar[j].pos[0];
          double dy = pGas[i].pos[1] - pStar[j].pos[1];
          double dz = pGas[i].pos[2] - pStar[j].pos[2];
          
          double r2 = dx*dx + dy*dy + dz*dz;
          double b;
          SPLINE_RAD(r2, twoh2, &b);
          dFlux += b*M_1_PI*0.25;
          }
      printf("%g\n", dFlux);
  }
}

              
