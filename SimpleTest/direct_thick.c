/*
 * Program to do the direct optically thick radiation calculation.
 * Usage: direct_thick dRadSoft < simple.tbin > direct.radFlux
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

/*
 * Return the integral of the kernel along a line passing through a
 * sphere a distance b from the center.
 * This is for the M4 kernel 
 */
double kernelIntegral(double h, double b){
    double norm = M_PI*h*h;
    b = b/h;
    double b2 = b*b;
    if (b2 > 4) return 0.0;
    if(b2 < 1e-14) return 1.5/norm;

    double c = 21./8.*b2*sqrt(4-b2) - 3*sqrt(4-b2) + pow(4-b2,1.5) + 3./16.*b2*b2*log(-sqrt(4-b2)+2) + 3*b2*log(-sqrt(4-b2)+2);
    assert(b > 0);
    if(b < 1){
        c = c - 21./4.*b2*sqrt(1-b2) + 3./2.*sqrt(1-b2) - 2*pow(1-b2,1.5) - 3./4.*b2*b2*log(-sqrt(1-b2)+1) - 3*b2*log(-sqrt(1-b2)+1) + 9./32.*b2*b2*log(b2);
        }
    else if(b < 2){
        c = fabs(c - 3./2.*b2*log(b2) - 3./32.*b2*b2*log(b2));
        }
    else
        c = 0;
    return c/norm;
    }

double dist2(const struct gas_particle *p1, const struct gas_particle *p2)
{
    double dx = p1->pos[0] - p2->pos[0];
    double dy = p1->pos[1] - p2->pos[1];
    double dz = p1->pos[2] - p2->pos[2];
    return dx*dx + dy*dy + dz*dz;
    }

struct gas_particle *target_gp = NULL;

/*
 * Sort function for stupidly determining kth nearest neighbor
 */
int compare_neighbors(const void *p, const void *q)
{
    const struct gas_particle *p1 = (const struct gas_particle *)p;
    const struct gas_particle *p2 = (const struct gas_particle *)q;
    double d1 = dist2(target_gp, p1);
    double d2 = dist2(target_gp, p2);
    if(d1 < d2)
        return -1;
    else if (d1 > d2)
        return 1;
    return 0;
    }

int
main(argc, argv)
     int argc;
     char **argv;
{
  XDR xdrs;
  struct dump h;
  int i, j, k;
  int nSmooth = 32;
  double dKappa = 1.0;

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
  
/* Slow but sure way of finding nSmooth nearest neighbors
 */
  struct gas_particle *pGasNeib = malloc(nSph*sizeof(*pGas));
  for(i = 0; i < nSph; i++) {
      pGasNeib[i] = pGas[i];
  }
  for(i = 0; i < nSph; i++) {
      target_gp = &pGas[i];
      qsort(pGasNeib, nSph, sizeof(*pGasNeib), &compare_neighbors);
      pGas[i].hsmooth = sqrt(dist2(&pGas[i], &pGasNeib[nSmooth-1]));
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
          double tmpFlux = b*M_1_PI*0.25;
          double rmag = sqrt(r2);
          double tau = 0.0;
          for(k = 0; k < nSph; k++) {
              if(i == k) continue; /* do not self absorb */
              double dxg = pGas[k].pos[0] - pStar[j].pos[0];
              double dyg = pGas[k].pos[1] - pStar[j].pos[1];
              double dzg = pGas[k].pos[2] - pStar[j].pos[2];
              /* distance projected long radius */
              double bpar = (dxg*dx + dyg*dy + dzg*dz)/rmag;
              /* either behind the star or behind the target */
              if(bpar < 0.0 || bpar > rmag) continue;
              /* perpendicular distance from Pythagoras */
              double bperp = dxg*dxg + dyg*dyg + dzg*dzg - bpar*bpar;
              assert(bperp >= 0.0);
              bperp = sqrt(bperp);
              double tau_i = dKappa*pGas[k].mass
                  * kernelIntegral(pGas[k].hsmooth, bperp);
              assert(tau_i >= 0.0);
              tau += tau_i;
          }
          dFlux += tmpFlux*exp(-tau);
      }
      printf("%g\n", dFlux);
  }
}
