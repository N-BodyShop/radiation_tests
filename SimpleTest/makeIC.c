/* program to create an initial distribution of particles for a
 * simple radiation test.
 * The test has gas particles at the 8 corners of the cube to enforce
 * a common division of the tree between Gasoline and ChaNGa.
 * The rest of the cube is filled with a Poisson distribution of gas
 * In the upper corner is a Poisson distribution of stars.
 *
 * The point of the test is to see how well Multipole expansions
 * represent the true radiation field from the stars.
 */
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <assert.h>

int
main(argc, argv)
     int argc;
     char **argv;
{
  double mass;
  double xmin;
  double xmax;
  double eps;
  struct dump h;
  double deltax;
  int i;
  int seed;
  XDR xdrs;

  if(argc != 4) {
      fprintf(stderr, "Usage: makeIC ngas nstar seed\n");
      return 1;
      }
  
  seed = atoi(argv[3]);

  srand(seed);
  
  h.nsph = atoi(argv[1]);
  assert(h.nsph > 8);
  h.ndark = 0;
  h.nstar = atoi(argv[2]);
  assert(h.nstar >= 0);
  h.nbodies = h.nsph + h.nstar;
  h.ndim = 3;
  h.time = 0.0;

  xmax = 0.5;
  xmin = -0.5;

  mass = 1.0/h.nsph; /* scale so the total mass is 1 */
  eps = pow(h.nsph, -1.0/3.0)*(xmax - xmin)/20.0;
  
  xdrstdio_create(&xdrs, stdout, XDR_ENCODE);

  xdr_header(&xdrs, &h);

  deltax = (xmax - xmin);
  for(i = 0; i < 8; i++) {
	  struct gas_particle gp;

	  gp.mass = mass;
          /* N.B. quick way to convert integar to 8 corners of the cube.
           */
	  gp.pos[0] = xmin + (i % 2) * deltax;
	  gp.pos[1] = xmin + ((i>>1) % 2) * deltax;
	  gp.pos[2] = xmin + ((i>>2) % 2) * deltax;
	  gp.vel[0] = 0.0;
	  gp.vel[1] = 0.0;
	  gp.vel[2] = 0.0;
	  gp.hsmooth = eps;
	  gp.rho = 0.0;
	  gp.temp = 100.0;
	  gp.phi = 0.0;

	  xdr_gas(&xdrs, &gp);
  }
  
  /* Now fill in the rest of the gas */
  for(i = 8; i < h.nsph; i++)
    {
	  struct gas_particle gp;

	  gp.mass = mass;
	  gp.pos[0] = xmin + rand()/((double) RAND_MAX)*deltax;
	  gp.pos[1] = xmin + rand()/((double) RAND_MAX)*deltax;
	  gp.pos[2] = xmin + rand()/((double) RAND_MAX)*deltax;
	  gp.vel[0] = 0.0;
	  gp.vel[1] = 0.0;
	  gp.vel[2] = 0.0;
	  gp.hsmooth = eps;
	  gp.rho = 0.0;
	  gp.temp = 100.0;
	  gp.phi = 0.0;

	  xdr_gas(&xdrs, &gp);
	}
  /* Now for the stars */
  xmax = 0.5;
  xmin = 0.25;
  deltax = (xmax - xmin);

  mass = 1.0/h.nstar; /* scale so the total mass is 1 */
  eps = pow(h.nstar, -1.0/3.0)*(xmax - xmin)/20.0;

  for(i = 0; i < h.nstar; i++)
    {
	  struct star_particle gp;

	  gp.mass = mass;
	  gp.pos[0] = xmin + rand()/((double) RAND_MAX)*deltax;
	  gp.pos[1] = xmin + rand()/((double) RAND_MAX)*deltax;
	  gp.pos[2] = xmin + rand()/((double) RAND_MAX)*deltax;
	  gp.vel[0] = 0.0;
	  gp.vel[1] = 0.0;
	  gp.vel[2] = 0.0;
	  gp.eps = eps;
	  gp.metals = 0.0;
	  gp.tform = 0.0;
	  gp.phi = 0.0;

	  xdr_star(&xdrs, &gp);
	}

  xdr_destroy(&xdrs);
  return(0);
}
