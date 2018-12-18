#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"
#if dimension == 3
# include "lambda2.h"
#endif
#include "view.h"
#include "maxruntime.h"
#define RHOR 1000.
#define MUR 100.
#if BUBBLE19
// Bubble 19 of Cano-Lozano et al, P.R.Fluids, 2016
# define Ga 100.8
# define Bo 4.
# define MAXTIME 82
#else
// Bubble 26 of Cano-Lozano et al, P.R.Fluids, 2016
# define Ga 100.25
# define Bo 10.
# define MAXTIME 110
#endif
#define WIDTH 120.0
#define Z0 3.5
int LEVEL = 12;
int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi (argv[1]);
size (WIDTH);
  origin (-L0/2, 0, -L0/2);
  init_grid (128);
 rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1./Ga;
  mu2 = 1./(MUR*Ga);
  f.σ = 1./Bo;
 TOLERANCE = 1e-4;
  run();
}
event init (t = 0) {
  if (!restore (file = "restart")) {
    refine (sq(x) + sq(y - Z0) + sq(z) - sq(0.75) < 0 && level < LEVEL);
    fraction (f, sq(x) + sq(y - Z0) + sq(z) - sq(.5));
  }
}
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1.;
}
event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);
}
event logfile (i += 10) {
  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb)
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
	  reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    yb += y*dv;
    zb += z*dv;
    vbx += u.x[]*dv;
    vby += u.y[]*dv;
    vbz += u.z[]*dv;
    sb += dv;
  }
  fprintf (ferr,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb,
	   xb/sb, yb/sb, zb/sb,
	   vbx/sb, vby/sb, vbz/sb);
  fflush (ferr);
}
event snapshot (t = 1; t <= MAXTIME; t++)
{
  scalar l2[], omegay[];
#if dimension == 3
  lambda2 (u, l2);
  foreach()
    omegay[] = (u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Δ);
  boundary ({omegay});
#endif
 char name[80];
  sprintf (name, "dump-%03d", (int) t);
  dump (file = name);
}
event movie (t = 21; t <= MAXTIME; t += 0.25)
{

#if BUBBLE19
  
  view (fov = 5.0278, quat = {-0.132839,0.513023,0.0748175,0.844727},
	tx = 0.00149469, ty = -0.355489, width = 300, height = 800);
 travelling (50, 60, fov = 2.07254, tx = 0.00524944, ty = -0.513744);
 travelling (60, 82, tx = 0.00898601, ty = -0.703841);

#else  // BUBBLE26
 view (fov = 5.0278, quat = {-0.0487657,0.654023,0.0506598,0.7532},
	tx = -0.00353642, ty = -0.302285, width = 300, height = 800);

  travelling (50, 60, fov = 2.07254, tx = -0.0022909, ty = -0.402237);

  travelling (60, 110, tx = -0.00644264, ty = -0.741016);
#endif
  clear();
  draw_vof ("f", fc = {0.13,0.47,0.77});
#if dimension == 3
  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -0.0002);
#endif

  save ("bubble.mp4");
}
local% qcc -source -grid=octree -D_MPI=1 bubble.c
local% scp _bubble.c occigen.cines.fr:
