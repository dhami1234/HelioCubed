#pragma once

#include "Proto_Accel.H"
#include <cmath>

namespace MHD_Riemann_Solvers
{
// A maps contravariant vectors to Cartesian components. The characteristic
// speed in coordinate direction d uses row d of A^{-1}, not column d of A.
// Use one pointwise geometry for the normal and its scale; ratios of averaged
// metrics introduce angle dependence even on a spherical radial face.
inline ACCEL_DECORATION double mappedFastSpeed(
    const double w[8], const double A[3][3], int d,
    double gamma, double densityFloor)
{
  const int t1 = (d+1)%3, t2 = (d+2)%3;
  double normal[3], bCart[3] = {}, normal2 = 0.0, det = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    const int j = (i+1)%3, k = (i+2)%3;
    normal[i] = A[j][t1]*A[k][t2] - A[k][t1]*A[j][t2];
    normal2 += normal[i]*normal[i];
    det += A[i][d]*normal[i];
    for (int c = 0; c < 3; ++c) bCart[i] += A[i][c]*w[5+c];
  }
  const double normalMag = sqrt(normal2);
  double b2 = 0.0, bn = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    b2 += bCart[i]*bCart[i];
    bn += bCart[i]*normal[i]/normalMag;
  }
  constexpr double inv4Pi = 0.079577471545947667884441881686257181;
  const double rho = fmax(w[0],densityFloor);
  const double sound2 = gamma*fmax(w[4],0.0)/rho;
  const double magnetic2 = inv4Pi*b2/rho;
  double tangent2 = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    const double tangent = bCart[i]-bn*normal[i]/normalMag;
    tangent2 += tangent*tangent;
  }
  tangent2 *= inv4Pi/rho;
  // This form avoids subtracting nearly equal numbers for parallel B when
  // sound and Alfven speeds coincide (the alternative loses sqrt(epsilon)).
  const double discriminant = hypot(sound2-magnetic2,2.0*sqrt(sound2*tangent2));
  const double fast = sqrt(0.5*(sound2+magnetic2+discriminant));
  return fabs(w[1+d]) + fast*normalMag/fabs(det);
}

// Positive diffusive coefficient: the caller adds +div(nu grad(U)). Shock
// viscosity acts only in compression. A cold compressing state saturates the
// limiter; expansion and zero divergence return zero without a 0/0.
inline ACCEL_DECORATION double compressionViscosity(
    double length, double divergence, double speedSquared,
    double alpha, double beta)
{
  const double lambda = length*fmax(-divergence,0.0);
  if (lambda == 0.0) return 0.0;
  const double fraction = speedSquared > 0.0
      ? fmin(1.0,beta*lambda*lambda/speedSquared) : 1.0;
  return alpha*lambda*fraction;
}
}
