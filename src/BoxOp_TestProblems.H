#pragma once

#include "MHD_Constants.H"

// Analytic initial and radial-boundary states for the idealized test problems.
// All magnetic fields are specified in Cartesian physical components and then
// transformed to the cubed-sphere basis used by the solver.

template <typename T, MemType MEM>
inline void testProblemBasisInverse(
    BoxData<T, DIM, MEM, DIM>& a_Ainv,
    const Box& a_box,
    T a_h,
    uint a_block)
{
  BoxData<T, DIM, MEM, DIM> A =
      CubedSphereShell::detA_A<T, MEM>(a_box, a_h, a_block, 0);
  BoxData<T, 1, MEM> detA =
      CubedSphereShell::detA<T, MEM>(a_box, a_h, 0);

  forallInPlace(
      [] PROTO_LAMBDA(
          Var<T, DIM, MEM, DIM>& a_A,
          const Var<T, 1, MEM>& a_detA)
      {
        const T invDetA = 1.0/a_detA(0);
        for (int ii = 0; ii < DIM; ++ii)
        {
          for (int jj = 0; jj < DIM; ++jj)
          {
            a_A(ii, jj) *= invDetA;
          }
        }
      }, A, detA);

  a_Ainv.setToZero();
  forallInPlace(f_matinv3by3, a_Ainv, A);
}

template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_smoothSphericalWaveInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM>& a_W,
    T a_gamma,
    int a_nradius)
{
  const T p0 = 1.0e-7;
  const T rho0 = 700*c_MP;
  const T eps = 0.1;
  const T arg = (a_pt[0] + 0.5 - 0.5*a_nradius)/(2.0*a_nradius);
  T amplitude = 0.0;
  if (abs(arg) < 0.25)
  {
    amplitude = eps*pow(cos(2.0*M_PI*arg), 6);
  }

  const T rho = rho0*(1.0 + amplitude);
  const T p = p0*pow(rho/rho0, a_gamma);
  const T ur = amplitude*sqrt(a_gamma*p0/rho0);
  a_W(iRHO) = rho;
  a_W(iVX) = 0.0;//ur;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}
PROTO_KERNEL_END(f_smoothSphericalWaveInit_F, f_smoothSphericalWaveInit)

template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_nonRadialSphericalWaveInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM>& a_W,
    const Var<T, DIM, MEM>& a_Xcart,
    T a_gamma)
{
  const T p0 = 1.0e-7;
  const T rho0 = 700*c_MP;
  const T eps = 0.1;

  // Put a compact spherical pulse halfway through the shell, in a direction
  // that is not aligned with a cubed-sphere coordinate axis.
  const T dirNorm = sqrt(0.9*0.9 + 0.4*0.4);
  const T centerRadius = 0.5*(CUBED_SPHERE_SHELL_R0 + CUBED_SPHERE_SHELL_R1);
  const T center[DIM] = {
      -0.9*centerRadius/dirNorm,
      -0.4*centerRadius/dirNorm,
      0.0};
  const T waveRadius = 0.2*(CUBED_SPHERE_SHELL_R1 - CUBED_SPHERE_SHELL_R0);
  T displacement[DIM];
  T distance2 = 0.0;
  for (int dir = 0; dir < DIM; ++dir)
  {
    displacement[dir] = a_Xcart(dir) - center[dir];
    distance2 += displacement[dir]*displacement[dir];
  }
  const T distance = sqrt(distance2);
  const T q = distance/waveRadius;
  T amplitude = 0.0;
  if (q < 1.0)
  {
    amplitude = eps*pow(cos(0.5*M_PI*q), 6);
  }

  const T rho = rho0*(1.0 + amplitude);
  const T p = p0*pow(rho/rho0, a_gamma);

  // Start the density/pressure pulse at rest, including at its center.
  a_W(iRHO) = rho;
  a_W(iVX) = 0.0;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}
PROTO_KERNEL_END(f_nonRadialSphericalWaveInit_F, f_nonRadialSphericalWaveInit)

template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_symmetricStrongShockInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM>& a_W,
    T a_gamma,
    T a_strength,
    int a_nradius)
{
  const T p0 = 1.0e-7;
  const T rho0 = 700*c_MP;
  const T arg = (a_pt[0] + 0.5 - 0.5*a_nradius)/(2.0*a_nradius);
  T amplitude = 0.0;
  if (arg < 0.0)
  {
    amplitude = a_strength*pow(cos(2.0*M_PI*arg), 6);
  }

  const T p = p0*(1.0 + amplitude);
  const T shockImpedance = sqrt(
      a_gamma*p0*rho0*
      (1.0 + ((a_gamma + 1.0)/(2.0*a_gamma))*(p - p0)/p0));
  const T rho = 1.0/(1.0/rho0 - (p - p0)/(shockImpedance*shockImpedance));
  const T ur = (p - p0)/shockImpedance;
  a_W(iRHO) = rho;
  a_W(iVX) = ur;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}
PROTO_KERNEL_END(f_symmetricStrongShockInit_F, f_symmetricStrongShockInit)

template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_setConstantCartesianBPrimitive_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM>& a_W,
    const Var<T, DIM, MEM, DIM>& a_Ainv,
    T a_Bx,
    T a_By,
    T a_Bz)
{
  const T Bcart[DIM] = {a_Bx, a_By, a_Bz};
  for (int ii = 0; ii < DIM; ++ii)
  {
    a_W(iBX + ii) = 0.0;
    for (int jj = 0; jj < DIM; ++jj)
    {
      a_W(iBX + ii) += a_Ainv(ii, jj)*Bcart[jj];
    }
  }
}
PROTO_KERNEL_END(
    f_setConstantCartesianBPrimitive_F,
    f_setConstantCartesianBPrimitive)

template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_testProblemHydroBoundary_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM>& a_U,
    T a_gamma)
{
  const T p0 = 1.0e-7;
  const T rho0 = 700*c_MP;
  a_U(iRHO) = rho0;
  a_U(iMOMX) = 0.0;
  a_U(iMOMY) = 0.0;
  a_U(iMOMZ) = 0.0;
  a_U(iE) = p0/(a_gamma - 1.0);
  a_U(iBX) = 0.0;
  a_U(iBY) = 0.0;
  a_U(iBZ) = 0.0;
}
PROTO_KERNEL_END(f_testProblemHydroBoundary_F, f_testProblemHydroBoundary)

template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_testProblemMagneticBoundary_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM>& a_U,
    const Var<T, DIM, MEM, DIM>& a_Ainv,
    T a_gamma)
{
  const T p0 = 1.0e-7;
  const T rho0 = 700*c_MP;
  const T Bcart[DIM] = {1.0e-2, 0.0, 0.0};

  a_U(iRHO) = rho0;
  a_U(iMOMX) = 0.0;
  a_U(iMOMY) = 0.0;
  a_U(iMOMZ) = 0.0;
  a_U(iE) = p0/(a_gamma - 1.0);
  for (int ii = 0; ii < DIM; ++ii)
  {
    a_U(iBX + ii) = 0.0;
    for (int jj = 0; jj < DIM; ++jj)
    {
      a_U(iBX + ii) += a_Ainv(ii, jj)*Bcart[jj];
    }
  }
  a_U(iE) += 1.0e-8/(8.0*M_PI);
}
PROTO_KERNEL_END(
    f_testProblemMagneticBoundary_F,
    f_testProblemMagneticBoundary)
