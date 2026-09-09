#pragma once

#define c_AU 14959787070000.0 // cm
#define c_SR 69600000000 // cm
#define c_MP 1.67262192e-24 // gm
#define c_PI 3.14159265358979323846
#define c_M0 1.98847e33 // Solar mass in grams
#define c_G 6.67430e-8 // Gravitational constant in cgs units

// Uniform thermodynamic background shared by the idealized pulse and
// strong-shock tests.
#define c_TEST_RHO0 (700.0*c_MP)
#define c_TEST_P0 1.0e-7

// Constant Cartesian Bx background for each magnetized test problem.
// These are separate so each problem can use a different field strength.
#define c_RADIAL_PULSE_BX0 1.0e-2       // Test problem 3
#define c_NONRADIAL_PULSE_BX0 1.0e-2    // Test problem 5
#define c_STRONG_SHOCK_BX0 1.0e-1       // Test problem 7

// FLASH-style off-center blasts (8: hydro, 9: MHD), in CGS units.
#define c_BLAST_PRESSURE_RATIO 1.0e4
#define c_BLAST_RADIUS_FRACTION 0.1
#define c_BLAST_BX0 1.0e-1

// Stationary pulse/blast backgrounds only. Case 11 evolves the full B(t).
inline bool isConstantMagneticFieldTest(int a_initConditionType)
{
  return a_initConditionType == 3
      || a_initConditionType == 5
      || a_initConditionType == 7
      || a_initConditionType == 9;
}

inline double testProblemBackgroundBx(int a_initConditionType)
{
  switch (a_initConditionType)
  {
    case 3: return c_RADIAL_PULSE_BX0;
    case 5: return c_NONRADIAL_PULSE_BX0;
    case 7: return c_STRONG_SHOCK_BX0;
    case 9: return c_BLAST_BX0;
    default: return 0.0;
  }
}
