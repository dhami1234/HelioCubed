#pragma once

#define c_AU 14959787070000.0 // cm
#define c_SR 69600000000 // cm
#define c_MP 1.67262192e-24 // gm
#define c_PI 3.14159265358979323846
#define c_M0 1.98847e33 // Solar mass in grams
#define c_G 6.67430e-8 // Gravitational constant in cgs units

// Uniform background state shared by the idealized pulse and strong-shock
// tests. Test problems 3, 5, and 7 evolve delta-B about this Cartesian field.
#define c_TEST_RHO0 (700.0*c_MP)
#define c_TEST_P0 1.0e-7
#define c_TEST_BX0 1.0e-1

inline bool isConstantMagneticFieldTest(int a_initConditionType)
{
  return a_initConditionType == 3
      || a_initConditionType == 5
      || a_initConditionType == 7;
}
