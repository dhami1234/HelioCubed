#pragma once

#include <array>
#include <cmath>

namespace HelioCubed
{
inline double probeInterpolationWeight(double coordinate, int center, int neighbor)
{
    // Both samples are cell-centered; their half-cell offsets cancel.
    return (coordinate - (center + 0.5)) / (neighbor - center);
}

// theta is colatitude and phi is longitude, both in radians. Using angles
// directly also defines the transverse basis at the poles for the given phi.
inline std::array<double,3> projectProbeVectorToRTN(
    double x, double y, double z, double theta, double phi)
{
    const double st=std::sin(theta), ct=std::cos(theta);
    const double sp=std::sin(phi), cp=std::cos(phi);
    const double radial=st*(cp*x+sp*y)+ct*z;
    const double polar=ct*(cp*x+sp*y)-st*z;
    const double azimuth=-sp*x+cp*y;
    return {radial,azimuth,-polar}; // RTN: T=e_phi, N=-e_theta.
}
}
