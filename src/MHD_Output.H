#pragma once

#include "Proto.H"
#include "Inputs_Parsing.H"
#include <cmath>
#include <string>
#include <vector>

namespace HelioCubed
{
inline std::vector<std::string> outputPrimitiveNames()
{
    if (ParseInputs::get_output_coordinates() == 1)
        return {"density", "Vx", "Vy", "Vz", "P", "Bx", "By", "Bz"};
    return {"density", "Vr", "Vt", "Vp", "P", "Br", "Bt", "Bp"};
}

// Values may be a Proto Var or an accessor returning a reference by component.
template<class Values>
inline ACCEL_DECORATION void sphericalOutputToCartesian(
    Values& values, double nx, double ny, double nz)
{
    const double phi = atan2(ny, nx), cp = cos(phi), sp = sin(phi);
    const double st = sqrt(nx*nx + ny*ny);
    for (int start : {1, 5})
    {
        const double radial = values(start), theta = values(start+1), azimuth = values(start+2);
        values(start) = nx*radial + nz*cp*theta - sp*azimuth;
        values(start+1) = ny*radial + nz*sp*theta + cp*azimuth;
        values(start+2) = nz*radial - st*theta;
    }
}

template<class Values>
inline ACCEL_DECORATION void formatCartesianOutput(Values& values,
    double x, double y, double z, double gamma, bool cartesian, bool scaled)
{
    if (cartesian && !scaled) return;
    const double radius = sqrt(x*x + y*y + z*z);
    const double nx = x/radius, ny = y/radius, nz = z/radius;
    const double phi = atan2(ny, nx), cp = cos(phi), sp = sin(phi);
    const double st = sqrt(nx*nx + ny*ny);
    for (int start : {1, 5})
    {
        const double vx = values(start), vy = values(start+1), vz = values(start+2);
        values(start) = nx*vx + ny*vy + nz*vz;
        values(start+1) = nz*(cp*vx + sp*vy) - st*vz;
        values(start+2) = -sp*vx + cp*vy;
    }
    if (scaled)
    {
        const double rOverAU = radius/c_AU;
        values(0) *= rOverAU*rOverAU/c_MP;
        values(4) *= pow(rOverAU, 2*gamma)*1e12;
        for (int d = 0; d < DIM; ++d) values(1+d) *= 1e-5;
        values(5) *= rOverAU*rOverAU*1e6;
        values(6) *= rOverAU*1e6;
        values(7) *= rOverAU*1e6;
    }
    // Apply the radial/transverse magnetic scaling before choosing the basis.
    if (cartesian) sphericalOutputToCartesian(values, nx, ny, nz);
}
}
