#ifndef _MHD_HANCOCK_H_
#define _MHD_HANCOCK_H_

#include "MHDOp.H"

namespace MHD_Hancock {

/// MUSCL--Hancock integrator for HelioCubed's second-order spherical operator.
class Hancock
{
public:
    void advance(double& a_time, double& a_dt, MHDLevelDataState& a_state);
};

inline void Hancock::advance(
        double& a_time,
        double& a_dt,
        MHDLevelDataState& a_state)
{
    MHDOp::step_spherical_2O_hancock(a_state, a_time, a_dt);
}

} // namespace MHD_Hancock

#endif // _MHD_HANCOCK_H_
