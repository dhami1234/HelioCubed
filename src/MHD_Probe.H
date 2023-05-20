#pragma once

#include "Proto.H"
#include "MHDOp.H"
#include "MHDLevelDataRK4.H"

namespace MHD_Probe
{
    void Probe(MHDLevelDataState &state,
               const double a_time,
               bool give_space);
    double getPhysTime(double a_time);
}