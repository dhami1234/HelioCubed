#pragma once

#include "Proto.H"
#include "MHDOp.H"
#include "MHDLevelDataRK4.H"

namespace MHD_Pre_Time_Step {
    void Define_CME(MHDLevelDataState& a_state);
    
    void Insert_CME(MHDLevelDataState& state,
                const double a_time);  
}