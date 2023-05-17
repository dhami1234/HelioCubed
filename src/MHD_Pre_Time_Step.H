#pragma once

#include "Proto.H"
#include "MHDOp.H"
#include "MHDLevelDataRK4.H"

namespace MHD_Pre_Time_Step {
    void Define_CME(MHDLevelDataState& a_state);
    
    void Insert_CME(AMRData<double,NUMCOMPS>& a_U,
					AMRGrid& a_grid,
					const double a_dx,
					const double a_dy,
					const double a_dz,
                    const int a_k,
                    const double a_time,
                    const double a_dt);  
}