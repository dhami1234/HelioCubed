#pragma once

#include "Proto.H"
#include "MHDOp.H"
#include "MHDLevelDataRK4.H"

namespace MHD_Output_Writer {
    
	void Write_data(MHDLevelDataState& state,
				const int k,
				const double time,
				const double dt);

	void Write_checkpoint(MHDLevelDataState& state,
				const int k,
				const double time,
				const double dt);											   
}
