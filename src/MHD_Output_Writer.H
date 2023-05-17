#pragma once

#include "Proto.H"
#include "MHDOp.H"
#include "MHDLevelDataRK4.H"

namespace MHD_Output_Writer {
    
	void Write_data(AMRData<double,NUMCOMPS>& U,
					AMRGrid& grid,
					const double dx,
					const double dy,
					const double dz,
					const int k,
					const double time,
					const double dt,
					bool CME_data);

	void Write_checkpoint(AMRData<double,NUMCOMPS>& U,
					const int k,
					const double time,
					const double dt,
					bool CME_checkpoint);											   
}
