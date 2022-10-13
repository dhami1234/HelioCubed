# HelioCubed
MHD code to study and forecast space weather 

This code uses proto library developed by the Applied Numerical Algorithms Group (ANAG) at Lawrence Berkeley National Lab (LBNL). The code can be downloaded from https://github.com/applied-numerical-algorithms-group-lbnl/proto and its ath included in the makefile. Proto handels the domain decomposition and CPU/GPU parallelization.  

This code also requires blis and mpi-hdf5, whose paths should be included in makefile
