#pragma once

#include "Proto.H"
#include <cstdio>
#include <cstdlib>

namespace MHD {

// Input and I/O failures must remain errors in optimized NDEBUG builds.
template<typename... Args>
inline void require(bool condition, const char* message, Args... args)
{
    if (condition) return;
    int rank = 0;
#ifdef PR_MPI
    int initialized = 0, finalized = 0;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Finalized(&finalized);
    const bool mpiActive = initialized && !finalized;
    if (mpiActive) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    std::fprintf(stderr, "\nHelioCubed ERROR [MPI rank %d]: ", rank);
    std::fprintf(stderr, message, args...);
    std::fputc('\n', stderr);
    std::fflush(stderr);
#ifdef PR_MPI
    if (mpiActive) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    std::abort();
}

} // namespace MHD
