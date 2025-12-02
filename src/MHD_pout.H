#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <streambuf>

#ifdef PR_MPI
#include <mpi.h>
#endif

// ======================================================================
//  Null stream buffer that discards everything (for muted ranks)
// ======================================================================
class HelioCubedNullBuffer : public std::streambuf
{
protected:
    int overflow(int c) override
    {
        // Pretend we successfully wrote the character
        return traits_type::not_eof(c);
    }

    int sync() override
    {
        return 0;
    }
};

// ======================================================================
//  TeeBuf: writes simultaneously to two underlying streambufs
// ======================================================================
class TeeBuf : public std::streambuf
{
public:
    TeeBuf(std::streambuf* sb1, std::streambuf* sb2)
        : m_sb1(sb1), m_sb2(sb2) {}

protected:
    int overflow(int c) override
    {
        if (c == EOF) return !EOF;

        const int r1 = m_sb1->sputc(c);
        const int r2 = m_sb2->sputc(c);

        if (r1 == EOF || r2 == EOF) return EOF;
        return c;
    }

    int sync() override
    {
        int r1 = m_sb1->pubsync();
        int r2 = m_sb2->pubsync();
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }

private:
    std::streambuf* m_sb1;
    std::streambuf* m_sb2;
};

// ======================================================================
//  Internal helpers to get rank / MPI-initialized status safely
// ======================================================================
inline int helioGetRank()
{
    int rank = 0;

#ifdef PR_MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
#endif

    return rank; // 0 in non-MPI or MPI not initialized
}

// ======================================================================
//  Base pout(): one file per rank; rank 0 also mirrors to std::cout
// ======================================================================
inline std::ostream& pout()
{
    // One per thread; constructed on first use
    static thread_local std::ostream* os = nullptr;

    if (os != nullptr)
    {
        return *os;
    }

    const int rank = helioGetRank();

    // Open per-rank file: pout.<rank>
    std::ostringstream fname;
    fname << "pout." << rank;

    auto fileStream = new std::ofstream(fname.str(),
                                        std::ios::out | std::ios::trunc);

    if (!fileStream->is_open())
    {
        // Failed to open file: fallback behavior differs for rank 0 vs others
        if (rank == 0)
        {
            std::cerr << "WARNING: could not open " << fname.str()
                      << ", falling back to std::cout\n";
            os = &std::cout;
        }
        else
        {
            static thread_local HelioCubedNullBuffer nullBuf;
            static thread_local std::ostream nullStream(&nullBuf);
            std::cerr << "WARNING: could not open " << fname.str()
                      << ", discarding pout() output on rank " << rank << "\n";
            os = &nullStream;
        }
        return *os;
    }

    // Rank 0: tee to both file and std::cout
    if (rank == 0)
    {
        auto teeBuf    = new TeeBuf(fileStream->rdbuf(), std::cout.rdbuf());
        auto teeStream = new std::ostream(teeBuf);
        teeStream->setf(std::ios::unitbuf); // flush after each insertion
        os = teeStream;
        return *os;
    }

    // Other ranks: file only
    fileStream->setf(std::ios::unitbuf);
    os = fileStream;
    return *os;
}

// ======================================================================
//  pout(rootRank): only rootRank writes; others discard
//  - On rootRank: same behavior as pout() (file + cout on rank 0)
//  - On other ranks: discard output
// ======================================================================
inline std::ostream& pout(int rootRank)
{
    const int rank = helioGetRank();

    if (rank != rootRank)
    {
        static thread_local HelioCubedNullBuffer nullBuf;
        static thread_local std::ostream nullStream(&nullBuf);
        return nullStream;
    }

    // This is the selected root rank
    return pout();
}