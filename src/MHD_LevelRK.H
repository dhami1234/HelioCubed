#pragma once

#include "Proto.H"
#include <stdexcept>

namespace HelioCubed
{
// Preserve Proto's existing Euler/RK3/RK4 path. SSPRK2 evaluates the complete
// RHS twice, including sources and the reconciled panel flux at each stage.
template<template<typename,typename,Proto::MemType> class OPType,
         class MAP, typename T, Proto::MemType MEM = Proto::MEMTYPE_DEFAULT>
class MBLevelRK : public Proto::MBLevelRK4<OPType,MAP,T,Proto::NullBC,MEM>
{
    using Base = Proto::MBLevelRK4<OPType,MAP,T,Proto::NullBC,MEM>;
    using OP = typename Base::OP;
public:
    using Base::Base;
    using Base::advance;
    using Level = typename Base::LevelStateData;
    using ScalarLevel = typename Base::LevelScalarData;

    void advance(Level& state, ScalarLevel& volume, T dt, T time, int order)
    {
        if (order != 2)
        {
            if (order != 1 && order != 3 && order != 4)
                throw std::invalid_argument("Unsupported temporal order");
            Base::advance(state,volume,dt,time,order);
            return;
        }
        if (!(dt > 0)) throw std::invalid_argument("RK2 requires a positive timestep");
        PROTO_ASSERT(this->m_stageInterp,
            "The MHD integrator requires a mutable interpolation operator");
        PR_TIMERS("HelioCubed::SSPRK2");
        if (!this->m_stageDefined)
        {
            this->m_stageState.define(state.layout(),OP::ghost());
            this->m_stageDefined = true;
        }
        for (int stage = 0; stage < 2; ++stage)
        {
            state.copyTo(this->m_stageState);
            if (stage == 1) this->m_stageState.increment(this->m_k1,dt);
            this->m_blockreg.clear();
            const T stageTime = time + stage*dt;
            OP::preStageLevel(this->m_stageState,this->m_k,volume,
                *this->m_stageInterp,this->m_blockreg,this->m_consSum,
                stageTime,dt,T(0.5),stage,true);
            OP::StageLevel(this->m_k,this->m_stageState,this->m_blockreg,
                this->m_consSum,stageTime,dt,T(0.5),stage);
            // The register stores time-integrated fluxes. Convert its
            // correction back to a rate before constructing the next stage.
            this->m_delta.setVal(0);
            this->m_blockreg.reflux(this->m_delta);
            this->m_k.increment(this->m_delta,T(2)/dt);
            if (stage == 0) this->m_k.copyTo(this->m_k1);
        }
        state.increment(this->m_k1,T(0.5)*dt);
        state.increment(this->m_k,T(0.5)*dt);
        ++this->m_counter;
    }
};
}
