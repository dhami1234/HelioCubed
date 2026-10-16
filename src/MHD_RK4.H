#pragma once

namespace MHD_RK4 {

// Classical RK4 for the shifted-state interface used by MHDLevelDataRK4Op.
// F replaces dX with dt*f(state + dX, time); dX::init initializes it to zero.
// Kept locally because Proto no longer provides examples/_common/RK4.H.
template <class X, class F, class dX>
class RK4
{
public:
    void advance(double& time, double& dt, X& state)
    {
        double sixth = 1.0/6.0, third = 1.0/3.0, half = 0.5;
        m_delta.init(state);
        m_k.init(state);
        m_f(m_k, time, dt, state);
        m_delta.increment(sixth, m_k);
        m_k *= half;
        m_f(m_k, time + half*dt, dt, state);
        m_delta.increment(third, m_k);
        m_k *= half;
        m_f(m_k, time + half*dt, dt, state);
        m_delta.increment(third, m_k);
        m_f(m_k, time + dt, dt, state);
        m_delta.increment(sixth, m_k);
        state.increment(m_delta);
    }

private:
    dX m_k;
    dX m_delta;
    F m_f;
};

} // namespace MHD_RK4
