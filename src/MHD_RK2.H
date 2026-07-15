#ifndef _MHD_RK2_H_
#define _MHD_RK2_H_

#include "Proto.H"

namespace MHD_RK2 {

/// Generic second-order explicit midpoint Runge--Kutta algorithm.
/** Given y' = f(y,t), advance from y_n to y_{n+1} with
 *
 *      y_{n+1/2} = y_n + (dt/2) f(y_n, t_n)
 *      y_{n+1}   = y_n + dt f(y_{n+1/2}, t_n + dt/2).
 *
 * F must interpret the dX passed to it as the current shift from a_state and
 * replace that dX with a_dt*f(a_state + shift, a_time). This is the same
 * interface used by Proto's generic RK4 integrator and MHDLevelDataRK4Op.
 */
template <class X, class F, class dX>
class RK2
{
public:
    void advance(double& a_time, double& a_dt, X& a_state);

protected:
    dX m_k;
    F m_f;
};

template <class X, class F, class dX>
void RK2<X, F, dX>::advance(double& a_time, double& a_dt, X& a_state)
{
    const double half_dt = 0.5*a_dt;

    m_k.init(a_state);
    m_f(m_k, a_time, half_dt, a_state);
    m_f(m_k, a_time + half_dt, a_dt, a_state);
    a_state.increment(m_k);
}

} // namespace MHD_RK2

#endif // _MHD_RK2_H_
