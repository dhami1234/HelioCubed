#pragma once

#include "Proto.H"
#include <cmath>

namespace HelioCubed
{
// MS-FLUKSS GETRELGRAD followed by MAGNITUDE. Only radial faces are
// physical boundaries; panel seams and patch edges use exchanged neighbors.
template<class Sample>
inline double amrRelativeGradientMagnitude(Proto::Point cell,const Proto::Box& domain,
    const Sample& sample)
{
    double squared=0;
    for(int d=0;d<DIM;++d)
    {
        const auto e=Proto::Point::Basis(d);
        const bool lower=d==0 && cell[d]==domain.low()[d];
        const bool upper=d==0 && cell[d]==domain.high()[d];
        const double low=sample(lower ? cell : cell-e);
        const double high=sample(upper ? cell : cell+e);
        const double difference=(lower || upper ? 1.0 : 0.5)*(high-low);
        const double mean=0.5*(high+low);
        const double gradient=std::abs(mean)<1e-32 ? 0 : difference/mean;
        squared+=gradient*gradient;
    }
    return std::sqrt(squared);
}
}
