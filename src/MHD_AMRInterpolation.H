#pragma once
#include "Proto.H"
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace HelioCubed
{
// Integrate the quartic matching five coarse cell averages over a child cell.
// The mean of the children is the central coarse average, for any integer
// ratio. An unrefined direction is the exact identity.
inline std::array<double,5> amrFourthOrderWeights(int ratio,int child)
{
    if(ratio==1) return {0,0,1,0,0};
    const double x=(child+0.5)/ratio-0.5,inverse=1.0/ratio;
    const std::array<double,5> moment={1,x,x*x+inverse*inverse/12,
        x*x*x+x*inverse*inverse/4,
        x*x*x*x+x*x*inverse*inverse/2+inverse*inverse*inverse*inverse/80};
    static constexpr double inverseMoments[5][5]={
        {3.0/640,-29.0/480,1067.0/960,-29.0/480,3.0/640},
        {5.0/48,-17.0/24,0.0/1,17.0/24,-5.0/48},
        {-1.0/16,3.0/4,-11.0/8,3.0/4,-1.0/16},
        {-1.0/12,1.0/6,0.0/1,-1.0/6,1.0/12},
        {1.0/24,-1.0/6,1.0/4,-1.0/6,1.0/24}};
    std::array<double,5> weight{};
    for(int i=0;i<5;++i) for(int n=0;n<5;++n) weight[i]+=moment[n]*inverseMoments[n][i];
    return weight;
}

// All components are Cartesian averages in the same block coordinates.
template<unsigned C>
inline bool amrInterpolateFourthOrder(const Proto::BoxData<double,C,Proto::HOST>& source,
    Proto::Point p,Proto::Point ratio,std::array<double,C>& value)
{
    using namespace Proto;
    Point q;std::array<std::array<double,5>,DIM> weights;
    for(int d=0;d<DIM;++d)
    {
        q[d]=int(floor(double(p[d])/ratio[d]));
        weights[d]=amrFourthOrderWeights(ratio[d],p[d]-q[d]*ratio[d]);
    }
    value.fill(0);
    for(auto offset : Box::Kernel(2))
    {
        double weight=1;
        for(int d=0;d<DIM;++d) weight*=weights[d][offset[d]+2];
        if(weight==0) continue;
        const Point from=q+offset;
        if(!source.box().containsPoint(from)) return false;
        for(unsigned c=0;c<C;++c)
        {
            if(!std::isfinite(source(from,c))) return false;
            value[c]+=weight*source(from,c);
        }
    }
    return true;
}

// Apply the tensor interpolant as three one-dimensional passes. This avoids
// 125 donor reads per component and ghost cell in every RK stage.
template<unsigned C>
inline Proto::BoxData<double,C,Proto::HOST> amrInterpolateFourthOrder(
    const Proto::BoxData<double,C,Proto::HOST>& source,const Proto::Box& target,Proto::Point ratio)
{
    using namespace Proto;
    Point qlo,qhi,radius;
    std::array<std::vector<std::array<double,5>>,DIM> weights;
    for(int d=0;d<DIM;++d)
    {
        qlo[d]=int(floor(double(target.low()[d])/ratio[d]));
        qhi[d]=int(floor(double(target.high()[d])/ratio[d]));
        radius[d]=ratio[d]>1 ? 2 : 0;
        for(int child=0;child<ratio[d];++child) weights[d].push_back(amrFourthOrderWeights(ratio[d],child));
    }
    BoxData<double,C,HOST> previous;
    for(int d=0;d<DIM;++d)
    {
        Point low=target.low(),high=target.high();
        for(int axis=d+1;axis<DIM;++axis) {low[axis]=qlo[axis]-radius[axis];high[axis]=qhi[axis]+radius[axis];}
        BoxData<double,C,HOST> next;
        next.defineForOverwrite(Box(low,high));
        const auto& input=d ? previous : source;
        for(auto p : next.box())
        {
            Point from=p;from[d]=int(floor(double(p[d])/ratio[d]));
            const auto& weight=weights[d][p[d]-from[d]*ratio[d]];
            for(unsigned c=0;c<C;++c)
            {
                double value=0;
                for(int k=-radius[d];k<=radius[d];++k)
                {
                    if(weight[k+2]==0) continue;
                    const Point point=from+Point::Basis(d,k);
                    if(!input.box().containsPoint(point)) {value=std::numeric_limits<double>::quiet_NaN();break;}
                    value+=weight[k+2]*input(point,c);
                }
                next(p,c)=value;
            }
        }
        previous=std::move(next);
    }
    return previous;
}
}
