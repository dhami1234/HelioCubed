#pragma once
#include "Proto.H"
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace HelioCubed
{
// Retain the original absolute-coordinate arithmetic, including at negative
// ghost indices and non-power-of-two ratios. Only the data change each stage.
struct AMRSecondOrderInterpolationPlan
{
    struct Coordinate {int parent;std::array<double,3> weights;};
    Proto::Box target;
    std::array<std::vector<Coordinate>,DIM> coordinates;
    std::vector<Proto::Point> offsets;

    AMRSecondOrderInterpolationPlan(const Proto::Box& box,Proto::Point ratio):target(box)
    {
        for(int d=0;d<DIM;++d)
            for(int p=box.low()[d];p<=box.high()[d];++p)
            {
                const int q=int(floor(double(p)/ratio[d]));
                const double x=(p+0.5)/ratio[d]-(q+0.5);
                const double quadratic=x*x+(1.0/(double(ratio[d])*ratio[d])-1)/12;
                coordinates[d].push_back({q,{0.5*(quadratic-x),1-quadratic,0.5*(quadratic+x)}});
            }
        for(auto offset : Proto::Box::Kernel(1))
        {
            bool used=true;
            for(int d=0;d<DIM;++d) if(ratio[d]==1 && offset[d]) used=false;
            if(used) offsets.push_back(offset);
        }
    }
};

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

struct AMRFourthOrderInterpolationPlan
{
    struct Coordinate {int parent;std::array<double,5> weights;};
    Proto::Box target;
    Proto::Point radius;
    std::array<Proto::Box,DIM> passes;
    std::array<std::vector<Coordinate>,DIM> coordinates;

    AMRFourthOrderInterpolationPlan(const Proto::Box& box,Proto::Point ratio):target(box)
    {
        Proto::Point qlo,qhi;
        for(int d=0;d<DIM;++d)
        {
            qlo[d]=int(floor(double(box.low()[d])/ratio[d]));
            qhi[d]=int(floor(double(box.high()[d])/ratio[d]));
            radius[d]=ratio[d]>1 ? 2 : 0;
            std::vector<std::array<double,5>> weights;
            for(int child=0;child<ratio[d];++child) weights.push_back(amrFourthOrderWeights(ratio[d],child));
            for(int p=box.low()[d];p<=box.high()[d];++p)
            {
                const int q=int(floor(double(p)/ratio[d]));
                coordinates[d].push_back({q,weights[p-q*ratio[d]]});
            }
        }
        for(int d=0;d<DIM;++d)
        {
            auto low=box.low(),high=box.high();
            for(int axis=d+1;axis<DIM;++axis) {low[axis]=qlo[axis]-radius[axis];high[axis]=qhi[axis]+radius[axis];}
            passes[d]=Proto::Box(low,high);
        }
    }
};

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
    const Proto::BoxData<double,C,Proto::HOST>& source,const AMRFourthOrderInterpolationPlan& plan)
{
    using namespace Proto;
    BoxData<double,C,HOST> previous;
    for(int d=0;d<DIM;++d)
    {
        BoxData<double,C,HOST> next;
        next.defineForOverwrite(plan.passes[d]);
        const auto& input=d ? previous : source;
        for(auto p : next.box())
        {
            const auto& coordinate=plan.coordinates[d][p[d]-plan.target.low()[d]];
            Point from=p;from[d]=coordinate.parent;
            const auto& weight=coordinate.weights;
            for(unsigned c=0;c<C;++c)
            {
                double value=0;
                for(int k=-plan.radius[d];k<=plan.radius[d];++k)
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

template<unsigned C>
inline Proto::BoxData<double,C,Proto::HOST> amrInterpolateFourthOrder(
    const Proto::BoxData<double,C,Proto::HOST>& source,const Proto::Box& target,Proto::Point ratio)
{
    return amrInterpolateFourthOrder(source,AMRFourthOrderInterpolationPlan(target,ratio));
}
}
