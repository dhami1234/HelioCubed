#pragma once

#include "Proto.H"
#include "Inputs_Parsing.H"
#include "MHD_CachedCartesianConversion.H"

namespace HelioCubed
{
// Recover Cartesian conserved cell averages before interblock interpolation.
// Unlike the semispherical path, this does not preserve radial symmetry exactly.
// Both divisions and the final basis change include fourth-order average
// corrections; pointwise JU/J and A^{-1}U would only be second-order accurate.
template<typename T, unsigned C, Proto::MemType MEM>
void consToSphInterpEulerCartesian(
    Proto::MBLevelBoxData<T,C,MEM>& a_JU,
    Proto::MBInterpOp& a_iop,
    Proto::MBLevelBoxData<T,1,MEM>& a_dVolr)
{
    using namespace Proto;
    PR_TIMERS("HelioCubed::consToSphInterpEulerCartesian()");
    static_assert(DIM == 3 && C == 2*DIM+2,
        "Cartesian state conversion requires three-dimensional MHD");
    const auto& layout = a_JU.layout();
    const T h = T(1)/layout.domain().getBlock(0).sizes()[1];
    const bool localCopies = MEM == HOST && ParseInputs::get_spatial_operator() == 2;

    // JU is valid only in each block's own coordinate system. Exchange only
    // within a block before dividing; one-sided differences handle block edges.
    Array<Point,DIM+1> localGhost;
    localGhost.fill(Point::Zeros());
    localGhost[0] = Point::Ones(2);
    MBLevelBoxData<T,C,MEM> ju;
    if (localCopies) ju.defineForOverwrite(layout,localGhost);
    else ju.define(layout,localGhost);
    // Level setVal also exchanges. Zero local storage now; exchange only
    // after the source data have been populated.
    for (auto index : layout) ju[index].setVal(T(0));
    if (localCopies)
        for (auto index : layout)
            a_JU[index].copyTo(ju[index],layout[index],Point::Zeros());
    else
        a_JU.copyTo(ju);
    ju.exchange0();

    // The final matrix quotient consumes one angular layer. Interpolation
    // sources require two more. No interblock buffers are needed on the output.
    const auto ghost = a_JU.ghost() + Point(0,1,1);
    localGhost[0] = ghost[0];
    MBLevelBoxData<T,C,MEM> input,cart;
    if (localCopies)
    {
        input.defineForOverwrite(layout,ghost + Point(0,2,2));
        cart.defineForOverwrite(layout,localGhost);
    }
    else
    {
        input.define(layout,ghost + Point(0,2,2));
        cart.define(layout,localGhost);
    }
    for (auto index : layout)
    {
        input[index].setVal(T(0));
        if (!localCopies) cart[index].setVal(T(0));
    }
    for (auto index : layout)
    {
        PR_TIMERS("Cartesian conversion: JU to U");
        ScopedHostStack<MEM> scratch;
        const Box valid = layout[index];
        const Box block = layout.getBlock(layout.block(index)).domain().box();
        // Retain the radial neighbors needed by the second division, including
        // two inward neighbors for the one-sided physical-boundary derivative.
        const Box radialBox = valid.grow(0,2) & block;
        const Box quotientBox = radialBox.grow(Point(0,2,2)) & block;
        const Box metricBox = quotientBox.grow(Point(0,1,1));
        auto det = CubedSphereShell::detA<T,MEM>(metricBox,h,0);
        auto quotient = forallOverwrite<T,C,MEM>([] PROTO_LAMBDA(
            Var<T,C,MEM>& out, const Var<T,C,MEM>& num,
            const Var<T,1,MEM>& den)
        {
            const T inverse = T(1)/den(0);
            for (unsigned c = 0; c < C; ++c) out(c) = num(c)*inverse;
        },quotientBox,ju[index],det);
        BoxData<T,C,MEM> correction(radialBox,T(0));
        for (int dir = 1; dir < DIM; ++dir)
        {
            auto dq = CubedSphereShell::deriv3case(quotient,radialBox,block,dir);
            BoxData<T,1,MEM> dd = Stencil<T>::Derivative(1,dir)(det,radialBox);
            forallInPlace([] PROTO_LAMBDA(
                Var<T,C,MEM>& corr, const Var<T,C,MEM>& dq,
                const Var<T,1,MEM>& dd)
            {
                for (unsigned c = 0; c < C; ++c) corr(c) += dq(c)*dd(0);
            },correction,dq,dd);
        }
        auto angular = forallOverwrite<T,C,MEM>([] PROTO_LAMBDA(
            Var<T,C,MEM>& out, const Var<T,C,MEM>& num,
            const Var<T,C,MEM>& corr, const Var<T,1,MEM>& det)
        {
            const T inverse = T(1)/det(0);
            for (unsigned c = 0; c < C; ++c)
                out(c) = (num(c) - corr(c)/T(12))*inverse;
        },radialBox,ju[index],correction,det);
        auto radialQuotient = forallOverwrite<T,C,MEM>([] PROTO_LAMBDA(
            Var<T,C,MEM>& out, const Var<T,C,MEM>& num,
            const Var<T,1,MEM>& volume)
        {
            const T inverse = T(1)/volume(0);
            for (unsigned c = 0; c < C; ++c) out(c) = num(c)*inverse;
        },radialBox,angular,a_dVolr[index]);
        auto dq = CubedSphereShell::deriv3case(radialQuotient,valid,block,0);
        BoxData<T,1,MEM> dv = Stencil<T>::Derivative(1,0)(a_dVolr[index],valid);
        forallInPlace([] PROTO_LAMBDA(
            Var<T,C,MEM>& out, const Var<T,C,MEM>& num,
            const Var<T,1,MEM>& volume, const Var<T,C,MEM>& dq,
            const Var<T,1,MEM>& dv)
        {
            const T inverse = T(1)/volume(0);
            for (unsigned c = 0; c < C; ++c)
                out(c) = (num(c) - dq(c)*dv(0)/T(12))*inverse;
        },valid,input[index],angular,a_dVolr[index],dq,dv);
    }

    // Every vector component has the same global Cartesian meaning on all
    // blocks. Interpolate these averages before rotating into the local basis.
    input.exchange();
    if (localCopies)
    {
        // exchange() has populated input's wider same-block halos. Copy those
        // complete patch values locally; interpolation changes only points
        // outside the block, leaving same-block halos ready for conversion.
        for (auto index : layout) input[index].copyTo(cart[index]);
    }
    else
        input.copyTo(cart);
    a_iop.apply(cart,input);
    if (!localCopies) cart.exchange0();
    for (auto index : layout)
    {
        PR_TIMERS("Cartesian conversion: U to Usph");
        ScopedHostStack<MEM> scratch;
        const Box outputBox = a_JU[index].box();
        if constexpr (MEM == HOST)
        {
            if (ParseInputs::get_spatial_operator() == 2)
            {
                cachedNativeConversion(a_JU[index],cart[index],layout,
                    a_JU.ghost()[0],index.global(),layout.block(index),h);
                continue;
            }
        }
        const Box metricBox = outputBox.grow(Point(0,1,1));
        auto matrix = CubedSphereShell::AMatrix<T,MEM>(metricBox,h,layout.block(index),0);
        // AMatrix averages derivatives with respect to angular radians. The
        // evolution operator uses normalized block coordinates: alpha,beta =
        // -pi/4 + (pi/2)*xi. Match f_Amatrix and detA_A's tangential columns.
        forallInPlace([] PROTO_LAMBDA(Var<T,DIM,MEM,DIM>& a)
        {
            for (int row = 0; row < DIM; ++row)
                for (int col = 1; col < DIM; ++col) a(row,col) *= T(M_PI/2);
        },matrix);
        auto inverse = forallOverwrite<T,DIM,MEM,DIM>(f_matinv3by3,matrix);
        cart[index].copyTo(a_JU[index]);
        for (unsigned component : {1u, unsigned(DIM+2)})
        {
            auto source = slice<T,C,DIM,MEM>(cart[index],component);
            auto target = slice<T,C,DIM,MEM>(a_JU[index],component);
            BoxData<T,DIM,MEM> quotient;
            quotient.defineForOverwrite(metricBox);
            CubedSphereShell::cellMatrixQuotientMatAndInv(
                target,quotient,matrix,inverse,source,false);
        }
    }
    // Fill same-block patch ghosts from the owner's converted result.
    a_JU.exchange0();
}

template<typename T, unsigned C, Proto::MemType MEM>
void consToSphInterpEuler(
    Proto::MBLevelBoxData<T,C,MEM>& a_JU,
    Proto::MBInterpOp& a_iop,
    Proto::MBLevelBoxData<T,1,MEM>& a_dVolr,
    unsigned a_order = 4, int a_counter = 0)
{
    if (ParseInputs::get_state_conversion() == 1)
        consToSphInterpEulerCartesian(a_JU,a_iop,a_dVolr);
    else
        Proto::CubedSphereShell::consToSphInterpEuler(
            a_JU,a_iop,a_dVolr,a_order,a_counter);
}
}
