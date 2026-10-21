#pragma once

#include "Proto.H"
#include <array>
#include <map>
#include <vector>

namespace HelioCubed
{
// Only angular coefficients are retained: radial planes share this geometry.
template<typename T>
struct NativeAngularCoefficients
{
    Proto::Box box;
    std::vector<std::array<T,27>> values;
    // Only AMR boundary extension needs the forward matrix. Keep its storage
    // separate so uniform conversion retains the original coefficient layout.
    std::vector<std::array<T,9>> matrices;

    NativeAngularCoefficients(const Proto::Box& a_output,T a_h,unsigned a_block,bool retainMatrix=false)
    {
        using namespace Proto;
        const Point lo = a_output.low(),hi = a_output.high();
        const Box output(Point(0,lo[1],lo[2]),Point(0,hi[1],hi[2]));
        box = output.grow(Point(0,1,1));
        auto matrix = CubedSphereShell::AMatrix<T,HOST>(box,a_h,a_block,0);
        forallInPlace([] PROTO_LAMBDA(Var<T,DIM,HOST,DIM>& a)
        {
            for (int row = 0; row < DIM; ++row)
                for (int col = 1; col < DIM; ++col) a(row,col) *= T(M_PI/2);
        },matrix);
        auto inverse = forallOverwrite<T,DIM,HOST,DIM>(f_matinv3by3,matrix);
        BoxData<T,DIM,HOST,DIM> first = Stencil<T>::Derivative(1,1)(matrix,output);
        BoxData<T,DIM,HOST,DIM> second = Stencil<T>::Derivative(1,2)(matrix,output);
        values.resize(box.size());
        if(retainMatrix) matrices.resize(box.size());
        for (auto point : box)
        {
            auto& coefficient = values[index(point[1],point[2])];
            for (unsigned row = 0; row < DIM; ++row)
                for (unsigned col = 0; col < DIM; ++col)
                {
                    coefficient[3*row+col] = inverse(point,row,col);
                    if(retainMatrix) matrices[index(point[1],point[2])][3*row+col]=matrix(point,row,col);
                    if (output.containsPoint(point))
                    {
                        coefficient[9+3*row+col] = first(point,row,col);
                        coefficient[18+3*row+col] = second(point,row,col);
                    }
                }
        }
    }

    size_t index(int a_j,int a_k) const
    {
        return (a_k-box.low()[2])*box.size(1) + a_j-box.low()[1];
    }
    const std::array<T,27>& at(int a_j,int a_k) const
    {
        return values[index(a_j,a_k)];
    }
    const std::array<T,9>& matrixAt(int a_j,int a_k) const
    {
        return matrices[index(a_j,a_k)];
    }
};

template<typename T>
inline void nativeQuotientLine(
    T* __restrict__ a_output,
    const T* a_first,const T* a_second,const T* a_third,
    const T* a_inverse,unsigned a_row,int a_count)
{
#ifdef __clang__
#pragma clang loop vectorize(enable)
#endif
    for (int cell = 0; cell < a_count; ++cell)
    {
        T value = T(0);
        value += a_inverse[3*a_row]*a_first[cell]*T(1);
        value += a_inverse[3*a_row+1]*a_second[cell]*T(1);
        value += a_inverse[3*a_row+2]*a_third[cell]*T(1);
        a_output[cell] = value;
    }
}

template<typename T>
inline void nativeCorrectedLine(
    T* __restrict__ a_out0,T* __restrict__ a_out1,T* __restrict__ a_out2,
    const T* const* a_source,const T* const* a_first,const T* const* a_second,
    const T* a_coefficient,int a_count)
{
#ifdef __clang__
#pragma clang loop vectorize(enable)
#endif
    for (int cell = 0; cell < a_count; ++cell)
    {
        T corrected[3];
        for (unsigned row = 0; row < DIM; ++row)
        {
            T first = T(0),second = T(0);
            for (unsigned col = 0; col < DIM; ++col)
            {
                first += a_coefficient[9+3*row+col]*a_first[col][cell]*T(-1.0/12.0);
                second += a_coefficient[18+3*row+col]*a_second[col][cell]*T(-1.0/12.0);
            }
            T correction = T(0);
            correction += first;
            correction += second;
            corrected[row] = a_source[row][cell] + correction;
        }
        T result[3] = {T(0),T(0),T(0)};
        for (unsigned row = 0; row < DIM; ++row)
            for (unsigned col = 0; col < DIM; ++col)
                result[row] += a_coefficient[3*row+col]*corrected[col]*T(1);
        a_out0[cell] = result[0];
        a_out1[cell] = result[1];
        a_out2[cell] = result[2];
    }
}

template<typename T>
struct NativeConversionCache
{
    Proto::MBDisjointBoxLayout layout;
    Proto::Point ghost;
    bool retainMatrix;
    std::map<int,NativeAngularCoefficients<T>> patches;
    NativeConversionCache(const Proto::MBDisjointBoxLayout& grid,Proto::Point g,bool keepMatrix=false)
        :layout(grid),ghost(g),retainMatrix(keepMatrix) {}
};

// The same fourth-order quotient as cellMatrixQuotientMatAndInv, with the
// two vector groups sharing derivatives and radius-independent coefficients.
template<typename T,unsigned C>
inline void cachedNativeConversion(
    Proto::BoxData<T,C,Proto::HOST>& a_output,
    const Proto::BoxData<T,C,Proto::HOST>& a_cart,
    const Proto::MBDisjointBoxLayout& a_layout,
    Proto::Point a_ghost,int a_index,unsigned a_block,T a_h,
    NativeConversionCache<T>* retained=nullptr)
{
    using namespace Proto;
    PR_TIMERS("Cartesian conversion: cached U to Usph");
    static_assert(DIM == 3 && C == 8);
    static std::unique_ptr<NativeConversionCache<T>> fallback;
    if(!retained)
    {
        if(!fallback || &fallback->layout.partition()!=&a_layout.partition() || fallback->ghost!=a_ghost)
            fallback=std::make_unique<NativeConversionCache<T>>(a_layout,a_ghost);
        retained=fallback.get();
    }
    PROTO_ASSERT(&retained->layout.partition()==&a_layout.partition() && retained->ghost==a_ghost,
        "Native conversion geometry must match its AMR level");
    auto entry = retained->patches.find(a_index);
    if (entry == retained->patches.end())
        entry = retained->patches.emplace(std::piecewise_construct,
            std::forward_as_tuple(a_index),
            std::forward_as_tuple(a_output.box(),a_h,a_block,retained->retainMatrix)).first;
    const auto& geometry = entry->second;
    const Box outputBox = a_output.box();
    const Box quotientBox = outputBox.grow(Point(0,1,1));
    BoxData<T,6,HOST> quotient;
    quotient.defineForOverwrite(quotientBox);
    Point lo = quotientBox.low(),hi = quotientBox.high();
    for (int k = lo[2]; k <= hi[2]; ++k)
        for (int j = lo[1]; j <= hi[1]; ++j)
        {
            const auto& coefficient = geometry.at(j,k);
            const Point start(lo[0],j,k);
            for (unsigned group = 0; group < 2; ++group)
            {
                const unsigned component = group == 0 ? 1 : DIM+2;
                for (unsigned row = 0; row < DIM; ++row)
                    nativeQuotientLine(quotient.data(start,group*3+row),
                        a_cart.data(start,component),a_cart.data(start,component+1),
                        a_cart.data(start,component+2),coefficient.data(),row,
                        quotientBox.size(0));
            }
        }
    BoxData<T,6,HOST> first = Stencil<T>::Derivative(1,1)(quotient,outputBox);
    BoxData<T,6,HOST> second = Stencil<T>::Derivative(1,2)(quotient,outputBox);
    a_cart.copyTo(a_output);
    lo = outputBox.low();
    hi = outputBox.high();
    for (int k = lo[2]; k <= hi[2]; ++k)
        for (int j = lo[1]; j <= hi[1]; ++j)
        {
            const auto& coefficient = geometry.at(j,k);
            const Point start(lo[0],j,k);
            for (unsigned group = 0; group < 2; ++group)
            {
                const unsigned component = group == 0 ? 1 : DIM+2;
                const T* source[3];
                const T* derivativeFirst[3];
                const T* derivativeSecond[3];
                for (unsigned col = 0; col < DIM; ++col)
                {
                    source[col] = a_cart.data(start,component+col);
                    derivativeFirst[col] = first.data(start,group*3+col);
                    derivativeSecond[col] = second.data(start,group*3+col);
                }
                nativeCorrectedLine(a_output.data(start,component),
                    a_output.data(start,component+1),a_output.data(start,component+2),
                    source,derivativeFirst,derivativeSecond,coefficient.data(),
                    outputBox.size(0));
            }
        }
}
}
