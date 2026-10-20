#pragma once

#include "Proto.H"
#include <type_traits>

namespace HelioCubed
{
// Scratch initialization must not exchange an array before its donors are filled.
template<unsigned C>
inline void amrSetValLocal(Proto::MBLevelBoxData<double,C,Proto::HOST>& data,double value)
{
    for(auto index : data.layout()) data[index].setVal(value);
}

// Stage arrays use the same patch ownership. Their valid-cell copies require
// neither communication nor construction of a general redistribution plan.
template<unsigned C>
inline void amrCopyValid(const Proto::MBLevelBoxData<double,C,Proto::HOST>& source,
    Proto::MBLevelBoxData<double,C,Proto::HOST>& destination)
{
    PROTO_ASSERT(source.layout().compatible(destination.layout()),
        "AMR local copy requires identical patch ownership");
    for(auto index : source.layout())
        source[index].copyTo(destination[index],source.layout()[index],Proto::Point::Zeros());
}

// Geometry and MPI packing remain fixed between regrids. Rebind only the data
// arrays when the same transfer is used for a different Runge--Kutta stage.
template<unsigned C,bool Full=false>
class AMRTransferCopier : public std::conditional_t<Full,
    Proto::MBLevelFullCopier<double,C,Proto::HOST,Proto::HOST,Proto::PR_CELL>,
    Proto::MBLevelCopier<double,C,Proto::HOST,Proto::HOST,Proto::PR_CELL>>
{
    using Data=Proto::MBLevelBoxData<double,C,Proto::HOST>;
    using Op=std::conditional_t<Full,
        Proto::MBLevelFullCopierOp<double,C,Proto::HOST,Proto::HOST,Proto::PR_CELL>,
        Proto::MBLevelCopierOp<double,C,Proto::HOST,Proto::HOST,Proto::PR_CELL>>;
    Proto::MBDisjointBoxLayout m_source,m_destination;
public:
    AMRTransferCopier(Data& source,Data& destination,bool unique=false)
        :m_source(source.layout()),m_destination(destination.layout())
    {
        Op op(source,destination);
        if constexpr(Full) op.m_uniqueDonors=unique;
        this->define(op);
    }

    void copy(const Data& source,Data& destination)
    {
        PROTO_ASSERT(m_source.compatible(source.layout()) && m_destination.compatible(destination.layout()),
            "AMR transfer must be rebuilt when its layouts change");
        this->m_op.m_src=const_cast<Data*>(&source);
        this->m_op.m_dst=&destination;
        this->execute();
    }
};
}
