#pragma once

#include "Proto.H"
#include <type_traits>
#include <map>
#include <vector>

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
    struct SourceRegion
    {
        Proto::MBIndex index;
        std::vector<Proto::Box> boxes;
    };
private:
    std::vector<SourceRegion> m_sources;
public:
    AMRTransferCopier(Data& source,Data& destination,bool unique=false,bool retainSources=false)
        :m_source(source.layout()),m_destination(destination.layout())
    {
        Op op(source,destination);
        if constexpr(Full) op.m_uniqueDonors=unique;
        this->define(op);
        // Local copies and outgoing messages are the only readers of this
        // rank's source arrays. Retain their union without changing the copy
        // plans or the donor selected for any destination cell.
        if constexpr(Full)
        {
            if(!retainSources) return;
            std::map<uint64_t,size_t> byIndex;
            for(auto type : {Proto::LOCAL,Proto::FROM})
                for(const auto& item : this->motionPlan(type))
                {
                    const auto key=item.fromIndex.global();
                    auto entry=byIndex.emplace(key,m_sources.size());
                    if(entry.second) m_sources.push_back({item.fromIndex,{}});
                    auto& boxes=m_sources[entry.first->second].boxes;
                    std::vector<Proto::Box> pieces{item.fromRegion};
                    for(const auto& existing : boxes)
                    {
                        std::vector<Proto::Box> remainder;
                        for(const auto& box : pieces)
                        {
                            const auto overlap=box&existing;
                            if(overlap.empty()) {remainder.push_back(box);continue;}
                            auto low=box.low(),high=box.high();
                            for(int d=0;d<DIM;++d)
                            {
                                if(low[d]<overlap.low()[d])
                                {
                                    auto end=high;end[d]=overlap.low()[d]-1;
                                    remainder.emplace_back(low,end);low[d]=overlap.low()[d];
                                }
                                if(high[d]>overlap.high()[d])
                                {
                                    auto begin=low;begin[d]=overlap.high()[d]+1;
                                    remainder.emplace_back(begin,high);high[d]=overlap.high()[d];
                                }
                            }
                        }
                        pieces.swap(remainder);
                        if(pieces.empty()) break;
                    }
                    boxes.insert(boxes.end(),pieces.begin(),pieces.end());
                }
            // Dense refinement often covers one rectangular source region.
            // Recover that single loop instead of traversing its copy-plan
            // fragments at every stage, without filling any additional cells.
            for(auto& source : m_sources)
            {
                auto low=source.boxes.front().low(),high=source.boxes.front().high();
                size_t cells=0;
                for(const auto& box : source.boxes)
                {
                    cells+=box.size();
                    for(int d=0;d<DIM;++d)
                    {low[d]=std::min(low[d],box.low()[d]);high[d]=std::max(high[d],box.high()[d]);}
                }
                const Proto::Box bounds(low,high);
                if(cells==bounds.size()) source.boxes={bounds};
            }
        }
    }

    const std::vector<SourceRegion>& sources() const {return m_sources;}

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
