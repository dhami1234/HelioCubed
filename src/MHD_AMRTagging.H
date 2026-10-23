#pragma once

#include "Proto.H"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace HelioCubed
{
// Physical tag bounds are independent of panel coordinates and AMR level.
// They select cell centers; buffering, whole patches and nesting may extend
// the allocated refinement beyond these bounds.
struct AMRRefinementBounds
{
    struct Range
    {
        bool enabled=false;
        double low=0,high=0;
        bool contains(double value) const {return !enabled || (value>=low && value<=high);}
    };
    Range radius,colatitude,longitude;
    bool fullLongitude=false;

    bool active() const {return radius.enabled || colatitude.enabled || longitude.enabled;}
    bool angular() const {return colatitude.enabled || (longitude.enabled && !fullLongitude);}

    static Range parseRange(std::string text,const char* name)
    {
        if(text.empty()) return {};
        text=text.substr(0,text.find('#'));
        std::replace(text.begin(),text.end(),',',' ');
        std::istringstream input(text);std::string low,high,extra;
        const std::string message=std::string(name)+" requires two finite numbers: min,max";
        if(!(input>>low>>high) || input>>extra) throw std::invalid_argument(message);
        auto number=[&](const std::string& word)
        {
            size_t end=0;double value;
            try {value=std::stod(word,&end);} catch(...) {throw std::invalid_argument(message);}
            if(end!=word.size() || !std::isfinite(value)) throw std::invalid_argument(message);
            return value;
        };
        return {true,number(low),number(high)};
    }

    static double normalizeLongitude(double degrees)
    {
        const double value=std::fmod(degrees,360.0);
        return value<0 ? value+360 : value;
    }

    static AMRRefinementBounds parse(const std::string& r,const std::string& theta,const std::string& phi,
        const std::string& suffix="")
    {
        AMRRefinementBounds result;
        const auto rName="amr_refine_r_range"+suffix;
        const auto thetaName="amr_refine_colat_range"+suffix;
        const auto phiName="amr_refine_lon_range"+suffix;
        result.radius=parseRange(r,rName.c_str());
        result.colatitude=parseRange(theta,thetaName.c_str());
        result.longitude=parseRange(phi,phiName.c_str());
        if(result.radius.enabled && (result.radius.low<0 || result.radius.high<result.radius.low))
            throw std::invalid_argument(rName+" requires 0 <= min <= max (cm)");
        if(result.colatitude.enabled && (result.colatitude.low<0 || result.colatitude.high>180
            || result.colatitude.high<result.colatitude.low))
            throw std::invalid_argument(thetaName+" requires 0 <= min <= max <= 180 (degrees)");
        if(result.longitude.enabled)
        {
            result.fullLongitude=std::abs(result.longitude.high-result.longitude.low)>=360;
            result.longitude.low=normalizeLongitude(result.longitude.low);
            result.longitude.high=normalizeLongitude(result.longitude.high);
        }
        return result;
    }

    bool containsDirection(double x,double y,double z) const
    {
        // atan2 is stable at the poles; longitude is immaterial at either pole.
        const double transverse=std::hypot(x,y);
        const double theta=std::atan2(transverse,z)*180/M_PI;
        constexpr double tolerance=64*std::numeric_limits<double>::epsilon()*360;
        if(colatitude.enabled && (theta<colatitude.low-tolerance || theta>colatitude.high+tolerance)) return false;
        if(!longitude.enabled || fullLongitude || transverse==0) return true;
        const double phi=normalizeLongitude(std::atan2(y,x)*180/M_PI);
        // A directed arc also handles negative longitudes and ranges across 0.
        const double width=normalizeLongitude(longitude.high-longitude.low);
        const double offset=normalizeLongitude(phi-longitude.low);
        return offset<=width+tolerance || offset>=360-tolerance;
    }
};

struct AMRRefinementRegions
{
    std::vector<AMRRefinementBounds> regions;
    bool active() const {return !regions.empty();}
};

// Cache one radial interval and one angular plane per local parent patch.
// Regrids then filter tags without trigonometry or a full 3D mask per cell.
class AMRRefinementRegionMask
{
    Proto::Box m_box;
    int m_low,m_high;
    bool m_empty=false;
    std::vector<unsigned char> m_angular;
public:
    AMRRefinementRegionMask(const Proto::Box& box,Proto::Point sizes,unsigned block,
        const AMRRefinementBounds& bounds,const Proto::CubedSphereShellGeometry& geometry)
        :m_box(box),m_low(box.low()[0]),m_high(box.high()[0])
    {
        if(bounds.radius.enabled)
        {
            m_low=box.high()[0]+1;m_high=box.low()[0]-1;
            for(int r=box.low()[0];r<=box.high()[0];++r)
            {
                Proto::Array<double,DIM> xi,position;double jacobian;
                xi[0]=(r+0.5)/sizes[0];xi[1]=xi[2]=0.5;
                Proto::f_cubedSphereMap(position,jacobian,xi,block,geometry);
                if(bounds.radius.contains(std::hypot(position[0],position[1],position[2])))
                {m_low=std::min(m_low,r);m_high=r;}
            }
            m_empty=m_low>m_high;
        }
        if(!m_empty && bounds.angular())
        {
            const int width=box.sizes()[1];
            m_angular.resize(size_t(width)*box.sizes()[2]);
            // Angles are independent of radius. Unit directions avoid losing
            // angular precision through the physical shell's length scale.
            Proto::CubedSphereShellGeometry unit;unit.r0=1;unit.r1=2;unit.radialStretch=0;
            bool any=false;
            for(int k=box.low()[2];k<=box.high()[2];++k)
                for(int j=box.low()[1];j<=box.high()[1];++j)
                {
                    Proto::Array<double,DIM> xi,position;double jacobian;
                    xi[0]=0;xi[1]=(j+0.5)/sizes[1];xi[2]=(k+0.5)/sizes[2];
                    Proto::f_cubedSphereMap(position,jacobian,xi,block,unit);
                    const bool allowed=bounds.containsDirection(position[0],position[1],position[2]);
                    m_angular[size_t(k-box.low()[2])*width+j-box.low()[1]]=allowed;
                    any|=allowed;
                }
            m_empty=!any;
        }
    }

    bool empty() const {return m_empty;}
    bool contains(Proto::Point cell) const
    {
        if(m_empty || cell[0]<m_low || cell[0]>m_high) return false;
        return m_angular.empty() || m_angular[size_t(cell[2]-m_box.low()[2])*m_box.sizes()[1]+cell[1]-m_box.low()[1]];
    }

    template<class Tag> void tagRegions(const Tag& tag) const
    {
        if(m_empty) return;
        if(m_angular.empty())
        {
            auto low=m_box.low(),high=m_box.high();low[0]=m_low;high[0]=m_high;
            tag(Proto::Box(low,high));return;
        }
        for(int k=m_box.low()[2];k<=m_box.high()[2];++k)
        {
            int j=m_box.low()[1];
            while(j<=m_box.high()[1])
            {
                while(j<=m_box.high()[1] && !contains(Proto::Point(m_low,j,k))) ++j;
                const int start=j;
                while(j<=m_box.high()[1] && contains(Proto::Point(m_low,j,k))) ++j;
                if(j>start) tag(Proto::Box(Proto::Point(m_low,start,k),Proto::Point(m_high,j-1,k)));
            }
        }
    }
};

// Keep each region's radial and angular bounds together. Combining radial and
// angular unions separately would incorrectly include space between regions.
class AMRRefinementMask
{
    std::vector<AMRRefinementRegionMask> m_regions;
public:
    AMRRefinementMask(const Proto::Box& box,Proto::Point sizes,unsigned block,
        const AMRRefinementBounds& bounds,const Proto::CubedSphereShellGeometry& geometry)
    {
        AMRRefinementRegionMask mask(box,sizes,block,bounds,geometry);
        if(!mask.empty()) m_regions.push_back(std::move(mask));
    }

    AMRRefinementMask(const Proto::Box& box,Proto::Point sizes,unsigned block,
        const AMRRefinementRegions& bounds,const Proto::CubedSphereShellGeometry& geometry)
    {
        m_regions.reserve(bounds.regions.size());
        if(!bounds.active()) m_regions.emplace_back(box,sizes,block,AMRRefinementBounds{},geometry);
        for(const auto& region : bounds.regions)
        {
            AMRRefinementRegionMask mask(box,sizes,block,region,geometry);
            if(!mask.empty()) m_regions.push_back(std::move(mask));
        }
    }

    bool empty() const {return m_regions.empty();}
    bool contains(Proto::Point cell) const
    {
        for(const auto& region : m_regions) if(region.contains(cell)) return true;
        return false;
    }

    template<class Tag> void tagRegions(const Tag& tag) const
    {
        // The caller inserts patches into a set, so overlapping regions do
        // not duplicate fine patches or increase reported coverage.
        for(const auto& region : m_regions) region.tagRegions(tag);
    }
};

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
