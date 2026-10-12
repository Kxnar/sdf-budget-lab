#pragma once
#include "sdf/precision.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace sdf {
struct MeasurementRay { Vec3 origin,direction; double low{},high{}; };
struct RaySegment { double low{},high{}; Box box; std::array<double,8> corners{}; };
struct RayTransitions { std::vector<double> crossings,tangencies,jumps; std::size_t zero_segments{}; bool ambiguous{}; };
using Cubic=std::array<long double,4>;
inline Cubic ray_polynomial(const RaySegment& segment,const MeasurementRay& ray) {
    Cubic result{};const auto width=segment.box.hi-segment.box.lo;
    for(unsigned c=0;c<8;++c){Cubic term{segment.corners[c],0,0,0};
        for(unsigned axis=0;axis<3;++axis){long double a=(static_cast<long double>(ray.origin[axis])+static_cast<long double>(ray.direction[axis])*segment.low-segment.box.lo[axis])/width[axis];long double b=static_cast<long double>(ray.direction[axis])*(segment.high-segment.low)/width[axis];if(!(c&(1u<<axis))){a=1-a;b=-b;}Cubic next{};for(unsigned k=0;k<4;++k){next[k]+=term[k]*a;if(k<3)next[k+1]+=term[k]*b;}term=next;}
        for(unsigned k=0;k<4;++k)result[k]+=term[k];
    }return result;
}
inline long double polynomial_value(const Cubic& c,long double u){return ((c[3]*u+c[2])*u+c[1])*u+c[0];}
inline std::vector<long double> critical_points(const Cubic& c){
    std::vector<long double> result;const auto a=3*c[3],b=2*c[2],d=c[1];
    if(a==0){if(b!=0)result.push_back(-d/b);}else{auto disc=b*b-4*a*d;if(disc>=0){auto r=std::sqrt(disc);auto q=-.5L*(b+std::copysign(r,b));if(q!=0){result.push_back(q/a);if(disc>0)result.push_back(d/q);}else result.push_back(-b/(2*a));}}
    result.erase(std::remove_if(result.begin(),result.end(),[](long double v){return !(v>0&&v<1);}),result.end());std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());return result;
}
inline int polynomial_side(const Cubic& c,long double u,int direction){
    auto v=polynomial_value(c,u);if(v!=0)return v>0?1:-1;
    auto d=(3*c[3]*u+2*c[2])*u+c[1];if(d!=0)return (d*direction)>0?1:-1;
    auto second=3*c[3]*u+c[2];if(second!=0)return second>0?1:-1;
    if(c[3]!=0)return c[3]*direction>0?1:-1;return 0;
}
inline RayTransitions segmented_transitions(const std::vector<RaySegment>& segments,const MeasurementRay& ray){
    RayTransitions result;std::vector<Cubic> polynomials;polynomials.reserve(segments.size());
    for(const auto& segment:segments){auto c=ray_polynomial(segment,ray);polynomials.push_back(c);if(std::all_of(c.begin(),c.end(),[](long double x){return x==0;})){++result.zero_segments;result.ambiguous=true;continue;}
        long double scale=0;for(auto x:c)scale+=std::abs(x);for(auto x:segment.corners)scale=std::max(scale,static_cast<long double>(std::abs(x)));auto uncertainty=128*std::numeric_limits<long double>::epsilon()*scale;
        if((segment.low+segment.high)/2==segment.low||(segment.low+segment.high)/2==segment.high)result.ambiguous=true;
        auto critical=critical_points(c);for(auto u:critical)if(std::abs(polynomial_value(c,u))<=uncertainty)result.ambiguous=true;
        for(long double u:{0.L,1.L}){auto v=polynomial_value(c,u);if(v!=0&&std::abs(v)<=uncertainty)result.ambiguous=true;}
        std::vector<long double> knots{0};knots.insert(knots.end(),critical.begin(),critical.end());knots.push_back(1);
        for(auto u:critical)if(polynomial_value(c,u)==0){if(polynomial_side(c,u,-1)!=polynomial_side(c,u,1))result.crossings.push_back(segment.low+double(u)*(segment.high-segment.low));else result.tangencies.push_back(segment.low+double(u)*(segment.high-segment.low));}
        for(std::size_t i=1;i<knots.size();++i){auto lo=knots[i-1],hi=knots[i];auto a=polynomial_value(c,lo),b=polynomial_value(c,hi);if(a==0||b==0||(a<0)==(b<0))continue;for(unsigned j=0;j<100;++j){auto mid=(lo+hi)/2;if(mid==lo||mid==hi)break;auto value=polynomial_value(c,mid);if(value==0){lo=hi=mid;break;}if((value<0)==(a<0))lo=mid;else hi=mid;}auto u=(lo+hi)/2;result.crossings.push_back(segment.low+double(u)*(segment.high-segment.low));}
    }
    for(std::size_t i=1;i<segments.size();++i){if(segments[i-1].high!=segments[i].low){result.ambiguous=true;continue;}const auto left=polynomial_side(polynomials[i-1],1,-1),right=polynomial_side(polynomials[i],0,1);if(!left||!right){result.ambiguous=true;continue;}if(left!=right){auto t=segments[i].low;result.crossings.push_back(t);if(polynomial_value(polynomials[i-1],1)!=0||polynomial_value(polynomials[i],0)!=0)result.jumps.push_back(t);}else if(polynomial_value(polynomials[i-1],1)==0&&polynomial_value(polynomials[i],0)==0)result.tangencies.push_back(segments[i].low);}
    std::sort(result.crossings.begin(),result.crossings.end());std::sort(result.tangencies.begin(),result.tangencies.end());return result;
}
inline bool ray_box_interval(const MeasurementRay& ray,const Box& box,double& lo,double& hi){lo=ray.low;hi=ray.high;for(unsigned axis=0;axis<3;++axis){auto d=ray.direction[axis],o=ray.origin[axis];if(d==0){if(o<box.lo[axis]||o>box.hi[axis]||(o==box.hi[axis]&&box.hi[axis]!=1))return false;}else{auto a=(box.lo[axis]-o)/d,b=(box.hi[axis]-o)/d;lo=std::max(lo,std::min(a,b));hi=std::min(hi,std::max(a,b));}}return hi>lo;}
inline std::vector<RaySegment> packed_ray_segments(const unsigned char* data,const MeasurementRay& ray){std::vector<RaySegment> result;std::function<void(std::size_t,Box)> visit=[&](std::size_t offset,Box box){double lo,hi;if(!ray_box_interval(ray,box,lo,hi))return;std::uint32_t descriptor;std::memcpy(&descriptor,data+offset,4);if(!(descriptor&3)){for(unsigned c=0;c<8;++c)visit((descriptor&~3u)+4*c,box.child(c));return;}EncodedCorners encoded;encoded.encoding=Encoding(descriptor&3);std::memcpy(encoded.data.data(),data+(descriptor&~3u),encoded.payload_bytes());RaySegment segment{lo,hi,box,{}};auto values=encoded.decode();std::copy(values.begin(),values.end(),segment.corners.begin());result.push_back(segment);};visit(0,{});std::sort(result.begin(),result.end(),[](const auto& a,const auto& b){return a.low<b.low;});return result;}
inline std::vector<double> grid_ray_boundaries(const MeasurementRay& ray,std::size_t cells){std::vector<double> times{ray.low,ray.high};for(unsigned axis=0;axis<3;++axis)if(ray.direction[axis]!=0)for(std::size_t i=1;i<cells;++i){auto t=(-1+2.*i/cells-ray.origin[axis])/ray.direction[axis];if(t>ray.low&&t<ray.high)times.push_back(t);}std::sort(times.begin(),times.end());times.erase(std::unique(times.begin(),times.end()),times.end());return times;}
}
