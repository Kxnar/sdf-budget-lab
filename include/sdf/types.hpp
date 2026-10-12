#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace sdf {
struct Vec3 {
    double x{}, y{}, z{};
    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
};
inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
inline Vec3 operator*(Vec3 a, double s) { return {a.x*s,a.y*s,a.z*s}; }
inline double dot(Vec3 a, Vec3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
inline double norm(Vec3 a) { return std::sqrt(dot(a,a)); }
struct Box {
    Vec3 lo{-1,-1,-1}, hi{1,1,1};
    bool contains(Vec3 p) const {
        return p.x>=lo.x && p.y>=lo.y && p.z>=lo.z && p.x<=hi.x && p.y<=hi.y && p.z<=hi.z;
    }
    Vec3 point(Vec3 t) const { return {lo.x+(hi.x-lo.x)*t.x,lo.y+(hi.y-lo.y)*t.y,lo.z+(hi.z-lo.z)*t.z}; }
    Vec3 corner(unsigned i) const { return point({double(i&1),double((i>>1)&1),double((i>>2)&1)}); }
    unsigned octant(Vec3 p) const {
        const auto mid=(lo+hi)*0.5;
        return (p.x>=mid.x?1u:0u)|(p.y>=mid.y?2u:0u)|(p.z>=mid.z?4u:0u);
    }
    Box child(unsigned i) const {
        const auto mid=(lo+hi)*0.5;
        return {{i&1?mid.x:lo.x,i&2?mid.y:lo.y,i&4?mid.z:lo.z},
                {i&1?hi.x:mid.x,i&2?hi.y:mid.y,i&4?hi.z:mid.z}};
    }
    double scale() const { return std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z}); }
    void validate() const {
        if (!(hi.x>lo.x && hi.y>lo.y && hi.z>lo.z) || !std::isfinite(norm(hi-lo)))
            throw std::invalid_argument("Invalid domain");
    }
};
using Corners=std::array<float,8>;
inline double interpolate(const Corners& values, const Box& b, Vec3 p) {
    const double u=(p.x-b.lo.x)/(b.hi.x-b.lo.x), v=(p.y-b.lo.y)/(b.hi.y-b.lo.y), w=(p.z-b.lo.z)/(b.hi.z-b.lo.z);
    double value=0;
    for(unsigned i=0;i<8;++i)
        value+=values[i]*(i&1?u:1-u)*(i&2?v:1-v)*(i&4?w:1-w);
    return value;
}
struct Sample { Vec3 point; double distance; };
using Samples=std::vector<Sample>;
} // namespace sdf
