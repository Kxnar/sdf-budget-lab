#include "sdf/oracle.hpp"
#include <fstream>
#include <sstream>

namespace sdf {
double unit(std::mt19937_64& rng) { return double(rng()>>11)*0x1.0p-53; }
namespace {
constexpr double pi=3.14159265358979323846;
Vec3 direction(std::mt19937_64& rng) {
    const double z=2*unit(rng)-1, a=2*pi*unit(rng), r=std::sqrt(1-z*z);
    return {r*std::cos(a),r*std::sin(a),z};
}
double box_distance(Vec3 p,Vec3 h) {
    const Vec3 q{std::abs(p.x)-h.x,std::abs(p.y)-h.y,std::abs(p.z)-h.z};
    return norm({std::max(q.x,0.0),std::max(q.y,0.0),std::max(q.z,0.0)})+std::min(std::max({q.x,q.y,q.z}),0.0);
}
class Analytic final:public Oracle {
    std::string shape_; double feature_;
public:
    Analytic(std::string shape,double feature):shape_(std::move(shape)),feature_(feature) {
        if(shape_!="sphere" && shape_!="torus" && shape_!="thin_box" && shape_!="double_sphere") throw std::invalid_argument("Unknown analytic shape");
        if(!std::isfinite(feature) || feature<0.005 || feature>0.25) throw std::invalid_argument("Feature must be in [0.005,0.25]");
    }
    double distance(Vec3 p) const override {
        if(shape_=="sphere") return norm(p)-0.55;
        if(shape_=="torus") return std::hypot(std::hypot(p.x,p.y)-0.55,p.z)-feature_;
        if(shape_=="thin_box") return box_distance(p,{0.6,0.6,feature_});
        const double c=0.3+feature_/2;
        return std::min(norm(p-Vec3{c,0,0}),norm(p+Vec3{c,0,0}))-0.3;
    }
    Vec3 surface(std::mt19937_64& rng) const override {
        if(shape_=="sphere") return direction(rng)*0.55;
        if(shape_=="double_sphere") {
            const double c=(unit(rng)<0.5?-1:1)*(0.3+feature_/2);
            return direction(rng)*0.3+Vec3{c,0,0};
        }
        if(shape_=="torus") {
            // Area-uniform rejection sampling, unlike the upstream angular sampler.
            const double a=2*pi*unit(rng); double b;
            do { b=2*pi*unit(rng); } while(unit(rng)*(0.55+feature_)>0.55+feature_*std::cos(b));
            const double r=0.55+feature_*std::cos(b);
            return {r*std::cos(a),r*std::sin(a),feature_*std::sin(b)};
        }
        Vec3 p{(2*unit(rng)-1)*0.6,(2*unit(rng)-1)*0.6,(2*unit(rng)-1)*feature_};
        const double a=unit(rng)*(0.36+1.2*feature_), sign=unit(rng)<0.5?-1:1;
        if(a<0.6*feature_) p.x=sign*0.6;
        else if(a<1.2*feature_) p.y=sign*0.6;
        else p.z=sign*feature_;
        return p;
    }
    std::string name() const override { return "analytic:"+shape_; }
};
}
std::unique_ptr<Oracle> analytic_oracle(const std::string& shape,double feature) { return std::make_unique<Analytic>(shape,feature); }
#ifndef SDF_WITH_CGAL
std::unique_ptr<Oracle> mesh_oracle(const std::string&) { throw std::runtime_error("Rebuild with SDF_WITH_CGAL=ON for meshes"); }
std::unique_ptr<Oracle> mesh_oracle_fast(const std::string&) { throw std::runtime_error("Rebuild with SDF_WITH_CGAL=ON for meshes"); }
#endif
Samples sample_queries(const Oracle& oracle,std::size_t n,std::uint64_t seed,double fraction,double noise,Box domain) {
    if(n==0 || fraction<0 || fraction>1 || !std::isfinite(fraction) || noise<=0 || !std::isfinite(noise)) throw std::invalid_argument("Invalid sampling settings");
    domain.validate();
    std::mt19937_64 rng(seed); std::normal_distribution<double> gaussian(0,noise);
    Samples result; result.reserve(n);
    const auto near_count=static_cast<std::size_t>(double(n)*fraction);
    for(std::size_t i=0;i<n;++i) {
        Vec3 p;
        if(i<near_count) {
            p=oracle.surface(rng)+Vec3{gaussian(rng),gaussian(rng),gaussian(rng)};
            p={std::clamp(p.x,domain.lo.x,domain.hi.x),std::clamp(p.y,domain.lo.y,domain.hi.y),std::clamp(p.z,domain.lo.z,domain.hi.z)};
        } else p=domain.point({unit(rng),unit(rng),unit(rng)});
        result.push_back({p,oracle.distance(p)});
    }
    std::shuffle(result.begin(),result.end(),rng);
    return result;
}
Samples feature_probes(const Oracle& oracle,const std::string& shape,double feature) {
    Samples result;
    auto add=[&](Vec3 p) { result.push_back({p,oracle.distance(p)}); };
    if(shape=="thin_box") {
        for(int i=0;i<16;++i) for(int j=0;j<16;++j) {
            const double x=-0.5+i/15.0,y=-0.5+j/15.0;
            add({x,y,0}); add({x,y,1.5*feature}); add({x,y,-1.5*feature});
        }
    } else if(shape=="torus") {
        for(int i=0;i<256;++i) {
            const double a=2*pi*i/256.0,c=std::cos(a),s=std::sin(a);
            add({0.55*c,0.55*s,0}); add({(0.55+1.5*feature)*c,(0.55+1.5*feature)*s,0});
        }
    } else if(shape=="double_sphere") {
        add({-0.3-feature/2,0,0}); add({0.3+feature/2,0,0});
        for(int i=0;i<128;++i) add({feature*(-0.4+0.8*i/127.0),0,0});
    } else { add({0,0,0}); add({0.7,0,0}); }
    return result;
}
Samples read_probes(const Oracle& oracle,const std::string& path,Box domain) {
    std::ifstream file(path); if(!file) throw std::runtime_error("Cannot open probes: "+path);
    Samples samples; std::string line; std::size_t line_no=0;
    while(std::getline(file,line)) {
        ++line_no; if(line.empty() || line[0]=='#' || line=="x,y,z") continue;
        std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line); Vec3 p; std::string extra;
        if(!(row>>p.x>>p.y>>p.z) || (row>>extra) || !domain.contains(p)) throw std::runtime_error("Invalid probe row "+std::to_string(line_no));
        samples.push_back({p,oracle.distance(p)});
    }
    if(samples.empty()) throw std::invalid_argument("Probe file is empty");
    return samples;
}
} // namespace sdf
