#pragma once
#include "sdf/oracle.hpp"
#include <fstream>
#include <iomanip>

namespace sdf {
// Proper orthogonal rotation plus translation; no distance-changing deformation.
struct Similarity {
    std::array<double,9> r{1,0,0,0,1,0,0,0,1}; Vec3 t{};
    Vec3 forward(Vec3 p) const { return {r[0]*p.x+r[1]*p.y+r[2]*p.z+t.x,r[3]*p.x+r[4]*p.y+r[5]*p.z+t.y,r[6]*p.x+r[7]*p.y+r[8]*p.z+t.z}; }
    Vec3 inverse(Vec3 p) const { p=p-t; return {r[0]*p.x+r[3]*p.y+r[6]*p.z,r[1]*p.x+r[4]*p.y+r[7]*p.z,r[2]*p.x+r[5]*p.y+r[8]*p.z}; }
};
struct GeometrySpec {
    std::string id, shape, mesh; double a{},b{},c{},noise{}; Similarity transform;
    static GeometrySpec read(const std::string& path) {
        std::ifstream f(path); GeometrySpec s; f>>s.id>>s.shape>>s.a>>s.b>>s.c>>s.noise;
        for(auto& x:s.transform.r) f>>x;
        f>>s.transform.t.x>>s.transform.t.y>>s.transform.t.z>>std::quoted(s.mesh);
        if(!f || !std::isfinite(s.a+s.b+s.c+s.noise+norm(s.transform.t)) || s.noise<=0 || s.a<=0 || s.a>.25 || s.b<0 || s.c<0) throw std::invalid_argument("Invalid geometry spec");
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) { double v=0; for(int k=0;k<3;++k) v+=s.transform.r[3*k+i]*s.transform.r[3*k+j]; if(std::abs(v-(i==j?1:0))>1e-12) throw std::invalid_argument("Nonorthogonal transform"); }
        const auto& r=s.transform.r; double det=r[0]*(r[4]*r[8]-r[5]*r[7])-r[1]*(r[3]*r[8]-r[5]*r[6])+r[2]*(r[3]*r[7]-r[4]*r[6]);
        if(std::abs(det-1)>1e-12) throw std::invalid_argument("Transform is not a proper rotation");
        if((s.shape=="shell" && s.b<=s.a) || ((s.shape=="capsule"||s.shape=="cylinder"||s.shape=="rounded_box") && s.b<=0) || (s.shape=="rounded_box"&&s.c<=0)) throw std::invalid_argument("Invalid novel shape parameters");
        return s;
    }
};
class GeneralizationOracle final:public Oracle {
    GeometrySpec s_; std::unique_ptr<Oracle> base_;
    static Vec3 direction(std::mt19937_64& r) { double z=2*unit(r)-1,a=6.283185307179586*unit(r),h=std::sqrt(1-z*z); return {h*std::cos(a),h*std::sin(a),z}; }
    static double box(Vec3 p,Vec3 h) { Vec3 q{std::abs(p.x)-h.x,std::abs(p.y)-h.y,std::abs(p.z)-h.z}; return norm({std::max(q.x,0.),std::max(q.y,0.),std::max(q.z,0.)})+std::min(std::max({q.x,q.y,q.z}),0.); }
public:
    explicit GeneralizationOracle(GeometrySpec s):s_(std::move(s)) {
        if(s_.shape=="mesh") base_=s_.b==1?mesh_oracle_fast(s_.mesh):mesh_oracle(s_.mesh);
        else if(s_.shape=="torus"||s_.shape=="thin_box"||s_.shape=="double_sphere") base_=analytic_oracle(s_.shape,s_.a);
        else if(s_.shape!="capsule"&&s_.shape!="cylinder"&&s_.shape!="rounded_box"&&s_.shape!="shell") throw std::invalid_argument("Unknown generalization shape");
    }
    double distance(Vec3 p) const override {
        if(s_.shape=="mesh") return base_->distance(p); // normalized mesh already transformed
        p=s_.transform.inverse(p);
        if(base_) return base_->distance(p);
        if(s_.shape=="capsule") return norm(p-Vec3{0,0,std::clamp(p.z,-s_.b,s_.b)})-s_.a;
        if(s_.shape=="cylinder") { double x=std::hypot(p.x,p.y)-s_.a,y=std::abs(p.z)-s_.b; return std::min(std::max(x,y),0.)+std::hypot(std::max(x,0.),std::max(y,0.)); }
        if(s_.shape=="rounded_box") return box(p,{s_.b,s_.c,.12})-s_.a;
        return std::abs(norm(p)-s_.b)-s_.a;
    }
    Vec3 surface(std::mt19937_64& r) const override {
        if(s_.shape=="mesh") return base_->surface(r);
        Vec3 p;
        if(base_) p=base_->surface(r);
        else if(s_.shape=="shell") { double o=s_.b+s_.a,i=s_.b-s_.a; p=direction(r)*(unit(r)<o*o/(o*o+i*i)?o:i); }
        else if(s_.shape=="capsule") { double ang=6.283185307179586*unit(r); if(unit(r)<s_.b/(s_.b+s_.a)) p={s_.a*std::cos(ang),s_.a*std::sin(ang),(2*unit(r)-1)*s_.b}; else { p=direction(r)*s_.a;p.z+=(p.z>=0?s_.b:-s_.b); } }
        else if(s_.shape=="cylinder") { double ang=6.283185307179586*unit(r); if(unit(r)<2*s_.b/(2*s_.b+s_.a)) p={s_.a*std::cos(ang),s_.a*std::sin(ang),(2*unit(r)-1)*s_.b}; else { double h=s_.a*std::sqrt(unit(r));p={h*std::cos(ang),h*std::sin(ang),(unit(r)<.5?-1:1)*s_.b}; } }
        else {
            // Area-uniform patches: six rectangles, twelve quarter cylinders,
            // eight sphere octants. Patch weights sum to the exact area.
            Vec3 h{s_.b,s_.c,.12}; double weights[3]={8*h.y*h.z,8*h.x*h.z,8*h.x*h.y};
            double face=weights[0]+weights[1]+weights[2],edge=12.566370614359172*s_.a*(h.x+h.y+h.z),sphere=12.566370614359172*s_.a*s_.a,q=unit(r)*(face+edge+sphere);
            if(q<face) { int axis=q<weights[0]?0:q<weights[0]+weights[1]?1:2; double v[3]={(2*unit(r)-1)*h.x,(2*unit(r)-1)*h.y,(2*unit(r)-1)*h.z};v[axis]=(unit(r)<.5?-1:1)*(h[axis]+s_.a);p={v[0],v[1],v[2]}; }
            else if(q<face+edge) { q=unit(r)*(h.x+h.y+h.z);int axis=q<h.x?0:q<h.x+h.y?1:2;double v[3]={0,0,0};v[axis]=(2*unit(r)-1)*h[axis];double angle=1.5707963267948966*unit(r);int u=(axis+1)%3,w=(axis+2)%3;v[u]=(unit(r)<.5?-1:1)*(h[u]+s_.a*std::cos(angle));v[w]=(unit(r)<.5?-1:1)*(h[w]+s_.a*std::sin(angle));p={v[0],v[1],v[2]}; }
            else { p=direction(r);p={p.x*s_.a+(p.x>=0?h.x:-h.x),p.y*s_.a+(p.y>=0?h.y:-h.y),p.z*s_.a+(p.z>=0?h.z:-h.z)}; }
        }
        return s_.transform.forward(p);
    }
    std::string name() const override { return s_.id; }
};
} // namespace sdf
