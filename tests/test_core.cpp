#include "sdf/metrics.hpp"
#include <filesystem>
#include <iostream>
#include <limits>

namespace {
using namespace sdf;
int checks=0;
void check(bool condition,const std::string& message) { ++checks; if(!condition) throw std::runtime_error(message); }
void close(double a,double b,double tolerance,const std::string& message) { check(std::abs(a-b)<=tolerance,message); }
template<class F> void throws(F f,const std::string& message) {
    bool caught=false; try { f(); } catch(const std::exception&) { caught=true; } check(caught,message);
}
void fields() {
    const Distance plane=[](Vec3 p){ return 0.25*p.x-0.5*p.y+0.125*p.z+0.375; };
    auto sphere=analytic_oracle("sphere",0.04);
    const auto samples=sample_queries(*sphere,1024,12,0,0.02);
    for(const auto budget:{512u,1024u,8192u,32768u}) {
        UniformGrid grid(plane,budget); AdaptiveField adaptive(plane,samples,budget);
        check(grid.bytes()<=budget,"Grid exceeded budget"); check(adaptive.bytes()<=budget,"Adaptive exceeded budget");
        check(adaptive.node_count()%8==1,"Broken octree node count");
        check(adaptive.leaves().size()==adaptive.leaf_count(),"Broken leaf count");
        for(const auto& s:samples) { close(grid.query(s.point),plane(s.point),1e-7,"Grid affine interpolation"); close(adaptive.query(s.point),plane(s.point),1e-7,"Octree affine interpolation"); }
        for(unsigned i=0;i<8;++i) { close(grid.query(Box{}.corner(i)),plane(Box{}.corner(i)),1e-7,"Grid upper/lower boundary"); close(adaptive.query(Box{}.corner(i)),plane(Box{}.corner(i)),1e-7,"Adaptive boundary"); }
        throws([&]{grid.query({1.01,0,0});},"Grid should reject out-of-domain query");
        throws([&]{adaptive.query({std::numeric_limits<double>::quiet_NaN(),0,0});},"Adaptive should reject NaN");
    }
    throws([&]{UniformGrid field(plane,1);},"Small grid budget");
    throws([&]{AdaptiveField field(plane,samples,1);},"Small octree budget");
    throws([&]{AdaptiveField field(plane,{},1024);},"Empty training pool");
    throws([&]{UniformGrid field(plane,1024,{{0,0,0},{0,1,1}});},"Invalid domain");
    throws([&]{UniformGrid field([](Vec3){return std::numeric_limits<double>::infinity();},1024);},"Nonfinite oracle");
    Box domain{{-2,-3,-4},{4,5,6}}; UniformGrid grid(plane,8192,domain);
    close(grid.query({4,5,6}),plane({4,5,6}),1e-6,"Rectangular domain upper boundary");
    AdaptiveField a([&](Vec3 p){return sphere->distance(p);},samples,8192);
    AdaptiveField b([&](Vec3 p){return sphere->distance(p);},samples,8192);
    for(const auto& s:samples) close(a.query(s.point),b.query(s.point),0,"Deterministic allocation");
}
void oracles() {
    auto sphere=analytic_oracle("sphere",0.04); auto torus=analytic_oracle("torus",0.04);
    auto plate=analytic_oracle("thin_box",0.04); auto gap=analytic_oracle("double_sphere",0.04);
    close(sphere->distance({0,0,0}),-0.55,1e-12,"Sphere SDF");
    close(torus->distance({0.55,0,0}),-0.04,1e-12,"Torus SDF");
    close(plate->distance({0,0,0}),-0.04,1e-12,"Plate SDF");
    close(gap->distance({0,0,0}),0.02,1e-12,"Gap SDF");
    std::mt19937_64 rng(1);
    for(const auto* oracle:{sphere.get(),torus.get(),plate.get(),gap.get()}) {
        for(int i=0;i<100;++i) close(oracle->distance(oracle->surface(rng)),0,1e-12,"Surface sampling");
        const auto probes=feature_probes(*oracle,oracle==torus.get()?"torus":oracle==plate.get()?"thin_box":oracle==gap.get()?"double_sphere":"sphere",0.04);
        const auto metric=evaluate([&](Vec3 p){return oracle->distance(p);},probes);
        close(metric.inside_recall,1,0,"Inside probes"); close(metric.outside_recall,1,0,"Outside probes");
    }
    const auto a=sample_queries(*torus,100,2,0.5,0.02),b=sample_queries(*torus,100,2,0.5,0.02);
    for(std::size_t i=0;i<a.size();++i) { close(norm(a[i].point-b[i].point),0,0,"Reproducible sampling"); check(Box{}.contains(a[i].point),"Sample in domain"); }
    throws([&]{sample_queries(*sphere,0,0,0.5,0.02);},"Zero samples");
    throws([&]{analytic_oracle("torus",0);},"Invalid feature");
}
void metrics() {
    Samples s; for(int i=1;i<=100;++i) s.push_back({{double(i),0,0},0});
    const auto m=evaluate([](Vec3 p){return p.x;},s);
    close(m.mae,50.5,1e-12,"MAE definition"); close(m.p99,99,0,"p99 nearest rank"); close(m.cvar99,100,0,"CVaR upper 1%");
    const auto signs=evaluate([](Vec3 p){return p.x;},{{{-1,0,0},-1},{{0,0,0},1},{{0,0,0},0}});
    check(signs.sign_count==2,"Exclude reference boundary"); close(signs.sign_error,0.5,0,"Zero prediction is a wrong sign");
    close(signs.inside_recall,1,0,"Inside recall"); close(signs.outside_recall,0,0,"Outside recall");
    throws([&]{evaluate([](Vec3){return 0.;},{});},"Empty metric input");
    const auto timing=time_queries([](Vec3 p){return p.x;},s,3);
    check(timing.median_ns>=0 && std::isfinite(timing.checksum),"Valid timing/checksum");
}
void meshes(const std::filesystem::path& fixtures) {
#ifdef SDF_WITH_CGAL
    auto cube=mesh_oracle((fixtures/"cube.off").string());
    close(cube->distance({0,0,0}),-0.5,1e-12,"Mesh interior sign/distance");
    close(cube->distance({0.7,0,0}),0.2,1e-12,"Mesh exterior distance");
    close(cube->distance({0.5,0,0}),0,1e-12,"Mesh boundary");
    close(cube->distance({0.7,0.7,0.7}),std::sqrt(0.12),1e-12,"Mesh vertex distance");
    std::mt19937_64 rng(5); for(int i=0;i<40;++i) close(cube->distance(cube->surface(rng)),0,1e-12,"Mesh area sampling");
    for(const auto* name:{"open.off","degenerate.off","intersecting.off","outside.off"})
        throws([&]{mesh_oracle((fixtures/name).string());},std::string("Reject invalid mesh ")+name);
    auto wall=mesh_oracle((fixtures/"thin_box.off").string());
    auto exact=analytic_oracle("thin_box",0.04);
    const auto samples=sample_queries(*exact,100,4,0.5,0.02);
    for(const auto& s:samples) close(wall->distance(s.point),s.distance,1e-12,"CGAL vs analytic box");
#else
    (void)fixtures;
    throws([]{mesh_oracle("unused");},"CGAL-disabled diagnostic");
#endif
}
}
int main(int argc,char** argv) {
    try { fields(); oracles(); metrics(); meshes(argc>1?argv[1]:"fixtures"); std::cout<<checks<<" checks passed\n"; return 0; }
    catch(const std::exception& e) { std::cerr<<"Failed after "<<checks<<" checks: "<<e.what()<<'\n'; return 1; }
}
