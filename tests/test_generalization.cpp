#include "sdf/generalization.hpp"
#include "sdf/metrics.hpp"
#include "sdf/precision.hpp"
#include <iostream>
using namespace sdf;
void check(bool condition){if(!condition)throw std::runtime_error("Generalization oracle check failed");}
int main(){
    auto test=[](const std::string& s,double a,double b,double c,Vec3 p,double expected){GeometrySpec g;g.id=s;g.shape=s;g.a=a;g.b=b;g.c=c;g.noise=.02;g.transform.r={0,-1,0,1,0,0,0,0,1};g.transform.t={.03,-.02,.01};GeneralizationOracle o(g);check(std::abs(o.distance(g.transform.forward(p))-expected)<1e-12);std::mt19937_64 r(999);for(int i=0;i<300;++i)check(std::abs(o.distance(o.surface(r)))<1e-12);};
    test("capsule",.1,.4,0,{0,0,.6},.1);test("capsule",.1,.4,0,{0,0,0},-.1);
    test("cylinder",.2,.3,0,{.3,0,.4},std::sqrt(.02));test("cylinder",.2,.3,0,{0,0,0},-.2);
    test("rounded_box",.05,.3,.2,{.3,.2,.12},-.05);test("rounded_box",.05,.3,.2,{.35,.2,.12},0);
    test("shell",.04,.5,0,{0,0,0},.46);test("shell",.04,.5,0,{.5,0,0},-.04);
    GeometrySpec g;g.id="refinement-test";g.shape="torus";g.a=.035;g.noise=.02;GeneralizationOracle o(g);Distance q=[&](Vec3 p){return o.distance(p);};auto fit=sample_queries(o,1000,777,.5,.02);PrecisionTree tree(q,fit,5);auto six=tree.refine_depth_six(q,fit,8);auto before=tree.nodes.size();auto seven=tree.refine_depth_seven(q,fit,4);check(!six.empty()&&!seven.empty());check(tree.nodes.size()==before+8*seven.size());for(auto r:seven){check(tree.spatial.nodes[r.node].field.first_child>=before);check(std::abs(r.box.scale()-1./32)<1e-12);}PrecisionDP dp(tree,32768,Encoding::q8);PackedField field(tree,dp.allocate(32768));check(field.bytes()<=32768);for(auto s:fit)check(std::isfinite(field.query(s.point)));
    std::cout<<"Exact new SDF points, rigid-transform distance invariance and surface samplers passed\n";
}
