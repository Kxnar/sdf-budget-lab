#include "sdf/ray.hpp"
#include <iostream>
using namespace sdf;
void require(bool value,const char* text){if(!value)throw std::runtime_error(text);}
RaySegment segment(double lo,double hi,const std::function<double(Vec3)>& field){RaySegment s{lo,hi,{{0,0,0},{1,1,1}},{}};for(unsigned c=0;c<8;++c)s.corners[c]=field(s.box.corner(c));return s;}
int main(){try{MeasurementRay diagonal{{0,0,0},{1,1,1},0,1};
auto three=segmented_transitions({segment(0,1,[](Vec3 p){return (p.x-.25)*(p.y-.5)*(p.z-.75);})},diagonal);require(three.crossings.size()==3,"Cubic three crossings");
auto narrow=segmented_transitions({segment(0,1,[](Vec3 p){return (p.x-(.5-1e-6))*(p.y-(.5+1e-6));})},diagonal);require(narrow.crossings.size()==2&&narrow.crossings[1]-narrow.crossings[0]<3e-6,"Narrow double crossing");
auto tangent=segmented_transitions({segment(0,1,[](Vec3 p){return (p.x-.5)*(p.y-.5);})},diagonal);require(tangent.crossings.empty()&&tangent.tangencies.size()==1,"Tangent separately recorded");
auto zero=segmented_transitions({segment(0,1,[](Vec3){return 0.;})},diagonal);require(zero.ambiguous&&zero.zero_segments==1,"Zero segment ambiguous");
auto jump=segmented_transitions({segment(0,.5,[](Vec3){return 1.;}),segment(.5,1,[](Vec3){return -1.;})},diagonal);require(jump.crossings==std::vector<double>{.5}&&jump.jumps.size()==1,"Interface signed jump");
auto near=segmented_transitions({segment(0,.5,[](Vec3 p){return p.x-(.5-1e-9);}),segment(.5,1,[](Vec3 p){return p.x-(.5-1e-9);})},diagonal);require(near.crossings.size()==1&&near.crossings[0]<.5,"Near-interface root no duplicate");
auto boundary=segmented_transitions({segment(0,.5,[](Vec3 p){return p.x-.5;}),segment(.5,1,[](Vec3 p){return p.x-.5;})},diagonal);require(boundary.crossings==std::vector<double>{.5},"Owned boundary crossing exactly once");
for(double minimum:{1e-8,-1e-8}){auto r=segmented_transitions({segment(0,1,[&](Vec3 p){return (p.x-.5)*(p.y-.5)+minimum;})},diagonal);require(r.crossings.size()==(minimum>0?0:2)&&!r.ambiguous,"Resolved near tangent minimum");}
for(double minimum:{1e-16,-1e-16}){auto r=segmented_transitions({segment(0,1,[&](Vec3 p){return (p.x-.5)*(p.y-.5)+minimum;})},diagonal);require(r.ambiguous,"Near tangent numeric uncertainty explicit");}
for(int exponent:{-20,0,20})for(double direction:{-1.,1.}){
    auto perturbed=segment(0,1,[](Vec3 p){return (p.x-.5)*(p.y-.5);});
    for(auto& corner:perturbed.corners)corner=std::nextafter(std::ldexp(corner,exponent),std::copysign(std::numeric_limits<double>::infinity(),direction));
    auto r=segmented_transitions({perturbed},diagonal);
    require(r.ambiguous,"One-ULP binary64 corner perturbation remains ambiguous across amplitudes");
}
auto repeated=segmented_transitions({segment(0,1,[](Vec3 p){return (p.x-.5)*(p.y-(.5+1e-8))*(p.z-.8);})},diagonal);require(repeated.ambiguous,"Nearly repeated cubic uncertainty");
double next=std::nextafter(.5,1.);auto tiny=segmented_transitions({segment(.5,next,[](Vec3 p){return p.x+1.;})},diagonal);require(tiny.ambiguous,"Unrepresentable segment midpoint");
auto nearplanes=grid_ray_boundaries({{0,0,0},{1,std::nextafter(1.,2.),0},0,1},4);require(std::find(nearplanes.begin(),nearplanes.end(),.5)!=nearplanes.end()&&std::find(nearplanes.begin(),nearplanes.end(),.5/std::nextafter(1.,2.))!=nearplanes.end(),"Distinct near-coincident planes retained");
std::cout<<"Segmented-ray regression suite passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
