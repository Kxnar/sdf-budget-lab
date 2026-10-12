#include "sdf/precision.hpp"
#include "sdf/surface.hpp"
#include <iostream>
#include <random>
#include <limits>

using namespace sdf;
namespace {
int checks=0;
void check(bool ok,const char* message) { ++checks; if(!ok) throw std::runtime_error(message); }
template<class F> void rejects(F f) { bool caught=false; try { f(); } catch(const std::exception&) { caught=true; } check(caught,"Expected exception"); }
void exhaustive() {
    PrecisionTree t([](Vec3 p){return p.x;},{{{0,0,0},0}},1);
    std::mt19937_64 rng(81);
    for(unsigned trial=0;trial<12;++trial) {
        for(auto& node:t.nodes) for(auto& loss:node.loss) loss=double(rng()%37);
        for(auto mode:{Encoding::branch,Encoding::f32,Encoding::q16,Encoding::q8}) {
            const auto max=sizeof(PackedField)+300;
            PrecisionDP dp(t,max,mode);
            for(std::size_t cap=sizeof(PackedField)+36;cap<=max;++cap) {
                double best=std::numeric_limits<double>::infinity(); std::size_t best_bytes=0;
                auto consider=[&](double loss,std::size_t bytes) {
                    bytes+=sizeof(PackedField);
                    if(bytes<=cap && (loss<best || (loss==best && bytes<best_bytes))) { best=loss; best_bytes=bytes; }
                };
                for(unsigned e=0;e<3;++e) if(mode==Encoding::branch || unsigned(mode)==e+1) consider(t.nodes[0].loss[e],36-8*e);
                // Enumerate all 3^8 child encoding assignments independently of DP.
                for(unsigned code=0;code<6561;++code) {
                    unsigned value=code; double loss=0; std::size_t bytes=4; bool valid=true;
                    for(unsigned c=0;c<8;++c) {
                        unsigned e=value%3; value/=3;
                        valid&=mode==Encoding::branch || unsigned(mode)==e+1;
                        loss+=t.nodes[1+c].loss[e]; bytes+=36-8*e;
                    }
                    if(valid) consider(loss,bytes);
                }
                const auto a=dp.allocate(cap);
                check(a.loss==best && a.bytes==best_bytes,"Mixed-cost DP differs from exhaustive optimum");
                check(PackedField(t,a).bytes()==a.bytes,"Actual byte cost differs");
            }
        }
    }
}
void codecs() {
    for(auto e:{Encoding::f32,Encoding::q16,Encoding::q8}) {
        for(auto v:{Corners{{-1,-.25f,0,.15f,.3f,.51f,.87f,1}},Corners{{-.01f,-.01f,-.01f,-.01f,-.01f,-.01f,-.01f,-.01f}}}) {
            const auto encoded=encode_corners(v,e); const auto decoded=encoded.decode();
            const double tolerance=e==Encoding::f32?0:e==Encoding::q16?2./65535:2./255;
            for(unsigned i=0;i<8;++i) check(std::abs(v[i]-decoded[i])<=tolerance,"Quantization error exceeds one step");
        }
    }
    Corners invalid{}; invalid[1]=std::numeric_limits<float>::infinity();
    rejects([&]{encode_corners(invalid,Encoding::q8);});
    rejects([&]{encode_corners(Corners{},Encoding::branch);});
}
void real_tree() {
    auto oracle=analytic_oracle("torus",.047); auto fit=sample_queries(*oracle,1024,27,.5,.02);
    unsigned calls=0;
    PrecisionTree t([&](Vec3 p){++calls;return oracle->distance(p);},fit,3);
    check(calls==8*t.nodes.size(),"Unexpected precision oracle calls");
    PrecisionDP joint(t,8192);
    for(auto mode:{Encoding::f32,Encoding::q16,Encoding::q8,Encoding::branch}) {
        PrecisionDP dp(t,8192,mode);
        for(auto cap:{100u,512u,2048u,8192u}) {
            const auto a=dp.allocate(cap); PackedField f(t,a); double loss=0;
            for(const auto& s:fit) { double e=(f.query(s.point)-s.distance)/2; loss+=e*e; }
            check(std::abs(loss-a.loss)<1e-11,"Decoded fitting loss mismatch");
            check(f.bytes()<=cap,"Exceeded byte budget");
            check(joint.allocate(cap).loss<=a.loss+1e-12,"Joint worse than fixed encoding optimum");
            rejects([&]{f.query({2,0,0});});
        }
    }
    rejects([&]{joint.allocate(83);}); rejects([&]{joint.allocate(8193);});
    check(calls==8*t.nodes.size(),"Allocation called oracle");
    // Exact fit constant prefers one cheap leaf despite a large budget.
    PrecisionTree constant([](Vec3){return -1.;},{{{0,0,0},-1}},2);
    auto a=PrecisionDP(constant,4096).allocate(4096);
    check(a.bytes==sizeof(PackedField)+20 && a.branches==0 && a.loss==0,"Constant cell/tie handling");
    check(PackedField(constant,a).query({1,1,1})==-1,"Domain boundary query");
}
void full_fields() {
    const Distance oracle=[](Vec3 p){return std::hypot(std::hypot(p.x,p.y)-.55,p.z)-.035;};
    for(unsigned depth:{0u,2u,5u}) {
        PrecisionTree tree(oracle,{{{0,0,0},oracle({0,0,0})}},depth);
        for(auto e:{Encoding::f32,Encoding::q8}) {
            PrecisionAllocation a; a.choice.resize(tree.nodes.size(),e);
            for(std::size_t v=0;v<tree.nodes.size();++v) {
                if(tree.spatial.nodes[v].field.first_child) { a.choice[v]=Encoding::branch; ++a.branches; }
                else ++a.leaves[unsigned(e)-1];
            }
            a.bytes=sizeof(PackedField)+4*a.branches+(e==Encoding::f32?36:20)*a.leaves[unsigned(e)-1];
            PackedField original(tree,a); std::size_t calls=0;
            PackedField direct(oracle,depth,e,calls);
            check(direct.bytes()==original.bytes(),"Full builder byte mismatch");
            check(std::equal(direct.data(),direct.data()+direct.payload_bytes(),original.data()),"Full builder differs from candidate packing");
            check(calls==8*a.leaves[unsigned(e)-1],"Full builder oracle count");
        }
    }
    std::size_t calls=0;
    PackedField affine([](Vec3 p){return p.x+2*p.y-p.z;},6,Encoding::f32,calls);
    check(calls==8*262144,"Depth-six oracle count");
    for(Vec3 p:std::vector<Vec3>{{-1,-1,-1},{1,1,1},{.017,-.023,.011},{0,0,0}})
        check(std::abs(affine.query(p)-(p.x+2*p.y-p.z))<1e-12,"Depth-six affine reproduction");
    rejects([&]{PackedField f(oracle,7,Encoding::q8,calls);});
    rejects([&]{PackedField f(oracle,1,Encoding::branch,calls);});
    rejects([&]{PackedField f([](Vec3){return std::numeric_limits<double>::infinity();},0,Encoding::f32,calls);});
}
void selective_refinement() {
    auto oracle=analytic_oracle("torus",.035);
    const auto fit=sample_queries(*oracle,10000,201,.5,.0175);
    const Distance truth=[&](Vec3 p){return oracle->distance(p);};
    PrecisionTree tree(truth,fit,5);
    const auto original=PrecisionDP(tree,32768,Encoding::q8).allocate(32768);
    const auto old_nodes=tree.nodes;
    const auto none=tree.refine_depth_six(truth,fit,0);
    check(none.empty() && tree.nodes.size()==37449,"Zero selection must preserve baseline");
    std::vector<std::pair<double,std::size_t>> expected;
    for(std::size_t v=0;v<tree.nodes.size();++v)
        if(!tree.spatial.nodes[v].field.first_child && tree.nodes[v].loss[2]>0)
            expected.emplace_back(-tree.nodes[v].loss[2],v);
    std::sort(expected.begin(),expected.end());
    const auto selected=tree.refine_depth_six(truth,fit,64);
    check(selected.size()==64 && tree.nodes.size()==37449+64*8,"Selective candidate count");
    check(tree.spatial.oracle_calls==8*37449+64*64,"Selective oracle count");
    for(std::size_t i=0;i<selected.size();++i) {
        check(selected[i].node==expected[i].second,"Residual ranking mismatch");
        check(selected[i].fitting_count>0,"Empty selected parent");
        const auto first=tree.spatial.nodes[selected[i].node].field.first_child;
        for(unsigned c=0;c<8;++c) for(unsigned e=0;e<3;++e) {
            const auto b=selected[i].box.child(c); double loss=0;
            const auto corners=tree.nodes[first+c].corners[e].decode();
            for(const auto& s:fit) {
                // Route from root to avoid inclusive-boundary ambiguity.
                std::size_t v=0; Box box;
                for(unsigned d=0;d<6;++d) {
                    const auto f=tree.spatial.nodes[v].field.first_child;
                    if(!f) break;
                    const auto child=box.octant(s.point); v=f+child; box=box.child(child);
                }
                if(v==first+c) { const auto err=(interpolate(corners,b,s.point)-s.distance)/2; loss+=err*err; }
            }
            check(std::abs(loss-tree.nodes[first+c].loss[e])<1e-14,"Child fitting loss mismatch");
        }
    }
    for(std::size_t v=0;v<old_nodes.size();++v)
        check(tree.nodes[v].loss==old_nodes[v].loss,"Refinement changed baseline option loss");
    const auto a=PrecisionDP(tree,32768,Encoding::q8).allocate(32768); PackedField field(tree,a);
    check(a.loss<=original.loss+1e-12 && field.bytes()<=32768,"Nested candidate optimum/cap");
    double actual=0;
    for(const auto& s:fit) { const auto e=(field.query(s.point)-s.distance)/2; actual+=e*e; }
    check(std::abs(actual-a.loss)<1e-12,"Selective decoded SSE mismatch");
    rejects([&]{tree.refine_depth_six(truth,fit,64);});
    PrecisionTree small(truth,fit,2);
    rejects([&]{small.refine_depth_six(truth,fit,64);});
    PrecisionTree constant([](Vec3){return 1.;},{{{0,0,0},1}},5);
    check(constant.refine_depth_six([](Vec3){return 1.;},{{{0,0,0},1}},64).empty(),"Zero residual cells selected");
}
void topology() {
    for(const auto& shape:{"sphere","torus","thin_box","double_sphere"}) {
        auto oracle=analytic_oracle(shape,.10);
        auto mesh=extract_surface([&](Vec3 p){return oracle->distance(p-Vec3{.017,-.023,.011});},32);
        auto t=surface_topology(mesh);
        check(t.closed_manifold(),"Oracle surface is not a closed manifold");
        check(t.components==(std::string(shape)=="double_sphere"?2u:1u),"Wrong oracle component count");
        check(t.euler==(std::string(shape)=="torus"?0:std::string(shape)=="double_sphere"?4:2),"Wrong oracle Euler characteristic");
        mesh.triangles.pop_back(); check(!surface_topology(mesh).closed_manifold(),"Missed punctured mesh");
    }
    check(extract_surface([](Vec3){return 1.;},8).triangles.empty(),"Empty surface");
    rejects([]{extract_surface([](Vec3){return 0.;},0);});
}
}
int main() {
    try { codecs(); exhaustive(); real_tree(); full_fields(); selective_refinement(); topology(); std::cout<<checks<<" precision/surface checks passed\n"; }
    catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
