#include "sdf/allocation.hpp"
#include <iostream>
#include <limits>
#include <random>

using namespace sdf;
namespace {
int checks=0;
void check(bool c,const char* message) { ++checks; if(!c) throw std::runtime_error(message); }
void close(double a,double b,const char* message) { check(std::abs(a-b)<1e-10,message); }
template<class F> void rejects(F f) { bool caught=false; try { f(); } catch(const std::exception&) { caught=true; } check(caught,"Expected rejection"); }
void exhaustive() {
    auto oracle=analytic_oracle("torus",0.04);
    auto samples=sample_queries(*oracle,500,63,0.5,0.02);
    CandidateTree tree([&](Vec3 p){return oracle->distance(p);},samples,2);
    std::mt19937_64 rng(879);
    // Independently enumerate all 257 valid prunings of a depth-two octree.
    // Random costs cover nonmonotone losses, ties, and unused budget.
    for(int trial=0;trial<30;++trial) {
        if(trial) for(auto& node:tree.nodes) node.loss=double(rng()%31);
        for(std::size_t cap=0;cap<=12;++cap) {
            double optimum=tree.nodes[0].loss;
            for(unsigned mask=0;mask<256;++mask) {
                std::size_t splits=1; double loss=0;
                for(unsigned c=0;c<8;++c) {
                    const auto v=tree.nodes[0].field.first_child+c;
                    if(mask&(1u<<c)) {
                        ++splits;
                        for(unsigned j=0;j<8;++j) loss+=tree.nodes[tree.nodes[v].field.first_child+j].loss;
                    } else loss+=tree.nodes[v].loss;
                }
                if(splits<=cap) optimum=std::min(optimum,loss);
            }
            const auto bytes=AllocatedField::budget_for_splits(cap)+17;
            const auto exact=exact_allocate(tree,bytes);
            close(exact.loss,optimum,"DP differs from exhaustive enumeration");
            check(AllocatedField(tree,exact).bytes()<=bytes,"Exact cap exceeded");
            for(auto h:{1u,2u,4u}) {
                const auto a=lookahead_allocate(tree,bytes,h);
                check(a.loss+1e-10>=optimum,"Heuristic below optimum");
                check(AllocatedField(tree,a).bytes()<=bytes,"Heuristic cap exceeded");
                check(a.loss<=tree.nodes[0].loss+1e-10,"Lost best prefix");
            }
        }
    }
}
void unlock() {
    CandidateTree t([](Vec3){return 0.;},{{{0,0,0},0}},3);
    for(auto& n:t.nodes) n.loss=0;
    t.nodes[0].loss=100;
    const auto a=t.nodes[0].field.first_child,b=a+1,deep=t.nodes[a].field.first_child;
    t.nodes[a].loss=10; t.nodes[b].loss=9; t.nodes[deep].loss=11;
    t.nodes[t.nodes[b].field.first_child].loss=5;
    t.nodes[t.nodes[t.nodes[b].field.first_child].field.first_child].loss=5;
    const auto cap=AllocatedField::budget_for_splits(3);
    close(exact_allocate(t,cap).loss,9,"Unlock exact loss");
    close(lookahead_allocate(t,cap,1).loss,15,"Greedy trap");
    close(lookahead_allocate(t,cap,2).loss,9,"Two-split lookahead crosses negative split");
}
void real_tree() {
    auto oracle=analytic_oracle("thin_box",0.037);
    auto s=sample_queries(*oracle,1024,82,0.5,0.018);
    std::size_t calls=0;
    CandidateTree t([&](Vec3 p){++calls; return oracle->distance(p);},s,3);
    check(t.nodes.size()==585 && calls==4680 && t.oracle_calls==calls,"Candidate accounting");
    for(auto k:{0u,1u,7u,19u,73u,100u}) {
        const auto cap=AllocatedField::budget_for_splits(k);
        for(auto a:{exact_allocate(t,cap),lookahead_allocate(t,cap,1,false),lookahead_allocate(t,cap,4)}) {
            AllocatedField f(t,a); double loss=0;
            for(const auto& sample:s) { const double e=(f.query(sample.point)-sample.distance)/2; loss+=e*e; }
            close(a.loss,loss,"Reconstructed field SSE differs");
            check(f.bytes()==AllocatedField::budget_for_splits(a.splits),"Storage differs from recurrence");
            for(unsigned c=0;c<8;++c) check(std::isfinite(f.query(Box{}.corner(c))),"Boundary query");
            rejects([&]{f.query({2,0,0});});
        }
    }
    check(calls==4680,"Allocation queried oracle");
    rejects([&]{exact_allocate(t,1);});
    rejects([&]{lookahead_allocate(t,1000,0);});
    rejects([&]{CandidateTree bad([](Vec3){return 0.;},{},2);});
    rejects([&]{CandidateTree bad([](Vec3){return 0.;},s,6);});
    rejects([&]{CandidateTree bad([](Vec3){return std::numeric_limits<double>::infinity();},s,1);});
    CandidateTree leaf([](Vec3 p){return p.x;},{{{0,0,0},0}},0);
    close(exact_allocate(leaf,10000).loss,0,"Depth-zero tree");
}
}
int main() {
    try { exhaustive(); unlock(); real_tree(); std::cout<<checks<<" allocation checks passed\n"; return 0; }
    catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
