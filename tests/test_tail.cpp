#include "sdf/tail.hpp"
#include "sdf/metrics.hpp"
#include <iostream>
#include <limits>
using namespace sdf;
namespace {
int checks=0;
void check(bool ok,const char* m) { ++checks; if(!ok) throw std::runtime_error(m); }
template<class F> void rejects(F f) { bool caught=false; try { f(); } catch(const std::exception&) { caught=true; } check(caught,"Expected rejection"); }
PrecisionAllocation pruning(const PrecisionTree& tree,int mask) {
    PrecisionAllocation a; a.choice.resize(tree.nodes.size(),Encoding::q8);
    if(mask>=0) {
        a.choice[0]=Encoding::branch;
        for(unsigned c=0;c<8;++c) if(mask&(1<<c)) a.choice[tree.spatial.nodes[0].field.first_child+c]=Encoding::branch;
    }
    std::function<void(std::size_t)> visit=[&](std::size_t v) {
        if(a.choice[v]==Encoding::branch) { ++a.branches; for(unsigned c=0;c<8;++c) visit(tree.spatial.nodes[v].field.first_child+c); }
        else { ++a.leaves[2]; a.loss+=tree.nodes[v].loss[2]; }
    }; visit(0); a.bytes=sizeof(PackedField)+4*a.branches+20*a.leaves[2]; return a;
}
double direct(const PackedField& f,const Samples& fit,Q8Objective objective,double parameter) {
    double result=0;
    for(const auto& s:fit) {
        const double pred=f.query(s.point),e=std::abs(pred-s.distance);
        if(objective==Q8Objective::excess) result+=std::max(e-parameter,0.);
        else result+=e*e/4+(objective==Q8Objective::sign && std::abs(s.distance)>2e-9 && (s.distance<0?pred>=0:pred<=0)?parameter:0);
    } return result;
}
void exhaustive() {
    auto oracle=analytic_oracle("torus",.06); const auto fit=sample_queries(*oracle,300,19,.5,.03);
    PrecisionTree tree([&](Vec3 p){return oracle->distance(p);},fit,2); Q8Losses losses(tree,fit);
    for(auto objective:{Q8Objective::mse,Q8Objective::sign,Q8Objective::excess}) for(double parameter:{0.,.0025,.05,.3}) {
        losses.apply(tree,objective,parameter);
        std::vector<std::pair<std::size_t,double>> brute;
        for(int mask=-1;mask<256;++mask) {
            const auto a=pruning(tree,mask); const PackedField f(tree,a); const double loss=direct(f,fit,objective,parameter);
            check(std::abs(loss-a.loss)<1e-11,"Additive loss differs from packed field"); brute.push_back({a.bytes,loss});
        }
        PrecisionDP dp(tree,1600,Encoding::q8);
        for(std::size_t budget=sizeof(PackedField)+20;budget<=1600;budget+=11) {
            double best=std::numeric_limits<double>::infinity();
            for(auto option:brute) if(option.first<=budget) best=std::min(best,option.second);
            const auto a=dp.allocate(budget);
            check(std::abs(best-a.loss)<1e-11,"DP differs from all 257 depth-two prunings");
        }
    }
    const auto tail=allocate_cvar99(tree,losses,1000,9);
    double brute=std::numeric_limits<double>::infinity();
    for(const auto& step:tail.search) for(int mask=-1;mask<256;++mask) {
        const auto a=pruning(tree,mask); if(a.bytes>1000) continue;
        brute=std::min(brute,step.threshold+direct(PackedField(tree,a),fit,Q8Objective::excess,step.threshold)/3);
    }
    check(std::abs(brute-tail.objective)<1e-11,"Finite search differs from brute force");
    check(tail.fit_cvar99<=tail.objective+1e-11 && tail.objective<=tail.upper+1e-11,"CVaR baseline upper bound failed");
    check(std::abs(tail.fit_cvar99-evaluate([&](Vec3 p){return PackedField(tree,tail.allocation).query(p);},fit).cvar99)<1e-11,"Independent sorted-tail mismatch");
    rejects([&]{allocate_cvar99(tree,losses,1000,1);});
    rejects([&]{losses.apply(tree,Q8Objective::excess,-1);});
    rejects([&]{losses.apply(tree,Q8Objective::sign,std::numeric_limits<double>::infinity());});
    PrecisionTree constant([](Vec3){return -1.;},{{{0,0,0},-1.}},2); Q8Losses zero(constant,{{{0,0,0},-1.}});
    const auto a=allocate_cvar99(constant,zero,1000,9);
    check(a.objective==0 && a.search.size()==1 && a.allocation.bytes==sizeof(PackedField)+20,"Zero-error threshold deduplication");
}
}
int main() {
    try { exhaustive(); std::cout<<checks<<" tail checks passed\n"; }
    catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
