#include "sdf/allocation.hpp"
#include <limits>
#include <numeric>
#include <queue>

namespace sdf {
namespace {
using Index=std::uint32_t;
constexpr double infinity=std::numeric_limits<double>::infinity();
struct State { double loss=infinity; std::array<std::size_t,8> children{}; };
using Tables=std::vector<std::vector<State>>;

std::size_t capacity(const CandidateTree& t,std::size_t bytes) {
    if(t.nodes.empty() || bytes<AllocatedField::budget_for_splits(0))
        throw std::invalid_argument("Empty candidate or insufficient budget");
    return std::min((bytes-AllocatedField::budget_for_splits(0))/(8*sizeof(AdaptiveField::Node)),(t.nodes.size()-1)/8);
}
Tables tables(const CandidateTree& t,std::size_t cap,Allocation& stats) {
    Tables dp(t.nodes.size());
    for(std::size_t v=t.nodes.size();v-->0;) {
        dp[v].push_back({t.nodes[v].loss,{}});
        const auto first=t.nodes[v].field.first_child;
        if(first && cap) {
            std::vector<State> previous(1); previous[0].loss=0;
            for(unsigned c=0;c<8;++c) {
                const auto& child=dp[first+c];
                std::vector<State> next(std::min(cap,previous.size()+child.size()-1));
                for(std::size_t a=0;a<previous.size();++a)
                    for(std::size_t b=0;b<child.size() && a+b<next.size();++b) {
                        ++stats.transitions;
                        const double loss=previous[a].loss+child[b].loss;
                        if(loss<next[a+b].loss) {
                            next[a+b]=previous[a]; next[a+b].loss=loss; next[a+b].children[c]=b;
                        }
                    }
                previous=std::move(next);
            }
            dp[v].insert(dp[v].end(),previous.begin(),previous.end());
        }
        stats.table_states+=dp[v].size();
    }
    return dp;
}
void recover(const CandidateTree& t,const Tables& dp,Index v,std::size_t k,Allocation& a) {
    if(!k) return;
    a.split[v]=1;
    const auto first=t.nodes[v].field.first_child;
    for(unsigned c=0;c<8;++c) recover(t,dp,first+c,dp[v][k].children[c],a);
}
double selected_loss(const CandidateTree& t,const Allocation& a,Index v=0) {
    if(!a.split[v]) return t.nodes[v].loss;
    double result=0;
    for(unsigned c=0;c<8;++c) result+=selected_loss(t,a,t.nodes[v].field.first_child+c);
    return result;
}
}
CandidateTree::CandidateTree(const Distance& oracle,const Samples& samples,unsigned depth,Box box):domain(box),sample_count(samples.size()) {
    box.validate();
    // Deliberate small-tree safeguard; depth 5 has 37,449 candidate nodes.
    if(depth>5 || samples.empty()) throw std::invalid_argument("Candidate depth must be <=5 and samples nonempty");
    for(const auto& s:samples) if(!box.contains(s.point) || !std::isfinite(s.distance)) throw std::invalid_argument("Invalid fitting sample");
    nodes.resize(1);
    std::function<void(Index,Box,unsigned,const std::vector<std::size_t>&)> build;
    build=[&](Index v,Box b,unsigned d,const std::vector<std::size_t>& ids) {
        for(unsigned c=0;c<8;++c) {
            const double value=oracle(b.corner(c)); ++oracle_calls;
            if(!std::isfinite(value) || std::abs(value)>std::numeric_limits<float>::max()) throw std::invalid_argument("Invalid candidate corner");
            nodes[v].field.values[c]=static_cast<float>(value);
        }
        for(auto id:ids) {
            const auto& s=samples[id];
            const double e=(interpolate(nodes[v].field.values,b,s.point)-s.distance)/domain.scale();
            nodes[v].loss+=e*e;
        }
        if(!std::isfinite(nodes[v].loss)) throw std::invalid_argument("Fitting loss overflow");
        if(d==depth) return;
        const auto first=static_cast<Index>(nodes.size()); nodes[v].field.first_child=first;
        nodes.resize(nodes.size()+8);
        std::array<std::vector<std::size_t>,8> groups;
        for(auto id:ids) groups[b.octant(samples[id].point)].push_back(id);
        for(unsigned c=0;c<8;++c) build(first+c,b.child(c),d+1,groups[c]);
    };
    std::vector<std::size_t> ids(samples.size()); std::iota(ids.begin(),ids.end(),0);
    build(0,box,0,ids);
}
std::size_t AllocatedField::budget_for_splits(std::size_t k) {
    if(k>(std::numeric_limits<std::size_t>::max()-sizeof(AllocatedField)-sizeof(AdaptiveField::Node))/(8*sizeof(AdaptiveField::Node)))
        throw std::invalid_argument("Budget overflow");
    return sizeof(AllocatedField)+(1+8*k)*sizeof(AdaptiveField::Node);
}
AllocatedField::AllocatedField(const CandidateTree& t,const Allocation& a):domain_(t.domain) {
    if(a.split.size()!=t.nodes.size()) throw std::invalid_argument("Allocation size mismatch");
    std::vector<AdaptiveField::Node> packed(1);
    std::function<void(Index,Index)> copy=[&](Index src,Index dst) {
        packed[dst].values=t.nodes[src].field.values;
        if(!a.split[src]) return;
        const auto first=t.nodes[src].field.first_child;
        if(!first) throw std::invalid_argument("Cannot split terminal candidate");
        const auto out=static_cast<Index>(packed.size()); packed[dst].first_child=out; packed.resize(packed.size()+8);
        for(unsigned c=0;c<8;++c) copy(first+c,out+c);
    };
    copy(0,0); count_=packed.size(); nodes_=std::make_unique<AdaptiveField::Node[]>(count_);
    std::copy(packed.begin(),packed.end(),nodes_.get());
}
double AllocatedField::query(Vec3 p) const {
    if(!domain_.contains(p)) throw std::out_of_range("Query outside candidate domain");
    Index v=0; Box b=domain_;
    while(nodes_[v].first_child) { const auto c=b.octant(p); v=nodes_[v].first_child+c; b=b.child(c); }
    return interpolate(nodes_[v].values,b,p);
}
Allocation exact_allocate(const CandidateTree& t,std::size_t bytes) {
    const auto cap=capacity(t,bytes);
    Allocation a; a.split.resize(t.nodes.size());
    const auto dp=tables(t,cap,a);
    for(std::size_t k=1;k<dp[0].size();++k) if(dp[0][k].loss<dp[0][a.splits].loss) a.splits=k;
    recover(t,dp,0,a.splits,a); a.loss=selected_loss(t,a); return a;
}
Allocation lookahead_allocate(const CandidateTree& t,std::size_t bytes,std::size_t horizon,bool best_prefix) {
    if(!horizon || horizon>64) throw std::invalid_argument("Lookahead horizon must be in [1,64]");
    const auto cap=capacity(t,bytes);
    Allocation a; a.split.resize(t.nodes.size()); a.loss=t.nodes[0].loss;
    const auto dp=tables(t,std::min(cap,horizon),a);
    Allocation best=a;
    struct Proposal { double score; Index node; std::size_t k; };
    struct Lower { bool operator()(const Proposal& x,const Proposal& y) const {
        if(x.score!=y.score) return x.score<y.score;
        if(x.node!=y.node) return x.node>y.node;
        return x.k>y.k;
    }};
    std::priority_queue<Proposal,std::vector<Proposal>,Lower> queue;
    auto propose=[&](Index v) {
        for(std::size_t k=1;k<dp[v].size();++k) queue.push({(t.nodes[v].loss-dp[v][k].loss)/double(k*8*sizeof(AdaptiveField::Node)),v,k});
    };
    propose(0);
    while(!queue.empty() && a.splits<cap) {
        const auto p=queue.top(); queue.pop();
        if(a.split[p.node] || p.k>cap-a.splits) continue;
        if(p.score<0) ++a.negative_steps;
        recover(t,dp,p.node,p.k,a); a.splits+=p.k;
        a.loss=selected_loss(t,a);
        if(a.loss<best.loss) best=a;
        std::function<void(Index)> frontier=[&](Index v) {
            if(!a.split[v]) { propose(v); return; }
            for(unsigned c=0;c<8;++c) frontier(t.nodes[v].field.first_child+c);
        };
        frontier(p.node);
    }
    if(best_prefix) { best.negative_steps=a.negative_steps; return best; }
    return a;
}
} // namespace sdf
