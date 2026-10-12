#include "sdf/tail.hpp"
#include <limits>
#include <numeric>

namespace sdf {
namespace {
double tail_mean(const std::vector<double>& sorted) {
    const double mass=.01*sorted.size();
    const auto whole=static_cast<std::size_t>(mass);
    double sum=std::accumulate(sorted.end()-whole,sorted.end(),0.);
    if(mass>whole) sum+=(mass-whole)*sorted[sorted.size()-whole-1];
    return sum/mass;
}
}
Q8Losses::Q8Losses(const PrecisionTree& t,const Samples& fit):nodes_(t.nodes.size()),count_(fit.size()) {
    if(fit.empty() || t.nodes.empty()) throw std::invalid_argument("Empty q8 loss input");
    for(const auto& s:fit) {
        if(!t.spatial.domain.contains(s.point) || !std::isfinite(s.distance)) throw std::invalid_argument("Invalid loss sample");
        std::size_t v=0; Box b=t.spatial.domain;
        for(;;) {
            const auto prediction=interpolate(t.nodes[v].corners[2].decode(),b,s.point);
            const auto error=std::abs(prediction-s.distance);
            auto& n=nodes_[v]; n.errors.push_back(error);
            n.mse+=std::pow(error/t.spatial.domain.scale(),2);
            if(std::abs(s.distance)>2e-9 && (s.distance<0?prediction>=0:prediction<=0)) ++n.wrong;
            if(!std::isfinite(n.mse)) throw std::invalid_argument("Q8 loss overflow");
            const auto first=t.spatial.nodes[v].field.first_child;
            if(!first) break;
            const auto c=b.octant(s.point); v=first+c; b=b.child(c);
        }
    }
    for(auto& n:nodes_) {
        std::sort(n.errors.begin(),n.errors.end()); n.prefix.push_back(0);
        for(double e:n.errors) n.prefix.push_back(n.prefix.back()+e);
    }
}
void Q8Losses::apply(PrecisionTree& t,Q8Objective objective,double parameter) const {
    if(t.nodes.size()!=nodes_.size() || !std::isfinite(parameter) || parameter<0)
        throw std::invalid_argument("Invalid loss parameter/tree");
    for(std::size_t v=0;v<nodes_.size();++v) {
        const auto& n=nodes_[v]; double loss=n.mse;
        if(objective==Q8Objective::sign) loss+=parameter*n.wrong;
        else if(objective==Q8Objective::excess) {
            const auto k=std::upper_bound(n.errors.begin(),n.errors.end(),parameter)-n.errors.begin();
            loss=std::max(0.,n.prefix.back()-n.prefix[k]-parameter*(n.errors.size()-k));
        } else if(objective!=Q8Objective::mse) throw std::invalid_argument("Unknown q8 objective");
        t.nodes[v].loss[2]=loss;
    }
}
std::vector<double> Q8Losses::errors(const PrecisionTree& t,const PrecisionAllocation& a) const {
    if(t.nodes.size()!=nodes_.size() || a.choice.size()!=nodes_.size()) throw std::invalid_argument("Loss allocation mismatch");
    std::vector<double> result; result.reserve(count_);
    std::function<void(std::size_t)> visit=[&](std::size_t v) {
        if(a.choice[v]==Encoding::q8) result.insert(result.end(),nodes_[v].errors.begin(),nodes_[v].errors.end());
        else if(a.choice[v]==Encoding::branch && t.spatial.nodes[v].field.first_child) {
            for(unsigned c=0;c<8;++c) visit(t.spatial.nodes[v].field.first_child+c);
        } else throw std::invalid_argument("Expected q8 allocation");
    }; visit(0);
    if(result.size()!=count_) throw std::invalid_argument("Loss sample coverage mismatch");
    std::sort(result.begin(),result.end()); return result;
}
double Q8Losses::cvar99(const PrecisionTree& t,const PrecisionAllocation& a) const { return tail_mean(errors(t,a)); }
double Q8Losses::p99(const PrecisionTree& t,const PrecisionAllocation& a) const {
    const auto e=errors(t,a); return e[static_cast<std::size_t>(std::ceil(.99*e.size()))-1];
}
TailAllocation allocate_cvar99(PrecisionTree& tree,const Q8Losses& losses,std::size_t budget,unsigned grid_points) {
    if(grid_points<2 || grid_points>4097) throw std::invalid_argument("Threshold grid must contain 2..4097 points");
    losses.apply(tree,Q8Objective::mse);
    const auto baseline=PrecisionDP(tree,budget,Encoding::q8).allocate(budget);
    TailAllocation result; result.upper=losses.cvar99(tree,baseline); result.objective=std::numeric_limits<double>::infinity();
    std::vector<double> thresholds{losses.p99(tree,baseline)};
    for(unsigned i=0;i<grid_points;++i) thresholds.push_back(result.upper*i/(grid_points-1));
    std::sort(thresholds.begin(),thresholds.end()); thresholds.erase(std::unique(thresholds.begin(),thresholds.end()),thresholds.end());
    for(double t:thresholds) {
        losses.apply(tree,Q8Objective::excess,t);
        PrecisionDP dp(tree,budget,Encoding::q8); auto a=dp.allocate(budget);
        const double objective=t+a.loss/(.01*losses.count());
        result.search.push_back({t,a.loss,objective,a.bytes,dp.states,dp.transitions});
        if(objective<result.objective || (objective==result.objective && a.bytes<result.allocation.bytes)) {
            result.threshold=t; result.objective=objective; result.allocation=std::move(a);
        }
    }
    result.fit_cvar99=losses.cvar99(tree,result.allocation);
    // Leave the tree consistent with the returned additive allocation loss.
    losses.apply(tree,Q8Objective::excess,result.threshold);
    return result;
}
} // namespace sdf
