#include "sdf/field.hpp"
#include <limits>
#include <numeric>
#include <queue>

namespace sdf {
namespace {
Corners corners(const Distance& f, const Box& b) {
    Corners result{};
    for(unsigned i=0;i<8;++i) {
        result[i]=static_cast<float>(f(b.corner(i)));
        if(!std::isfinite(result[i])) throw std::runtime_error("Non-finite reference distance");
    }
    return result;
}
void check_point(Vec3 p, const Box& b) {
    if(!b.contains(p)) throw std::out_of_range("Query outside field domain or non-finite");
}
}
UniformGrid::UniformGrid(const Distance& f, std::size_t budget, Box domain):domain_(domain) {
    domain.validate();
    if(budget<sizeof(*this)+8*sizeof(float)) throw std::invalid_argument("Grid budget too small");
    const auto available=(budget-sizeof(*this))/sizeof(float);
    n_=static_cast<std::size_t>(std::cbrt(static_cast<double>(available)));
    while((n_+1)<=available/(n_+1)/(n_+1)) ++n_;
    while(n_>available/n_/n_) --n_;
    values_=std::make_unique<float[]>(n_*n_*n_);
    for(std::size_t z=0;z<n_;++z) for(std::size_t y=0;y<n_;++y) for(std::size_t x=0;x<n_;++x) {
        const double d=f(domain.point({double(x)/(n_-1),double(y)/(n_-1),double(z)/(n_-1)}));
        if(!std::isfinite(d) || std::abs(d)>std::numeric_limits<float>::max()) throw std::runtime_error("Invalid distance");
        values_[x+n_*(y+n_*z)]=static_cast<float>(d);
    }
}
double UniformGrid::query(Vec3 p) const {
    check_point(p,domain_);
    std::array<std::size_t,3> base{};
    Vec3 t;
    double* coords[3]={&t.x,&t.y,&t.z};
    for(int axis=0;axis<3;++axis) {
        const double q=(p[axis]-domain_.lo[axis])/(domain_.hi[axis]-domain_.lo[axis])*double(n_-1);
        base[axis]=std::min(static_cast<std::size_t>(q),n_-2);
        *coords[axis]=q-double(base[axis]);
    }
    Corners c{};
    for(unsigned i=0;i<8;++i) c[i]=values_[base[0]+(i&1)+n_*(base[1]+((i>>1)&1)+n_*(base[2]+((i>>2)&1)))];
    return interpolate(c,{{0,0,0},{1,1,1}},t);
}
AdaptiveField::AdaptiveField(const Distance& f, const Samples& training, std::size_t budget,
                             AdaptiveOptions opt, Box domain):domain_(domain) {
    domain.validate();
    if(budget<sizeof(*this)+sizeof(Node)) throw std::invalid_argument("Adaptive budget too small");
    if(training.empty() || opt.max_depth>20 || opt.max_depth==0 || opt.sign_penalty<0 || !std::isfinite(opt.sign_penalty))
        throw std::invalid_argument("Invalid adaptive training/options");
    for(const auto& s:training) if(!domain.contains(s.point) || !std::isfinite(s.distance)) throw std::invalid_argument("Invalid training sample");
    const std::size_t max_nodes=std::min((budget-sizeof(*this))/sizeof(Node),std::size_t(std::numeric_limits<std::uint32_t>::max()));
    // Build scratch memory and oracle are discarded, leaving an exact-size query allocation.
    std::vector<Node> nodes(1);
    nodes[0].values=corners(f,domain);
    struct Candidate {
        double gain;
        std::uint32_t index;
        unsigned depth;
        Box box;
        std::array<Node,8> children;
        std::array<std::vector<std::size_t>,8> groups;
    };
    struct LowerPriority {
        bool operator()(const Candidate& a,const Candidate& b) const {
            if(a.gain!=b.gain) return a.gain<b.gain;
            return a.index>b.index;
        }
    };
    std::priority_queue<Candidate,std::vector<Candidate>,LowerPriority> queue;
    const double scale=domain.scale(), epsilon=scale*1e-9;
    auto loss=[&](double pred,double truth) {
        const double e=(pred-truth)/scale;
        const bool wrong=std::abs(truth)>epsilon && (truth<0 ? pred>=0 : pred<=0);
        return e*e+(wrong?opt.sign_penalty:0);
    };
    auto candidate=[&](std::uint32_t index,unsigned depth,Box box,const std::vector<std::size_t>& ids) {
        if(depth>=opt.max_depth || ids.empty()) return;
        Candidate c{}; c.index=index; c.depth=depth; c.box=box;
        for(unsigned i=0;i<8;++i) c.children[i].values=corners(f,box.child(i));
        for(const auto id:ids) {
            const auto& s=training[id]; const auto child=box.octant(s.point);
            c.groups[child].push_back(id);
            c.gain+=loss(interpolate(nodes[index].values,box,s.point),s.distance)
                   -loss(interpolate(c.children[child].values,box.child(child),s.point),s.distance);
        }
        // Each split adds eight 36-byte nodes. Negative gains remain eligible to
        // traverse refinement plateaus and spend the same cap; no optimality claim.
        c.gain/=double(8*sizeof(Node));
        queue.push(std::move(c));
    };
    std::vector<std::size_t> ids(training.size()); std::iota(ids.begin(),ids.end(),0);
    if(max_nodes>=9) candidate(0,0,domain,ids);
    while(!queue.empty() && nodes.size()+8<=max_nodes) {
        auto c=queue.top(); queue.pop();
        const auto first=static_cast<std::uint32_t>(nodes.size());
        nodes[c.index].first_child=first;
        for(const auto& child:c.children) nodes.push_back(child);
        if(nodes.size()+8<=max_nodes)
            for(unsigned i=0;i<8;++i) candidate(first+i,c.depth+1,c.box.child(i),c.groups[i]);
    }
    count_=nodes.size();
    nodes_=std::make_unique<Node[]>(count_);
    std::copy(nodes.begin(),nodes.end(),nodes_.get());
}
double AdaptiveField::query(Vec3 p) const {
    check_point(p,domain_);
    std::uint32_t index=0; auto box=domain_;
    while(nodes_[index].first_child) {
        const auto child=box.octant(p); index=nodes_[index].first_child+child; box=box.child(child);
    }
    return interpolate(nodes_[index].values,box,p);
}
std::vector<Box> AdaptiveField::leaves() const {
    std::vector<Box> result;
    std::function<void(std::uint32_t,Box)> visit=[&](std::uint32_t i,Box b) {
        if(!nodes_[i].first_child) result.push_back(b);
        else for(unsigned child=0;child<8;++child) visit(nodes_[i].first_child+child,b.child(child));
    };
    visit(0,domain_); return result;
}
} // namespace sdf
