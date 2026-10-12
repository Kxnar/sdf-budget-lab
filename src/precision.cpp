#include "sdf/precision.hpp"
#include <cstring>
#include <limits>

namespace sdf {
namespace {
template<class T> T read(const unsigned char* p) { T x; std::memcpy(&x,p,sizeof(x)); return x; }
template<class T> void write(unsigned char* p,T x) { std::memcpy(p,&x,sizeof(x)); }
std::size_t payload(Encoding e) {
    if(e==Encoding::f32) return 32;
    if(e==Encoding::q16) return 24;
    if(e==Encoding::q8) return 16;
    throw std::invalid_argument("Invalid leaf encoding");
}
Corners decode(const unsigned char* p,Encoding e) {
    Corners v;
    if(e==Encoding::f32) std::memcpy(v.data(),p,32);
    else {
        const double offset=read<float>(p),scale=read<float>(p+4);
        for(unsigned i=0;i<8;++i) {
            const auto q=e==Encoding::q16?read<std::uint16_t>(p+8+2*i):p[8+i];
            v[i]=static_cast<float>(offset+scale*q);
        }
    }
    return v;
}
}
std::size_t EncodedCorners::payload_bytes() const { return payload(encoding); }
Corners EncodedCorners::decode() const { return sdf::decode(data.data(),encoding); }
EncodedCorners encode_corners(const Corners& v,Encoding e) {
    EncodedCorners out; out.encoding=e; payload(e);
    for(float x:v) if(!std::isfinite(x)) throw std::invalid_argument("Nonfinite corner");
    if(e==Encoding::f32) { std::memcpy(out.data.data(),v.data(),32); return out; }
    const auto limits=std::minmax_element(v.begin(),v.end());
    const float offset=*limits.first;
    const double maxq=e==Encoding::q16?65535.:255.;
    const float scale=static_cast<float>((double(*limits.second)-offset)/maxq);
    write(out.data.data(),offset); write(out.data.data()+4,scale);
    for(unsigned i=0;i<8;++i) {
        const auto q=static_cast<unsigned>(scale>0?std::clamp(std::round((double(v[i])-offset)/scale),0.,maxq):0.);
        if(e==Encoding::q16) write(out.data.data()+8+2*i,static_cast<std::uint16_t>(q));
        else out.data[8+i]=static_cast<unsigned char>(q);
    }
    for(float x:out.decode()) if(!std::isfinite(x)) throw std::invalid_argument("Quantized corner overflow");
    return out;
}
PrecisionTree::PrecisionTree(const Distance& oracle,const Samples& fit,unsigned depth,Box box)
    :spatial(oracle,fit,depth,box),nodes(spatial.nodes.size()) {
    for(std::size_t v=0;v<nodes.size();++v)
        for(unsigned e=0;e<3;++e) nodes[v].corners[e]=encode_corners(spatial.nodes[v].field.values,static_cast<Encoding>(e+1));
    // Same points, routing, decoded arithmetic and normalization for every option.
    std::vector<std::array<Corners,3>> decoded(nodes.size());
    for(std::size_t v=0;v<nodes.size();++v) for(unsigned e=0;e<3;++e) decoded[v][e]=nodes[v].corners[e].decode();
    for(const auto& s:fit) {
        std::size_t v=0; Box b=box;
        for(;;) {
            for(unsigned e=0;e<3;++e) {
                const double error=(interpolate(decoded[v][e],b,s.point)-s.distance)/box.scale();
                nodes[v].loss[e]+=error*error;
                if(!std::isfinite(nodes[v].loss[e])) throw std::invalid_argument("Precision loss overflow");
            }
            const auto first=spatial.nodes[v].field.first_child;
            if(!first) break;
            const auto c=b.octant(s.point); v=first+c; b=b.child(c);
        }
    }
}
std::vector<RefinementCandidate> PrecisionTree::refine_depth_seven(const Distance& oracle,const Samples& ranking,std::size_t max_parents) {
    if(max_parents>32768 || ranking.empty())throw std::invalid_argument("Invalid depth-seven ranking request");
    std::vector<RefinementCandidate> terminal;std::vector<std::size_t> index(nodes.size(),std::numeric_limits<std::size_t>::max());
    std::function<void(std::size_t,Box,unsigned)> collect=[&](std::size_t v,Box b,unsigned d){auto first=spatial.nodes[v].field.first_child;if(first){if(d>=6)throw std::invalid_argument("Expected maximum depth six before depth-seven refinement");for(unsigned c=0;c<8;++c)collect(first+c,b.child(c),d+1);}else if(d==6){index[v]=terminal.size();terminal.push_back({v,0,b,0});}};collect(0,spatial.domain,0);
    for(auto s:ranking){if(!spatial.domain.contains(s.point)||!std::isfinite(s.distance))throw std::invalid_argument("Invalid depth-seven ranking sample");std::size_t v=0;Box b=spatial.domain;while(spatial.nodes[v].field.first_child){auto c=b.octant(s.point);v=spatial.nodes[v].field.first_child+c;b=b.child(c);}auto ix=index[v];if(ix!=std::numeric_limits<std::size_t>::max()){auto& r=terminal[ix];++r.fitting_count;double e=(interpolate(nodes[v].corners[2].decode(),b,s.point)-s.distance)/spatial.domain.scale();r.q8_loss+=e*e;}}
    terminal.erase(std::remove_if(terminal.begin(),terminal.end(),[](auto r){return !r.fitting_count||r.q8_loss<=0;}),terminal.end());std::sort(terminal.begin(),terminal.end(),[](auto a,auto b){return a.q8_loss>b.q8_loss||(a.q8_loss==b.q8_loss&&a.node<b.node);});if(terminal.size()>max_parents)terminal.resize(max_parents);
    for(auto r:terminal){auto first=static_cast<std::uint32_t>(nodes.size());spatial.nodes[r.node].field.first_child=first;spatial.nodes.resize(nodes.size()+8);nodes.resize(nodes.size()+8);for(unsigned c=0;c<8;++c){auto b=r.box.child(c);auto& sn=spatial.nodes[first+c];for(unsigned i=0;i<8;++i){double v=oracle(b.corner(i));++spatial.oracle_calls;if(!std::isfinite(v))throw std::invalid_argument("Invalid depth-seven corner");sn.field.values[i]=float(v);}for(unsigned e=0;e<3;++e)nodes[first+c].corners[e]=encode_corners(sn.field.values,static_cast<Encoding>(e+1));}}
    return terminal;
}
std::vector<RefinementCandidate> PrecisionTree::refine_depth_six(
    const Distance& oracle,const Samples& fit,std::size_t max_parents) {
    if(nodes.size()!=37449 || spatial.nodes.size()!=nodes.size() ||
       fit.size()!=spatial.sample_count || max_parents>32768)
        throw std::invalid_argument("Selective refinement requires full depth five and original fitting pool");
    std::vector<std::vector<std::size_t>> ids(nodes.size());
    for(std::size_t i=0;i<fit.size();++i) {
        if(!spatial.domain.contains(fit[i].point) || !std::isfinite(fit[i].distance))
            throw std::invalid_argument("Invalid refinement fitting sample");
        std::size_t v=0; Box b=spatial.domain;
        for(unsigned d=0;d<5;++d) {
            const auto first=spatial.nodes[v].field.first_child;
            if(!first || first+7>=nodes.size()) throw std::invalid_argument("Expected full depth-five tree");
            const auto c=b.octant(fit[i].point); v=first+c; b=b.child(c);
        }
        if(spatial.nodes[v].field.first_child) throw std::invalid_argument("Already refined");
        ids[v].push_back(i);
    }
    std::vector<RefinementCandidate> ranked;
    std::function<void(std::size_t,Box,unsigned)> collect=[&](std::size_t v,Box b,unsigned d) {
        const auto first=spatial.nodes[v].field.first_child;
        if(d==5) {
            if(first) throw std::invalid_argument("Already refined");
            if(!ids[v].empty() && nodes[v].loss[2]>0)
                ranked.push_back({v,ids[v].size(),b,nodes[v].loss[2]});
        } else {
            if(!first || first+7>=nodes.size()) throw std::invalid_argument("Expected full depth-five tree");
            for(unsigned c=0;c<8;++c) collect(first+c,b.child(c),d+1);
        }
    };
    collect(0,spatial.domain,0);
    std::sort(ranked.begin(),ranked.end(),[](const auto& a,const auto& b) {
        return a.q8_loss>b.q8_loss || (a.q8_loss==b.q8_loss && a.node<b.node);
    });
    if(ranked.size()>max_parents) ranked.resize(max_parents);
    for(const auto& r:ranked) {
        const auto first=static_cast<std::uint32_t>(nodes.size());
        spatial.nodes[r.node].field.first_child=first;
        spatial.nodes.resize(nodes.size()+8); nodes.resize(nodes.size()+8);
        for(unsigned c=0;c<8;++c) {
            const Box b=r.box.child(c); auto& sn=spatial.nodes[first+c]; auto& pn=nodes[first+c];
            for(unsigned corner=0;corner<8;++corner) {
                const double value=oracle(b.corner(corner)); ++spatial.oracle_calls;
                if(!std::isfinite(value) || std::abs(value)>std::numeric_limits<float>::max())
                    throw std::invalid_argument("Invalid refinement corner");
                sn.field.values[corner]=static_cast<float>(value);
            }
            for(unsigned e=0;e<3;++e) {
                pn.corners[e]=encode_corners(sn.field.values,static_cast<Encoding>(e+1));
                const auto decoded=pn.corners[e].decode();
                for(auto i:ids[r.node]) if(r.box.octant(fit[i].point)==c) {
                    const double err=(interpolate(decoded,b,fit[i].point)-fit[i].distance)/spatial.domain.scale();
                    pn.loss[e]+=err*err;
                }
                if(!std::isfinite(pn.loss[e])) throw std::invalid_argument("Refinement loss overflow");
            }
            sn.loss=pn.loss[0];
        }
    }
    return ranked;
}
PrecisionDP::PrecisionDP(const PrecisionTree& t,std::size_t max_bytes,Encoding fixed):tree_(t),budget_(max_bytes),tables_(t.nodes.size()) {
    if(t.nodes.empty() || t.nodes.size()!=t.spatial.nodes.size() || unsigned(fixed)>3 || max_bytes<sizeof(PackedField)+20 || max_bytes>=std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument("Invalid precision DP input/budget");
    const auto cap=(max_bytes-sizeof(PackedField))/4;
    auto prune=[](std::vector<State>& dense) {
        std::vector<State> result; double best=std::numeric_limits<double>::infinity();
        for(auto& s:dense) if(s.loss<best) { best=s.loss; result.push_back(std::move(s)); }
        return result;
    };
    for(std::size_t v=t.nodes.size();v-->0;) {
        std::vector<State> options;
        for(unsigned e=1;e<=3;++e) if(fixed==Encoding::branch || unsigned(fixed)==e) {
            const auto cost=(4+payload(static_cast<Encoding>(e)))/4;
            const double loss=t.nodes[v].loss[e-1];
            if(!std::isfinite(loss) || loss<0) throw std::invalid_argument("Invalid option loss");
            if(cost<=cap) options.push_back({cost,loss,static_cast<Encoding>(e),{}});
        }
        const auto first=t.spatial.nodes[v].field.first_child;
        if(first) {
            if(first<=v || first+7>=t.nodes.size()) throw std::invalid_argument("Invalid candidate ordering");
            std::vector<State> previous{{1,0,Encoding::branch,{}}};
            for(unsigned c=0;c<8 && !previous.empty();++c) {
                const auto& child=tables_[first+c];
                if(child.empty()) { previous.clear(); break; }
                const auto limit=std::min(cap,previous.back().units+child.back().units);
                std::vector<State> dense(limit+1);
                for(std::size_t k=0;k<dense.size();++k) { dense[k].units=k; dense[k].loss=std::numeric_limits<double>::infinity(); }
                for(const auto& a:previous) for(std::size_t j=0;j<child.size();++j) {
                    const auto& b=child[j]; const auto k=a.units+b.units;
                    if(k>limit) break;
                    ++transitions;
                    const double loss=a.loss+b.loss;
                    if(loss<dense[k].loss) { dense[k]=a; dense[k].units=k; dense[k].loss=loss; dense[k].children[c]=j; }
                }
                previous=prune(dense);
            }
            options.insert(options.end(),previous.begin(),previous.end());
        }
        std::stable_sort(options.begin(),options.end(),[](const State& a,const State& b) {
            return a.units<b.units || (a.units==b.units && a.loss<b.loss);
        });
        // Same-cost states after the first cannot survive strict loss pruning.
        tables_[v]=prune(options); states+=tables_[v].size();
    }
    if(tables_[0].empty()) throw std::invalid_argument("Budget cannot fit requested encoding");
}
PrecisionAllocation PrecisionDP::allocate(std::size_t bytes) const {
    if(bytes>budget_ || bytes<sizeof(PackedField)) throw std::invalid_argument("Budget outside DP range");
    const auto cap=(bytes-sizeof(PackedField))/4;
    const auto& root=tables_[0];
    auto it=std::upper_bound(root.begin(),root.end(),cap,[](std::size_t x,const State& s){return x<s.units;});
    if(it==root.begin()) throw std::invalid_argument("Insufficient encoding budget");
    --it; PrecisionAllocation result; result.choice.resize(tree_.nodes.size(),Encoding::f32);
    result.loss=it->loss; result.bytes=sizeof(PackedField)+it->units*4;
    std::function<void(std::size_t,const State&)> recover=[&](std::size_t v,const State& s) {
        result.choice[v]=s.encoding;
        if(s.encoding!=Encoding::branch) { ++result.leaves[unsigned(s.encoding)-1]; return; }
        ++result.branches; const auto first=tree_.spatial.nodes[v].field.first_child;
        for(unsigned c=0;c<8;++c) recover(first+c,tables_[first+c][s.children[c]]);
    };
    recover(0,*it); return result;
}
PackedField::PackedField(const PrecisionTree& t,const PrecisionAllocation& a):domain_(t.spatial.domain) {
    if(a.choice.size()!=t.nodes.size()) throw std::invalid_argument("Precision allocation size mismatch");
    std::vector<unsigned char> bytes(4);
    auto reserve=[&](std::size_t n) {
        const auto pos=bytes.size();
        if(pos+n>=std::numeric_limits<std::uint32_t>::max()) throw std::invalid_argument("Packed offset overflow");
        bytes.resize(pos+n); return static_cast<std::uint32_t>(pos);
    };
    std::function<void(std::size_t,std::size_t)> pack=[&](std::size_t v,std::size_t dst) {
        const auto e=a.choice[v];
        if(e==Encoding::branch) {
            const auto first=t.spatial.nodes[v].field.first_child;
            if(!first) throw std::invalid_argument("Split terminal node");
            const auto pos=reserve(32); write(bytes.data()+dst,pos);
            for(unsigned c=0;c<8;++c) pack(first+c,pos+4*c);
        } else {
            payload(e); const auto& corners=t.nodes[v].corners[unsigned(e)-1];
            const auto pos=reserve(corners.payload_bytes()); write(bytes.data()+dst,pos|unsigned(e));
            std::memcpy(bytes.data()+pos,corners.data.data(),corners.payload_bytes());
        }
    };
    pack(0,0); size_=bytes.size(); data_=std::make_unique<unsigned char[]>(size_);
    std::copy(bytes.begin(),bytes.end(),data_.get());
    if(this->bytes()!=a.bytes) throw std::invalid_argument("Allocation byte accounting mismatch");
}
PackedField::PackedField(const Distance& oracle,unsigned depth,Encoding encoding,
                         std::size_t& oracle_calls,Box domain):domain_(domain) {
    domain.validate(); payload(encoding);
    if(depth>6) throw std::invalid_argument("Full diagnostic depth must be <=6");
    const std::size_t leaves=std::size_t(1)<<(3*depth),branches=(leaves-1)/7;
    size_=4*branches+(4+payload(encoding))*leaves;
    data_=std::make_unique<unsigned char[]>(size_);
    std::size_t next=4; oracle_calls=0;
    std::function<void(std::size_t,Box,unsigned)> build=[&](std::size_t dst,Box b,unsigned d) {
        const auto pos=next;
        if(d<depth) {
            next+=32; write(data_.get()+dst,static_cast<std::uint32_t>(pos));
            for(unsigned c=0;c<8;++c) build(pos+4*c,b.child(c),d+1);
        } else {
            Corners values;
            for(unsigned c=0;c<8;++c) {
                const double value=oracle(b.corner(c)); ++oracle_calls;
                if(!std::isfinite(value) || std::abs(value)>std::numeric_limits<float>::max())
                    throw std::invalid_argument("Invalid full-field corner");
                values[c]=static_cast<float>(value);
            }
            const auto encoded=encode_corners(values,encoding);
            next+=encoded.payload_bytes();
            write(data_.get()+dst,static_cast<std::uint32_t>(pos)|unsigned(encoding));
            std::memcpy(data_.get()+pos,encoded.data.data(),encoded.payload_bytes());
        }
    };
    build(0,domain,0);
    if(next!=size_) throw std::logic_error("Full-field size mismatch");
}
double PackedField::query(Vec3 p) const {
    if(!domain_.contains(p)) throw std::out_of_range("Query outside packed domain");
    std::size_t v=0; Box b=domain_;
    for(;;) {
        const auto descriptor=read<std::uint32_t>(data_.get()+v);
        const auto e=static_cast<Encoding>(descriptor&3u); const auto offset=descriptor&~3u;
        if(e!=Encoding::branch) return interpolate(decode(data_.get()+offset,e),b,p);
        const auto c=b.octant(p); v=offset+4*c; b=b.child(c);
    }
}
} // namespace sdf
