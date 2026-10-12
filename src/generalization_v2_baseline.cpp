#include "sdf/generalization.hpp"
#include "sdf/precision.hpp"
#include "sdf/tail.hpp"
#include "sdf/metrics.hpp"
#include "sdf/surface.hpp"
#include "sdf/ray.hpp"
#include <chrono>
#include <filesystem>
#include <sstream>
#include <iostream>
#include <numeric>
#include <map>
#include <limits>
#include <cstring>
#include <unordered_map>

namespace {
using namespace sdf; namespace fs=std::filesystem; using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
std::ofstream out(const fs::path& p) { std::ofstream f(p);f.exceptions(std::ios::badbit|std::ios::failbit);f<<std::setprecision(17);return f; }
template<class T> void binary(const fs::path& p,const T* v,std::size_t n) { std::ofstream f(p,std::ios::binary);f.exceptions(std::ios::badbit|std::ios::failbit);f.write(reinterpret_cast<const char*>(v),std::streamsize(n*sizeof(T))); }
struct Key {std::array<std::uint64_t,3> bits;bool operator==(const Key& other)const{return bits==other.bits;}};
struct KeyHash {std::size_t operator()(const Key& k)const{std::size_t h=0;for(auto v:k.bits)h^=std::hash<std::uint64_t>{}(v)+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}};
class Memo final:public Oracle {
    GeneralizationOracle base_;fs::path cache_;std::string identity_;mutable std::unordered_map<Key,double,KeyHash> values_;mutable std::vector<Sample> added_;
    static Key key(Vec3 p){Key k;double v[3]={p.x,p.y,p.z};std::memcpy(k.bits.data(),v,24);return k;}
public:
    mutable std::size_t logical{},physical{},hits{};std::size_t loaded{};double load_seconds{};
    Memo(GeometrySpec spec,fs::path cache,std::string identity):base_(std::move(spec)),cache_(std::move(cache)),identity_(std::move(identity)){auto start=Clock::now();if(!cache_.empty()){if(identity_.size()!=64)throw std::invalid_argument("Shared cache needs a64-character spec/mesh/kernel identity");if(fs::exists(cache_)){if(fs::file_size(cache_)<65||(fs::file_size(cache_)-65)%sizeof(Sample))throw std::runtime_error("Truncated shared cache; preserve and rebuild separately");std::ifstream f(cache_,std::ios::binary);std::string header;std::getline(f,header);if(header!=identity_)throw std::runtime_error("Shared cache spec/mesh/kernel identity mismatch");Sample s;while(f.read(reinterpret_cast<char*>(&s),sizeof(s))){if(!Box{}.contains(s.point)||!std::isfinite(s.distance))throw std::runtime_error("Invalid shared cache sample");auto result=values_.emplace(key(s.point),s.distance);if(!result.second&&result.first->second!=s.distance)throw std::runtime_error("Contradictory shared cache labels");++loaded;}}}load_seconds=elapsed(start);}
    double distance(Vec3 p)const override{++logical;if(!cache_.empty()){auto k=key(p);auto it=values_.find(k);if(it!=values_.end()){++hits;return it->second;}double d=base_.distance(p);++physical;values_.emplace(k,d);added_.push_back({p,d});return d;}++physical;return base_.distance(p);}
    Vec3 surface(std::mt19937_64& rng)const override{return base_.surface(rng);}
    std::string name()const override{return base_.name();}
    void save(){if(cache_.empty())return;fs::create_directories(cache_.parent_path());bool exists=fs::exists(cache_);std::ofstream f(cache_,std::ios::binary|std::ios::app);f.exceptions(std::ios::badbit|std::ios::failbit);if(!exists)f<<identity_<<'\n';f.write(reinterpret_cast<const char*>(added_.data()),std::streamsize(added_.size()*sizeof(Sample)));}
};
Samples load(const fs::path& p) { std::ifstream f(p,std::ios::binary|std::ios::ate);if(!f || std::size_t(f.tellg())%sizeof(Sample)) throw std::runtime_error("Bad samples");Samples s(std::size_t(f.tellg())/sizeof(Sample));f.seekg(0);f.read(reinterpret_cast<char*>(s.data()),std::streamsize(s.size()*sizeof(Sample)));return s; }
struct Ray { int group;Vec3 o,d;double low,high; };
std::vector<Ray> rays(const fs::path& p) { std::ifstream f(p);std::vector<Ray> r;Ray v;while(f>>v.group>>v.o.x>>v.o.y>>v.o.z>>v.d.x>>v.d.y>>v.d.z>>v.low>>v.high) r.push_back(v);if(r.empty()) throw std::runtime_error("No rays");return r; }
std::vector<double> roots(const Distance& q,const Ray& r,bool& continuous) {
    std::vector<double> v;const unsigned n=1024;double t=r.low,a=q(r.o+r.d*t);continuous=true;
    for(unsigned i=1;i<=n;++i) { double u=r.low+(r.high-r.low)*i/n,b=q(r.o+r.d*u);if((a<0)!=(b<0)) { double l=t,h=u;bool sign=a<0;for(unsigned j=0;j<36;++j) { double m=(l+h)*.5;if((q(r.o+r.d*m)<0)==sign) l=m;else h=m; } double x=(l+h)*.5;v.push_back(x);continuous&=std::abs(q(r.o+r.d*x))<=1e-7; }t=u;a=b; }
    return v;
}
RayTransitions complete_roots(const PackedField* packed,const UniformGrid* grid,const Ray& r) {
    MeasurementRay ray{r.o,r.d,r.low,r.high};std::vector<RaySegment> segments;bool tiny=false;
    if(packed)segments=packed_ray_segments(packed->data(),ray);
    else {
        auto n=grid->resolution();auto times=grid_ray_boundaries(ray,n-1);
        for(std::size_t k=1;k<times.size();++k){auto low=times[k-1],high=times[k],mid=(low+high)*.5;
            tiny |= mid==low||mid==high;auto p=r.o+r.d*mid;std::array<std::size_t,3> index;
            for(unsigned a=0;a<3;++a)index[a]=std::min<std::size_t>(std::size_t(std::max(0.,(p[a]+1)*(n-1)/2)),n-2);
            Box box{Box{}.point({double(index[0])/(n-1),double(index[1])/(n-1),double(index[2])/(n-1)}),Box{}.point({double(index[0]+1)/(n-1),double(index[1]+1)/(n-1),double(index[2]+1)/(n-1)})};
            RaySegment segment{low,high,box,{}};for(unsigned c=0;c<8;++c)segment.corners[c]=grid->data()[index[0]+bool(c&1)+n*(index[1]+bool(c&2)+n*(index[2]+bool(c&4)))];segments.push_back(segment);
        }
    }
    auto result=segmented_transitions(segments,ray);result.ambiguous |= tiny;return result;
}
bool ray_signs(const Distance& q,const Ray& r,const std::vector<double>& roots) {
    if(!(q(r.o+r.d*r.low)>0&&q(r.o+r.d*r.high)>0))return false;
    std::vector<double> edges{r.low};edges.insert(edges.end(),roots.begin(),roots.end());edges.push_back(r.high);
    for(std::size_t j=1;j<edges.size();++j){auto mid=(edges[j-1]+edges[j])*.5;auto value=q(r.o+r.d*mid);if(!std::isfinite(value)||!(value*(j%2?-1.:1.)<0))return false;}
    return true;
}
Samples boundary_audit(const PackedField* packed,std::size_t grid_n,const Distance& truth) {
    std::vector<Vec3> points;
    auto faces=[&](Box b) {auto c=(b.lo+b.hi)*.5;for(int axis=0;axis<3;++axis)for(int side=0;side<2;++side){auto p=c;double* coord[3]={&p.x,&p.y,&p.z};*coord[axis]=side?b.hi[axis]:b.lo[axis];if(*coord[axis]<=-1||*coord[axis]>=1)continue;for(double e:{-1e-9,1e-9}){auto v=p;double* u[3]={&v.x,&v.y,&v.z};*u[axis]+=e;points.push_back(v);}}};
    if(packed) { std::function<void(std::size_t,Box)> visit=[&](std::size_t offset,Box b){std::uint32_t d;std::memcpy(&d,packed->data()+offset,4);if(d&3)faces(b);else for(unsigned c=0;c<8;++c)visit((d&~3u)+4*c,b.child(c));};visit(0,{}); }
    else {for(std::size_t k=1;k+1<grid_n;++k)for(int axis=0;axis<3;++axis)for(unsigned i=0;i<8;++i)for(unsigned j=0;j<8;++j){Vec3 p{};double* u[3]={&p.x,&p.y,&p.z};*u[axis]=-1+2.*k/(grid_n-1);*u[(axis+1)%3]=-1+2.*(i+.5)/8;*u[(axis+2)%3]=-1+2.*(j+.5)/8;for(double e:{-1e-9,1e-9}){auto v=p;double* q[3]={&v.x,&v.y,&v.z};*q[axis]+=e;points.push_back(v);}}}
    Samples s;const auto n=std::min<std::size_t>(4096,points.size());for(std::size_t i=0;i<n;++i){auto p=points[i*points.size()/n];s.push_back({p,truth(p)});}return s;
}
SurfaceTopology topology(std::ostream& f,const std::string& method,std::size_t cap,const Distance& q,unsigned n,const fs::path& dir) {
    auto start=Clock::now();auto mesh=extract_surface(q,n);auto t=surface_topology(mesh);
    f<<method<<','<<cap<<','<<n<<','<<t.closed_manifold()<<','<<t.components<<','<<t.euler<<','<<t.boundary_edges<<','<<t.nonmanifold_edges<<','<<t.nonmanifold_vertices<<',';
    for(std::size_t i=0;i<t.component_euler.size();++i) f<<(i?";":"")<<t.component_euler[i];f<<','<<elapsed(start)<<'\n';f.flush();
    mesh.write_obj((dir/(method+"_"+std::to_string(cap)+"_"+std::to_string(n)+".obj")).string());
    return t;
}
}
int main(int argc,char** argv) {
    try {
        if(argc<6) throw std::invalid_argument("Usage: qualify SPEC OUT COUNT SEED | fit SPEC REF OUT SEED METHOD [CAP ...]");
        const std::string mode=argv[1];const auto spec=GeometrySpec::read(argv[2]);fs::path cache_path;std::string cache_identity;for(int i=3;i+1<argc;++i){if(std::string(argv[i])=="--cache")throw std::invalid_argument("V2 forbids the persistent v1 oracle cache");}Memo oracle(spec,cache_path,cache_identity);std::size_t calls=0;Distance truth=[&](Vec3 p){++calls;return oracle.distance(p);};
        if(mode=="compare"){
            if(spec.shape!="mesh"||spec.b!=0)throw std::invalid_argument("Comparison must use EPECK geometry spec");fs::path dir=argv[3];fs::create_directories(dir);unsigned count=std::stoul(argv[4]),seed=std::stoul(argv[5]);auto fast=mesh_oracle_fast(spec.mesh);Samples samples;auto start=Clock::now();for(double near:{0.,1.}){auto s=sample_queries(oracle,count,seed+(near?100000:0),near,spec.noise);samples.insert(samples.end(),s.begin(),s.end());}std::ifstream points(dir/"comparison.points");Vec3 p;while(points>>p.x>>p.y>>p.z)samples.push_back({p,truth(p)});double exact_seconds=elapsed(start);std::vector<double> approximate;start=Clock::now();for(auto s:samples)approximate.push_back(fast->distance(s.point));double fast_seconds=elapsed(start);binary(dir/"epeck.bin",samples.data(),samples.size());binary(dir/"epick.bin",approximate.data(),approximate.size());std::mt19937_64 r1(seed+1),r2(seed+1);double sampler_delta=0;for(unsigned i=0;i<1024;++i)sampler_delta=std::max(sampler_delta,norm(oracle.surface(r1)-fast->surface(r2)));out(dir/"comparison.json")<<"{\"stratified_queries\":"<<2*count<<",\"total_queries\":"<<samples.size()<<",\"epeck_seconds\":"<<exact_seconds<<",\"epick_seconds\":"<<fast_seconds<<",\"sampler_max_delta\":"<<sampler_delta<<",\"epeck_input_validation\":true}";out(dir/"COMPLETE")<<"compared\n";return 0;
        }
        if(mode=="qualify") {
            const fs::path dir=argv[3];fs::create_directories(dir);unsigned count=std::stoul(argv[4]);unsigned seed=std::stoul(argv[5]);
            for(const auto& x:std::vector<std::pair<std::string,double>>{{"uniform",0},{"near",1}}) { auto s=sample_queries(oracle,count,seed+(x.second?100000:0),x.second,spec.noise);binary(dir/(x.first+".bin"),s.data(),s.size()); }
            for(const std::string name:{"probes","audit"}) { std::ifstream f(dir/(name+".points"));Samples s;Vec3 p;while(f>>p.x>>p.y>>p.z) { if(!Box{}.contains(p)) throw std::invalid_argument("Probe out of domain");s.push_back({p,truth(p)}); }binary(dir/(name+".bin"),s.data(),s.size()); }
            auto ray=out(dir/"reference_roots.csv");ray<<"ray,group,continuous,roots\n";auto rs=rays(dir/"rays.txt");
            std::vector<std::vector<double>> expected_roots;if(fs::exists(dir/"expected_ray_roots.csv")){std::ifstream f(dir/"expected_ray_roots.csv");std::string line;std::getline(f,line);while(std::getline(f,line)){auto comma=line.find(',');line=line.substr(comma+1);std::replace(line.begin(),line.end(),';',' ');std::istringstream row(line);double x;expected_roots.emplace_back();while(row>>x)expected_roots.back().push_back(x);}if(expected_roots.size()!=rs.size())throw std::invalid_argument("Bad independently intersected mesh ray roots");}
            for(std::size_t i=0;i<rs.size();++i) { bool c=true;auto v=expected_roots.empty()?roots(truth,rs[i],c):expected_roots[i];if(!expected_roots.empty()){for(auto t:v)c&=std::abs(truth(rs[i].o+rs[i].d*t))<1e-7;std::vector<double> edges{rs[i].low};edges.insert(edges.end(),v.begin(),v.end());edges.push_back(rs[i].high);for(std::size_t j=0;j+1<edges.size();++j){auto d=truth(rs[i].o+rs[i].d*((edges[j]+edges[j+1])*.5));if((d<0)!=(j%2==1)||std::abs(d)<=1e-9)throw std::invalid_argument("Triangle intersections disagree with reference signed intervals");}}
                ray<<i<<','<<rs[i].group<<','<<c<<',';for(std::size_t j=0;j<v.size();++j) ray<<(j?";":"")<<v[j];ray<<'\n'; }
            auto surfaces=out(dir/"reference_topology.csv");surfaces<<"method,cap,cells,closed,components,euler,boundary_edges,nonmanifold_edges,nonmanifold_vertices,component_euler,seconds\n";
            std::vector<long long> expected;{std::ifstream f(dir/"expected_euler.txt");long long e;while(f>>e)expected.push_back(e);}
            std::ifstream nf(dir/"levels.txt");unsigned n;bool prior=false;while(nf>>n) {
                const unsigned m=n+1;std::vector<double> lattice(std::size_t(m)*m*m);for(std::size_t i=0;i<lattice.size();++i) lattice[i]=truth(Box{}.point({double(i%m)/n,double((i/m)%m)/n,double(i/(m*m))/n}));
                binary(dir/("lattice_"+std::to_string(n)+".bin"),lattice.data(),lattice.size());
                Distance cached=[&](Vec3 p) {auto x=std::size_t(std::llround((p.x+1)*n/2)),y=std::size_t(std::llround((p.y+1)*n/2)),z=std::size_t(std::llround((p.z+1)*n/2));return lattice[x+m*(y+m*z)];};
                auto t=topology(surfaces,"oracle",0,cached,n,dir);bool valid=t.closed_manifold()&&t.component_euler==expected;
                std::cout<<spec.id<<" reference "<<n<<std::endl;
                if(valid&&prior)break;prior=valid;
            }
            out(dir/"oracle.json")<<"{\"sample_calls\":"<<2*count<<",\"other_calls\":"<<calls<<",\"oracle\":"<<std::quoted(oracle.name())<<"}";
            out(dir/"QUALIFIED")<<"reference only\n";return 0;
        }
        if(mode!="fit"||argc<8) throw std::invalid_argument("Invalid fit command");
        const fs::path ref=argv[3],dir=argv[4];if(fs::exists(dir)) throw std::invalid_argument("Use new job directory");fs::create_directories(dir);unsigned seed=std::stoul(argv[5]);std::string method=argv[6];std::vector<std::size_t> caps;std::size_t policy_k=128,policy_k7=0;double policy_lambda=.001;unsigned global_quadrature=0,covered_ranking=0;bool probe_max=false;
        for(int i=7;i<argc;++i){if(std::string(argv[i])=="--policy"&&i+1<argc){std::ifstream f(argv[++i]);if(!(f>>policy_k>>policy_lambda>>global_quadrature>>covered_ranking))throw std::invalid_argument("Incomplete policy");if(!(f>>policy_k7))policy_k7=0;if(policy_k7>32768||!policy_k||policy_k>32768||policy_lambda<0||!std::isfinite(policy_lambda)||global_quadrature>2||covered_ranking>1||(covered_ranking&&!global_quadrature))throw std::invalid_argument("Invalid policy");}else if((std::string(argv[i])=="--cache"||std::string(argv[i])=="--cache-key")&&i+1<argc)++i;else if(std::string(argv[i])=="--probe-max")probe_max=true;else caps.push_back(std::stoull(argv[i]));}
        if(caps.empty())throw std::invalid_argument("No caps");const bool is_candidate=method=="candidate"||method=="candidate_original"||method=="coverage_q8";
        if(method!="grid"&&method!="mse"&&method!="sign"&&method!="cvar99"&&method!="selective32"&&!is_candidate) throw std::invalid_argument("Unknown comparator");
        if(method=="candidate_original"){policy_k=128;policy_lambda=.001;global_quadrature=0;covered_ranking=0;policy_k7=0;}
        if(method=="coverage_q8"){policy_k=4096;policy_lambda=0;global_quadrature=1;covered_ranking=1;policy_k7=0;}
        out(dir/"policy.json")<<"{\"k\":"<<policy_k<<",\"k7\":"<<policy_k7<<",\"lambda\":"<<policy_lambda<<",\"quadrature\":4,\"global_quadrature\":"<<global_quadrature<<",\"covered_ranking\":"<<covered_ranking<<",\"probe_max_gate\":"<<(probe_max?"true":"false")<<"}";
        out(dir/"abi.json")<<"{\"packed_object_bytes\":"<<sizeof(PackedField)<<",\"grid_object_bytes\":"<<sizeof(UniformGrid)<<",\"sample_bytes\":"<<sizeof(Sample)<<",\"configuration\":"<<std::quoted(SDF_BUILD_CONFIGURATION)<<"}";
        std::vector<std::pair<std::string,Samples>> datasets;for(const std::string s:{"uniform","near","audit","probes"}) datasets.push_back({s,load(ref/(s+".bin"))});
        std::vector<int> groups; { std::ifstream f(ref/"probe_groups.txt");int g;while(f>>g)groups.push_back(g); }if(groups.size()!=datasets.back().second.size())throw std::runtime_error("Bad groups");
        auto rs=rays(ref/"rays.txt");std::vector<std::vector<double>> rr; {std::ifstream f(ref/"reference_roots.csv");std::string l;std::getline(f,l);while(std::getline(f,l)){auto a=l.find_last_of(',');std::replace(l.begin()+a+1,l.end(),';',' ');std::istringstream v(l.substr(a+1));double x;rr.emplace_back();while(v>>x)rr.back().push_back(x);} }
        std::vector<unsigned> levels;{ std::ifstream f(ref/"qualified_levels.txt");unsigned n;while(f>>n)levels.push_back(n); }
        auto summary=out(dir/"fields.csv");summary<<"method,seed,cap,bytes,grid_n,scalar_pass,probe_pass,width_pass,topology_status,build_tree_seconds,loss_seconds,dp_seconds,pack_seconds,oracle_calls,physical_build_oracle_calls,boundary_oracle_calls,physical_boundary_oracle_calls,boundary_label_seconds,fit_count,nodes,selected_parents,selected_d6_parents,dp_states,dp_transitions,objective\n";
        auto metric=out(dir/"metrics.csv");metric<<"cap,dataset,count,mae,p99,cvar99,max_error,inside_recall,outside_recall,sign_error\n";
        auto pg=out(dir/"groups.csv");pg<<"cap,group,inside_count,outside_count,inside_recall,outside_recall,pass\n";
        auto widths=out(dir/"widths.csv");widths<<"cap,ray,group,continuous,interval_count_ok,reference_roots,predicted_roots,mean_relative_error,p95_relative_error,pass,numerical_unresolved,signed_intervals_ok,tangencies,jumps,zero_segments\n";
        auto surface=out(dir/"topology.csv");surface<<"method,cap,cells,closed,components,euler,boundary_edges,nonmanifold_edges,nonmanifold_vertices,component_euler,seconds\n";
        auto timing=out(dir/"timing.csv");timing<<"cap,repeat,ns_per_query,checksum\n";
        auto start=Clock::now();Samples fit;if(method!="grid") { fit=sample_queries(oracle,10000,seed,.5,spec.noise);binary(dir/"fit.bin",fit.data(),fit.size()); }const double label_time=elapsed(start);out(dir/"label_seconds.txt")<<label_time;
        std::unique_ptr<PrecisionTree> tree;std::unique_ptr<Q8Losses> losses;std::unique_ptr<PrecisionDP> dp;double tree_time=0,loss_time=0,dp_time=0;std::size_t selected=0,selected7=0;
        if(method!="grid") { start=Clock::now();Samples global;
            if(is_candidate&&global_quadrature){for(unsigned iz=0;iz<32;++iz)for(unsigned iy=0;iy<32;++iy)for(unsigned ix=0;ix<32;++ix)for(unsigned z=0;z<global_quadrature;++z)for(unsigned y=0;y<global_quadrature;++y)for(unsigned x=0;x<global_quadrature;++x){auto p=Vec3{-1+(ix+(x+.5)/global_quadrature)/16,-1+(iy+(y+.5)/global_quadrature)/16,-1+(iz+(z+.5)/global_quadrature)/16};global.push_back({p,truth(p)});}binary(dir/"global.bin",global.data(),global.size());if(covered_ranking)fit.insert(fit.end(),global.begin(),global.end());}
            tree=std::make_unique<PrecisionTree>(truth,fit,5);
            {std::vector<float> corners;std::vector<std::uint32_t> children;for(auto v:tree->spatial.nodes){corners.insert(corners.end(),v.field.values.begin(),v.field.values.end());children.push_back(v.field.first_child);}binary(dir/"depth5_corners.bin",corners.data(),corners.size());binary(dir/"depth5_children.bin",children.data(),children.size());}
            if(method=="selective32"||is_candidate) {auto r=tree->refine_depth_six(truth,fit,is_candidate?policy_k:32);selected=r.size();auto selection=out(dir/"selection.csv");selection<<"rank,node,fit_count,q8_loss,lo_x,lo_y,lo_z,hi_x,hi_y,hi_z\n";unsigned rank=0;for(auto v:r)selection<<++rank<<','<<v.node<<','<<v.fitting_count<<','<<v.q8_loss<<','<<v.box.lo.x<<','<<v.box.lo.y<<','<<v.box.lo.z<<','<<v.box.hi.x<<','<<v.box.hi.y<<','<<v.box.hi.z<<'\n';
                if(is_candidate) { Samples extra;for(const auto& v:r) for(unsigned z=0;z<4;++z)for(unsigned y=0;y<4;++y)for(unsigned x=0;x<4;++x){auto p=v.box.point({(x+.5)/4,(y+.5)/4,(z+.5)/4});extra.push_back({p,truth(p)});} binary(dir/"extra.bin",extra.data(),extra.size());fit.insert(fit.end(),extra.begin(),extra.end());if(!covered_ranking)fit.insert(fit.end(),global.begin(),global.end()); }
            }
            if(is_candidate&&policy_k7){
                {std::vector<float> corners;std::vector<std::uint32_t> children;for(auto v:tree->spatial.nodes){corners.insert(corners.end(),v.field.values.begin(),v.field.values.end());children.push_back(v.field.first_child);}binary(dir/"depth6_corners.bin",corners.data(),corners.size());binary(dir/"depth6_children.bin",children.data(),children.size());}
                auto selected_d6=tree->refine_depth_seven(truth,fit,policy_k7);selected7=selected_d6.size();auto selection=out(dir/"selection7.csv");selection<<"rank,node,fit_count,q8_loss,lo_x,lo_y,lo_z,hi_x,hi_y,hi_z\n";unsigned rank=0;Samples extra;for(auto v:selected_d6){selection<<++rank<<','<<v.node<<','<<v.fitting_count<<','<<v.q8_loss<<','<<v.box.lo.x<<','<<v.box.lo.y<<','<<v.box.lo.z<<','<<v.box.hi.x<<','<<v.box.hi.y<<','<<v.box.hi.z<<'\n';for(unsigned z=0;z<4;++z)for(unsigned y=0;y<4;++y)for(unsigned x=0;x<4;++x){auto p=v.box.point({(x+.5)/4,(y+.5)/4,(z+.5)/4});extra.push_back({p,truth(p)});}}binary(dir/"extra7.bin",extra.data(),extra.size());fit.insert(fit.end(),extra.begin(),extra.end());
            }tree_time=elapsed(start);
            start=Clock::now();losses=std::make_unique<Q8Losses>(*tree,fit);losses->apply(*tree,method=="sign"||is_candidate?Q8Objective::sign:Q8Objective::mse,is_candidate?policy_lambda:method=="sign"?.0025:0);loss_time=elapsed(start);
            if(method!="cvar99") { start=Clock::now();dp=std::make_unique<PrecisionDP>(*tree,*std::max_element(caps.begin(),caps.end()),Encoding::q8);dp_time=elapsed(start); }
        }
        const auto amortized_calls=calls+10000,amortized_physical=oracle.physical;
        for(auto cap:caps) {
            std::size_t grid_n=0,bytes=0,states=0,transitions=0;double objective=0,allocation_time=0,pack_time=0;std::unique_ptr<PackedField> packed;std::unique_ptr<UniformGrid> grid;Distance q;auto before=calls,physical_before=oracle.physical;
            if(method=="grid") {start=Clock::now();grid=std::make_unique<UniformGrid>(truth,cap);pack_time=elapsed(start);grid_n=grid->resolution();bytes=grid->bytes();binary(dir/(std::to_string(cap)+".bin"),grid->data(),grid_n*grid_n*grid_n);q=[&](Vec3 p){return grid->query(p);};}
            else {start=Clock::now();PrecisionAllocation a;if(method=="cvar99") {auto tail=allocate_cvar99(*tree,*losses,cap,9);a=tail.allocation;auto f=out(dir/(std::to_string(cap)+"_thresholds.csv"));f<<"threshold,objective,bytes,states,transitions\n";for(auto s:tail.search){f<<s.threshold<<','<<s.objective<<','<<s.bytes<<','<<s.states<<','<<s.transitions<<'\n';states+=s.states;transitions+=s.transitions;}}else {a=dp->allocate(cap);states=dp->states;transitions=dp->transitions;}allocation_time=elapsed(start);objective=a.loss;start=Clock::now();packed=std::make_unique<PackedField>(*tree,a);pack_time=elapsed(start);bytes=packed->bytes();binary(dir/(std::to_string(cap)+".bin"),packed->data(),packed->payload_bytes());q=[&](Vec3 p){return packed->query(p);}; }
            const auto build_calls=method=="grid"?calls-before:amortized_calls,physical_build_calls=method=="grid"?oracle.physical-physical_before:amortized_physical;
            if(bytes>cap)throw std::logic_error("Byte overflow");bool scalar=true;for(const auto& d:datasets) {auto m=evaluate(q,d.second);metric<<cap<<','<<d.first<<','<<m.count<<','<<m.mae<<','<<m.p99<<','<<m.cvar99<<','<<m.max_error<<','<<m.inside_recall<<','<<m.outside_recall<<','<<m.sign_error<<'\n';if(d.first!="probes")scalar&=m.max_error<=.020&&(d.first!="uniform"||m.cvar99<=.015)&&(d.first!="near"||m.cvar99<=.010);std::vector<double> preds;preds.reserve(d.second.size());for(auto s:d.second)preds.push_back(q(s.point));binary(dir/(std::to_string(cap)+"_"+d.first+".pred"),preds.data(),preds.size()); }
            if(probe_max)scalar&=evaluate(q,datasets.back().second).max_error<=.020;
            const auto audit_before=calls,audit_physical_before=oracle.physical;start=Clock::now();auto boundary=boundary_audit(packed.get(),grid_n,truth);const double audit_label_time=elapsed(start);const auto audit_calls=calls-audit_before,audit_physical=oracle.physical-audit_physical_before;
            binary(dir/(std::to_string(cap)+"_boundary.bin"),boundary.data(),boundary.size());auto bm=evaluate(q,boundary);scalar&=bm.max_error<=.020;metric<<cap<<",boundary,"<<bm.count<<','<<bm.mae<<','<<bm.p99<<','<<bm.cvar99<<','<<bm.max_error<<','<<bm.inside_recall<<','<<bm.outside_recall<<','<<bm.sign_error<<'\n';std::vector<double> bp;for(auto s:boundary)bp.push_back(q(s.point));binary(dir/(std::to_string(cap)+"_boundary.pred"),bp.data(),bp.size());
            bool probes=true;for(int g=0;g<=*std::max_element(groups.begin(),groups.end());++g) {Samples s;for(std::size_t i=0;i<groups.size();++i)if(groups[i]==g)s.push_back(datasets.back().second[i]);if(s.empty())continue;auto m=evaluate(q,s);bool pass=m.inside_count&&m.outside_count&&m.inside_recall>=.95&&m.outside_recall>=.95;probes&=pass;pg<<cap<<','<<g<<','<<m.inside_count<<','<<m.outside_count<<','<<m.inside_recall<<','<<m.outside_recall<<','<<pass<<'\n'; }
            bool width=true;std::map<int,std::vector<double>> group_errors;
            for(std::size_t i=0;i<rs.size();++i) {
                auto events=complete_roots(packed.get(),grid.get(),rs[i]);auto v=events.crossings;
                bool cont=true;for(auto root:v)cont &= std::abs(q(rs[i].o+rs[i].d*root))<=1e-7;
                bool ok=v.size()==rr[i].size()&&!v.empty()&&v.size()%2==0;
                bool signed_ok=ray_signs(q,rs[i],rr[i])&&ray_signs(q,rs[i],v);
                std::vector<double> err;
                if(ok)for(std::size_t j=0;j+1<v.size();++j){double w=rr[i][j+1]-rr[i][j];if(w<=0)ok=false;else err.push_back(std::abs((v[j+1]-v[j])/w-1));}
                double mean=std::numeric_limits<double>::infinity(),p95=mean;
                if(ok){std::sort(err.begin(),err.end());mean=std::accumulate(err.begin(),err.end(),0.)/err.size();p95=err[std::size_t(std::ceil(.95*err.size()))-1];auto& errors=group_errors[rs[i].group];errors.insert(errors.end(),err.begin(),err.end());}
                bool pass=ok&&signed_ok&&!events.ambiguous;width &= pass;
                widths<<cap<<','<<i<<','<<rs[i].group<<','<<cont<<','<<ok<<',';
                for(std::size_t j=0;j<rr[i].size();++j)widths<<(j?";":"")<<rr[i][j];widths<<',';
                for(std::size_t j=0;j<v.size();++j)widths<<(j?";":"")<<v[j];
                widths<<','<<mean<<','<<p95<<','<<pass<<','<<events.ambiguous<<','<<signed_ok<<','<<events.tangencies.size()<<','<<events.jumps.size()<<','<<events.zero_segments<<'\n';
            }
            for(const auto& r:rs)if(group_errors[r.group].empty())width=false;
            for(auto& entry:group_errors){auto& errors=entry.second;if(errors.empty()){width=false;continue;}std::sort(errors.begin(),errors.end());width&=std::accumulate(errors.begin(),errors.end(),0.)/errors.size()<=.10&&errors[std::size_t(std::ceil(.95*errors.size()))-1]<=.15;}
            std::string status="not_evaluated_due_to_prior_failure";if(scalar&&probes&&width) {status=levels.size()==2?"evaluated":"oracle_resolution_unresolved";for(auto n:levels)topology(surface,method,cap,q,n,dir);}
            if(!(scalar&&probes&&width)&&method=="candidate"&&seed==240&&cap==32768){auto diagnostic=extract_surface(q,64);diagnostic.write_obj((dir/"candidate_32768_diagnostic64.obj").string());}
            volatile double warm=0;for(auto s:datasets[0].second)warm=warm+q(s.point);(void)warm;
            for(unsigned rep=0;rep<9;++rep){double sum=0;start=Clock::now();for(auto s:datasets[0].second)sum+=q(s.point);timing<<cap<<','<<rep<<','<<elapsed(start)*1e9/datasets[0].second.size()<<','<<sum<<'\n';}
            summary<<method<<','<<seed<<','<<cap<<','<<bytes<<','<<grid_n<<','<<scalar<<','<<probes<<','<<width<<','<<status<<','<<tree_time<<','<<loss_time<<','<<(allocation_time+(method=="cvar99"?0:dp_time))<<','<<pack_time<<','<<build_calls<<','<<physical_build_calls<<','<<audit_calls<<','<<audit_physical<<','<<audit_label_time<<','<<fit.size()<<','<<(tree?tree->nodes.size():0)<<','<<selected<<','<<selected7<<','<<states<<','<<transitions<<','<<objective<<'\n';summary.flush();metric.flush();pg.flush();widths.flush();timing.flush();std::cout<<spec.id<<' '<<seed<<' '<<method<<' '<<cap<<' '<<bytes<<' '<<scalar<<probes<<width<<std::endl;
        }
        oracle.save();out(dir/"cache.json")<<"{\"path\":"<<std::quoted(cache_path.string())<<",\"loaded_coordinates\":"<<oracle.loaded<<",\"load_seconds\":"<<oracle.load_seconds<<",\"logical_requests\":"<<oracle.logical<<",\"cache_hits\":"<<oracle.hits<<",\"physical_oracle_calls\":"<<oracle.physical<<"}";
        out(dir/"COMPLETE.tmp")<<"complete\n";fs::rename(dir/"COMPLETE.tmp",dir/"COMPLETE");
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
