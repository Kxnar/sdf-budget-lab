#include "sdf/metrics.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <sstream>
#ifdef SDF_WITH_CGAL
#include <CGAL/version.h>
#endif

namespace {
using namespace sdf;
namespace fs=std::filesystem;
struct Config {
    std::string shape="torus",mesh,probes,out="runs/first",plane;
    double feature=0.04,penalty=0.0025;
    std::size_t train=8192,test=8192,slice_size=257;
    unsigned depth=10,repeats=11;
    std::vector<std::size_t> budgets{8192,32768,131072},seeds{0,1,2};
};
std::size_t integer(const std::string& s) {
    if(s.empty() || s.find_first_not_of("0123456789")!=std::string::npos) throw std::invalid_argument("Expected nonnegative integer: "+s);
    std::size_t consumed=0; const auto v=std::stoull(s,&consumed);
    if(consumed!=s.size() || v>std::numeric_limits<std::size_t>::max()) throw std::invalid_argument("Integer out of range");
    return static_cast<std::size_t>(v);
}
double real(const std::string& s) {
    std::size_t n=0; const double v=std::stod(s,&n);
    if(n!=s.size() || !std::isfinite(v)) throw std::invalid_argument("Expected finite number: "+s);
    return v;
}
std::vector<std::size_t> list(const std::string& s) {
    std::vector<std::size_t> result; std::istringstream input(s); std::string part;
    while(std::getline(input,part,',')) result.push_back(integer(part));
    if(result.empty() || s.back()==',') throw std::invalid_argument("Empty list entry");
    auto sorted=result; std::sort(sorted.begin(),sorted.end());
    if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end()) throw std::invalid_argument("Duplicate list entry");
    return result;
}
Config parse(int argc,char** argv) {
    Config c;
    for(int i=1;i<argc;++i) {
        const std::string key=argv[i];
        if(key=="--help") {
            std::cout<<"sdf_bench [--shape sphere|torus|thin_box|double_sphere] [--feature 0.04]\n"
                "  [--mesh closed.off] [--probes xyz.csv] [--budgets 8192,32768,131072]\n"
                "  [--seeds 0,1,2] [--train 8192] [--test 8192] [--max-depth 10]\n"
                "  [--sign-penalty 0.0025] [--repeats 11] [--slice-size 257] [--slice-plane xy|xz] [--out NEW_DIR]\n";
            std::exit(0);
        }
        if(i+1==argc) throw std::invalid_argument("Missing value for "+key);
        const std::string value=argv[++i];
        if(key=="--shape") c.shape=value;
        else if(key=="--feature") c.feature=real(value);
        else if(key=="--mesh") c.mesh=value;
        else if(key=="--probes") c.probes=value;
        else if(key=="--out") c.out=value;
        else if(key=="--slice-plane") c.plane=value;
        else if(key=="--budgets") c.budgets=list(value);
        else if(key=="--seeds") c.seeds=list(value);
        else if(key=="--train") c.train=integer(value);
        else if(key=="--test") c.test=integer(value);
        else if(key=="--slice-size") c.slice_size=integer(value);
        else if(key=="--max-depth") { const auto v=integer(value); if(v>20) throw std::invalid_argument("Max depth > 20"); c.depth=static_cast<unsigned>(v); }
        else if(key=="--repeats") { const auto v=integer(value); if(v>1000) throw std::invalid_argument("Too many repeats"); c.repeats=static_cast<unsigned>(v); }
        else if(key=="--sign-penalty") c.penalty=real(value);
        else throw std::invalid_argument("Unknown option: "+key);
    }
    if(!c.train || !c.test || !c.depth || c.repeats<3 || c.slice_size<3 || c.slice_size>1025 || c.penalty<0 || c.feature<0.005 || c.feature>0.25)
        throw std::invalid_argument("Invalid benchmark configuration");
    for(auto seed:c.seeds) if(seed>=10000) throw std::invalid_argument("Seeds must be <10000; evaluation uses disjoint seed ranges");
    for(auto budget:c.budgets) if(budget<512) throw std::invalid_argument("Minimum budget is 512 bytes");
    if(c.plane.empty()) c.plane=c.mesh.empty() && c.shape=="torus"?"xy":"xz";
    if(c.plane!="xy" && c.plane!="xz") throw std::invalid_argument("Slice plane must be xy or xz");
    return c;
}
std::ofstream output(const fs::path& path) {
    std::ofstream file(path); file.exceptions(std::ios::badbit|std::ios::failbit); file<<std::setprecision(17); return file;
}
double elapsed(std::chrono::steady_clock::time_point start) { return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); }
void array_json(std::ostream& f,const std::vector<std::size_t>& values) {
    f<<'['; for(std::size_t i=0;i<values.size();++i) { if(i) f<<','; f<<values[i]; } f<<']';
}
void config(const Config& c,const Oracle& oracle,const fs::path& dir) {
    auto f=output(dir/"config.json");
    f<<"{\n  \"schema\": 1,\n  \"shape\": "<<std::quoted(c.mesh.empty()?c.shape:"mesh")
     <<",\n  \"oracle\": "<<std::quoted(oracle.name())<<",\n  \"mesh\": "<<std::quoted(c.mesh)
     <<",\n  \"probes\": "<<std::quoted(c.probes)<<",\n  \"feature\": "<<c.feature
     <<",\n  \"training_samples\": "<<c.train<<",\n  \"test_samples_per_distribution\": "<<c.test
     <<",\n  \"train_surface_fraction\": 0.5,\n  \"noise_sigma\": "<<c.feature/2
     <<",\n  \"sign_penalty\": "<<c.penalty<<",\n  \"max_depth\": "<<c.depth
     <<",\n  \"timing_repeats\": "<<c.repeats<<",\n  \"slice_size\": "<<c.slice_size
     <<",\n  \"slice_plane\": "<<std::quoted(c.plane)
     <<",\n  \"budgets\": "; array_json(f,c.budgets); f<<",\n  \"seeds\": "; array_json(f,c.seeds);
    f<<",\n  \"seed_streams\": {\"train\":0,\"uniform\":10000,\"near\":20000,\"shift\":30000}"
     <<",\n  \"node_bytes\": "<<sizeof(AdaptiveField::Node)<<",\n  \"grid_object_bytes\": "<<sizeof(UniformGrid)
     <<",\n  \"adaptive_object_bytes\": "<<sizeof(AdaptiveField)<<",\n  \"build_configuration\": "<<std::quoted(SDF_BUILD_CONFIGURATION);
#ifdef _MSC_VER
    f<<",\n  \"compiler\": \"MSVC "<<_MSC_FULL_VER<<'"';
#else
    f<<",\n  \"compiler\": "<<std::quoted(__VERSION__);
#endif
#ifdef SDF_WITH_CGAL
    f<<",\n  \"cgal_version\": "<<std::quoted(CGAL_VERSION_STR);
#endif
    f<<"\n}\n";
}
void row(std::ostream& f,const std::string& method,const std::string& dataset,std::size_t seed,std::size_t cap,
         std::size_t bytes,std::size_t nodes,double build,std::size_t calls,double label_seconds,
         const Metrics& m,const Timing& t) {
    f<<method<<','<<dataset<<','<<seed<<','<<cap<<','<<bytes<<','<<nodes<<','<<build<<','<<calls<<','<<label_seconds
     <<','<<m.count<<','<<m.rmse<<','<<m.mae<<','<<m.p99<<','<<m.cvar99<<','<<m.max_error<<','<<m.sign_error
     <<','<<m.inside_recall<<','<<m.outside_recall<<','<<m.inside_count<<','<<m.outside_count
     <<','<<t.median_ns<<','<<t.p95_batch_ns<<','<<t.checksum<<'\n';
}
void slice(const fs::path& dir,std::size_t cap,const Distance& truth,const UniformGrid& grid,
           const AdaptiveField& mse,const AdaptiveField& risk,std::size_t n,const std::string& plane) {
    const bool xy=plane=="xy";
    auto f=output(dir/("slice_"+std::to_string(cap)+".csv")); f<<(xy?"x,y":"x,z")<<",reference,uniform,adaptive_mse,adaptive_risk\n";
    for(std::size_t z=0;z<n;++z) for(std::size_t x=0;x<n;++x) {
        const double v=-1+2.0*z/(n-1);
        const Vec3 p{-1+2.0*x/(n-1),xy?v:0,xy?0:v};
        f<<p.x<<','<<v<<','<<truth(p)<<','<<grid.query(p)<<','<<mse.query(p)<<','<<risk.query(p)<<'\n';
    }
    auto leaves=output(dir/("leaves_"+std::to_string(cap)+".csv")); leaves<<"x0,z0,x1,z1\n";
    // Query uses the positive half at the fixed zero coordinate.
    for(const auto& b:risk.leaves()) if((xy?b.lo.z:b.lo.y)<=0 && (xy?b.hi.z:b.hi.y)>0)
        leaves<<b.lo.x<<','<<(xy?b.lo.y:b.lo.z)<<','<<b.hi.x<<','<<(xy?b.hi.y:b.hi.z)<<'\n';
}
int run(const Config& c) {
    auto oracle=c.mesh.empty()?analytic_oracle(c.shape,c.feature):mesh_oracle(c.mesh);
    const fs::path dir=c.out;
    if(fs::exists(dir)) throw std::runtime_error("Output already exists; choose a new directory");
    if(!fs::create_directories(dir)) throw std::runtime_error("Could not create output directory");
    config(c,*oracle,dir);
    auto csv=output(dir/"metrics.csv");
    csv<<"method,dataset,seed,budget_bytes,representation_bytes,nodes,build_seconds,build_oracle_queries,training_label_seconds,count,rmse,mae,p99,cvar99,max_error,sign_error,inside_recall,outside_recall,inside_count,outside_count,median_ns_per_query,p95_batch_ns_per_query,checksum\n";
    const Distance reference=[&](Vec3 p){ return oracle->distance(p); };
    for(const auto seed:c.seeds) {
        const auto start=std::chrono::steady_clock::now();
        const auto training=sample_queries(*oracle,c.train,seed,0.5,c.feature/2);
        const double label_seconds=elapsed(start);
        std::vector<std::pair<std::string,Samples>> test;
        test.push_back({"uniform",sample_queries(*oracle,c.test,10000+seed,0,c.feature/2)});
        test.push_back({"near",sample_queries(*oracle,c.test,20000+seed,1,c.feature/2)});
        test.push_back({"shift",sample_queries(*oracle,c.test,30000+seed,0.9,2*c.feature)});
        if(!c.probes.empty()) test.push_back({"probes",read_probes(*oracle,c.probes)});
        else if(c.mesh.empty()) test.push_back({"probes",feature_probes(*oracle,c.shape,c.feature)});
        const auto oracle_timing=time_queries(reference,test[0].second,c.repeats);
        for(const auto& data:test) row(csv,"oracle",data.first,seed,0,0,0,0,0,0,evaluate(reference,data.second),oracle_timing);
        for(const auto cap:c.budgets) {
            std::size_t calls=0;
            const Distance counted=[&](Vec3 p){ ++calls; return reference(p); };
            auto begin=std::chrono::steady_clock::now(); UniformGrid grid(counted,cap);
            const double grid_build=elapsed(begin); const auto grid_calls=calls;
            calls=0; begin=std::chrono::steady_clock::now(); AdaptiveField mse(counted,training,cap,{c.depth,0});
            const double mse_build=elapsed(begin); const auto mse_calls=calls;
            calls=0; begin=std::chrono::steady_clock::now(); AdaptiveField risk(counted,training,cap,{c.depth,c.penalty});
            const double risk_build=elapsed(begin); const auto risk_calls=calls;
            auto record=[&](const std::string& method,const Distance& f,std::size_t bytes,std::size_t nodes,double build,std::size_t build_calls) {
                const auto timing=time_queries(f,test[0].second,c.repeats);
                for(const auto& data:test) row(csv,method,data.first,seed,cap,bytes,nodes,build,build_calls,
                    method=="uniform"?0:label_seconds,evaluate(f,data.second),timing);
            };
            record("uniform",[&](Vec3 p){return grid.query(p);},grid.bytes(),grid.resolution()*grid.resolution()*grid.resolution(),grid_build,grid_calls);
            record("adaptive_mse",[&](Vec3 p){return mse.query(p);},mse.bytes(),mse.node_count(),mse_build,mse_calls);
            record("adaptive_risk",[&](Vec3 p){return risk.query(p);},risk.bytes(),risk.node_count(),risk_build,risk_calls);
            if(seed==c.seeds.front()) slice(dir,cap,reference,grid,mse,risk,c.slice_size,c.plane);
            csv.flush(); std::cout<<oracle->name()<<" seed="<<seed<<" budget="<<cap<<" complete\n";
        }
    }
    auto done=output(dir/"complete.json"); done<<"{\"status\":\"complete\"}\n";
    return 0;
}
}
int main(int argc,char** argv) {
    try { return run(parse(argc,argv)); }
    catch(const std::exception& e) { std::cerr<<"sdf_bench: "<<e.what()<<'\n'; return 1; }
}
