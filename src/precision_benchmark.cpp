#include "sdf/precision.hpp"
#include "sdf/surface.hpp"
#include "sdf/metrics.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace {
using namespace sdf;
namespace fs=std::filesystem;
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
std::ofstream output(const fs::path& p) { std::ofstream f(p); f.exceptions(std::ios::badbit|std::ios::failbit); f<<std::setprecision(17); return f; }
const Vec3 shift{.017,-.023,.011};
class Shifted final:public Oracle {
    std::unique_ptr<Oracle> base_;
public:
    Shifted(const std::string& shape,double feature):base_(analytic_oracle(shape,feature)) {}
    double distance(Vec3 p) const override { return base_->distance(p-shift); }
    Vec3 surface(std::mt19937_64& rng) const override { return base_->surface(rng)+shift; }
    std::string name() const override { return base_->name(); }
};
void save(const fs::path& path,const Samples& s) {
    auto f=output(path); f<<"x,y,z,distance\n";
    for(const auto& p:s) f<<p.point.x<<','<<p.point.y<<','<<p.point.z<<','<<p.distance<<'\n';
}
bool topology_pass(const std::string& shape,const SurfaceTopology& t) {
    const std::vector<long long> expected=shape=="torus"?std::vector<long long>{0}:shape=="double_sphere"?std::vector<long long>{2,2}:std::vector<long long>{2};
    return t.closed_manifold() && t.component_euler==expected;
}
}
int main(int argc,char** argv) {
    try {
        std::string out; bool smoke=false;
        for(int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--smoke") smoke=true;
            else if(arg=="--out" && i+1<argc) out=argv[++i];
            else throw std::invalid_argument("Usage: sdf_precision_bench --out NEW_DIR [--smoke]");
        }
        if(out.empty() || fs::exists(out)) throw std::invalid_argument("Output must be a new directory");
        const fs::path dir(out);
        for(const auto* p:{"samples","candidates","fields","meshes","choices"}) fs::create_directories(dir/p);
        const unsigned depth=smoke?2:4,nfit=smoke?256:4096,ntest=smoke?1024:32768,nseeds=smoke?1:2;
        const std::vector<std::string> shapes=smoke?std::vector<std::string>{"torus"}:std::vector<std::string>{"thin_box","torus","double_sphere"};
        const std::vector<double> features=smoke?std::vector<double>{.06}:std::vector<double>{.02,.035,.06,.10};
        const std::vector<unsigned> budgets=smoke?std::vector<unsigned>{2048}:std::vector<unsigned>{8192,16384,32768};
        const std::vector<unsigned> resolutions=smoke?std::vector<unsigned>{16,24}:std::vector<unsigned>{32,64,96};
        {
            auto f=output(dir/"config.json");
            f<<"{\"protocol\":\"precision-v1\",\"smoke\":"<<(smoke?"true":"false")<<",\"depth\":"<<depth<<",\"fit\":"<<nfit<<",\"test\":"<<ntest<<",\"seeds\":"<<nseeds<<",\"field_object_bytes\":"<<sizeof(PackedField)<<",\"configuration\":"<<std::quoted(SDF_BUILD_CONFIGURATION)<<",\"compiler\":";
#ifdef _MSC_VER
            f<<std::quoted("MSVC "+std::to_string(_MSC_FULL_VER));
#else
            f<<std::quoted(__VERSION__);
#endif
            f<<",\"budgets\":[";
            for(std::size_t i=0;i<budgets.size();++i) f<<(i?",":"")<<budgets[i];
            f<<"],\"resolutions\":[";
            for(std::size_t i=0;i<resolutions.size();++i) f<<(i?",":"")<<resolutions[i];
            f<<"],\"features\":[";
            for(std::size_t i=0;i<features.size();++i) f<<(i?",":"")<<features[i];
            f<<"]}\n";
        }
        auto metrics=output(dir/"metrics.csv");
        metrics<<"case,shape,feature,seed,budget,method,bytes,branches,f32_leaves,q16_leaves,q8_leaves,fit_sse,dataset,rmse,mae,p99,cvar99,inside_recall,outside_recall,sign_error\n";
        auto construction=output(dir/"construction.csv");
        construction<<"case,method,label_seconds,candidate_seconds,oracle_calls,dp_seconds,states,transitions\n";
        auto timings=output(dir/"timings.csv");
        timings<<"case,budget,method,repeat,packed_ns,topology_f32_ns,packed_checksum,topology_f32_checksum,backtrack_pack_seconds\n";
        auto surfaces=output(dir/"surfaces.csv");
        surfaces<<"case,shape,feature,seed,budget,method,cells,vertices,triangles,components,euler,boundary_edges,nonmanifold_edges,nonmanifold_vertices,topology_pass,probe_pass,success,seconds\n";
        const std::vector<std::string> methods{"f32","q16","q8","joint"};
        for(const auto& shape:shapes) for(double feature:features) for(unsigned seed=0;seed<nseeds;++seed) {
            const auto prefix=shape+"_"+std::to_string(feature)+"_"+std::to_string(seed);
            Shifted oracle(shape,feature); auto start=Clock::now();
            const auto fit=sample_queries(oracle,nfit,seed,.5,feature/2); const double label_time=seconds(start);
            start=Clock::now(); PrecisionTree tree([&](Vec3 p){return oracle.distance(p);},fit,depth); const double candidate_time=seconds(start);
            save(dir/"samples"/(prefix+"_fit.csv"),fit);
            auto uniform=sample_queries(oracle,ntest,seed+10000,0,feature/2),near=sample_queries(oracle,ntest,seed+20000,1,feature/2);
            auto base=analytic_oracle(shape,feature); auto probes=feature_probes(*base,shape,feature);
            for(auto& p:probes) { p.point=p.point+shift; p.distance=oracle.distance(p.point); }
            save(dir/"samples"/(prefix+"_uniform.csv"),uniform); save(dir/"samples"/(prefix+"_near.csv"),near); save(dir/"samples"/(prefix+"_probes.csv"),probes);
            {
                auto f=output(dir/"candidates"/(prefix+".csv")); f<<"node,first_child,f32_loss,q16_loss,q8_loss\n";
                for(std::size_t v=0;v<tree.nodes.size();++v) f<<v<<','<<tree.spatial.nodes[v].field.first_child<<','<<tree.nodes[v].loss[0]<<','<<tree.nodes[v].loss[1]<<','<<tree.nodes[v].loss[2]<<'\n';
            }
            auto extract=[&](const Distance& query,unsigned budget,const std::string& method,bool probe_pass) {
                for(auto cells:resolutions) {
                    const auto ts=Clock::now(); auto mesh=extract_surface(query,cells); const auto topology=surface_topology(mesh); const double duration=seconds(ts);
                    const bool pass=topology_pass(shape,topology);
                    surfaces<<prefix<<','<<shape<<','<<feature<<','<<seed<<','<<budget<<','<<method<<','<<cells<<','<<mesh.vertices.size()<<','<<mesh.triangles.size()<<','<<topology.components<<','<<topology.euler<<','<<topology.boundary_edges<<','<<topology.nonmanifold_edges<<','<<topology.nonmanifold_vertices<<','<<pass<<','<<probe_pass<<','<<(pass && probe_pass)<<','<<duration<<'\n'; surfaces.flush();
                    if(cells==resolutions.back()) mesh.write_obj((dir/"meshes"/(prefix+"_"+std::to_string(budget)+"_"+method+".obj")).string());
                }
            };
            extract([&](Vec3 p){return oracle.distance(p);},0,"oracle",true);
            for(unsigned m=0;m<methods.size();++m) {
                const auto mode=m==3?Encoding::branch:static_cast<Encoding>(m+1);
                start=Clock::now(); PrecisionDP dp(tree,budgets.back(),mode); const double dp_time=seconds(start);
                construction<<prefix<<','<<methods[m]<<','<<label_time<<','<<candidate_time<<','<<tree.spatial.oracle_calls<<','<<dp_time<<','<<dp.states<<','<<dp.transitions<<'\n'; construction.flush();
                for(auto budget:budgets) {
                    start=Clock::now(); auto allocation=dp.allocate(budget); PackedField field(tree,allocation); const double pack_time=seconds(start);
                    if(field.bytes()>budget) throw std::runtime_error("Byte cap exceeded");
                    double fit_loss=0; for(const auto& s:fit) { const double e=(field.query(s.point)-s.distance)/2; fit_loss+=e*e; }
                    if(std::abs(fit_loss-allocation.loss)>1e-9*(1+allocation.loss)) throw std::runtime_error("Packed fitting SSE mismatch");
                    const auto name=prefix+"_"+std::to_string(budget)+"_"+methods[m];
                    {
                        std::ofstream f(dir/"fields"/(name+".bin"),std::ios::binary); f.exceptions(std::ios::badbit|std::ios::failbit);
                        f.write(reinterpret_cast<const char*>(field.data()),static_cast<std::streamsize>(field.payload_bytes()));
                        auto choices=output(dir/"choices"/(name+".csv")); choices<<"node,encoding\n";
                        std::function<void(std::size_t)> save_choice=[&](std::size_t v) {
                            choices<<v<<','<<unsigned(allocation.choice[v])<<'\n';
                            if(allocation.choice[v]==Encoding::branch) for(unsigned c=0;c<8;++c) save_choice(tree.spatial.nodes[v].field.first_child+c);
                        }; save_choice(0);
                    }
                    PrecisionAllocation copy=allocation;
                    for(auto& choice:copy.choice) if(choice!=Encoding::branch) choice=Encoding::f32;
                    copy.bytes=sizeof(PackedField)+4*copy.branches+36*(copy.leaves[0]+copy.leaves[1]+copy.leaves[2]);
                    PackedField float_copy(tree,copy);
                    volatile double sink=0;
                    auto time=[&](const PackedField& f) {
                        auto ts=Clock::now(); double sum=0; for(const auto& s:uniform) sum+=f.query(s.point); sink=sum;
                        return std::pair<double,double>{seconds(ts)*1e9/uniform.size(),sum};
                    };
                    time(field); time(float_copy);
                    for(unsigned repeat=0;repeat<3;++repeat) {
                        std::pair<double,double> compressed,reference;
                        if(repeat%2) { reference=time(float_copy); compressed=time(field); }
                        else { compressed=time(field); reference=time(float_copy); }
                        timings<<prefix<<','<<budget<<','<<methods[m]<<','<<repeat<<','<<compressed.first<<','<<reference.first<<','<<compressed.second<<','<<reference.second<<','<<pack_time<<'\n';
                    }
                    (void)sink; timings.flush();
                    for(const auto& set:std::vector<std::pair<std::string,const Samples*>>{{"uniform",&uniform},{"near",&near},{"probes",&probes}}) {
                        const auto metric=evaluate([&](Vec3 p){return field.query(p);},*set.second);
                        metrics<<prefix<<','<<shape<<','<<feature<<','<<seed<<','<<budget<<','<<methods[m]<<','<<field.bytes()<<','<<allocation.branches<<','<<allocation.leaves[0]<<','<<allocation.leaves[1]<<','<<allocation.leaves[2]<<','<<allocation.loss<<','<<set.first<<','<<metric.rmse<<','<<metric.mae<<','<<metric.p99<<','<<metric.cvar99<<','<<metric.inside_recall<<','<<metric.outside_recall<<','<<metric.sign_error<<'\n';
                    }
                    metrics.flush(); const auto retention=evaluate([&](Vec3 p){return field.query(p);},probes);
                    extract([&](Vec3 p){return field.query(p);},budget,methods[m],retention.inside_recall>=.95 && retention.outside_recall>=.95);
                }
            }
            std::cout<<prefix<<" complete"<<std::endl;
        }
        output(dir/"COMPLETE")<<"precision-v1 completed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
