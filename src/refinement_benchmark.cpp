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
    return t.closed_manifold() && t.component_euler==(shape=="torus"?std::vector<long long>{0}:std::vector<long long>{2,2});
}
}
int main(int argc,char** argv) {
    try {
        std::string out; bool smoke=false,development=false; std::size_t selected_k=0;
        for(int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--smoke") smoke=true;
            else if(arg=="--development") development=true;
            else if(arg=="--k" && i+1<argc) selected_k=std::stoul(argv[++i]);
            else if(arg=="--out" && i+1<argc) out=argv[++i];
            else throw std::invalid_argument("Usage: sdf_refinement_bench --out NEW_DIR [--development | --k K] [--smoke]");
        }
        if(out.empty() || fs::exists(out)) throw std::invalid_argument("Output must be a new directory");
        const fs::path dir(out);
        if(development) {
            fs::create_directories(dir);
            auto f=output(dir/"development.csv");
            f<<"shape,feature,seed,k,selected_parents,fit_sse,bytes,oracle_calls,total_seconds,retained_parents\n";
            for(const std::string shape:{"torus","double_sphere"}) for(double feature:{.03,.08}) for(unsigned seed:{201u,202u}) {
                Shifted oracle(shape,feature); Distance truth=[&](Vec3 p){return oracle.distance(p);};
                const auto ls=Clock::now(); const auto fit=sample_queries(oracle,10000,seed,.5,feature/2); const double label_time=seconds(ls);
                for(std::size_t k:{0u,32u,64u,128u,256u,512u,1024u}) {
                    const auto start=Clock::now(); PrecisionTree tree(truth,fit,5);
                    const auto selected=tree.refine_depth_six(truth,fit,k);
                    const auto a=PrecisionDP(tree,32768,Encoding::q8).allocate(32768); PackedField field(tree,a);
                    const double total=label_time+seconds(start);
                    std::size_t retained=0; for(const auto& r:selected) retained+=a.choice[r.node]==Encoding::branch;
                    f<<shape<<','<<feature<<','<<seed<<','<<k<<','<<selected.size()<<','<<a.loss<<','<<field.bytes()<<','<<(tree.spatial.oracle_calls+fit.size())<<','<<total<<','<<retained<<'\n'; f.flush();
                }
                std::cout<<shape<<' '<<feature<<' '<<seed<<" development complete"<<std::endl;
            }
            output(dir/"COMPLETE")<<"development completed\n"; return 0;
        }
        if(!selected_k || selected_k>32768) throw std::invalid_argument("Final run needs frozen positive K");
        for(const auto* p:{"samples","fields","meshes","predictions"}) fs::create_directories(dir/p);
        const unsigned nfit=smoke?200:10000,ntest=smoke?1000:100000;
        const std::vector<unsigned> seeds=smoke?std::vector<unsigned>{40}:std::vector<unsigned>{40,41};
        const std::vector<unsigned> resolutions=smoke?std::vector<unsigned>{16,24}:std::vector<unsigned>{32,64,96};
        output(dir/"config.json")<<"{\"protocol\":\"refinement-v1\",\"smoke\":"<<(smoke?"true":"false")<<",\"fit\":"<<nfit<<",\"test\":"<<ntest<<",\"field_object_bytes\":"<<sizeof(PackedField)<<",\"configuration\":"<<std::quoted(SDF_BUILD_CONFIGURATION)<<"}\n";
        output(dir/"selected_k.txt")<<selected_k;
        auto selection=output(dir/"selection.csv");
        selection<<"case,rank,node,fit_count,q8_loss,lo_x,lo_y,lo_z,hi_x,hi_y,hi_z,retained\n";
        auto allocation=output(dir/"allocation.csv");
        allocation<<"case,method,fit_sse,candidate_nodes,selected_parents,retained_parents,dp_states,dp_transitions\n";
        auto metrics=output(dir/"metrics.csv");
        metrics<<"case,shape,feature,seed,depth,method,bytes,dataset,count,mae,cvar99,inside_recall,outside_recall,sign_error\n";
        auto construction=output(dir/"construction.csv");
        construction<<"case,depth,method,label_seconds,candidate_seconds,refinement_seconds,allocation_seconds,pack_seconds,total_seconds,oracle_calls\n";
        auto surfaces=output(dir/"surfaces.csv");
        surfaces<<"case,shape,feature,seed,depth,method,cells,components,euler,boundary_edges,nonmanifold_edges,nonmanifold_vertices,topology_pass,probe_pass,success\n";
        for(const std::string shape:{"torus","double_sphere"}) {
        const auto features=smoke?std::vector<double>{.035}:shape=="torus"?std::vector<double>{.035,.06,.10}:std::vector<double>{.06,.10};
        for(double feature:features) for(auto seed:seeds) {
            const auto prefix=shape+"_"+std::to_string(feature)+"_"+std::to_string(seed);
            Shifted oracle(shape,feature); Distance truth=[&](Vec3 p){return oracle.distance(p);};
            const auto ls=Clock::now(); const auto fit=sample_queries(oracle,nfit,seed,.5,feature/2); const double label_time=seconds(ls);
            const auto uniform=sample_queries(oracle,ntest,seed+10000,0,feature/2),near=sample_queries(oracle,ntest,seed+20000,1,feature/2);
            auto base=analytic_oracle(shape,feature); auto probes=feature_probes(*base,shape,feature);
            for(auto& p:probes) { p.point=p.point+shift; p.distance=truth(p.point); }
            for(const auto& set:std::vector<std::pair<std::string,const Samples*>>{{"fit",&fit},{"uniform",&uniform},{"near",&near},{"probes",&probes}}) save(dir/"samples"/(prefix+"_"+set.first+".csv"),*set.second);
            auto extract=[&](const Distance& query,unsigned depth,const std::string& method,bool probe_pass) {
                for(auto cells:resolutions) {
                    auto mesh=extract_surface(query,cells); const auto topology=surface_topology(mesh); const bool pass=topology_pass(shape,topology);
                    surfaces<<prefix<<','<<shape<<','<<feature<<','<<seed<<','<<depth<<','<<method<<','<<cells<<','<<topology.components<<','<<topology.euler<<','<<topology.boundary_edges<<','<<topology.nonmanifold_edges<<','<<topology.nonmanifold_vertices<<','<<pass<<','<<probe_pass<<','<<(pass && probe_pass)<<'\n'; surfaces.flush();
                    if(cells==resolutions.back()) mesh.write_obj((dir/"meshes"/(prefix+"_"+std::to_string(depth)+"_"+method+".obj")).string());
                }
            };
            extract(truth,0,"oracle",true);
            for(const std::string method:{"mse_q8","selective_q8"}) {
                const unsigned depth=method=="mse_q8"?5:6;
                auto start=Clock::now(); PrecisionTree tree(truth,fit,5); const double candidate_time=seconds(start);
                std::vector<RefinementCandidate> selected; start=Clock::now();
                if(method=="selective_q8") selected=tree.refine_depth_six(truth,fit,selected_k);
                const double refinement_time=seconds(start);
                start=Clock::now(); PrecisionDP dp(tree,32768,Encoding::q8); const auto a=dp.allocate(32768); const double allocation_time=seconds(start);
                start=Clock::now(); auto field=std::make_unique<PackedField>(tree,a); const double pack_time=seconds(start);
                if(field->bytes()>32768) throw std::logic_error("Exceeded byte cap");
                std::size_t retained=0,rank=0;
                for(const auto& r:selected) {
                    const bool keep=a.choice[r.node]==Encoding::branch; retained+=keep;
                    selection<<prefix<<','<<++rank<<','<<r.node<<','<<r.fitting_count<<','<<r.q8_loss<<','<<r.box.lo.x<<','<<r.box.lo.y<<','<<r.box.lo.z<<','<<r.box.hi.x<<','<<r.box.hi.y<<','<<r.box.hi.z<<','<<keep<<'\n';
                }
                selection.flush();
                allocation<<prefix<<','<<method<<','<<a.loss<<','<<tree.nodes.size()<<','<<selected.size()<<','<<retained<<','<<dp.states<<','<<dp.transitions<<'\n'; allocation.flush();
                construction<<prefix<<','<<depth<<','<<method<<','<<label_time<<','<<candidate_time<<','<<refinement_time<<','<<allocation_time<<','<<pack_time<<','<<(label_time+candidate_time+refinement_time+allocation_time+pack_time)<<','<<(tree.spatial.oracle_calls+nfit)<<'\n'; construction.flush();
                const auto name=prefix+"_"+std::to_string(depth)+"_"+method;
                { std::ofstream f(dir/"fields"/(name+".bin"),std::ios::binary); f.exceptions(std::ios::badbit|std::ios::failbit);
                  f.write(reinterpret_cast<const char*>(field->data()),static_cast<std::streamsize>(field->payload_bytes())); }
                for(const auto& set:std::vector<std::pair<std::string,const Samples*>>{{"fit",&fit},{"uniform",&uniform},{"near",&near},{"probes",&probes}}) {
                    const auto m=evaluate([&](Vec3 p){return field->query(p);},*set.second);
                    metrics<<prefix<<','<<shape<<','<<feature<<','<<seed<<','<<depth<<','<<method<<','<<field->bytes()<<','<<set.first<<','<<m.count<<','<<m.mae<<','<<m.cvar99<<','<<m.inside_recall<<','<<m.outside_recall<<','<<m.sign_error<<'\n';
                    std::ofstream f(dir/"predictions"/(name+"_"+set.first+".bin"),std::ios::binary); f.exceptions(std::ios::badbit|std::ios::failbit);
                    for(const auto& p:*set.second) { const double v=field->query(p.point); f.write(reinterpret_cast<const char*>(&v),sizeof(v)); }
                }
                metrics.flush(); const auto retention=evaluate([&](Vec3 p){return field->query(p);},probes);
                extract([&](Vec3 p){return field->query(p);},depth,method,retention.inside_recall>=.95 && retention.outside_recall>=.95);
            }
            std::cout<<prefix<<" complete"<<std::endl;
        }}
        output(dir/"COMPLETE")<<"refinement completed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
