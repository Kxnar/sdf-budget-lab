#include "sdf/precision.hpp"
#include "sdf/tail.hpp"
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
PrecisionAllocation full(const PrecisionTree& tree) {
    PrecisionAllocation a; a.choice.resize(tree.nodes.size(),Encoding::q8);
    for(std::size_t v=0;v<a.choice.size();++v) {
        if(tree.spatial.nodes[v].field.first_child) { a.choice[v]=Encoding::branch; ++a.branches; }
        else { ++a.leaves[2]; a.loss+=tree.nodes[v].loss[2]; }
    }
    a.bytes=sizeof(PackedField)+4*a.branches+20*a.leaves[2]; return a;
}
void development(const fs::path& dir,bool smoke) {
    auto results=output(dir/"development.csv");
    results<<"case,budget,grid,objective,fit_cvar99,baseline_cvar99,threshold,bytes,seconds,thresholds\n";
    auto trace=output(dir/"thresholds.csv");
    trace<<"case,budget,grid,threshold,excess,objective,bytes,states,transitions\n";
    const std::vector<std::string> shapes=smoke?std::vector<std::string>{"torus"}:std::vector<std::string>{"thin_box","torus","double_sphere"};
    for(const auto& shape:shapes) for(double feature:{.03,.08}) {
        const auto name=shape+"_"+std::to_string(feature)+"_101";
        Shifted oracle(shape,feature); const auto fit=sample_queries(oracle,smoke?200:10000,101,.5,feature/2);
        save(dir/"samples"/(name+"_fit.csv"),fit);
        PrecisionTree tree([&](Vec3 p){return oracle.distance(p);},fit,smoke?2:5); Q8Losses losses(tree,fit);
        for(auto budget:{8192u,16384u,32768u}) for(auto grid:{9u,17u,33u}) {
            const auto start=Clock::now(); const auto a=allocate_cvar99(tree,losses,budget,grid); const auto duration=seconds(start);
            results<<name<<','<<budget<<','<<grid<<','<<a.objective<<','<<a.fit_cvar99<<','<<a.upper<<','<<a.threshold<<','<<a.allocation.bytes<<','<<duration<<','<<a.search.size()<<'\n'; results.flush();
            for(const auto& t:a.search) trace<<name<<','<<budget<<','<<grid<<','<<t.threshold<<','<<t.excess<<','<<t.objective<<','<<t.bytes<<','<<t.states<<','<<t.transitions<<'\n'; trace.flush();
        }
        std::cout<<name<<" development complete"<<std::endl;
    }
}
}
int main(int argc,char** argv) {
    try {
        std::string out,stage="resolution"; bool smoke=false; unsigned grid=0;
        for(int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--smoke") smoke=true;
            else if(arg=="--out" && i+1<argc) out=argv[++i];
            else if(arg=="--stage" && i+1<argc) stage=argv[++i];
            else if(arg=="--grid" && i+1<argc) grid=static_cast<unsigned>(std::stoul(argv[++i]));
            else throw std::invalid_argument("Usage: sdf_tail_bench --out NEW_DIR --stage resolution|development|evaluation [--grid N] [--smoke]");
        }
        if(stage!="resolution" && stage!="development" && stage!="evaluation") throw std::invalid_argument("Unknown stage");
        if(stage=="evaluation" && (grid<2 || grid>4097)) throw std::invalid_argument("Evaluation requires frozen --grid");
        if(out.empty() || fs::exists(out)) throw std::invalid_argument("Output must be a new directory");
        const fs::path dir(out);
        for(const auto* p:{"samples","fields","meshes","choices","candidates"}) fs::create_directories(dir/p);
        if(stage=="development") {
            output(dir/"config.json")<<"{\"stage\":\"development\",\"smoke\":"<<(smoke?"true":"false")<<"}\n";
            development(dir,smoke); output(dir/"COMPLETE")<<"development completed\n"; return 0;
        }
        const unsigned nfit=smoke?200:10000,ntest=smoke?1000:100000;
        const std::vector<std::string> shapes=smoke?std::vector<std::string>{"torus"}:std::vector<std::string>{"thin_box","torus","double_sphere"};
        const std::vector<double> features=smoke?std::vector<double>{.06}:std::vector<double>{.02,.035,.06,.10};
        const std::vector<unsigned> depths=stage=="resolution"?(smoke?std::vector<unsigned>{2,3}:std::vector<unsigned>{4,5}):std::vector<unsigned>{smoke?2u:5u};
        const std::vector<unsigned> seeds=stage=="resolution"?std::vector<unsigned>{20}:smoke?std::vector<unsigned>{30}:std::vector<unsigned>{30,31};
        const std::vector<unsigned> resolutions=smoke?std::vector<unsigned>{16,24}:std::vector<unsigned>{32,64,96};
        output(dir/"config.json")<<"{\"protocol\":\"tail-v1\",\"stage\":"<<std::quoted(stage)<<",\"grid\":"<<grid<<",\"smoke\":"<<(smoke?"true":"false")<<",\"fit\":"<<nfit<<",\"test\":"<<ntest<<",\"field_object_bytes\":"<<sizeof(PackedField)<<",\"configuration\":"<<std::quoted(SDF_BUILD_CONFIGURATION)<<"}\n";
        auto metrics=output(dir/"metrics.csv");
        metrics<<"case,shape,feature,seed,depth,budget,method,bytes,branches,f32_leaves,q16_leaves,q8_leaves,fit_loss,threshold,fit_objective,fit_cvar99,dataset,count,rmse,mae,p99,cvar99,inside_recall,outside_recall,sign_error\n";
        auto construction=output(dir/"construction.csv");
        construction<<"case,depth,budget,method,label_seconds,candidate_seconds,loss_cache_seconds,allocation_seconds,pack_seconds,oracle_calls\n";
        auto trace=output(dir/"thresholds.csv"); trace<<"case,budget,grid,threshold,excess,objective,bytes,states,transitions\n";
        auto surfaces=output(dir/"surfaces.csv");
        surfaces<<"case,shape,feature,seed,depth,budget,method,cells,vertices,triangles,components,euler,boundary_edges,nonmanifold_edges,nonmanifold_vertices,topology_pass,probe_pass,success,seconds\n";
        for(const auto& shape:shapes) for(double feature:features) for(auto seed:seeds) {
            const auto prefix=shape+"_"+std::to_string(feature)+"_"+std::to_string(seed);
            Shifted oracle(shape,feature);
            const auto label_start=Clock::now(); const auto fit=sample_queries(oracle,nfit,seed,.5,feature/2); const double label_time=seconds(label_start);
            auto uniform=sample_queries(oracle,ntest,seed+10000,0,feature/2),near=sample_queries(oracle,ntest,seed+20000,1,feature/2);
            auto base=analytic_oracle(shape,feature); auto probes=feature_probes(*base,shape,feature);
            for(auto& p:probes) { p.point=p.point+shift; p.distance=oracle.distance(p.point); }
            save(dir/"samples"/(prefix+"_fit.csv"),fit);
            save(dir/"samples"/(prefix+"_uniform.csv"),uniform); save(dir/"samples"/(prefix+"_near.csv"),near); save(dir/"samples"/(prefix+"_probes.csv"),probes);
            auto extract=[&](const Distance& query,unsigned depth,std::size_t budget,const std::string& method,bool probe_pass) {
                for(auto cells:resolutions) {
                    const auto ts=Clock::now(); auto mesh=extract_surface(query,cells); const auto topology=surface_topology(mesh); const double duration=seconds(ts);
                    const bool pass=topology_pass(shape,topology);
                    surfaces<<prefix<<','<<shape<<','<<feature<<','<<seed<<','<<depth<<','<<budget<<','<<method<<','<<cells<<','<<mesh.vertices.size()<<','<<mesh.triangles.size()<<','<<topology.components<<','<<topology.euler<<','<<topology.boundary_edges<<','<<topology.nonmanifold_edges<<','<<topology.nonmanifold_vertices<<','<<pass<<','<<probe_pass<<','<<(pass && probe_pass)<<','<<duration<<'\n'; surfaces.flush();
                    if(cells==resolutions.back()) mesh.write_obj((dir/"meshes"/(prefix+"_"+std::to_string(depth)+"_"+std::to_string(budget)+"_"+method+".obj")).string());
                }
            };
            extract([&](Vec3 p){return oracle.distance(p);},0,0,"oracle",true);
            for(auto depth:depths) {
                auto start=Clock::now(); PrecisionTree tree([&](Vec3 p){return oracle.distance(p);},fit,depth); const double candidate_time=seconds(start);
                start=Clock::now(); Q8Losses losses(tree,fit); const double cache_time=seconds(start);
                { auto f=output(dir/"candidates"/(prefix+"_"+std::to_string(depth)+".csv")); f<<"node,first_child\n";
                  for(std::size_t v=0;v<tree.nodes.size();++v) f<<v<<','<<tree.spatial.nodes[v].field.first_child<<'\n'; }
                const std::vector<std::string> methods=stage=="resolution"?std::vector<std::string>{"full_q8"}:std::vector<std::string>{"mse","sign","cvar99"};
                const std::vector<std::size_t> budgets=stage=="resolution"?std::vector<std::size_t>{full(tree).bytes}:std::vector<std::size_t>{8192,16384,32768};
                for(auto budget:budgets) for(const auto& method:methods) {
                start=Clock::now(); PrecisionAllocation allocation; double threshold=0,fit_objective=0;
                if(method=="full_q8") allocation=full(tree);
                else if(method=="cvar99") {
                    auto a=allocate_cvar99(tree,losses,budget,grid); allocation=std::move(a.allocation); threshold=a.threshold; fit_objective=a.objective;
                    for(const auto& t:a.search) trace<<prefix<<','<<budget<<','<<grid<<','<<t.threshold<<','<<t.excess<<','<<t.objective<<','<<t.bytes<<','<<t.states<<','<<t.transitions<<'\n'; trace.flush();
                } else {
                    losses.apply(tree,method=="mse"?Q8Objective::mse:Q8Objective::sign,.0025);
                    allocation=PrecisionDP(tree,budget,Encoding::q8).allocate(budget);
                }
                const double allocation_time=seconds(start); const double fit_cvar=losses.cvar99(tree,allocation);
                start=Clock::now(); PackedField field(tree,allocation); const double pack_time=seconds(start);
                if(field.bytes()>budget) throw std::runtime_error("Byte cap exceeded");
                construction<<prefix<<','<<depth<<','<<budget<<','<<method<<','<<label_time<<','<<candidate_time<<','<<cache_time<<','<<allocation_time<<','<<pack_time<<','<<tree.spatial.oracle_calls<<'\n'; construction.flush();
                const auto name=prefix+"_"+std::to_string(depth)+"_"+std::to_string(budget)+"_"+method;
                {
                    std::ofstream f(dir/"fields"/(name+".bin"),std::ios::binary); f.exceptions(std::ios::badbit|std::ios::failbit);
                    f.write(reinterpret_cast<const char*>(field.data()),static_cast<std::streamsize>(field.payload_bytes()));
                    auto choices=output(dir/"choices"/(name+".csv")); choices<<"node,encoding,loss\n";
                    std::function<void(std::size_t)> visit=[&](std::size_t v) {
                        choices<<v<<','<<unsigned(allocation.choice[v])<<','<<tree.nodes[v].loss[2]<<'\n';
                        if(allocation.choice[v]==Encoding::branch) for(unsigned c=0;c<8;++c) visit(tree.spatial.nodes[v].field.first_child+c);
                    }; visit(0);
                }
                for(const auto& set:std::vector<std::pair<std::string,const Samples*>>{{"uniform",&uniform},{"near",&near},{"probes",&probes}}) {
                    const auto m=evaluate([&](Vec3 p){return field.query(p);},*set.second);
                    metrics<<prefix<<','<<shape<<','<<feature<<','<<seed<<','<<depth<<','<<budget<<','<<method<<','<<field.bytes()<<','<<allocation.branches<<",0,0,"<<allocation.leaves[2]<<','<<allocation.loss<<','<<threshold<<','<<fit_objective<<','<<fit_cvar<<','<<set.first<<','<<m.count<<','<<m.rmse<<','<<m.mae<<','<<m.p99<<','<<m.cvar99<<','<<m.inside_recall<<','<<m.outside_recall<<','<<m.sign_error<<'\n';
                }
                metrics.flush(); const auto retention=evaluate([&](Vec3 p){return field.query(p);},probes);
                extract([&](Vec3 p){return field.query(p);},depth,budget,method,retention.inside_recall>=.95 && retention.outside_recall>=.95);
                }
            }
            std::cout<<prefix<<" complete"<<std::endl;
        }
        output(dir/"COMPLETE")<<stage<<" completed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
