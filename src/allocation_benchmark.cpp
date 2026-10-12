#include "sdf/allocation.hpp"
#include "sdf/metrics.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <set>

namespace {
using namespace sdf;
namespace fs=std::filesystem;
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }
std::ofstream output(const fs::path& p) { std::ofstream f(p); f.exceptions(std::ios::badbit|std::ios::failbit); f<<std::setprecision(17); return f; }
Vec3 rotate(Vec3 p,Vec3 angles,bool inverse=false) {
    const double rad=3.14159265358979323846/180;
    auto x=[](Vec3 v,double a) { return Vec3{v.x,std::cos(a)*v.y-std::sin(a)*v.z,std::sin(a)*v.y+std::cos(a)*v.z}; };
    auto y=[](Vec3 v,double a) { return Vec3{std::cos(a)*v.x+std::sin(a)*v.z,v.y,-std::sin(a)*v.x+std::cos(a)*v.z}; };
    auto z=[](Vec3 v,double a) { return Vec3{std::cos(a)*v.x-std::sin(a)*v.y,std::sin(a)*v.x+std::cos(a)*v.y,v.z}; };
    if(inverse) return x(y(z(p,-angles.z*rad),-angles.y*rad),-angles.x*rad);
    return z(y(x(p,angles.x*rad),angles.y*rad),angles.z*rad);
}
struct Case { std::string id,split,shape; double feature; Vec3 translation,angles; };
class Transformed final:public Oracle {
    std::unique_ptr<Oracle> base_; Case c_;
public:
    explicit Transformed(const Case& c):base_(analytic_oracle(c.shape,c.feature)),c_(c) {}
    Vec3 world(Vec3 p) const { return rotate(p,c_.angles)+c_.translation; }
    double distance(Vec3 p) const override { return base_->distance(rotate(p-c_.translation,c_.angles,true)); }
    Vec3 surface(std::mt19937_64& rng) const override { return world(base_->surface(rng)); }
    std::string name() const override { return c_.id; }
};
std::vector<Case> cases(const std::string& path) {
    std::ifstream f(path); if(!f) throw std::runtime_error("Cannot read case manifest");
    std::string line; std::getline(f,line); std::vector<Case> result; std::set<std::string> ids;
    while(std::getline(f,line)) {
        std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line); Case c; std::string extra;
        if(!(row>>c.id>>c.split>>c.shape>>c.feature>>c.translation.x>>c.translation.y>>c.translation.z>>c.angles.x>>c.angles.y>>c.angles.z) || row>>extra)
            throw std::invalid_argument("Invalid case row");
        if(!ids.insert(c.id).second || (c.split!="development" && c.split!="heldout") || !std::isfinite(norm(c.translation)) || !std::isfinite(norm(c.angles)))
            throw std::invalid_argument("Invalid case values");
        result.push_back(c);
    }
    if(result.empty()) throw std::invalid_argument("Empty case manifest");
    return result;
}
void save_samples(const fs::path& path,const Samples& samples) {
    auto f=output(path); f<<"x,y,z,distance\n";
    for(const auto& s:samples) f<<s.point.x<<','<<s.point.y<<','<<s.point.z<<','<<s.distance<<'\n';
}
}
int main(int argc,char** argv) {
    try {
        std::string case_path="fixtures/allocation_cases.csv",out; unsigned depth=4;
        for(int i=1;i<argc;++i) {
            const std::string key=argv[i];
            if(i+1==argc) throw std::invalid_argument("Missing argument");
            const std::string value=argv[++i];
            if(key=="--cases") case_path=value;
            else if(key=="--out") out=value;
            else if(key=="--depth") {
                if(value.size()!=1 || value[0]<'0' || value[0]>'5') throw std::invalid_argument("Depth must be 0..5");
                depth=unsigned(value[0]-'0');
            } else throw std::invalid_argument("Unknown option: "+key);
        }
        if(out.empty()) throw std::invalid_argument("Usage: sdf_allocation_bench --out NEW_DIR [--cases CSV] [--depth 4]");
        const auto manifest=cases(case_path); const fs::path dir(out);
        if(fs::exists(dir)) throw std::invalid_argument("Output must be a new directory");
        fs::create_directories(dir/"samples"); fs::create_directories(dir/"candidates");
        fs::copy_file(case_path,dir/"cases.csv");
        {
            auto f=output(dir/"config.json");
            f<<"{\"schema\":1,\"depth\":"<<depth<<",\"train\":4096,\"test_per_distribution\":4096,\"seeds\":[0,1,2],\"budgets\":[4096,8192,16384,32768],\"horizons\":[1,2,4,8],\"timing_repeats\":3,\"node_bytes\":"<<sizeof(AdaptiveField::Node)<<",\"field_object_bytes\":"<<sizeof(AllocatedField)<<",\"loss\":\"sum squared distance error / domain_scale^2\",\"rotation\":\"Rz Ry Rx; degrees\",\"seed_streams\":{\"fit\":0,\"uniform\":10000,\"near\":20000,\"surface\":30000},\"near_sigma\":\"feature/2\",\"candidate_corner_cache\":false}\n";
        }
        {
            auto f=output(dir/"build.json");
            f<<"{\"configuration\":"<<std::quoted(SDF_BUILD_CONFIGURATION)<<",\"compiler\":";
#ifdef _MSC_VER
            f<<std::quoted("MSVC "+std::to_string(_MSC_FULL_VER));
#else
            f<<std::quoted(__VERSION__);
#endif
            f<<"}\n";
        }
        auto costs=output(dir/"candidate_costs.csv");
        costs<<"case,split,seed,nodes,fit_label_seconds,fit_label_calls,candidate_seconds,candidate_oracle_calls,evaluation_label_calls\n";
        auto results=output(dir/"metrics.csv");
        results<<"case,split,seed,budget,method,bytes,splits,fit_sse,allocation_seconds,packing_seconds,total_construction_seconds,allocation_oracle_calls,dp_transitions,dp_states,negative_steps,dataset,rmse,mae,p99,sign_error,inside_recall,outside_recall\n";
        auto times=output(dir/"timings.csv");
        times<<"case,seed,budget,method,repeat,allocation_seconds,packing_seconds\n";
        auto selections=output(dir/"selections.csv");
        selections<<"case,seed,budget,method,split_node_ids\n";
        const std::vector<std::string> methods{"exact","greedy_spend","greedy_best","lookahead_2","lookahead_4","lookahead_8"};
        for(const auto& c:manifest) for(unsigned seed=0;seed<3;++seed) {
            Transformed oracle(c); const auto start=Clock::now();
            auto fit=sample_queries(oracle,4096,seed,0.5,c.feature/2);
            const double labels=elapsed(start);
            const auto build_start=Clock::now();
            CandidateTree tree([&](Vec3 p){return oracle.distance(p);},fit,depth);
            const double candidate_time=elapsed(build_start);
            const auto prefix=c.id+"_"+std::to_string(seed);
            {
                auto f=output(dir/"candidates"/(prefix+".csv")); f<<"node,first_child,loss\n";
                for(std::size_t v=0;v<tree.nodes.size();++v) f<<v<<','<<tree.nodes[v].field.first_child<<','<<tree.nodes[v].loss<<'\n';
            }
            save_samples(dir/"samples"/(prefix+"_fit.csv"),fit);
            std::vector<std::pair<std::string,Samples>> eval;
            eval.push_back({"uniform",sample_queries(oracle,4096,10000+seed,0,c.feature/2)});
            eval.push_back({"near",sample_queries(oracle,4096,20000+seed,1,c.feature/2)});
            Samples surface; std::mt19937_64 rng(30000+seed);
            for(unsigned i=0;i<4096;++i) { const auto p=oracle.surface(rng); if(!Box{}.contains(p)) throw std::runtime_error("Surface outside domain"); surface.push_back({p,oracle.distance(p)}); }
            eval.push_back({"surface",std::move(surface)});
            auto base=analytic_oracle(c.shape,c.feature); auto probes=feature_probes(*base,c.shape,c.feature);
            for(auto& p:probes) { p.point=oracle.world(p.point); p.distance=oracle.distance(p.point); }
            eval.push_back({"probes",std::move(probes)});
            std::size_t evaluation_calls=0;
            for(const auto& set:eval) { evaluation_calls+=set.second.size(); save_samples(dir/"samples"/(prefix+"_"+set.first+".csv"),set.second); }
            // feature_probes also labels its base-frame points before transformation.
            evaluation_calls+=eval.back().second.size();
            costs<<c.id<<','<<c.split<<','<<seed<<','<<tree.nodes.size()<<','<<labels<<','<<fit.size()<<','<<candidate_time<<','<<tree.oracle_calls<<','<<evaluation_calls<<'\n'; costs.flush();
            for(const auto budget:{4096u,8192u,16384u,32768u}) for(std::size_t m=0;m<methods.size();++m) {
                Allocation a; std::vector<double> alloc_times,pack_times;
                std::unique_ptr<AllocatedField> field;
                for(unsigned repeat=0;repeat<3;++repeat) {
                    const auto ts=Clock::now();
                    if(m==0) a=exact_allocate(tree,budget);
                    else a=lookahead_allocate(tree,budget,m<=2?1:std::size_t(1)<<(m-2),m!=1);
                    alloc_times.push_back(elapsed(ts)); const auto ps=Clock::now();
                    field=std::make_unique<AllocatedField>(tree,a); pack_times.push_back(elapsed(ps));
                    times<<c.id<<','<<seed<<','<<budget<<','<<methods[m]<<','<<repeat<<','<<alloc_times.back()<<','<<pack_times.back()<<'\n';
                }
                if(field->bytes()>budget) throw std::runtime_error("Storage cap exceeded");
                selections<<c.id<<','<<seed<<','<<budget<<','<<methods[m]<<',';
                for(std::size_t v=0;v<a.split.size();++v) if(a.split[v]) selections<<v<<';';
                selections<<'\n'; selections.flush();
                double measured_loss=0;
                for(const auto& s:fit) { const double e=(field->query(s.point)-s.distance)/tree.domain.scale(); measured_loss+=e*e; }
                if(std::abs(measured_loss-a.loss)>1e-9*(1+a.loss)) throw std::runtime_error("Recovered fitting loss mismatch");
                std::sort(alloc_times.begin(),alloc_times.end()); std::sort(pack_times.begin(),pack_times.end());
                for(const auto& set:eval) {
                    const auto metric=evaluate([&](Vec3 p){return field->query(p);},set.second);
                    results<<c.id<<','<<c.split<<','<<seed<<','<<budget<<','<<methods[m]<<','<<field->bytes()<<','<<a.splits<<','<<a.loss<<','<<alloc_times[1]<<','<<pack_times[1]<<','<<labels+candidate_time+alloc_times[1]+pack_times[1]<<",0,"<<a.transitions<<','<<a.table_states<<','<<a.negative_steps<<','<<set.first<<','<<metric.rmse<<','<<metric.mae<<','<<metric.p99<<','<<metric.sign_error<<','<<metric.inside_recall<<','<<metric.outside_recall<<'\n';
                }
                results.flush(); times.flush();
            }
            std::cout<<prefix<<" complete\n";
        }
        output(dir/"COMPLETE")<<"All cases completed. Generate provenance manifest with allocation_report.py.\n";
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
