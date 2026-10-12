#include "sdf/surface.hpp"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <unordered_map>

namespace sdf {
namespace {
std::uint64_t edge(std::uint32_t a,std::uint32_t b) {
    if(a>b) std::swap(a,b);
    return (std::uint64_t(a)<<32)|b;
}
struct DSU {
    std::vector<std::size_t> p;
    explicit DSU(std::size_t n):p(n) { std::iota(p.begin(),p.end(),0); }
    std::size_t root(std::size_t a) { while(a!=p[a]) { p[a]=p[p[a]]; a=p[a]; } return a; }
    void join(std::size_t a,std::size_t b) { p[root(a)]=root(b); }
};
}
Surface extract_surface(const Distance& field,unsigned cells,Box domain) {
    domain.validate(); if(cells<2 || cells>512) throw std::invalid_argument("Extraction cells must be 2..512");
    const unsigned n=cells+1;
    auto id=[&](unsigned x,unsigned y,unsigned z){return (z*n+y)*n+x;};
    auto point=[&](unsigned i) { return domain.point({double(i%n)/cells,double((i/n)%n)/cells,double(i/(n*n))/cells}); };
    std::vector<double> values(std::size_t(n)*n*n);
    for(unsigned i=0;i<values.size();++i) {
        values[i]=field(point(i));
        if(!std::isfinite(values[i])) throw std::invalid_argument("Nonfinite extraction value");
        // Deterministic infinitesimal outside perturbation avoids zero-level vertices.
        if(values[i]==0) values[i]=1e-14*domain.scale();
    }
    Surface mesh; std::unordered_map<std::uint64_t,std::uint32_t> edges;
    auto vertex=[&](unsigned a,unsigned b) {
        const auto key=edge(a,b); const auto it=edges.find(key); if(it!=edges.end()) return it->second;
        const auto index=static_cast<std::uint32_t>(mesh.vertices.size());
        const double u=values[a]/(values[a]-values[b]);
        mesh.vertices.push_back(point(a)*(1-u)+point(b)*u); edges.emplace(key,index); return index;
    };
    constexpr unsigned tetra[6][4]={{0,1,3,7},{0,3,2,7},{0,2,6,7},{0,6,4,7},{0,4,5,7},{0,5,1,7}};
    for(unsigned z=0;z<cells;++z) for(unsigned y=0;y<cells;++y) for(unsigned x=0;x<cells;++x) {
        unsigned cube[8]; bool negative=false,positive=false;
        for(unsigned c=0;c<8;++c) { cube[c]=id(x+(c&1),y+((c>>1)&1),z+((c>>2)&1)); negative|=values[cube[c]]<0; positive|=values[cube[c]]>0; }
        if(!negative || !positive) continue;
        for(const auto& tet:tetra) {
            unsigned inside[4],outside[4],ni=0,no=0;
            for(auto c:tet) { const auto v=cube[c]; if(values[v]<0) inside[ni++]=v; else outside[no++]=v; }
            if(ni==0 || ni==4) continue;
            auto triangle=[&](unsigned a,unsigned b,unsigned c) {
                // Outward winding from the negative centroid toward the positive side.
                Vec3 in{},out{};
                for(unsigned j=0;j<ni;++j) in=in+point(inside[j])*(1./ni);
                for(unsigned j=0;j<no;++j) out=out+point(outside[j])*(1./no);
                if(dot(cross(mesh.vertices[b]-mesh.vertices[a],mesh.vertices[c]-mesh.vertices[a]),out-in)<0) std::swap(b,c);
                mesh.triangles.push_back({a,b,c});
            };
            if(ni==1) triangle(vertex(inside[0],outside[0]),vertex(inside[0],outside[1]),vertex(inside[0],outside[2]));
            else if(ni==3) triangle(vertex(outside[0],inside[0]),vertex(outside[0],inside[1]),vertex(outside[0],inside[2]));
            else {
                const auto a=vertex(inside[0],outside[0]),b=vertex(inside[0],outside[1]),c=vertex(inside[1],outside[0]),d=vertex(inside[1],outside[1]);
                triangle(a,b,c); triangle(b,d,c);
            }
        }
    }
    return mesh;
}
SurfaceTopology surface_topology(const Surface& mesh) {
    SurfaceTopology result; DSU components(mesh.vertices.size());
    std::unordered_map<std::uint64_t,unsigned> edges;
    std::vector<std::vector<std::array<unsigned,2>>> links(mesh.vertices.size());
    for(const auto& t:mesh.triangles) for(unsigned c=0;c<3;++c) {
        if(t[c]>=mesh.vertices.size()) throw std::invalid_argument("Invalid triangle index");
        ++edges[edge(t[c],t[(c+1)%3])]; components.join(t[c],t[(c+1)%3]);
        links[t[c]].push_back({t[(c+1)%3],t[(c+2)%3]});
    }
    std::unordered_map<std::size_t,long long> chi;
    for(std::size_t v=0;v<mesh.vertices.size();++v) if(!links[v].empty()) {
        ++chi[components.root(v)];
        // A closed manifold vertex has a single cycle as its link.
        std::unordered_map<unsigned,std::vector<unsigned>> adj;
        for(auto p:links[v]) { adj[p[0]].push_back(p[1]); adj[p[1]].push_back(p[0]); }
        bool valid=true; for(const auto& p:adj) valid&=p.second.size()==2;
        if(valid) {
            std::vector<unsigned> stack{adj.begin()->first}; std::unordered_map<unsigned,bool> seen;
            while(!stack.empty()) { const auto a=stack.back(); stack.pop_back(); if(seen[a]) continue; seen[a]=true; for(auto b:adj[a]) if(!seen[b]) stack.push_back(b); }
            valid=seen.size()==adj.size();
        }
        if(!valid) ++result.nonmanifold_vertices;
    }
    for(const auto& p:edges) { --chi[components.root(unsigned(p.first>>32))]; if(p.second==1) ++result.boundary_edges; else if(p.second!=2) ++result.nonmanifold_edges; }
    for(const auto& t:mesh.triangles) ++chi[components.root(t[0])];
    result.components=chi.size();
    for(const auto& p:chi) { result.euler+=p.second; result.component_euler.push_back(p.second); }
    std::sort(result.component_euler.begin(),result.component_euler.end()); return result;
}
void Surface::write_obj(const std::string& path) const {
    std::ofstream f(path); f.exceptions(std::ios::badbit|std::ios::failbit); f<<std::setprecision(17);
    for(auto p:vertices) f<<"v "<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    for(auto t:triangles) f<<"f "<<t[0]+1<<' '<<t[1]+1<<' '<<t[2]+1<<'\n';
}
} // namespace sdf
