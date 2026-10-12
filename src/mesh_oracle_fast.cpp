#include "sdf/oracle.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
namespace sdf { namespace {
using K=CGAL::Exact_predicates_inexact_constructions_kernel;using M=CGAL::Surface_mesh<K::Point_3>;using P=CGAL::AABB_face_graph_triangle_primitive<M>;using T=CGAL::AABB_tree<CGAL::AABB_traits_3<K,P>>;
Vec3 vec(const K::Point_3& p){return {p.x(),p.y(),p.z()};}
class Fast final:public Oracle {
    M mesh_;T tree_;std::unique_ptr<CGAL::Side_of_triangle_mesh<M,K>> side_;std::vector<std::array<Vec3,3>> triangles_;std::vector<double> areas_;
public:
    explicit Fast(const std::string& path){
        if(!CGAL::IO::read_polygon_mesh(path,mesh_)||mesh_.is_empty()||!CGAL::is_triangle_mesh(mesh_)||!CGAL::is_closed(mesh_))throw std::invalid_argument("Fast oracle input must be closed triangle mesh");
        double area=0;for(auto v:mesh_.vertices())if(!Box{}.contains(vec(mesh_.point(v))))throw std::invalid_argument("Fast oracle vertex out of domain");
        for(auto face:mesh_.faces()){auto h=mesh_.halfedge(face);auto a=mesh_.point(mesh_.source(h)),b=mesh_.point(mesh_.target(h)),c=mesh_.point(mesh_.target(mesh_.next(h)));if(CGAL::collinear(a,b,c))throw std::invalid_argument("Fast oracle degenerate triangle");auto av=vec(a),bv=vec(b),cv=vec(c);triangles_.push_back({av,bv,cv});area+=norm(cross(bv-av,cv-av))*.5;areas_.push_back(area);}
        // Exact EPECK self-intersection/degeneracy validation is mandatory in
        // the separate pre-freeze comparison. This path never repairs input.
        tree_.insert(faces(mesh_).first,faces(mesh_).second,mesh_);tree_.build();tree_.accelerate_distance_queries();side_=std::make_unique<CGAL::Side_of_triangle_mesh<M,K>>(tree_);
    }
    double distance(Vec3 p)const override{K::Point_3 q(p.x,p.y,p.z);auto side=(*side_)(q);double d=std::sqrt(tree_.squared_distance(q));return side==CGAL::ON_BOUNDARY?0:side==CGAL::ON_BOUNDED_SIDE?-d:d;}
    Vec3 surface(std::mt19937_64& rng)const override{auto ix=std::size_t(std::upper_bound(areas_.begin(),areas_.end(),unit(rng)*areas_.back())-areas_.begin());auto t=triangles_[std::min(ix,triangles_.size()-1)];double u=std::sqrt(unit(rng)),v=unit(rng);return t[0]*(1-u)+t[1]*(u*(1-v))+t[2]*(u*v);}
    std::string name()const override{return "CGAL:EPICK-audited";}
};}
std::unique_ptr<Oracle> mesh_oracle_fast(const std::string& path){return std::make_unique<Fast>(path);}
}
