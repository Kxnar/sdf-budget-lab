#include "sdf/oracle.hpp"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#if __has_include(<CGAL/AABB_traits_3.h>)
#include <CGAL/AABB_traits_3.h>
#else
#include <CGAL/AABB_traits.h>
#endif
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

namespace sdf {
namespace {
using Kernel=CGAL::Exact_predicates_exact_constructions_kernel;
using Mesh=CGAL::Surface_mesh<Kernel::Point_3>;
using Primitive=CGAL::AABB_face_graph_triangle_primitive<Mesh>;
#if __has_include(<CGAL/AABB_traits_3.h>)
using Tree=CGAL::AABB_tree<CGAL::AABB_traits_3<Kernel,Primitive>>;
#else
using Tree=CGAL::AABB_tree<CGAL::AABB_traits<Kernel,Primitive>>;
#endif
Vec3 vec(const Kernel::Point_3& p) { return {CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z())}; }
class MeshOracle final:public Oracle {
    Mesh mesh_;
    Tree tree_;
    std::unique_ptr<CGAL::Side_of_triangle_mesh<Mesh,Kernel>> side_;
    std::vector<std::array<Vec3,3>> triangles_;
    std::vector<double> areas_;
public:
    explicit MeshOracle(const std::string& path) {
        if(!CGAL::IO::read_polygon_mesh(path,mesh_) || mesh_.is_empty() || !CGAL::is_triangle_mesh(mesh_) || !CGAL::is_closed(mesh_))
            throw std::invalid_argument("Mesh must be nonempty, closed, and triangulated");
        for(const auto v:mesh_.vertices()) if(!Box{}.contains(vec(mesh_.point(v)))) throw std::invalid_argument("Normalize mesh vertices into [-1,1]^3 before loading");
        double total=0;
        for(const auto face:mesh_.faces()) {
            auto h=mesh_.halfedge(face);
            const auto a=mesh_.point(mesh_.source(h)), b=mesh_.point(mesh_.target(h)), c=mesh_.point(mesh_.target(mesh_.next(h)));
            if(CGAL::collinear(a,b,c)) throw std::invalid_argument("Degenerate triangle in mesh");
            const auto av=vec(a),bv=vec(b),cv=vec(c);
            const double area=norm(cross(bv-av,cv-av))*0.5;
            if(!(area>0) || !std::isfinite(area)) throw std::invalid_argument("Invalid triangle area");
            triangles_.push_back({av,bv,cv}); total+=area; areas_.push_back(total);
        }
        if(CGAL::Polygon_mesh_processing::does_self_intersect(mesh_)) throw std::invalid_argument("Self-intersecting mesh");
        tree_.insert(faces(mesh_).first,faces(mesh_).second,mesh_);
        tree_.build(); tree_.accelerate_distance_queries();
        side_=std::make_unique<CGAL::Side_of_triangle_mesh<Mesh,Kernel>>(tree_);
    }
    double distance(Vec3 p) const override {
        const Kernel::Point_3 point(p.x,p.y,p.z);
        // Exact squared distance to input triangles; final sqrt is rounded to double.
        const double d=std::sqrt(CGAL::to_double(tree_.squared_distance(point)));
        const auto side=(*side_)(point);
        return side==CGAL::ON_BOUNDARY?0.0:(side==CGAL::ON_BOUNDED_SIDE?-d:d);
    }
    Vec3 surface(std::mt19937_64& rng) const override {
        const auto index=static_cast<std::size_t>(std::upper_bound(areas_.begin(),areas_.end(),unit(rng)*areas_.back())-areas_.begin());
        const auto& t=triangles_[std::min(index,triangles_.size()-1)];
        const double u=std::sqrt(unit(rng)),v=unit(rng);
        return t[0]*(1-u)+t[1]*(u*(1-v))+t[2]*(u*v);
    }
    std::string name() const override { return "CGAL:EPECK"; }
};
}
std::unique_ptr<Oracle> mesh_oracle(const std::string& path) { return std::make_unique<MeshOracle>(path); }
} // namespace sdf
