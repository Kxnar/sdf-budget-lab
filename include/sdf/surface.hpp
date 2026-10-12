#pragma once
#include "sdf/field.hpp"
#include <string>

namespace sdf {
struct Surface {
    std::vector<Vec3> vertices;
    std::vector<std::array<std::uint32_t,3>> triangles;
    void write_obj(const std::string& path) const;
};
struct SurfaceTopology {
    std::size_t components{}, boundary_edges{}, nonmanifold_edges{}, nonmanifold_vertices{};
    long long euler{};
    std::vector<long long> component_euler;
    bool closed_manifold() const { return components && !boundary_edges && !nonmanifold_edges && !nonmanifold_vertices; }
};
// Six conforming tetrahedra per cube. Shared lattice-edge vertices are indexed.
Surface extract_surface(const Distance& field,unsigned cells,Box domain={});
SurfaceTopology surface_topology(const Surface& surface);
} // namespace sdf
