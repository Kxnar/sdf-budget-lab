#include "sdf/compact_mesh.hpp"
#include <map>
#ifdef SDF_WITH_CGAL
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Surface_mesh.h>
#endif
namespace sdf {
CompactDecodedValidation validate_decoded_mesh(const CompactMeshField &field) {
  field.validate_layout();
  auto decoded = field.decoded_geometry();
  Surface surface{decoded.vertices, decoded.triangles};
  auto topology = surface_topology(surface);
  if (!topology.closed_manifold())
    throw std::invalid_argument(
        "Decoded mesh is not a closed edge/vertex manifold");
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::pair<unsigned, int>>
      edges;
  for (auto t : decoded.triangles)
    for (unsigned j = 0; j < 3; ++j) {
      auto a = t[j], b = t[(j + 1) % 3];
      auto &e = edges[std::minmax(a, b)];
      ++e.first;
      e.second += (a < b ? 1 : -1);
    }
  for (auto e : edges)
    if (e.second.first != 2 || e.second.second != 0)
      throw std::invalid_argument(
          "Decoded mesh has inconsistent oriented edge pairs");
  CompactDecodedValidation result{topology, false};
#ifdef SDF_WITH_CGAL
  using K = CGAL::Exact_predicates_exact_constructions_kernel;
  using M = CGAL::Surface_mesh<K::Point_3>;
  M mesh;
  std::vector<M::Vertex_index> vertices;
  for (auto p : decoded.vertices)
    vertices.push_back(mesh.add_vertex({p.x, p.y, p.z}));
  for (auto t : decoded.triangles) {
    if (CGAL::collinear(mesh.point(vertices[t[0]]), mesh.point(vertices[t[1]]),
                        mesh.point(vertices[t[2]])))
      throw std::invalid_argument("Decoded mesh has exact degenerate triangle");
    if (mesh.add_face(vertices[t[0]], vertices[t[1]], vertices[t[2]]) ==
        M::null_face())
      throw std::invalid_argument(
          "Decoded topology cannot form a surface mesh");
  }
  if (CGAL::Polygon_mesh_processing::does_self_intersect(mesh))
    throw std::invalid_argument("Decoded mesh self-intersects");
  result.self_intersection_checked = true;
#endif
  return result;
}
} // namespace sdf
