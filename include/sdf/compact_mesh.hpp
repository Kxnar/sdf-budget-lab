#pragma once
#include "sdf/surface.hpp"
#include "sdf/types.hpp"
#include <memory>
#include <string>

namespace sdf {
enum class MeshCodec : std::uint8_t {
  f64 = 1,
  f32 = 2,
  q24 = 3,
  q28 = 4,
  xor32 = 5
};
struct TriangleGeometry {
  std::vector<Vec3> vertices;
  std::vector<std::array<std::uint32_t, 3>> triangles;
};
TriangleGeometry read_triangle_obj(const std::string &path);
struct CompactQueryStats {
  std::uint64_t queries{}, classification_queries{},nearest_nodes{}, distance_triangles{},
      parity_nodes{}, parity_triangles{}, exact_fallbacks{},
      direction_retries{}, unresolved{};
};
struct MeshQuery {
  double distance{};
  bool resolved{}, boundary{};
};
struct MeshRayEvents {
  std::vector<double> transitions, tangencies;
  bool coplanar{}, unresolved{};
  std::size_t exact_shared_events{}, triangle_tests{};
};
struct CompactBuildPhases {
  double reorder_seconds{},encoding_seconds{},bvh_seconds{},packing_seconds{},layout_validation_seconds{};
};
// Exactly one owning contiguous buffer and counted counters. No source/oracle,
// decoded mesh, BVH object, adjacency or external permutation is retained.
class CompactMeshField {
  std::unique_ptr<unsigned char[]> data_;
  std::size_t size_{};
  mutable CompactQueryStats statistics_{};
  Vec3 vertex(std::uint32_t index) const;
  std::array<std::uint32_t, 3> triangle(std::uint32_t index) const;
  MeshQuery parity(Vec3 point) const;

public:
  CompactMeshField(const TriangleGeometry &mesh, MeshCodec encoding,CompactBuildPhases *phases=nullptr,
                   std::vector<unsigned char> *serialized_evidence=nullptr);
  explicit CompactMeshField(const std::vector<unsigned char> &serialized);
  CompactMeshField(CompactMeshField &&) = default;
  CompactMeshField &operator=(CompactMeshField &&) = default;
  const unsigned char *data() const { return data_.get(); }
  std::size_t payload_bytes() const { return size_; }
  std::size_t bytes() const { return sizeof(*this) + size_; }
  std::uint32_t vertex_count() const;
  std::uint32_t triangle_count() const;
  std::uint32_t node_count() const;
  MeshCodec codec() const;
  TriangleGeometry decoded_geometry() const;
  MeshQuery evaluate(Vec3 point) const;
  MeshQuery classify(Vec3 point) const;
  double query(Vec3 point) const;
  MeshRayEvents ray(Vec3 origin, Vec3 direction, double low, double high) const;
  CompactQueryStats statistics() const { return statistics_; }
  void reset_statistics() const { statistics_ = {}; }
  void validate_layout() const;
};
struct CompactDecodedValidation {
  SurfaceTopology topology;
  bool self_intersection_checked{};
};
// Mandatory acceptance gate outside structural buffer validation. Throws on
// winding/closure/link/degeneracy/domain defects and (with CGAL) intersections.
// OFF builds remain queryable but cannot claim geometric validation of inputs.
CompactDecodedValidation validate_decoded_mesh(const CompactMeshField &field);
} // namespace sdf
