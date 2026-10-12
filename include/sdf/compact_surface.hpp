#pragma once
#include "sdf/compact_mesh.hpp"
namespace sdf {
struct CompactExtraction {
  Surface surface;
  std::vector<double> lattice;
  std::size_t classification_sites{},active_distance_sites{},unresolved_sites{};
};
// Full exact-filtered signs at all lattice sites, distances at every vertex of
// every mixed-sign cube; the original extractor consumes this same lattice.
CompactExtraction extract_compact_surface(const CompactMeshField &field,unsigned cells);
}
