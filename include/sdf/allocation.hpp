#pragma once
#include "sdf/field.hpp"

namespace sdf {
// Build-time data only. Each node's SSE uses the SAME routed fitting samples.
struct CandidateNode {
    AdaptiveField::Node field;
    double loss{};
};
struct CandidateTree {
    Box domain;
    std::vector<CandidateNode> nodes;
    std::size_t sample_count{}, oracle_calls{};
    CandidateTree(const Distance& oracle, const Samples& samples, unsigned depth, Box box={});
};
struct Allocation {
    std::vector<unsigned char> split;
    double loss{};
    std::size_t splits{}, transitions{}, table_states{}, negative_steps{};
};
// A cap includes the query object and all nodes, including retained parents.
class AllocatedField {
    Box domain_;
    std::size_t count_{};
    std::unique_ptr<AdaptiveField::Node[]> nodes_;
public:
    AllocatedField(const CandidateTree& tree, const Allocation& allocation);
    double query(Vec3 point) const;
    std::size_t bytes() const { return sizeof(*this)+count_*sizeof(AdaptiveField::Node); }
    static std::size_t budget_for_splits(std::size_t splits);
};
Allocation exact_allocate(const CandidateTree& tree, std::size_t bytes);
// horizon counts the root split AND descendant splits in a committed bundle.
// best_prefix retains the lowest-loss feasible intermediate result.
Allocation lookahead_allocate(const CandidateTree& tree, std::size_t bytes,
                              std::size_t horizon=1, bool best_prefix=true);
} // namespace sdf
