#pragma once
#include "sdf/allocation.hpp"

namespace sdf {
enum class Encoding : unsigned char { branch=0, f32=1, q16=2, q8=3 };
struct EncodedCorners {
    std::array<unsigned char,32> data{};
    Encoding encoding=Encoding::f32;
    std::size_t payload_bytes() const;
    Corners decode() const;
};
EncodedCorners encode_corners(const Corners& values, Encoding encoding);
struct PrecisionNode {
    std::array<EncodedCorners,3> corners;
    std::array<double,3> loss{};
};
struct RefinementCandidate {
    std::size_t node{}, fitting_count{};
    Box box;
    double q8_loss{};
};
struct PrecisionTree {
    CandidateTree spatial;
    std::vector<PrecisionNode> nodes;
    PrecisionTree(const Distance& oracle,const Samples& fit,unsigned depth,Box box={});
    // Append children of the top positive-SSE depth-five cells only. Ranking
    // uses the original fitting pool and decoded q8 residuals, never probes.
    // Call once, on a full depth-five tree, with its original fitting samples.
    std::vector<RefinementCandidate> refine_depth_six(const Distance& oracle,
        const Samples& fit,std::size_t max_parents);
    // Generalization experiment only: rank existing depth-six terminal cells
    // on an explicit allocation/ranking pool, then append depth-seven children.
    std::vector<RefinementCandidate> refine_depth_seven(const Distance& oracle,
        const Samples& ranking,std::size_t max_parents);
};
struct PrecisionAllocation {
    std::vector<Encoding> choice;
    double loss{};
    std::size_t bytes{}, branches{};
    std::array<std::size_t,3> leaves{};
};
// One contiguous allocation: 4-byte tagged descriptors and aligned leaf payloads.
// The tag occupies the low two bits of a byte offset; branches point to eight
// contiguous child descriptors. No internal corner values are retained.
class PackedField {
    Box domain_;
    std::size_t size_{};
    std::unique_ptr<unsigned char[]> data_;
public:
    PackedField(const PrecisionTree& tree,const PrecisionAllocation& allocation);
    // Diagnostic full subdivision, independent of fitting samples and tree DP.
    // Corner queries are made only at terminal cells; duplicates are counted.
    PackedField(const Distance& oracle,unsigned depth,Encoding encoding,
                std::size_t& oracle_calls,Box domain={});
    double query(Vec3 p) const;
    std::size_t bytes() const { return sizeof(*this)+size_; }
    const unsigned char* data() const { return data_.get(); }
    std::size_t payload_bytes() const { return size_; }
};
// Exact Pareto-frontier tree knapsack; no cost rounding beyond exact 4-byte units.
// A fixed encoding restricts leaf choices; branch permits all three encodings.
class PrecisionDP {
    struct State {
        std::size_t units{};
        double loss{};
        Encoding encoding=Encoding::branch;
        std::array<std::size_t,8> children{};
    };
    const PrecisionTree& tree_;
    std::size_t budget_;
    std::vector<std::vector<State>> tables_;
public:
    std::size_t transitions{}, states{};
    PrecisionDP(const PrecisionTree& tree,std::size_t max_bytes,Encoding fixed=Encoding::branch);
    PrecisionAllocation allocate(std::size_t bytes) const;
};
} // namespace sdf
