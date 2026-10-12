#pragma once
#include "sdf/oracle.hpp"
#include <functional>
#include <memory>

namespace sdf {
using Distance=std::function<double(Vec3)>;
class UniformGrid {
    Box domain_;
    std::size_t n_{};
    std::unique_ptr<float[]> values_;
public:
    UniformGrid(const Distance& reference, std::size_t budget, Box domain={});
    double query(Vec3 p) const;
    std::size_t bytes() const { return sizeof(*this)+n_*n_*n_*sizeof(float); }
    std::size_t resolution() const { return n_; }
    const float* data() const { return values_.get(); }
};
struct AdaptiveOptions { unsigned max_depth=10; double sign_penalty=0.0025; };
class AdaptiveField {
public:
    // Internal nodes retain their corners; ALL nodes count toward the byte cap.
    struct Node { Corners values{}; std::uint32_t first_child=0; };
private:
    Box domain_;
    std::size_t count_{};
    std::unique_ptr<Node[]> nodes_;
public:
    AdaptiveField(const Distance& reference, const Samples& training, std::size_t budget,
                  AdaptiveOptions options={}, Box domain={});
    double query(Vec3 p) const;
    std::size_t bytes() const { return sizeof(*this)+count_*sizeof(Node); }
    std::size_t node_count() const { return count_; }
    std::size_t leaf_count() const { return (count_-1)/8*7+1; }
    std::vector<Box> leaves() const;
};
} // namespace sdf
