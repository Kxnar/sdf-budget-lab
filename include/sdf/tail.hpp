#pragma once
#include "sdf/precision.hpp"

namespace sdf {
enum class Q8Objective { mse, sign, excess };
// Cache fitting errors for every q8 candidate, using the actual packed decoder.
// Sorted errors and prefix sums make repeated excess-loss calculation cheap.
class Q8Losses {
    struct Node {
        std::vector<double> errors, prefix;
        double mse{}, wrong{};
    };
    std::vector<Node> nodes_;
    std::size_t count_;
public:
    Q8Losses(const PrecisionTree& tree,const Samples& fit);
    void apply(PrecisionTree& tree,Q8Objective objective,double parameter=0) const;
    double cvar99(const PrecisionTree& tree,const PrecisionAllocation& a) const;
    double p99(const PrecisionTree& tree,const PrecisionAllocation& a) const;
    std::vector<double> errors(const PrecisionTree& tree,const PrecisionAllocation& a) const;
    std::size_t count() const { return count_; }
};
struct ThresholdResult {
    double threshold{}, excess{}, objective{};
    std::size_t bytes{}, states{}, transitions{};
};
struct TailAllocation {
    PrecisionAllocation allocation;
    double threshold{}, objective{}, fit_cvar99{}, upper{};
    std::vector<ThresholdResult> search;
};
// Exact for each threshold and over this finite threshold set, NOT over all t.
// Uses a uniform grid [0, MSE-field CVaR99] plus the MSE-field p99 threshold.
TailAllocation allocate_cvar99(PrecisionTree& tree,const Q8Losses& losses,
                               std::size_t budget,unsigned grid_points);
} // namespace sdf
