#pragma once
#include "sdf/field.hpp"

namespace sdf {
struct Metrics {
    std::size_t count{}, sign_count{}, inside_count{}, outside_count{};
    double rmse{}, mae{}, p99{}, cvar99{}, max_error{}, sign_error{}, inside_recall{}, outside_recall{};
};
Metrics evaluate(const Distance& f,const Samples& samples,double boundary_epsilon=2e-9);
struct Timing { double median_ns{}, p95_batch_ns{}, checksum{}; };
Timing time_queries(const Distance& f,const Samples& samples,unsigned repeats=9);
} // namespace sdf
