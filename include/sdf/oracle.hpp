#pragma once
#include "sdf/types.hpp"
#include <memory>
#include <random>
#include <string>

namespace sdf {
class Oracle {
public:
    virtual ~Oracle()=default;
    virtual double distance(Vec3 p) const=0;
    virtual Vec3 surface(std::mt19937_64& rng) const=0;
    virtual std::string name() const=0;
};
double unit(std::mt19937_64& rng);
std::unique_ptr<Oracle> analytic_oracle(const std::string& shape, double feature);
std::unique_ptr<Oracle> mesh_oracle(const std::string& path);
std::unique_ptr<Oracle> mesh_oracle_fast(const std::string& path);
Samples sample_queries(const Oracle& oracle, std::size_t n, std::uint64_t seed,
                       double surface_fraction, double noise, Box domain={});
Samples feature_probes(const Oracle& oracle, const std::string& shape, double feature);
Samples read_probes(const Oracle& oracle, const std::string& path, Box domain={});
} // namespace sdf
