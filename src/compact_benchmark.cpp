#include "sdf/compact_mesh.hpp"
#include "sdf/compact_surface.hpp"
#include "sdf/generalization.hpp"
#include "sdf/metrics.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#endif
using namespace sdf;
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point t) {
  return std::chrono::duration<double>(Clock::now() - t).count();
}
std::ofstream output(const fs::path &p) {
  std::ofstream f(p);
  f.exceptions(std::ios::badbit | std::ios::failbit);
  f << std::setprecision(17);
  return f;
}
template <class T> void binary(const fs::path &p, const std::vector<T> &v) {
  std::ofstream f(p, std::ios::binary);
  f.exceptions(std::ios::badbit | std::ios::failbit);
  f.write(reinterpret_cast<const char *>(v.data()),
          std::streamsize(v.size() * sizeof(T)));
}
Samples samples(const fs::path &p) {
  std::ifstream f(p, std::ios::binary | std::ios::ate);
  if (!f || std::size_t(f.tellg()) % sizeof(Sample))
    throw std::runtime_error("Invalid query samples");
  Samples v(std::size_t(f.tellg()) / sizeof(Sample));
  f.seekg(0);
  f.read(reinterpret_cast<char *>(v.data()),
         std::streamsize(v.size() * sizeof(Sample)));
  for (auto s : v)
    if (!Box{}.contains(s.point) || !std::isfinite(s.distance))
      throw std::runtime_error("Invalid query coordinate/label");
  return v;
}
std::size_t peak() {
#ifdef _WIN32
  PROCESS_MEMORY_COUNTERS p{};
  if (!GetProcessMemoryInfo(GetCurrentProcess(), &p, sizeof(p)))
    throw std::runtime_error("Memory query failed");
  return p.PeakWorkingSetSize;
#else
  return 0;
#endif
}
void guard() {
  if (peak() > 1610612736ull)
    throw std::runtime_error(
        "resource_unresolved: exact 1.5GiB process guard exceeded");
}
struct Ray {
  int group;
  Vec3 o, d;
  double low, high;
};
std::vector<Ray> rays(const fs::path &p) {
  std::ifstream f(p);
  std::vector<Ray> v;
  Ray r;
  while (f >> r.group >> r.o.x >> r.o.y >> r.o.z >> r.d.x >> r.d.y >> r.d.z >>
         r.low >> r.high)
    v.push_back(r);
  if (v.empty())
    throw std::runtime_error("No rays");
  return v;
}
int main(int argc, char **argv) {
  try {
    if (argc < 5)
      throw std::invalid_argument(
          "Usage: sdf_compact_bench SPEC REF OUT CODEC [--pilot count]");
    auto spec = GeometrySpec::read(argv[1]);
    fs::path ref = argv[2], dir = argv[3];
    auto codec = MeshCodec(std::stoul(argv[4]));
    if (fs::exists(dir))
      throw std::invalid_argument(
          "Output exists; preserve incomplete attempts");
    fs::create_directories(dir);
    std::size_t limit =
        argc > 6 && std::string(argv[5]) == "--pilot" ? std::stoul(argv[6]) : 0;
    auto start = Clock::now();
    auto original = read_triangle_obj(spec.mesh);
    double source_seconds = seconds(start);
    guard();
    const auto source_peak=peak();
    start = Clock::now();
    CompactBuildPhases phases;
    std::vector<unsigned char> encoded_evidence;std::unique_ptr<CompactMeshField> owned;
    try {owned=std::make_unique<CompactMeshField>(original,codec,&phases,&encoded_evidence);}
    catch(const std::invalid_argument &e){if(!encoded_evidence.empty())binary(dir/"field.bin",encoded_evidence);
      output(dir/"STATUS")<<"decoded_invalid\n"<<e.what()<<'\n';
      output(dir/"construction.json")<<"{\"object_bytes\":"<<sizeof(CompactMeshField)<<",\"buffer_bytes\":"<<encoded_evidence.size()<<",\"actual_bytes\":"<<sizeof(CompactMeshField)+encoded_evidence.size()<<"}";
      output(dir/"COMPLETE")<<"decoded_invalid\n";return 0;}
    auto &field=*owned;encoded_evidence={};
    double build_seconds = seconds(start);
    guard();
    const auto encoded_peak=peak();
    {
      std::ofstream f(dir / "field.bin", std::ios::binary);
      f.write(reinterpret_cast<const char *>(field.data()),
              std::streamsize(field.payload_bytes()));
    }
    auto decoded = field.decoded_geometry();
    output(dir/"construction.json")<<"{\"object_bytes\":"<<sizeof(CompactMeshField)<<",\"buffer_bytes\":"<<field.payload_bytes()<<",\"actual_bytes\":"<<field.bytes()<<",\"source_seconds\":"<<source_seconds<<",\"encoding_bvh_packing_seconds\":"<<build_seconds<<",\"peak_working_set_bytes\":"<<peak()<<"}";
    Surface{decoded.vertices, decoded.triangles}.write_obj(
        (dir / "decoded.obj").string());
    start = Clock::now();
    CompactDecodedValidation validation;
    try {
      validation = validate_decoded_mesh(field);
      std::ifstream expected_file(ref/"expected_euler.txt");std::vector<long long> expected;long long value;
      while(expected_file>>value)expected.push_back(value);
      if(expected.empty()||validation.topology.component_euler!=expected)throw std::invalid_argument("Decoded component Euler differs from frozen original topology");
      if (!validation.self_intersection_checked)
        throw std::runtime_error("CGAL exact decoded validation unavailable");
    } catch (const std::invalid_argument &e) {
      output(dir/"construction.json")<<"{\"object_bytes\":"<<sizeof(CompactMeshField)<<",\"buffer_bytes\":"<<field.payload_bytes()<<",\"actual_bytes\":"<<field.bytes()<<",\"source_seconds\":"<<source_seconds<<",\"encoding_bvh_packing_seconds\":"<<build_seconds<<",\"decoded_validation_seconds\":"<<seconds(start)<<",\"source_cumulative_peak_working_set_bytes\":"<<source_peak<<",\"encoded_cumulative_peak_working_set_bytes\":"<<encoded_peak<<",\"validated_construction_cumulative_peak_working_set_bytes\":"<<peak()<<",\"total_benchmark_peak_working_set_bytes\":"<<peak()<<"}";
      output(dir / "STATUS") << "decoded_invalid\n" << e.what() << '\n';
      output(dir / "COMPLETE") << "decoded_invalid\n";
      return 0;
    }
    double validation_seconds = seconds(start);
    guard();
    const auto validated_peak=peak();
    // The original oracle is evaluation-only and never retained by the field.
    start = Clock::now();
    GeneralizationOracle oracle(spec);
    double oracle_setup_seconds = seconds(start);
    original = {};
    decoded = {};
    Distance q = [&](Vec3 p) { return field.query(p); };
    auto metrics = output(dir / "metrics.csv");
    metrics << "set,count,unresolved,nonfinite,rmse,cvar99,max_error,seconds\n";
    bool scalar = true, unresolved = false;
    for (std::string name : {"uniform", "near", "probes", "audit"}) {
      auto s = samples(ref / (name + ".bin"));
      if (limit && (name == "uniform" || name == "near") && s.size() > limit)
        s.resize(limit);
      std::vector<double> values;
      std::vector<unsigned char> flags;
      std::size_t bad = 0, nonfinite = 0;
      start = Clock::now();
      for (auto x : s) {
        auto r = field.evaluate(x.point);
        values.push_back(r.distance);
        flags.push_back(static_cast<unsigned char>((r.resolved ? 1 : 0) |
                                                   (r.boundary ? 2 : 0)));
        bad += !r.resolved;
        nonfinite += !std::isfinite(r.distance);
      }
      auto elapsed = seconds(start);
      guard();
      binary(dir / (name + "_pred.bin"), values);
      binary(dir / (name + "_flags.bin"), flags);
      if (bad || nonfinite) {
        metrics << name << ',' << s.size() << ',' << bad << ',' << nonfinite
                << ",nan,nan,nan," << elapsed << '\n';
        scalar = false;
        unresolved = true;
      } else {
        std::size_t i = 0;
        auto m = evaluate([&](Vec3) { return values[i++]; }, s);
        metrics << name << ',' << s.size() << ",0,0," << m.rmse << ','
                << m.cvar99 << ',' << m.max_error << ',' << elapsed << '\n';
        scalar &= m.max_error <= .020 &&
                  (name != "uniform" || m.cvar99 <= .015) &&
                  (name != "near" || m.cvar99 <= .010);
      }
      std::cout << name << " evaluated " << s.size() << " unresolved " << bad
                << '\n'
                << std::flush;
    }
    // Both sides of decoded faces/edges/vertices. Keep source labels and raw
    // flags.
    start = Clock::now();
    auto geometry = field.decoded_geometry();
    Samples boundary;
    const auto count = std::min<std::size_t>(512, geometry.triangles.size());
    for (std::size_t k = 0; k < count; ++k) {
      auto t = geometry.triangles[k * geometry.triangles.size() / count];
      auto a = geometry.vertices[t[0]], b = geometry.vertices[t[1]],
           c = geometry.vertices[t[2]], n = cross(b - a, c - a);
      n = n * (1 / norm(n));
      for (auto p : {(a + b + c) * (1. / 3), a, (a + b) * .5})
        for (double side : {-1e-9, 1e-9}) {
          auto x = p + n * side;
          if (Box{}.contains(x))
            boundary.push_back({x, oracle.distance(x)});
        }
    }
    auto audit_label_seconds = seconds(start);
    binary(dir / "boundary.bin", boundary);
    std::vector<double> bv;
    std::vector<unsigned char> bf;
    double maximum = 0;
    std::size_t bad = 0;
    for (auto s : boundary) {
      auto r = field.evaluate(s.point);
      bv.push_back(r.distance);
      bf.push_back(static_cast<unsigned char>((r.resolved ? 1 : 0) |
                                              (r.boundary ? 2 : 0)));
      if (!r.resolved || !std::isfinite(r.distance)) {
        ++bad;
        unresolved = true;
      } else
        maximum = std::max(maximum, std::abs(r.distance - s.distance));
    }
    binary(dir / "boundary_pred.bin", bv);
    binary(dir / "boundary_flags.bin", bf);
    scalar &= bad == 0 && maximum <= .020;
    metrics << "boundary," << boundary.size() << ',' << bad << ',' << bad
            << ",nan,nan," << maximum << ",nan\n";
    auto ray_rows = output(dir / "roots.csv");
    ray_rows << "ray,group,numerical_unresolved,coplanar,shared_events,"
                "triangle_tests,roots,tangencies\n";
    auto witnesses = output(dir / "ray_witnesses.csv");
    witnesses << "ray,t,prediction,resolved,boundary\n";
    auto rs = rays(ref / "rays.txt");
    std::ifstream rr(ref / "reference_roots.csv");
    std::string line;
    std::getline(rr, line);
    std::size_t ri = 0;
    start = Clock::now();
    for (auto r : rs) {
      auto events = field.ray(r.o, r.d, r.low, r.high);
      ray_rows << ri << ',' << r.group << ',' << events.unresolved << ','
               << events.coplanar << ',' << events.exact_shared_events << ','
               << events.triangle_tests << ',';
      for (std::size_t j = 0; j < events.transitions.size(); ++j)
        ray_rows << (j ? ";" : "") << events.transitions[j];
      ray_rows << ',';
      for (std::size_t j = 0; j < events.tangencies.size(); ++j)
        ray_rows << (j ? ";" : "") << events.tangencies[j];
      ray_rows << '\n';
      std::getline(rr, line);
      auto p = line.find_last_of(',');
      std::string values = line.substr(p + 1);
      std::replace(values.begin(), values.end(), ';', ' ');
      std::istringstream v(values);
      double t;
      std::vector<double> original_roots;
      while (v >> t)
        original_roots.push_back(t);
      for (auto roots : {original_roots, events.transitions}) {
        std::vector<double> cuts{r.low};
        cuts.insert(cuts.end(), roots.begin(), roots.end());
        cuts.push_back(r.high);
        for (std::size_t j = 1; j < cuts.size(); ++j) {
          auto m = (cuts[j - 1] + cuts[j]) * .5;
          auto answer = field.evaluate(r.o + r.d * m);
          witnesses << ri << ',' << m << ',' << answer.distance << ','
                    << answer.resolved << ',' << answer.boundary << '\n';
        }
      }
      for (auto t : {r.low, r.high}) {
        auto answer = field.evaluate(r.o + r.d * t);
        witnesses << ri << ',' << t << ',' << answer.distance << ','
                  << answer.resolved << ',' << answer.boundary << '\n';
      }
      unresolved |= events.unresolved;
      ++ri;
      guard();
    }
    double ray_seconds = seconds(start);
    auto topology = output(dir / "topology.csv");
    topology << "cells,closed,components,euler,boundary_edges,nonmanifold_"
                "edges,nonmanifold_vertices,component_euler,seconds\n";
    bool surface_check=false;for(int i=5;i<argc;++i)surface_check |= std::string(argv[i])=="--surface-check";
    auto extraction=output(dir/"extraction.csv");extraction<<"cells,classification_sites,active_distance_sites,unresolved_sites,accelerated_seconds,full_query_seconds,byte_identical\n";
    if ((!limit||surface_check) && !unresolved) {
      // Only independently qualified reference extraction levels are measured.
      // levels.txt retains the entire qualification ladder, including unneeded
      // or unresolved levels, so it is not an acceptance-level manifest.
      std::ifstream reference_topology(ref/"qualified_levels.txt");
      std::vector<unsigned> qualified;
      unsigned cells;while(reference_topology>>cells)qualified.push_back(cells);
      if(!qualified.empty()&&qualified.size()!=2)throw std::runtime_error("Invalid frozen qualified extraction levels");
      unsigned n;
      for(auto level:qualified) { n=level;
        start = Clock::now();
        auto accelerated=extract_compact_surface(field,n);auto accelerated_seconds=seconds(start);
        if(accelerated.unresolved_sites){unresolved=true;extraction<<n<<','<<accelerated.classification_sites<<','<<accelerated.active_distance_sites<<','<<accelerated.unresolved_sites<<','<<accelerated_seconds<<",nan,0\n";break;}
        auto mesh=std::move(accelerated.surface);double full_seconds=std::numeric_limits<double>::quiet_NaN();bool identical=false;
        if(surface_check){start=Clock::now();auto full=extract_surface(q,n);full_seconds=seconds(start);
          identical=full.triangles==mesh.triangles && full.vertices.size()==mesh.vertices.size();
          if(identical)for(std::size_t i=0;i<full.vertices.size();++i)if(std::memcmp(&full.vertices[i],&mesh.vertices[i],sizeof(Vec3))){identical=false;break;}
          full.write_obj((dir/("full_surface_"+std::to_string(n)+".obj")).string());
          if(!identical)throw std::runtime_error("Accelerated/full-query surfaces differ");}
        extraction<<n<<','<<accelerated.classification_sites<<','<<accelerated.active_distance_sites<<','<<accelerated.unresolved_sites<<','<<accelerated_seconds<<','<<full_seconds<<','<<(surface_check?(identical?"1":"0"):"not_checked")<<'\n';
        auto result = surface_topology(mesh);
        topology << n << ',' << result.closed_manifold() << ','
                 << result.components << ',' << result.euler << ','
                 << result.boundary_edges << ',' << result.nonmanifold_edges
                 << ',' << result.nonmanifold_vertices << ',';
        for (std::size_t j = 0; j < result.component_euler.size(); ++j)
          topology << (j ? ";" : "") << result.component_euler[j];
        topology << ',' << seconds(start) << '\n';
        mesh.write_obj(
            (dir / ("surface_" + std::to_string(n) + ".obj")).string());
        guard();
        std::cout << "extraction " << n << " completed\n" << std::flush;
      }
    }
    auto timing = output(dir / "timing.csv");
    timing << "repeat,count,ns,checksum\n";
    auto warm = samples(ref / "uniform.bin");
    warm.resize(std::min<std::size_t>(4096, warm.size()));
    if (!unresolved) {
      for (auto s : warm)
        q(s.point);
      for (unsigned i = 0; i < 9; ++i) {
        start = Clock::now();
        double checksum = 0;
        for (auto s : warm)
          checksum += q(s.point);
        timing << i << ',' << warm.size() << ','
               << seconds(start) * 1e9 / warm.size() << ',' << checksum << '\n';
      }
    }
    auto st = field.statistics();
    output(dir / "construction.json")
        << "{\"source_seconds\":" << source_seconds
        << ",\"encoding_bvh_packing_seconds\":" << build_seconds
        << ",\"decoded_validation_seconds\":" << validation_seconds
        << ",\"oracle_setup_seconds\":" << oracle_setup_seconds
        << ",\"boundary_label_seconds\":" << audit_label_seconds
        << ",\"ray_seconds\":" << ray_seconds
        << ",\"oracle_audit_labels\":" << boundary.size()
        << ",\"object_bytes\":" << sizeof(CompactMeshField)
        << ",\"buffer_bytes\":" << field.payload_bytes()
        << ",\"reorder_seconds\":"<<phases.reorder_seconds<<",\"encoding_seconds\":"<<phases.encoding_seconds
        << ",\"bvh_seconds\":"<<phases.bvh_seconds<<",\"packing_seconds\":"<<phases.packing_seconds
        << ",\"layout_validation_seconds\":"<<phases.layout_validation_seconds
        << ",\"actual_bytes\":" << field.bytes()
        << ",\"peak_working_set_bytes\":" << peak()
        << ",\"total_benchmark_peak_working_set_bytes\":"<<peak()
        << ",\"source_cumulative_peak_working_set_bytes\":"<<source_peak
        << ",\"encoded_cumulative_peak_working_set_bytes\":"<<encoded_peak
        << ",\"validated_construction_cumulative_peak_working_set_bytes\":"<<validated_peak
        << ",\"vertices\":" << field.vertex_count()
        << ",\"faces\":" << field.triangle_count()
        << ",\"nodes\":" << field.node_count() << ",\"queries\":" << st.queries
        << ",\"classification_queries\":"<<st.classification_queries
        << ",\"nearest_nodes\":" << st.nearest_nodes
        << ",\"distance_triangles\":" << st.distance_triangles
        << ",\"parity_nodes\":" << st.parity_nodes
        << ",\"parity_triangles\":" << st.parity_triangles
        << ",\"exact_fallbacks\":" << st.exact_fallbacks
        << ",\"direction_retries\":" << st.direction_retries
        << ",\"unresolved\":" << st.unresolved
        << ",\"pilot\":" << (limit ? "true" : "false")
        << ",\"scalar_pass\":" << (scalar ? "true" : "false") << "}";
    output(dir / "STATUS") << (unresolved ? "numerical_unresolved"
                                          : "evaluated")
                           << '\n';
    metrics.flush();ray_rows.flush();witnesses.flush();topology.flush();timing.flush();extraction.flush();
    {auto done=output(dir/"COMPLETE.tmp");done<<"deterministic compact field evaluated once\n";}
    fs::rename(dir/"COMPLETE.tmp",dir/"COMPLETE");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
