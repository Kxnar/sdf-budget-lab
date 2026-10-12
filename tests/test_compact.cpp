#include "sdf/compact_mesh.hpp"
#include "sdf/compact_surface.hpp"
#include <cstring>
#include <iostream>
#include <filesystem>
#include <fstream>
using namespace sdf;
void check(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}
TriangleGeometry cube(double half = .5) {
  TriangleGeometry mesh;
  for (unsigned c = 0; c < 8; ++c)
    mesh.vertices.push_back(
        {c & 1 ? half : -half, c & 2 ? half : -half, c & 4 ? half : -half});
  mesh.triangles = {{{0, 2, 3}}, {{0, 3, 1}}, {{4, 5, 7}}, {{4, 7, 6}},
                    {{0, 1, 5}}, {{0, 5, 4}}, {{2, 6, 7}}, {{2, 7, 3}},
                    {{0, 4, 6}}, {{0, 6, 2}}, {{1, 3, 7}}, {{1, 7, 5}}};
  return mesh;
}
int main(int argc,char **argv) {
  try {
    TriangleGeometry skinny;
    skinny.vertices = {
        {0, 0, 0},
        {.40720951557159424, .741518440656364, .8091662498190999},
        {.40720951184630394, .7415184332057834, .809166251681745},
        {-.5, 0, 0}};
    skinny.triangles = {{{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}}, {{1, 2, 3}}};
    Vec3 boundary{.20360475685447454, .37075921846553683, .40458312537521124};
    for (auto codec : {MeshCodec::f64, MeshCodec::xor32}) {
      CompactMeshField field(skinny, codec);
      auto answer = field.evaluate(boundary);
      if (!answer.boundary)
        std::cerr << "Skinny answer " << answer.distance << " resolved "
                  << answer.resolved << " fallback "
                  << field.statistics().exact_fallbacks << " codec "
                  << int(codec) << '\n';
      check(answer.resolved && answer.boundary && answer.distance == 0,
            "Skinny exact dyadic boundary");
    }
    // Two distinct binary64 crossings with no representable binary64 midpoint.
    auto thin = cube();
    double x0 = .5, x1 = std::nextafter(.5, 1.);
    for (auto &p : thin.vertices)
      p.x = p.x < 0 ? x0 : x1;
    CompactMeshField wedge(thin, MeshCodec::f64);
    auto narrow = wedge.ray({0, 0, 0}, {1, 0, 0}, .4, .9);
    check(!narrow.unresolved && narrow.transitions.size() == 2 &&
              narrow.transitions[0] == x0 && narrow.transitions[1] == x1,
          "Exact midpoint preserves adjacent representable transitions");
    auto nested = cube();
    auto inner = cube(.25);
    auto offset = nested.vertices.size();
    nested.vertices.insert(nested.vertices.end(), inner.vertices.begin(),
                           inner.vertices.end());
    for (auto t : inner.triangles) {
      for (auto &i : t)
        i += std::uint32_t(offset);
      std::swap(t[1], t[2]);
      nested.triangles.push_back(t);
    }
    CompactMeshField hollow(nested, MeshCodec::f64);
    check(hollow.query({0, 0, 0}) > 0 && hollow.query({.3, 0, 0}) < 0,
          "Nested component parity");
    check(validate_decoded_mesh(hollow).topology.components == 2,
          "Decoded nested topology gate");
    auto all = hollow.ray({0, 0, 0}, {1, 0, 0}, -.9, .9);
    check(!all.unresolved && all.transitions.size() == 4,
          "Nested ray intervals");
    for (auto encoding : {MeshCodec::f64, MeshCodec::f32, MeshCodec::q24,
                          MeshCodec::q28, MeshCodec::xor32}) {
      auto mesh = cube();
      CompactMeshField built(mesh, encoding);
      std::vector<unsigned char> buffer(built.data(),
                                        built.data() + built.payload_bytes());
      mesh = {};
      CompactMeshField loaded(buffer);
      if(argc>1){std::filesystem::create_directories(argv[1]);std::ofstream f(std::filesystem::path(argv[1])/(std::to_string(int(encoding))+".bin"),std::ios::binary);f.write(reinterpret_cast<const char*>(loaded.data()),std::streamsize(loaded.payload_bytes()));}
      check(validate_decoded_mesh(loaded).topology.closed_manifold(),
            "Decoded closure/vertex-link/oriented-edge gate");
      check(loaded.bytes() == sizeof(CompactMeshField) + buffer.size(),
            "Exact retained bytes");
      check(std::abs(loaded.query({0, 0, 0}) + .5) < 2e-7,
            "Inside cube signed distance");
      check(std::abs(loaded.query({.75, .2, .1}) - .25) < 2e-7,
            "Outside cube distance");
      for(unsigned i=0;i<64;++i){Vec3 p{-.9+1.8*(i%4)/3,-.9+1.8*((i/4)%4)/3,-.9+1.8*(i/16)/3};
        Vec3 delta{std::abs(p.x)-.5,std::abs(p.y)-.5,std::abs(p.z)-.5};
        auto expected=norm({std::max(0.,delta.x),std::max(0.,delta.y),std::max(0.,delta.z)})+std::min(0.,std::max({delta.x,delta.y,delta.z}));
        check(std::abs(loaded.query(p)-expected)<2e-7,"Unpruned analytic cube distances and tied nearest faces");}
      check(loaded.evaluate({.5, .5, .5}).boundary, "Exact vertex boundary");
      check(loaded.evaluate({.5, 0, 0}).boundary, "Exact face boundary");
      check(loaded.classify({.5,0,0}).boundary && loaded.classify({0,0,0}).distance<0,"Classification exact boundary/inside");
      for(unsigned resolution:{8u,9u}) {
        auto accelerated=extract_compact_surface(loaded,resolution);auto full=extract_surface([&](Vec3 p){return loaded.query(p);},resolution);
        check(!accelerated.unresolved_sites && full.triangles==accelerated.surface.triangles && full.vertices.size()==accelerated.surface.vertices.size(),"Sign-first extraction topology equivalence");
        for(std::size_t i=0;i<full.vertices.size();++i)check(!std::memcmp(&full.vertices[i],&accelerated.surface.vertices[i],sizeof(Vec3)),"Sign-first extraction exact vertex bits including zero lattice sites");
      }
      auto events = loaded.ray({0, 0, 0}, {1, 0, 0}, -.9, .9);
      check(!events.unresolved && events.transitions.size() == 2,
            "Shared face diagonal exact ray events");
      check(std::abs(events.transitions[0] + .5) < 2e-7, "Entering event");
      auto coplanar = loaded.ray({0, .5, 0}, {1, 0, 0}, -.9, .9);
      check(coplanar.coplanar && coplanar.unresolved, "Coplanar ray explicit");
      auto disjoint = loaded.ray({0, .5, 0}, {1, 0, 0}, .6, .9);
      check(!disjoint.unresolved && disjoint.transitions.empty(),
            "Disjoint coplanar supporting plane");
      bool finite_rejected = false;
      try {
        loaded.ray({std::numeric_limits<double>::infinity(), 0, 0}, {1, 0, 0},
                   -.9, .9);
      } catch (const std::invalid_argument &) {
        finite_rejected = true;
      }
      check(finite_rejected, "Nonfinite ray origin rejected");
      auto tangent = loaded.ray({0, 0, 1}, {1, 1, -1}, .1, .9);
      check(!tangent.unresolved && tangent.transitions.empty() &&
                tangent.tangencies.size() == 1,
            "Isolated exact vertex tangent");
      if (encoding == MeshCodec::xor32) {
        auto malformed = buffer;
        auto set32 = [&](std::size_t p, std::uint32_t x) {
          for (unsigned a = 0; a < 4; ++a)
            malformed[p + a] = static_cast<unsigned char>(x >> (8 * a));
        };
        set32(8, 32);
        set32(84, 37);
        malformed[88] = 1;
        malformed[89] = 32;
        malformed[114] = 63;
        malformed[115] = 64;
        malformed[116] = 64;
        bool rejected = false;
        try {
          CompactMeshField invalid(malformed);
        } catch (const std::invalid_argument &) {
          rejected = true;
        }
        check(rejected, "Short XOR block with large advertised residual "
                        "rejected before tail access");
      }
      buffer.pop_back();
      bool rejected = false;
      try {
        CompactMeshField invalid(buffer);
      } catch (const std::exception &) {
        rejected = true;
      }
      check(rejected, "Truncated buffer rejected");
      // Erase the serialized input as well as the source geometry; the query
      // object must still own everything needed for a standalone answer.
      buffer.clear();buffer.shrink_to_fit();check(loaded.query({0,0,0})<0,"Loaded field independent of destroyed input buffer");
    }
    auto open = cube();
    open.triangles.pop_back();
    CompactMeshField open_field(open, MeshCodec::f64);
    bool open_rejected = false;
    try {
      validate_decoded_mesh(open_field);
    } catch (const std::invalid_argument &) {
      open_rejected = true;
    }
    check(open_rejected, "Open decoded mesh rejected before acceptance");
    TriangleGeometry planar{{{.2, .2, .2}, {.3, .2, .2}, {.2, .3, .2}},
                            {{{0, 1, 2}}}};
    for (auto codec : {MeshCodec::q24, MeshCodec::q28}) {
      CompactMeshField plane(planar, codec);
      auto decoded = plane.decoded_geometry();
      for (auto p : decoded.vertices)
        check(p.z == .2, "Constant quantizer axis");
      auto raw = std::vector<unsigned char>(
          plane.data(), plane.data() + plane.payload_bytes());
      if (codec == MeshCodec::q28) {
        raw[80 + 48 + 31] |= 128;
        bool bad = false;
        try {
          CompactMeshField invalid(raw);
        } catch (const std::invalid_argument &) {
          bad = true;
        }
        check(bad, "Nonzero final Q28 word tail rejected");
      }
    }
    auto crowded = cube();
    crowded.vertices.resize(65536, {.1, .1, .1});
    CompactMeshField short_index(crowded, MeshCodec::f32);
    check(short_index.data()[7] == 2, "65536 index-width boundary");
    crowded.vertices.push_back({.1, .1, .1});
    CompactMeshField long_index(crowded, MeshCodec::f32);
    check(long_index.data()[7] == 4, "65537 index-width boundary");
    TriangleGeometry compressible;
    for(unsigned i=0;i<64;++i)compressible.vertices.push_back({.4+(i%4)*std::ldexp(1.,-32),.4+((i/4)%4)*std::ldexp(1.,-32),.4+(i/16)*std::ldexp(1.,-32)});
    compressible.triangles={{{0,1,4}}};CompactMeshField compressed(compressible,MeshCodec::xor32);
    auto data=compressed.data();auto block_offset=std::uint32_t(data[80])|(std::uint32_t(data[81])<<8)|(std::uint32_t(data[82])<<16)|(std::uint32_t(data[83])<<24);
    check(data[80+block_offset]==1,"Compressible XOR block chooses residual encoding");
    auto decoded=compressed.decoded_geometry();
    auto order=[](Vec3 a,Vec3 b){if(a.x!=b.x)return a.x<b.x;if(a.y!=b.y)return a.y<b.y;return a.z<b.z;};
    std::sort(decoded.vertices.begin(),decoded.vertices.end(),order);std::sort(compressible.vertices.begin(),compressible.vertices.end(),order);
    for(unsigned i=0;i<64;++i)check(!std::memcmp(&decoded.vertices[i],&compressible.vertices[i],sizeof(Vec3)),"Compressed XOR binary64 coordinate bits");
    if(argc>1){std::ofstream f(std::filesystem::path(argv[1])/"xor-compressed.bin",std::ios::binary);f.write(reinterpret_cast<const char*>(data),std::streamsize(compressed.payload_bytes()));}
    std::cout << "Compact codec roundtrip, standalone distance/sign, exact "
                 "boundary/event, and truncated-buffer checks passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
