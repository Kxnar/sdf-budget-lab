#include "sdf/compact_mesh.hpp"
#include <algorithm>
#include <chrono>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <sstream>

namespace sdf {
namespace {
using Buffer = std::vector<unsigned char>;
using Rational = boost::multiprecision::cpp_rational;
using Integer = boost::multiprecision::cpp_int;
constexpr std::size_t header = 80;
std::uint64_t integer(const unsigned char *data, std::size_t bytes) {
  std::uint64_t v = 0;
  for (std::size_t i = 0; i < bytes; ++i)
    v |= std::uint64_t(data[i]) << (8 * i);
  return v;
}
void put(Buffer &b, std::size_t offset, std::uint64_t value,
         std::size_t bytes) {
  for (std::size_t i = 0; i < bytes; ++i)
    b.at(offset + i) = static_cast<unsigned char>(value >> (8 * i));
}
void append(Buffer &b, std::uint64_t value, std::size_t bytes) {
  auto offset = b.size();
  b.resize(offset + bytes);
  put(b, offset, value, bytes);
}
std::uint64_t bits(double value) {
  std::uint64_t r;
  std::memcpy(&r, &value, 8);
  return r;
}
double floating(const unsigned char *p) {
  auto v = integer(p, 8);
  double d;
  std::memcpy(&d, &v, 8);
  return d;
}
float single(const unsigned char *p) {
  auto v = static_cast<std::uint32_t>(integer(p, 4));
  float f;
  std::memcpy(&f, &v, 4);
  return f;
}
void append_double(Buffer &b, double value) { append(b, bits(value), 8); }
void append_single(Buffer &b, double value) {
  float f = static_cast<float>(value);
  std::uint32_t v;
  std::memcpy(&v, &f, 4);
  append(b, v, 4);
}
void align(Buffer &b, std::size_t alignment) {
  while (b.size() % alignment)
    b.push_back(0);
}
void bit_write(Buffer &b, std::size_t start, std::size_t &cursor,
               std::uint64_t value, unsigned width) {
  for (unsigned j = 0; j < width; ++j, ++cursor)
    if ((value >> j) & 1)
      b.at(start + cursor / 8) |=
          static_cast<unsigned char>(1u << (cursor % 8));
}
std::uint64_t bit_read(const unsigned char *p, std::size_t &cursor,
                       unsigned width) {
  std::uint64_t v = 0;
  unsigned written = 0;
  while (written < width) {
    unsigned offset = unsigned(cursor % 8),
             take = std::min(8 - offset, width - written);
    auto chunk = (p[cursor / 8] >> offset) & ((1u << take) - 1);
    v |= std::uint64_t(chunk) << written;
    cursor += take;
    written += take;
  }
  return v;
}
unsigned bit_width(std::uint64_t value) {
  unsigned w = 0;
  while (value) {
    ++w;
    value >>= 1;
  }
  return w;
}
Box domain(const unsigned char *p) {
  return {{floating(p + 32), floating(p + 40), floating(p + 48)},
          {floating(p + 56), floating(p + 64), floating(p + 72)}};
}
std::uint32_t count(const unsigned char *p, unsigned i) {
  return static_cast<std::uint32_t>(integer(p + 8 + 4 * i, 4));
}
std::uint32_t section(const unsigned char *p, unsigned i) {
  return static_cast<std::uint32_t>(integer(p + 20 + 4 * i, 4));
}
Vec3 decode_vertex(const unsigned char *p, std::uint32_t index) {
  auto codec = MeshCodec(p[6]);
  auto v = count(p, 0);
  if (index >= v)
    throw std::invalid_argument("Vertex index out of range");
  auto start = section(p, 0);
  Vec3 result;
  if (codec == MeshCodec::f64 || codec == MeshCodec::f32) {
    auto stride = codec == MeshCodec::f64 ? 24 : 12;
    for (unsigned a = 0; a < 3; ++a)
      (&result.x)[a] = stride == 24
                           ? floating(p + start + stride * index + 8 * a)
                           : single(p + start + stride * index + 4 * a);
  } else if (codec == MeshCodec::q24 || codec == MeshCodec::q28) {
    unsigned width = codec == MeshCodec::q24 ? 24 : 28;
    std::size_t cursor = std::size_t(index) * 3 * width;
    for (unsigned a = 0; a < 3; ++a)
      (&result.x)[a] = floating(p + start + 8 * a) +
                       floating(p + start + 24 + 8 * a) *
                           double(bit_read(p + start + 48, cursor, width));
  } else if (codec == MeshCodec::xor32) {
    auto block = index / 32, local = index % 32;
    auto offset = integer(p + start + 4 * block, 4);
    auto q = p + start + offset;
    auto mode = q[0], n = q[1];
    if (local >= n)
      throw std::invalid_argument("Invalid XOR block count");
    if (mode == 0) {
      for (unsigned a = 0; a < 3; ++a)
        (&result.x)[a] = floating(q + 2 + 24 * local + 8 * a);
    } else {
      std::array<std::uint64_t, 3> base{integer(q + 2, 8), integer(q + 10, 8),
                                        integer(q + 18, 8)};
      std::size_t cursor =
          local ? (local - 1) * std::size_t(q[26] + q[27] + q[28]) : 0;
      if (local)
        for (unsigned a = 0; a < 3; ++a) {
          auto residual = bit_read(q + 29, cursor, q[26 + a]);
          auto value = base[a] ^ residual;
          std::memcpy((&result.x) + a, &value, 8);
        }
      if (local == 0)
        for (unsigned a = 0; a < 3; ++a)
          std::memcpy((&result.x) + a, &base[a], 8);
    }
  } else
    throw std::invalid_argument("Unknown vertex codec");
  return result;
}
std::uint64_t morton(Vec3 p, Box box) {
  std::uint32_t q[3]{};
  for (unsigned a = 0; a < 3; ++a) {
    auto w = box.hi[a] - box.lo[a];
    if (w)
      q[a] = static_cast<std::uint32_t>(
          std::clamp((p[a] - box.lo[a]) / w, 0., 1.) * 2097151.);
  }
  std::uint64_t code = 0;
  for (unsigned bit = 0; bit < 21; ++bit)
    for (unsigned a = 0; a < 3; ++a)
      code |= std::uint64_t((q[a] >> bit) & 1) << (3 * bit + a);
  return code;
}
Box bounds(const std::vector<Vec3> &vertices) {
  Box b{vertices[0], vertices[0]};
  for (auto p : vertices)
    for (unsigned a = 0; a < 3; ++a) {
      (&b.lo.x)[a] = std::min(b.lo[a], p[a]);
      (&b.hi.x)[a] = std::max(b.hi[a], p[a]);
    }
  return b;
}
Buffer vertex_buffer(const std::vector<Vec3> &vertices, MeshCodec encoding) {
  Buffer result;
  auto range = bounds(vertices);
  if (encoding == MeshCodec::f64 || encoding == MeshCodec::f32) {
    for (auto p : vertices)
      for (unsigned a = 0; a < 3; ++a) {
        if (encoding == MeshCodec::f64)
          append_double(result, p[a]);
        else
          append_single(result, p[a]);
      }
  } else if (encoding == MeshCodec::q24 || encoding == MeshCodec::q28) {
    unsigned width = encoding == MeshCodec::q24 ? 24 : 28;
    double scale[3];
    auto maximum = (std::uint64_t(1) << width) - 1;
    for (unsigned a = 0; a < 3; ++a)
      append_double(result, range.lo[a]);
    for (unsigned a = 0; a < 3; ++a) {
      scale[a] = (range.hi[a] - range.lo[a]) / double(maximum);
      append_double(result, scale[a]);
    }
    result.resize(48 + (std::size_t(vertices.size()) * 3 * width + 7) / 8);
    std::size_t cursor = 0;
    for (auto p : vertices)
      for (unsigned a = 0; a < 3; ++a) {
        auto q = scale[a] > 0 ? static_cast<std::uint64_t>(std::clamp(
                                    std::round((p[a] - range.lo[a]) / scale[a]),
                                    0., double(maximum)))
                              : 0;
        bit_write(result, 48, cursor, q, width);
      }
  } else if (encoding == MeshCodec::xor32) {
    auto blocks = (vertices.size() + 31) / 32;
    result.resize(4 * (blocks + 1));
    for (std::size_t block = 0; block < blocks; ++block) {
      put(result, 4 * block, result.size(), 4);
      auto first = 32 * block,
           n = std::min<std::size_t>(32, vertices.size() - first);
      std::array<std::uint64_t, 3> base{bits(vertices[first].x),
                                        bits(vertices[first].y),
                                        bits(vertices[first].z)};
      unsigned widths[3]{};
      for (std::size_t i = 1; i < n; ++i)
        for (unsigned a = 0; a < 3; ++a)
          widths[a] = std::max(
              widths[a], bit_width(base[a] ^ bits(vertices[first + i][a])));
      auto compressed =
          29 + ((n - 1) * (widths[0] + widths[1] + widths[2]) + 7) / 8;
      auto raw = 2 + 24 * n;
      bool pack = compressed < raw;
      result.push_back(pack ? 1 : 0);
      result.push_back(static_cast<unsigned char>(n));
      if (!pack) {
        for (std::size_t i = 0; i < n; ++i)
          for (unsigned a = 0; a < 3; ++a)
            append_double(result, vertices[first + i][a]);
      } else {
        for (auto value : base)
          append(result, value, 8);
        for (auto width : widths)
          result.push_back(static_cast<unsigned char>(width));
        auto start = result.size();
        result.resize(start + compressed - 29);
        std::size_t cursor = 0;
        for (std::size_t i = 1; i < n; ++i)
          for (unsigned a = 0; a < 3; ++a)
            bit_write(result, start, cursor,
                      base[a] ^ bits(vertices[first + i][a]), widths[a]);
      }
    }
    put(result, 4 * blocks, result.size(), 4);
  } else
    throw std::invalid_argument("Unknown vertex encoding");
  return result;
}
struct BuildNode {
  Box box;
  std::uint32_t first{}, second{};
};
Rational exact(double d);
double down(double x);
double up(double x);
std::array<std::uint16_t, 6> outward(Box box, Box scene) {
  std::array<std::uint16_t, 6> q{};
  for (unsigned a = 0; a < 3; ++a) {
    auto width = scene.hi[a] - scene.lo[a];
    auto lo = static_cast<long long>(
             std::floor((box.lo[a] - scene.lo[a]) / width * 65535)),
         hi = static_cast<long long>(
             std::ceil((box.hi[a] - scene.lo[a]) / width * 65535));
    lo = std::clamp<long long>(lo, 0, 65535);
    hi = std::clamp<long long>(hi, 0, 65535);
    while (lo > 0 &&
           exact(scene.lo[a]) + exact(width) * lo / 65535 > exact(box.lo[a]))
      --lo;
    while (hi < 65535 &&
           exact(scene.lo[a]) + exact(width) * hi / 65535 < exact(box.hi[a]))
      ++hi;
    q[a] = static_cast<std::uint16_t>(lo);
    q[a + 3] = static_cast<std::uint16_t>(hi);
  }
  return q;
}
Box node_box(const unsigned char *p, std::uint32_t index) {
  auto b = domain(p);
  auto n = p + section(p, 2) + 20 * index;
  for (unsigned a = 0; a < 3; ++a) {
    auto lo = b.lo[a], width = b.hi[a] - lo;
    // The declared domain is exactly [-1,1], so the product 2*q is exact.
    auto lower = down(width * integer(n + 2 * a, 2) / 65535.);
    auto upper = up(width * integer(n + 6 + 2 * a, 2) / 65535.);
    (&b.lo.x)[a] = down(lo + lower);
    (&b.hi.x)[a] = up(lo + upper);
  }
  return b;
}
Rational exact(double d) {
  if (!std::isfinite(d))
    throw std::invalid_argument("Nonfinite exact predicate input");
  if (d == 0)
    return Rational(0);
  auto value = bits(d);
  auto mantissa = value & ((std::uint64_t(1) << 52) - 1);
  int exponent = int((value >> 52) & 2047);
  if (exponent) {
    mantissa |= std::uint64_t(1) << 52;
    exponent -= 1075;
  } else
    exponent = -1074;
  Rational r = Integer(mantissa);
  if (exponent < 0)
    r /= Integer(1) << (-exponent);
  else
    r *= Integer(1) << exponent;
  if (value >> 63)
    r = -r;
  return r;
}
struct RVec {
  Rational x, y, z;
};
RVec rational(Vec3 p) { return {exact(p.x), exact(p.y), exact(p.z)}; }
// Binary64 successor/predecessor, exactly std::nextafter(x,+/-infinity).
// Bit stepping avoids a CRT call for each interval operation. Signed zero and
// infinities are explicit; memcpy preserves aliasing rules.
double down(double x) {
  if(std::isnan(x)||x==-std::numeric_limits<double>::infinity())return x;
  if(x==0)return -std::numeric_limits<double>::denorm_min();
  auto word=bits(x);word=x>0?word-1:word+1;std::memcpy(&x,&word,8);return x;
}
double up(double x) {
  if(std::isnan(x)||x==std::numeric_limits<double>::infinity())return x;
  if(x==0)return std::numeric_limits<double>::denorm_min();
  auto word=bits(x);word=x>0?word+1:word-1;std::memcpy(&x,&word,8);return x;
}
// Each elementary operation encloses its exact real result by directed one-step
// widening under IEEE round-to-nearest. No empirical epsilon decides a sign.
struct Interval {
  double lo, hi;
};
Interval operator+(Interval a, Interval b) {
  return {down(a.lo + b.lo), up(a.hi + b.hi)};
}
Interval operator-(Interval a, Interval b) {
  return {down(a.lo - b.hi), up(a.hi - b.lo)};
}
Interval operator*(Interval a, Interval b) {
  double p[4] = {a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
  return {down(*std::min_element(p, p + 4)), up(*std::max_element(p, p + 4))};
}
bool uncertain(Interval x) { return !(x.lo > 0 || x.hi < 0); }
struct IVec {
  Interval x, y, z;
};
IVec interval(Vec3 p) { return {{p.x, p.x}, {p.y, p.y}, {p.z, p.z}}; }
IVec operator-(IVec a, IVec b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
IVec cross_interval(IVec a, IVec b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
Interval dot_interval(IVec a, IVec b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
RVec subtract(const RVec &a, const RVec &b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
RVec cross_exact(const RVec &a, const RVec &b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
Rational dot_exact(const RVec &a, const RVec &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
bool boundary_exact(Vec3 p, Vec3 a, Vec3 b, Vec3 c) {
  auto ar = rational(a), br = rational(b), cr = rational(c), pr = rational(p);
  auto ab = subtract(br, ar), ac = subtract(cr, ar), ap = subtract(pr, ar),
       n = cross_exact(ab, ac);
  Rational height = dot_exact(ap, n);
  if (height != 0)
    return false;
  return dot_exact(cross_exact(ab, ap), n) >= 0 &&
         dot_exact(cross_exact(subtract(cr, br), subtract(pr, br)), n) >= 0 &&
         dot_exact(cross_exact(subtract(ar, cr), subtract(pr, cr)), n) >= 0;
}
double exact_squared(Vec3 p, Vec3 a, Vec3 b, Vec3 c,double &upper) {
  auto pr = rational(p), ar = rational(a), br = rational(b), cr = rational(c);
  auto ab = subtract(br, ar), ac = subtract(cr, ar), ap = subtract(pr, ar),
       n = cross_exact(ab, ac);
  Rational nn = dot_exact(n, n), height = dot_exact(ap, n);
  Rational best = -1;
  if (nn != 0) {
    RVec projection{pr.x - n.x * height / nn, pr.y - n.y * height / nn,
                    pr.z - n.z * height / nn};
    if (dot_exact(cross_exact(ab, subtract(projection, ar)), n) >= 0 &&
        dot_exact(cross_exact(subtract(cr, br), subtract(projection, br)), n) >=
            0 &&
        dot_exact(cross_exact(subtract(ar, cr), subtract(projection, cr)), n) >=
            0)
      best = height * height / nn;
  }
  for (auto edge :
       std::array<std::pair<RVec, RVec>, 3>{{{ar, br}, {br, cr}, {cr, ar}}}) {
    auto e = subtract(edge.second, edge.first), s = subtract(pr, edge.first);
    Rational length = dot_exact(e, e),
             t = length != 0 ? Rational(dot_exact(s, e) / length) : Rational(0);
    if (t < 0)
      t = 0;
    if (t > 1)
      t = 1;
    RVec delta{s.x - t * e.x, s.y - t * e.y, s.z - t * e.z};
    Rational d = dot_exact(delta, delta);
    if (best < 0 || d < best)
      best = d;
  }
  auto rounded=best.convert_to<double>();upper=up(rounded);
  while(exact(upper)<best)upper=up(upper);
  return rounded;
}
Interval square(Interval x) {
  auto lo=uncertain(x)?0.:std::min(std::abs(x.lo),std::abs(x.hi));
  auto hi=std::max(std::abs(x.lo),std::abs(x.hi));
  return {std::max(0.,down(lo*lo)),up(hi*hi)};
}
double squared_triangle(Vec3 p, Vec3 a, Vec3 b, Vec3 c,
                        std::uint64_t &fallbacks,double &upper) {
  auto ab = b - a, ac = c - a, n = cross(ab, ac);
  double nn = dot(n, n), height = dot(p - a, n);
  auto ai = interval(a), bi = interval(b), ci = interval(c), pi = interval(p);
  auto ni = cross_interval(bi - ai, ci - ai);
  auto hi = dot_interval(pi - ai, ni);
  auto scale = std::max({dot(ab, ab), dot(ac, ac), dot(c - b, c - b)});
  // Conditioning is a performance choice only: uncertain signs always use exact
  // arithmetic, and skinny triangles bypass the unstable projection formula.
  if (nn <= 1e-12 * scale * scale || uncertain(hi)) {
    ++fallbacks;
    return exact_squared(p, a, b, c,upper);
  }
  Interval side[3] = {dot_interval(cross_interval(bi - ai, pi - ai), ni),
                      dot_interval(cross_interval(ci - bi, pi - bi), ni),
                      dot_interval(cross_interval(ai - ci, pi - ci), ni)};
  for (auto value : side)
    if (uncertain(value)) {
      ++fallbacks;
      return exact_squared(p, a, b, c,upper);
    }
  double best = std::numeric_limits<double>::infinity();
  upper=std::numeric_limits<double>::infinity();
  if (side[0].lo > 0 && side[1].lo > 0 && side[2].lo > 0) {
    best = height * height / nn;
    auto normal_squared=square(ni.x)+square(ni.y)+square(ni.z);
    if(normal_squared.lo<=0){++fallbacks;return exact_squared(p,a,b,c,upper);}
    upper=up(square(hi).hi/normal_squared.lo);
    return best;
  }
  Vec3 chosen_a{},chosen_b{};double chosen_t=0;
  for (auto edge :
       std::array<std::pair<Vec3, Vec3>, 3>{{{a, b}, {b, c}, {c, a}}}) {
    auto e = edge.second - edge.first;
    auto t = std::clamp(dot(p - edge.first, e) / dot(e, e), 0., 1.);
    auto delta = p - edge.first - e * t;
    auto distance=dot(delta,delta);
    if(distance<best){best=distance;chosen_a=edge.first;chosen_b=edge.second;chosen_t=t;}
  }
  auto ei=interval(chosen_b)-interval(chosen_a);auto si=interval(p)-interval(chosen_a);
  Interval ti{chosen_t,chosen_t};IVec di{si.x-ei.x*ti,si.y-ei.y*ti,si.z-ei.z*ti};
  upper=(square(di.x)+square(di.y)+square(di.z)).hi;
  return best;
}
double box_squared(Box b, Vec3 p) {
  double result = 0;
  for (unsigned a = 0; a < 3; ++a) {
    auto d = std::max({down(b.lo[a] - p[a]), down(p[a] - b.hi[a]), 0.});
    result = std::max(0., down(result + std::max(0., down(d * d))));
  }
  return result;
}
bool ray_box(Box b, Vec3 origin, Vec3 direction, double low, double high) {
  for (unsigned a = 0; a < 3; ++a) {
    if (direction[a] == 0) {
      if (origin[a] < b.lo[a] || origin[a] > b.hi[a])
        return false;
    } else {
      double lo, hi;
      if (direction[a] > 0) {
        lo = down(down(b.lo[a] - origin[a]) / direction[a]);
        hi = up(up(b.hi[a] - origin[a]) / direction[a]);
      } else {
        lo = down(up(b.hi[a] - origin[a]) / direction[a]);
        hi = up(down(b.lo[a] - origin[a]) / direction[a]);
      }
      low = std::max(low, lo);
      high = std::min(high, hi);
    }
  }
  return high >= low;
}
struct ExactHit {
  bool hit{}, edge{}, coplanar{};
  Rational t{};
};
ExactHit triangle_hit_rational(const RVec &origin, const RVec &direction,
                               Vec3 a, Vec3 b, Vec3 c) {
  auto ar = rational(a), e1 = subtract(rational(b), ar),
       e2 = subtract(rational(c), ar), d = direction, s = subtract(origin, ar),
       h = cross_exact(d, e2);
  Rational det = dot_exact(e1, h);
  if (det == 0)
    return {false, false, dot_exact(s, cross_exact(e1, e2)) == 0, {}};
  Rational u = dot_exact(s, h) / det;
  auto q = cross_exact(s, e1);
  Rational v = dot_exact(d, q) / det;
  if (u < 0 || v < 0 || u + v > 1)
    return {};
  Rational t = dot_exact(e2, q) / det;
  return {true, u == 0 || v == 0 || u + v == 1, false, t};
}
ExactHit triangle_hit_exact(Vec3 origin, Vec3 direction, Vec3 a, Vec3 b,
                            Vec3 c) {
  return triangle_hit_rational(rational(origin), rational(direction), a, b, c);
}
// In a coplanar case clip the parameter interval against all three exact
// triangle halfplanes. An unrelated supporting plane cannot trigger a retry.
bool coplanar_overlap(const RVec &o, const RVec &d, Vec3 a, Vec3 b, Vec3 c,
                      Rational low, const Rational *high,
                      Rational *touch = nullptr) {
  auto ar = rational(a), br = rational(b), cr = rational(c);
  auto n = cross_exact(subtract(br, ar), subtract(cr, ar));
  Rational upper = high ? *high : Rational(0);
  bool bounded = high != nullptr;
  for (auto edge :
       std::array<std::pair<RVec, RVec>, 3>{{{ar, br}, {br, cr}, {cr, ar}}}) {
    auto e = subtract(edge.second, edge.first);
    Rational value = dot_exact(cross_exact(e, subtract(o, edge.first)), n);
    Rational slope = dot_exact(cross_exact(e, d), n);
    if (slope == 0) {
      if (value < 0)
        return false;
    } else {
      Rational t = -value / slope;
      if (slope > 0) {
        if (t > low)
          low = t;
      } else {
        if (!bounded || t < upper)
          upper = t;
        bounded = true;
      }
    }
    if (bounded && upper < low)
      return false;
  }
  if (touch && bounded && upper == low)
    *touch = low;
  return !bounded || upper >= low;
}
struct FilterHit {
  bool hit{}, degenerate{},boundary{};
};
FilterHit parity_hit(Vec3 origin, Vec3 direction, Vec3 a, Vec3 b, Vec3 c,
                     std::uint64_t &fallbacks) {
  auto e1 = b - a, e2 = c - a, s = origin - a, h = cross(direction, e2),
       q = cross(s, e1);
  auto det = dot(e1, h), u = dot(s, h), v = dot(direction, q), t = dot(e2, q),
       w = det - u - v;
  auto e1i = interval(b) - interval(a), e2i = interval(c) - interval(a),
       si = interval(origin) - interval(a);
  auto di = interval(direction), hi = cross_interval(di, e2i),
       qi = cross_interval(si, e1i);
  auto deti = dot_interval(e1i, hi), ui = dot_interval(si, hi),
       vi = dot_interval(di, qi), ti = dot_interval(e2i, qi);
  auto wi = deti - ui - vi;
  if (!uncertain(deti) && !uncertain(ui) && !uncertain(vi) && !uncertain(wi) &&
      !uncertain(ti)) {
    if (det < 0) {
      u = -u;
      v = -v;
      w = -w;
      t = -t;
    }
    return {u > 0 && v > 0 && w > 0 && t > 0, false};
  }
  ++fallbacks;
  auto hit = triangle_hit_exact(origin, direction, a, b, c);
  if (hit.coplanar)
    return {false, coplanar_overlap(rational(origin), rational(direction), a, b,
                                    c, 0, nullptr)};
  if(hit.hit && hit.t==0)return {false,false,true};
  if (!hit.hit || hit.t < 0)
    return {};
  return {!hit.edge, hit.edge};
}
constexpr std::array<Vec3, 12> directions{{{.738179, .521087, .429771},
                                           {.193337, .834731, .515219},
                                           {-.623911, .270113, .733879},
                                           {.449213, -.781137, .434129},
                                           {.873119, .337871, -.351731},
                                           {-.271913, -.576187, .770731},
                                           {.311729, .922177, .229731},
                                           {-.912131, .179119, -.368177},
                                           {.537113, -.281831, -.794719},
                                           {-.151373, .638291, -.754111},
                                           {.713179, -.649771, .262331},
                                           {-.491719, -.337117, -.803129}}};
} // namespace

TriangleGeometry read_triangle_obj(const std::string &path) {
  std::ifstream stream(path);
  if (!stream)
    throw std::invalid_argument("Mesh input unavailable");
  TriangleGeometry mesh;
  std::string line;
  while (std::getline(stream, line)) {
    std::istringstream row(line);
    std::string type;
    row >> type;
    if (type == "v") {
      Vec3 p;
      if (!(row >> p.x >> p.y >> p.z))
        throw std::invalid_argument("Invalid OBJ vertex");
      mesh.vertices.push_back(p);
    } else if (type == "f") {
      std::array<std::uint32_t, 3> face;
      std::string token;
      for (unsigned i = 0; i < 3; ++i) {
        if (!(row >> token))
          throw std::invalid_argument("Nontriangle OBJ");
        auto id = std::stoll(token.substr(0, token.find('/')));
        auto index =
            id > 0 ? id - 1 : static_cast<long long>(mesh.vertices.size()) + id;
        if (index < 0 || index >= static_cast<long long>(mesh.vertices.size()))
          throw std::invalid_argument("OBJ index out of range");
        face[i] = static_cast<std::uint32_t>(index);
      }
      if (row >> token)
        throw std::invalid_argument(
            "Nontriangle OBJ; no triangulation allowed");
      mesh.triangles.push_back(face);
    }
  }
  return mesh;
}

CompactMeshField::CompactMeshField(const TriangleGeometry &source,
                                  MeshCodec encoding,CompactBuildPhases *phases,std::vector<unsigned char> *serialized_evidence) {
  using Clock=std::chrono::steady_clock;auto phase=Clock::now();
  auto lap=[&](){auto end=Clock::now();auto seconds=std::chrono::duration<double>(end-phase).count();phase=end;return seconds;};
  if (source.vertices.empty() || source.triangles.empty() ||
      source.vertices.size() > std::numeric_limits<std::uint32_t>::max() ||
      source.triangles.size() > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument("Invalid mesh size");
  Box scene;
  for (auto p : source.vertices)
    if (!scene.contains(p) || !std::isfinite(p.x) || !std::isfinite(p.y) ||
        !std::isfinite(p.z))
      throw std::invalid_argument("Mesh outside finite domain");
  auto range = bounds(source.vertices);
  std::vector<std::uint32_t> order(source.vertices.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) {
    auto x = morton(source.vertices[a], range),
         y = morton(source.vertices[b], range);
    return x < y || (x == y && a < b);
  });
  std::vector<Vec3> vertices;
  std::vector<std::uint32_t> inverse(order.size());
  for (auto id : order) {
    inverse[id] = static_cast<std::uint32_t>(vertices.size());
    vertices.push_back(source.vertices[id]);
  }
  auto faces = source.triangles;
  for (auto &f : faces)
    for (auto &id : f) {
      if (id >= inverse.size())
        throw std::invalid_argument("Triangle index out of range");
      id = inverse[id];
    }
  Buffer data(header);
  data[0] = 'S';
  data[1] = 'D';
  data[2] = 'F';
  data[3] = 'M';
  data[4] = 1;
  data[6] = static_cast<unsigned char>(encoding);
  data[7] = vertices.size() <= 65536 ? 2 : 4;
  put(data, 8, vertices.size(), 4);
  put(data, 12, faces.size(), 4);
  put(data, 20, header, 4);
  for (unsigned a = 0; a < 3; ++a) {
    put(data, 32 + 8 * a, bits(scene.lo[a]), 8);
    put(data, 56 + 8 * a, bits(scene.hi[a]), 8);
  }
  auto reorder_time=lap();
  auto vertex_data = vertex_buffer(vertices, encoding);
  data.insert(data.end(), vertex_data.begin(), vertex_data.end());
  for (std::uint32_t i = 0; i < vertices.size(); ++i)
    vertices[i] = decode_vertex(data.data(), i);
  auto encoding_time=lap();
  std::vector<std::uint32_t> triangle_order(faces.size());
  std::iota(triangle_order.begin(), triangle_order.end(), 0);
  std::vector<BuildNode> nodes;
  std::function<void(std::size_t, std::size_t)> build = [&](std::size_t first,
                                                            std::size_t last) {
    auto index = nodes.size();
    nodes.emplace_back();
    Box box{vertices[faces[triangle_order[first]][0]],
            vertices[faces[triangle_order[first]][0]]};
    std::vector<Vec3> centroids;
    centroids.reserve(last - first);
    for (auto j = first; j < last; ++j) {
      auto f = faces[triangle_order[j]];
      for (auto id : f)
        for (unsigned a = 0; a < 3; ++a) {
          (&box.lo.x)[a] = std::min(box.lo[a], vertices[id][a]);
          (&box.hi.x)[a] = std::max(box.hi[a], vertices[id][a]);
        }
      centroids.push_back((vertices[f[0]] + vertices[f[1]] + vertices[f[2]]) *
                          (1. / 3));
    }
    nodes[index].box = box;
    if (last - first <= 16) {
      nodes[index].first = static_cast<std::uint32_t>(first);
      nodes[index].second = static_cast<std::uint32_t>(last - first);
      return;
    }
    auto extent = bounds(centroids);
    unsigned axis = 0;
    for (unsigned a = 1; a < 3; ++a)
      if (extent.hi[a] - extent.lo[a] > extent.hi[axis] - extent.lo[axis])
        axis = a;
    std::stable_sort(triangle_order.begin() + first,
                     triangle_order.begin() + last, [&](auto a, auto b) {
                       auto fa = faces[a], fb = faces[b];
                       auto ca =
                                (vertices[fa[0]][axis] + vertices[fa[1]][axis] +
                                 vertices[fa[2]][axis]) /
                                3,
                            cb =
                                (vertices[fb[0]][axis] + vertices[fb[1]][axis] +
                                 vertices[fb[2]][axis]) /
                                3;
                       return ca < cb || (ca == cb && a < b);
                     });
    auto middle = (first + last) / 2;
    build(first, middle);
    nodes[index].first = static_cast<std::uint32_t>(nodes.size());
    build(middle, last);
  };
  build(0, faces.size());
  auto bvh_time=lap();
  align(data, 4);
  put(data, 24, data.size(), 4);
  for (auto id : triangle_order)
    for (auto vertex : faces[id])
      append(data, vertex, data[7]);
  align(data, 4);
  put(data, 28, data.size(), 4);
  put(data, 16, nodes.size(), 4);
  for (auto node : nodes) {
    auto q = outward(node.box, scene);
    for (auto value : q)
      append(data, value, 2);
    append(data, node.first, 4);
    append(data, node.second, 4);
  }
  if (data.size() > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument("Serialized mesh exceeds offset range");
  size_ = data.size();
  data_ = std::make_unique<unsigned char[]>(size_);
  std::copy(data.begin(), data.end(), data_.get());
  if(serialized_evidence)*serialized_evidence=data;
  auto packing_time=lap();
  validate_layout();
  auto validation_time=lap();
  if(phases)*phases={reorder_time,encoding_time,bvh_time,packing_time,validation_time};
}
CompactMeshField::CompactMeshField(const Buffer &serialized)
    : data_(std::make_unique<unsigned char[]>(serialized.size())),
      size_(serialized.size()) {
  std::copy(serialized.begin(), serialized.end(), data_.get());
  validate_layout();
}
std::uint32_t CompactMeshField::vertex_count() const {
  return count(data_.get(), 0);
}
std::uint32_t CompactMeshField::triangle_count() const {
  return count(data_.get(), 1);
}
std::uint32_t CompactMeshField::node_count() const {
  return count(data_.get(), 2);
}
MeshCodec CompactMeshField::codec() const { return MeshCodec(data_[6]); }
Vec3 CompactMeshField::vertex(std::uint32_t index) const {
  return decode_vertex(data_.get(), index);
}
std::array<std::uint32_t, 3>
CompactMeshField::triangle(std::uint32_t index) const {
  if (index >= triangle_count())
    throw std::invalid_argument("Triangle index out of range");
  auto p =
      data_.get() + section(data_.get(), 1) + std::size_t(index) * 3 * data_[7];
  return {static_cast<std::uint32_t>(integer(p, data_[7])),
          static_cast<std::uint32_t>(integer(p + data_[7], data_[7])),
          static_cast<std::uint32_t>(integer(p + 2 * data_[7], data_[7]))};
}
TriangleGeometry CompactMeshField::decoded_geometry() const {
  TriangleGeometry result;
  for (std::uint32_t i = 0; i < vertex_count(); ++i)
    result.vertices.push_back(vertex(i));
  for (std::uint32_t i = 0; i < triangle_count(); ++i)
    result.triangles.push_back(triangle(i));
  return result;
}
void CompactMeshField::validate_layout() const {
  if (size_ < header || std::memcmp(data_.get(), "SDFM", 4) || data_[4] != 1 ||
      data_[5] != 0 || data_[6] < 1 || data_[6] > 5 ||
      (data_[7] != 2 && data_[7] != 4))
    throw std::invalid_argument("Invalid compact header");
  auto v = vertex_count(), f = triangle_count(), n = node_count();
  if (!v || !f || !n || (data_[7] == 2) != (v <= 65536))
    throw std::invalid_argument("Invalid compact counts/index width");
  auto vp = section(data_.get(), 0), tp = section(data_.get(), 1),
       bp = section(data_.get(), 2);
  if (vp != header || tp < vp || bp < tp || tp % 4 || bp % 4 ||
      std::uint64_t(bp) + 20ull * n != size_ ||
      std::uint64_t(tp) + std::uint64_t(f) * 3 * data_[7] > bp)
    throw std::invalid_argument("Invalid section layout/truncated buffer");
  auto scene = domain(data_.get());
  scene.validate();
  for (unsigned a = 0; a < 3; ++a)
    if (scene.lo[a] != -1 || scene.hi[a] != 1)
      throw std::invalid_argument("Unexpected compact domain");
  std::size_t end = vp;
  if (codec() == MeshCodec::f64 || codec() == MeshCodec::f32)
    end += std::size_t(v) * (codec() == MeshCodec::f64 ? 24 : 12);
  else if (codec() == MeshCodec::q24 || codec() == MeshCodec::q28) {
    unsigned width = codec() == MeshCodec::q24 ? 24 : 28;
    auto bitcount = std::size_t(v) * 3 * width;
    end += 48 + (bitcount + 7) / 8;
    if (end > tp)
      throw std::invalid_argument("Truncated quantized vertices");
    for (unsigned a = 0; a < 3; ++a)
      if (!std::isfinite(floating(data_.get() + vp + 8 * a)) ||
          !std::isfinite(floating(data_.get() + vp + 24 + 8 * a)) ||
          floating(data_.get() + vp + 24 + 8 * a) < 0)
        throw std::invalid_argument("Invalid quantizer parameters");
    if (bitcount % 8 && (data_[end - 1] >> (bitcount % 8)))
      throw std::invalid_argument("Nonzero quantizer tail bits");
  } else {
    auto blocks = (std::size_t(v) + 31) / 32, dir = 4 * (blocks + 1);
    if (vp + dir > tp || integer(data_.get() + vp, 4) != dir)
      throw std::invalid_argument("Invalid XOR directory");
    for (std::size_t block = 0; block < blocks; ++block) {
      auto begin = integer(data_.get() + vp + 4 * block, 4),
           finish = integer(data_.get() + vp + 4 * (block + 1), 4);
      if (finish <= begin || vp + finish > tp || finish - begin < 2)
        throw std::invalid_argument("Invalid XOR block extent");
      auto p = data_.get() + vp + begin;
      auto expected = std::min<std::size_t>(32, v - 32 * block);
      if (p[1] != expected || p[0] > 1)
        throw std::invalid_argument("Invalid XOR mode/count");
      std::size_t length = 2 + 24 * expected;
      if (p[0]) {
        if (finish - begin < 29 || p[26] > 64 || p[27] > 64 || p[28] > 64)
          throw std::invalid_argument("Invalid XOR residual widths");
        auto bits = (expected - 1) * (p[26] + p[27] + p[28]);
        length = 29 + (bits + 7) / 8;
        if (finish - begin != length)
          throw std::invalid_argument("XOR block length mismatch");
        if (bits % 8 && (p[length - 1] >> (bits % 8)))
          throw std::invalid_argument("Nonzero XOR tail bits");
      }
      if (finish - begin != length)
        throw std::invalid_argument("XOR block length mismatch");
    }
    end = vp + integer(data_.get() + vp + 4 * blocks, 4);
  }
  if (end > tp || tp - end > 3)
    throw std::invalid_argument("Vertex section padding invalid");
  for (auto i = end; i < tp; ++i)
    if (data_[i])
      throw std::invalid_argument("Nonzero vertex padding");
  end = tp + std::size_t(f) * 3 * data_[7];
  if (bp - end > 3)
    throw std::invalid_argument("Triangle section padding invalid");
  for (auto i = end; i < bp; ++i)
    if (data_[i])
      throw std::invalid_argument("Nonzero triangle padding");
  for (std::uint32_t i = 0; i < v; ++i) {
    auto p = vertex(i);
    if (!scene.contains(p) || !std::isfinite(p.x) || !std::isfinite(p.y) ||
        !std::isfinite(p.z))
      throw std::invalid_argument("Invalid decoded vertex");
  }
  for (std::uint32_t i = 0; i < f; ++i) {
    auto indices = triangle(i);
    if (indices[0] == indices[1] || indices[1] == indices[2] ||
        indices[0] == indices[2])
      throw std::invalid_argument("Repeated triangle index");
    auto a = rational(vertex(indices[0])), b = rational(vertex(indices[1])),
         c = rational(vertex(indices[2]));
    auto normal = cross_exact(subtract(b, a), subtract(c, a));
    if (dot_exact(normal, normal) == 0)
      throw std::invalid_argument("Degenerate decoded triangle");
  }
  std::vector<unsigned char> seen(n), covered(f);
  std::function<void(std::uint32_t)> visit = [&](std::uint32_t index) {
    if (index >= n || seen[index])
      throw std::invalid_argument("BVH cycle/shared node");
    seen[index] = 1;
    auto p = data_.get() + bp + 20 * index;
    for (unsigned a = 0; a < 3; ++a)
      if (integer(p + 2 * a, 2) > integer(p + 6 + 2 * a, 2))
        throw std::invalid_argument("Inverted BVH bounds");
    auto first = static_cast<std::uint32_t>(integer(p + 12, 4)),
         second = static_cast<std::uint32_t>(integer(p + 16, 4));
    if (second) {
      if (second > 16 || std::uint64_t(first) + second > f)
        throw std::invalid_argument("Invalid BVH leaf");
      for (auto i = first; i < first + second; ++i) {
        if (covered[i]++)
          throw std::invalid_argument("Duplicate BVH triangle");
        for (auto id : triangle(i)) {
          auto v = vertex(id);
          for (unsigned a = 0; a < 3; ++a) {
            Rational lo =
                exact(scene.lo[a]) + (exact(scene.hi[a]) - exact(scene.lo[a])) *
                                         integer(p + 2 * a, 2) / 65535;
            Rational hi =
                exact(scene.lo[a]) + (exact(scene.hi[a]) - exact(scene.lo[a])) *
                                         integer(p + 6 + 2 * a, 2) / 65535;
            if (exact(v[a]) < lo || exact(v[a]) > hi)
              throw std::invalid_argument("Nonconservative BVH bounds");
          }
        }
      }
    } else {
      if (first <= index + 1 || first >= n)
        throw std::invalid_argument("Invalid BVH right child");
      for (auto child : {index + 1, first}) {
        auto b = data_.get() + bp + 20 * child;
        for (unsigned a = 0; a < 3; ++a)
          if (integer(b + 2 * a, 2) < integer(p + 2 * a, 2) ||
              integer(b + 6 + 2 * a, 2) > integer(p + 6 + 2 * a, 2))
            throw std::invalid_argument(
                "BVH parent bounds do not contain child");
        visit(child);
      }
    }
  };
  visit(0);
  if (std::find(seen.begin(), seen.end(), 0) != seen.end() ||
      std::find(covered.begin(), covered.end(), 0) != covered.end())
    throw std::invalid_argument("Uncovered BVH node/triangle");
}
MeshQuery CompactMeshField::evaluate(Vec3 point) const {
  if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
      !std::isfinite(point.z))
    throw std::invalid_argument("Nonfinite query");
  ++statistics_.queries;
  double best = std::numeric_limits<double>::infinity();
  double prune_upper = best;
  std::vector<std::uint32_t> pending{0};
  while (!pending.empty()) {
    auto index = pending.back();
    pending.pop_back();
    ++statistics_.nearest_nodes;
    // A decoded vertex is a genuine mesh point. This upper bound is independent
    // of the rounded triangle-distance result and cannot prune the nearest
    // face.
    if (box_squared(node_box(data_.get(), index), point) > prune_upper)
      continue;
    auto p = data_.get() + section(data_.get(), 2) + 20 * index;
    auto first = static_cast<std::uint32_t>(integer(p + 12, 4)),
         second = static_cast<std::uint32_t>(integer(p + 16, 4));
    if (!second) {
      auto left = index + 1, right = first;
      auto a = box_squared(node_box(data_.get(), left), point),
           b = box_squared(node_box(data_.get(), right), point);
      if (a < b) {
        pending.push_back(right);
        pending.push_back(left);
      } else {
        pending.push_back(left);
        pending.push_back(right);
      }
    } else
      for (auto i = first; i < first + second; ++i) {
        auto indices = triangle(i);
        auto a = vertex(indices[0]), b = vertex(indices[1]),
             c = vertex(indices[2]);
        ++statistics_.distance_triangles;
        double triangle_upper;
        auto d = squared_triangle(point, a, b, c, statistics_.exact_fallbacks,triangle_upper);
        prune_upper=std::min(prune_upper,triangle_upper);
        if (d == 0) {
          ++statistics_.exact_fallbacks;
          if (boundary_exact(point, a, b, c))
            return {0, true, true};
        }
        best = std::min(best, d);
      }
  }
  auto side=parity(point);if(!side.resolved||side.boundary)return side;
  return {std::sqrt(best)*side.distance,true,false};
}
MeshQuery CompactMeshField::classify(Vec3 point) const {
  if(!std::isfinite(point.x)||!std::isfinite(point.y)||!std::isfinite(point.z))throw std::invalid_argument("Nonfinite classification");
  ++statistics_.classification_queries;return parity(point);
}
MeshQuery CompactMeshField::parity(Vec3 point) const {
  std::vector<std::uint32_t> pending;
  for (std::size_t attempt = 0; attempt < directions.size(); ++attempt) {
    if (attempt)
      ++statistics_.direction_retries;
    std::size_t intersections = 0;
    bool degeneracy = false;
    pending = {0};
    while (!pending.empty() && !degeneracy) {
      auto index = pending.back();
      pending.pop_back();
      ++statistics_.parity_nodes;
      if (!ray_box(node_box(data_.get(), index), point, directions[attempt], 0,
                   std::numeric_limits<double>::infinity()))
        continue;
      auto p = data_.get() + section(data_.get(), 2) + 20 * index;
      auto first = static_cast<std::uint32_t>(integer(p + 12, 4)),
           second = static_cast<std::uint32_t>(integer(p + 16, 4));
      if (!second) {
        pending.push_back(first);
        pending.push_back(index + 1);
      } else
        for (auto i = first; i < first + second; ++i) {
          auto indices = triangle(i);
          ++statistics_.parity_triangles;
          auto hit = parity_hit(point, directions[attempt], vertex(indices[0]),
                                vertex(indices[1]), vertex(indices[2]),
                                statistics_.exact_fallbacks);
          if(hit.boundary)return {0,true,true};
          if (hit.degenerate) {
            degeneracy = true;
            break;
          }
          intersections += hit.hit;
        }
    }
    if (!degeneracy)
      return {intersections % 2 ? -1. : 1., true, false};
  }
  ++statistics_.unresolved;
  return {std::numeric_limits<double>::quiet_NaN(), false, false};
}
double CompactMeshField::query(Vec3 point) const {
  auto result = evaluate(point);
  if (!result.resolved)
    throw std::runtime_error(
        "Compact parity direction retries exhausted: numerical_unresolved");
  return result.distance;
}
MeshRayEvents CompactMeshField::ray(Vec3 origin, Vec3 direction, double low,
                                    double high) const {
  if (!(high > low) || !std::isfinite(low) || !std::isfinite(high) ||
      !std::isfinite(origin.x) || !std::isfinite(origin.y) ||
      !std::isfinite(origin.z) || !std::isfinite(direction.x) ||
      !std::isfinite(direction.y) || !std::isfinite(direction.z) ||
      !std::isfinite(norm(direction)) || norm(direction) == 0)
    throw std::invalid_argument("Invalid mesh measurement ray");
  MeshRayEvents result;
  std::vector<Rational> events;
  auto lower = exact(low), upper = exact(high);
  std::vector<std::uint32_t> candidates, pending{0};
  while (!pending.empty()) {
    auto node = pending.back();
    pending.pop_back();
    if (!ray_box(node_box(data_.get(), node), origin, direction, low, high))
      continue;
    auto p = data_.get() + section(data_.get(), 2) + 20 * node;
    auto first = std::uint32_t(integer(p + 12, 4)),
         count = std::uint32_t(integer(p + 16, 4));
    if (count)
      for (auto i = first; i < first + count; ++i)
        candidates.push_back(i);
    else {
      pending.push_back(first);
      pending.push_back(node + 1);
    }
  }
  for (auto i : candidates) {
    auto indices = triangle(i);
    ++result.triangle_tests;
    auto hit = triangle_hit_exact(origin, direction, vertex(indices[0]),
                                  vertex(indices[1]), vertex(indices[2]));
    if (hit.coplanar) {
      Rational touch = lower - 1;
      if (coplanar_overlap(rational(origin), rational(direction),
                           vertex(indices[0]), vertex(indices[1]),
                           vertex(indices[2]), lower, &upper, &touch)) {
        result.coplanar = true;
        if (touch >= lower && touch <= upper) {
          if (touch > lower && touch < upper)
            events.push_back(touch);
        } else
          result.unresolved = true;
      }
      continue;
    }
    if (hit.hit && hit.t > lower && hit.t < upper)
      events.push_back(hit.t);
  }
  if (result.unresolved)
    return result;
  std::sort(events.begin(), events.end());
  auto end = std::unique(events.begin(), events.end());
  result.exact_shared_events = events.end() - end;
  events.erase(end, events.end());
  std::vector<double> times{low};
  for (const auto &t : events) {
    auto value = t.convert_to<double>();
    if (value <= times.back() || value >= high) {
      result.unresolved = true;
      return result;
    }
    times.push_back(value);
  }
  times.push_back(high);
  std::vector<int> signs;
  auto ro = rational(origin), rd = rational(direction);
  std::vector<Rational> cuts{lower};
  cuts.insert(cuts.end(), events.begin(), events.end());
  cuts.push_back(upper);
  for (std::size_t i = 1; i < cuts.size(); ++i) {
    Rational mid = (cuts[i - 1] + cuts[i]) / 2;
    RVec point{ro.x + rd.x * mid, ro.y + rd.y * mid, ro.z + rd.z * mid};
    Vec3 rounded{point.x.convert_to<double>(), point.y.convert_to<double>(),
                 point.z.convert_to<double>()};
    auto rp = rational(rounded);
    Vec3 error{};
    Rational delta[3] = {point.x - rp.x, point.y - rp.y, point.z - rp.z};
    for (unsigned a = 0; a < 3; ++a) {
      if (delta[a] < 0)
        delta[a] = -delta[a];
      (&error.x)[a] = up(delta[a].convert_to<double>());
      while (exact(error[a]) < delta[a])
        (&error.x)[a] = up(error[a]);
    }
    bool resolved = false;
    for (auto parity_direction : directions) {
      std::size_t crossings = 0;
      bool degenerate = false;
      auto d = rational(parity_direction);
      pending = {0};
      candidates.clear();
      while (!pending.empty()) {
        auto node = pending.back();
        pending.pop_back();
        auto box = node_box(data_.get(), node);
        for (unsigned a = 0; a < 3; ++a) {
          (&box.lo.x)[a] = down(box.lo[a] - error[a]);
          (&box.hi.x)[a] = up(box.hi[a] + error[a]);
        }
        if (!ray_box(box, rounded, parity_direction, 0,
                     std::numeric_limits<double>::infinity()))
          continue;
        auto p = data_.get() + section(data_.get(), 2) + 20 * node;
        auto first = std::uint32_t(integer(p + 12, 4)),
             count = std::uint32_t(integer(p + 16, 4));
        if (count)
          for (auto j = first; j < first + count; ++j)
            candidates.push_back(j);
        else {
          pending.push_back(first);
          pending.push_back(node + 1);
        }
      }
      for (auto j : candidates) {
        auto ids = triangle(j);
        auto a = vertex(ids[0]), b = vertex(ids[1]), c = vertex(ids[2]);
        auto hit = triangle_hit_rational(point, d, a, b, c);
        if (hit.coplanar) {
          if (coplanar_overlap(point, d, a, b, c, 0, nullptr)) {
            degenerate = true;
            break;
          }
        } else if (hit.hit && hit.t == 0) {
          result.unresolved = true;
          return result;
        } else if (hit.hit && hit.t > 0) {
          if (hit.edge) {
            degenerate = true;
            break;
          }
          ++crossings;
        }
      }
      if (!degenerate) {
        signs.push_back(crossings % 2 ? -1 : 1);
        resolved = true;
        break;
      }
    }
    if (!resolved) {
      result.unresolved = true;
      return result;
    }
  }
  for (std::size_t i = 1; i < times.size() - 1; ++i) {
    if (signs[i - 1] != signs[i])
      result.transitions.push_back(times[i]);
    else
      result.tangencies.push_back(times[i]);
  }
  return result;
}
} // namespace sdf
