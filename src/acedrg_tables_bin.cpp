// Copyright 2026 Global Phasing Ltd.
//
// Binary cache of AcedrgTables — the heavy bond/angle index structures.
// See include/gemmi/acedrg_tables.hpp::save_binary / load_binary.
//
// On-disk format:
//
//   magic      "GMMIACDB" (8 bytes)
//   version    uint32_t = 1
//   reserved   uint32_t = 0
//   then a sequence of typed records, recursively serialised by
//   write_value() / read_value() below. The order matches save_binary().
//
// The format is intentionally simple — no string interning, no global
// index. Reading is one big fread() followed by linear pointer walks.

#include "gemmi/acedrg_tables.hpp"
#include "gemmi/fail.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace gemmi {

namespace {

constexpr char kMagic[8]    = {'G','M','M','I','A','C','D','B'};
constexpr uint32_t kVersion = 1;

// ─── BinaryWriter ────────────────────────────────────────────────────────────

struct BW {
  FILE* f;
  void w(const void* p, size_t n) {
    if (std::fwrite(p, 1, n, f) != n)
      throw std::runtime_error("acedrg binary: short write");
  }
  void u32(uint32_t v) { w(&v, 4); }
  void i32(int32_t v)  { w(&v, 4); }
  void f64(double v)   { w(&v, 8); }
};

// ─── BinaryReader (buffer-walking; loads the whole file with one fread) ─────

struct BR {
  const char* p;
  const char* end;
  void check(size_t n) {
    if (size_t(end - p) < n)
      throw std::runtime_error("acedrg binary: short read");
  }
  void r(void* dst, size_t n) { check(n); std::memcpy(dst, p, n); p += n; }
  uint32_t u32() { uint32_t v; r(&v, 4); return v; }
  int32_t  i32() { int32_t v;  r(&v, 4); return v; }
  double   f64() { double v;   r(&v, 8); return v; }
};

// ─── write_value overloads ──────────────────────────────────────────────────

inline void write_value(BW& bw, int v) { bw.i32(v); }
inline void write_value(BW& bw, double v) { bw.f64(v); }

inline void write_value(BW& bw, const std::string& s) {
  bw.u32(static_cast<uint32_t>(s.size()));
  if (!s.empty()) bw.w(s.data(), s.size());
}

inline void write_value(BW& bw, const CodStats& v) {
  bw.f64(v.value); bw.f64(v.sigma); bw.i32(v.count); bw.i32(v.level);
}

inline void write_value(BW& bw, const AcedrgTables::Ccp4BondEntry& v) {
  bw.f64(v.length); bw.f64(v.sigma);
}

template<typename T>
inline void write_value(BW& bw, const std::vector<T>& v) {
  bw.u32(static_cast<uint32_t>(v.size()));
  for (const auto& x : v) write_value(bw, x);
}

template<typename K, typename V>
inline void write_value(BW& bw, const std::map<K, V>& m) {
  bw.u32(static_cast<uint32_t>(m.size()));
  for (const auto& kv : m) {
    write_value(bw, kv.first);
    write_value(bw, kv.second);
  }
}

template<typename K, typename V>
inline void write_value(BW& bw, const std::unordered_map<K, V>& m) {
  bw.u32(static_cast<uint32_t>(m.size()));
  for (const auto& kv : m) {
    write_value(bw, kv.first);
    write_value(bw, kv.second);
  }
}

template<typename T>
inline void write_value(BW& bw, const std::unordered_set<T>& s) {
  bw.u32(static_cast<uint32_t>(s.size()));
  for (const auto& x : s) write_value(bw, x);
}

// ─── read_value overloads ───────────────────────────────────────────────────

inline void read_value(BR& br, int& v) { v = br.i32(); }
inline void read_value(BR& br, double& v) { v = br.f64(); }

inline void read_value(BR& br, std::string& s) {
  uint32_t n = br.u32();
  s.resize(n);
  if (n) br.r(&s[0], n);
}

inline void read_value(BR& br, CodStats& v) {
  v.value = br.f64();
  v.sigma = br.f64();
  v.count = br.i32();
  v.level = br.i32();
}

inline void read_value(BR& br, AcedrgTables::Ccp4BondEntry& v) {
  v.length = br.f64();
  v.sigma  = br.f64();
}

template<typename T>
inline void read_value(BR& br, std::vector<T>& v) {
  uint32_t n = br.u32();
  v.resize(n);
  for (uint32_t i = 0; i < n; ++i)
    read_value(br, v[i]);
}

// std::map is serialised in iteration order (sorted ascending), so on
// read we can hint-insert at end() — amortized O(1) per insertion instead
// of O(log n).
template<typename K, typename V>
inline void read_value(BR& br, std::map<K, V>& m) {
  m.clear();
  uint32_t n = br.u32();
  for (uint32_t i = 0; i < n; ++i) {
    K k; read_value(br, k);
    V vv; read_value(br, vv);
    m.emplace_hint(m.end(), std::move(k), std::move(vv));
  }
}

template<typename K, typename V>
inline void read_value(BR& br, std::unordered_map<K, V>& m) {
  m.clear();
  uint32_t n = br.u32();
  m.reserve(n);
  for (uint32_t i = 0; i < n; ++i) {
    K k; read_value(br, k);
    V vv; read_value(br, vv);
    m.emplace(std::move(k), std::move(vv));
  }
}

template<typename T>
inline void read_value(BR& br, std::unordered_set<T>& s) {
  s.clear();
  uint32_t n = br.u32();
  s.reserve(n);
  for (uint32_t i = 0; i < n; ++i) {
    T x; read_value(br, x);
    s.insert(std::move(x));
  }
}

}  // anonymous namespace

// ─── save_binary ────────────────────────────────────────────────────────────

void AcedrgTables::save_binary(const std::string& path) const {
  FILE* f = std::fopen(path.c_str(), "wb");
  if (!f)
    throw std::runtime_error("acedrg binary: cannot open for write: " + path);
  BW bw{f};
  try {
    bw.w(kMagic, sizeof(kMagic));
    bw.u32(kVersion);
    bw.u32(0);  // reserved

    // Heavy bond indices.
    write_value(bw, bond_idx_1d_);
    write_value(bw, bond_idx_full_);
    write_value(bw, bond_idx_2d_);
    write_value(bw, bond_hasp_2d_);
    write_value(bw, bond_hasp_1d_);
    write_value(bw, bond_hasp_0d_);
    write_value(bw, bond_2d_hybr_keys_);
    write_value(bw, bond_full_4prefix_keys_);

    // Heavy angle indices.
    write_value(bw, angle_idx_1d_);
    write_value(bw, angle_idx_2d_);
    write_value(bw, angle_idx_3d_);
    write_value(bw, angle_idx_4d_);
    write_value(bw, angle_idx_5d_);
    write_value(bw, angle_idx_6d_);
  } catch (...) {
    std::fclose(f);
    throw;
  }
  std::fclose(f);
}

// ─── load_binary ────────────────────────────────────────────────────────────

void AcedrgTables::load_binary(const std::string& path) {
  FILE* f = std::fopen(path.c_str(), "rb");
  if (!f)
    throw std::runtime_error("acedrg binary: cannot open for read: " + path);

  // Read the whole file into one buffer; the format is dense enough that
  // sequential pointer-walking beats per-record fread overhead.
  std::fseek(f, 0, SEEK_END);
  size_t size = std::ftell(f);
  std::fseek(f, 0, SEEK_SET);
  std::vector<char> buf(size);
  if (size && std::fread(buf.data(), 1, size, f) != size) {
    std::fclose(f);
    throw std::runtime_error("acedrg binary: short read of " + path);
  }
  std::fclose(f);

  BR br{buf.data(), buf.data() + buf.size()};

  char magic[8];
  br.r(magic, sizeof(magic));
  if (std::memcmp(magic, kMagic, sizeof(magic)) != 0)
    throw std::runtime_error("acedrg binary: bad magic in " + path);
  uint32_t ver = br.u32();
  if (ver != kVersion)
    throw std::runtime_error("acedrg binary: version mismatch in " + path +
                             " (file=" + std::to_string(ver) +
                             ", code=" + std::to_string(kVersion) + ")");
  br.u32();  // reserved

  using Clock = std::chrono::steady_clock;
  auto t0 = Clock::now();
  auto step = [&](const char* label, size_t entries) {
    if (verbose >= 2) {
      auto t1 = Clock::now();
      double dt = std::chrono::duration<double>(t1 - t0).count();
      std::fprintf(stderr,
                   "    [bin] %-22s %6.3f s  (%zu top-level entries)\n",
                   label, dt, entries);
      t0 = t1;
    }
  };

  read_value(br, bond_idx_1d_);              step("bond_idx_1d_",     bond_idx_1d_.size());
  read_value(br, bond_idx_full_);            step("bond_idx_full_",   bond_idx_full_.size());
  read_value(br, bond_idx_2d_);              step("bond_idx_2d_",     bond_idx_2d_.size());
  read_value(br, bond_hasp_2d_);             step("bond_hasp_2d_",    bond_hasp_2d_.size());
  read_value(br, bond_hasp_1d_);             step("bond_hasp_1d_",    bond_hasp_1d_.size());
  read_value(br, bond_hasp_0d_);             step("bond_hasp_0d_",    bond_hasp_0d_.size());
  read_value(br, bond_2d_hybr_keys_);        step("bond_2d_hybr_keys",bond_2d_hybr_keys_.size());
  read_value(br, bond_full_4prefix_keys_);   step("bond_4prefix_keys",bond_full_4prefix_keys_.size());

  read_value(br, angle_idx_1d_);             step("angle_idx_1d_",    angle_idx_1d_.size());
  read_value(br, angle_idx_2d_);             step("angle_idx_2d_",    angle_idx_2d_.size());
  read_value(br, angle_idx_3d_);             step("angle_idx_3d_",    angle_idx_3d_.size());
  read_value(br, angle_idx_4d_);             step("angle_idx_4d_",    angle_idx_4d_.size());
  read_value(br, angle_idx_5d_);             step("angle_idx_5d_",    angle_idx_5d_.size());
  read_value(br, angle_idx_6d_);             step("angle_idx_6d_",    angle_idx_6d_.size());
}

}  // namespace gemmi
