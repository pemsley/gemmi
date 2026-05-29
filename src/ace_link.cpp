// Copyright 2026 Global Phasing Ltd.
//
// Implementation of make_joined_chemcomp and prepare_chemlink — see
// include/gemmi/ace_link.hpp for the design rationale.

#include "gemmi/ace_link.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

#include "gemmi/ace_cc.hpp"          // prepare_chemcomp, PrepareChemcompOptions
#include "gemmi/ace_graph.hpp"       // expected_valence_for_nonmetal, build_bond_adjacency
#include "gemmi/acedrg_tables.hpp"   // AcedrgTables
#include "gemmi/util.hpp"            // vector_remove_if

namespace gemmi {

// ─── atom-id prefix helpers ──────────────────────────────────────────────────

std::string make_joined_atom_id(int side, const std::string& bare) {
  return std::to_string(side) + "/" + bare;
}

std::pair<int, std::string> split_joined_atom_id(const std::string& joined_id) {
  if (joined_id.size() >= 3 && joined_id[1] == '/' &&
      (joined_id[0] == '1' || joined_id[0] == '2'))
    return {joined_id[0] - '0', joined_id.substr(2)};
  return {0, joined_id};
}

// ─── joined builder helpers (file-local) ────────────────────────────────────

namespace {

// Bond-order multiplier matching what ace_graph.cpp uses internally.
float order_of(BondType t) {
  switch (t) {
    case BondType::Single:   return 1.0f;
    case BondType::Double:   return 2.0f;
    case BondType::Triple:   return 3.0f;
    case BondType::Aromatic: return 1.5f;
    case BondType::Deloc:    return 1.5f;
    case BondType::Metal:    return 0.0f;
    case BondType::Unspec:   return 1.0f;
  }
  return 1.0f;
}

// Copy a ChemComp::Atom into the joined CC, prefixing the id.
ChemComp::Atom prefix_atom(const ChemComp::Atom& src, int side) {
  ChemComp::Atom out = src;
  out.id = make_joined_atom_id(side, src.id);
  if (!src.old_id.empty())
    out.old_id = make_joined_atom_id(side, src.old_id);
  return out;
}

// Append every atom of `src` (prefixed) to `dst`.
void append_atoms_prefixed(ChemComp& dst, const ChemComp& src, int side) {
  dst.atoms.reserve(dst.atoms.size() + src.atoms.size());
  for (const auto& a : src.atoms)
    dst.atoms.push_back(prefix_atom(a, side));
}

// Append every bond of src.rt into dst.rt, prefixing both ids.
void append_bonds_prefixed(ChemComp& dst, const ChemComp& src, int side) {
  dst.rt.bonds.reserve(dst.rt.bonds.size() + src.rt.bonds.size());
  for (const auto& b : src.rt.bonds) {
    Restraints::Bond out = b;
    out.id1 = {0, make_joined_atom_id(side, b.id1.atom)};
    out.id2 = {0, make_joined_atom_id(side, b.id2.atom)};
    dst.rt.bonds.push_back(out);
  }
}

// Delete the named atom by id, plus every bond/angle/torsion/chirality/plane
// touching it. Planes whose atom-list drops below 3 atoms are removed.
void purge_atom(ChemComp& cc, const std::string& atom_id) {
  vector_remove_if(cc.atoms,
                   [&](const ChemComp::Atom& a) { return a.id == atom_id; });
  vector_remove_if(cc.rt.bonds, [&](const Restraints::Bond& b) {
    return b.id1.atom == atom_id || b.id2.atom == atom_id;
  });
  vector_remove_if(cc.rt.angles, [&](const Restraints::Angle& a) {
    return a.id1.atom == atom_id || a.id2.atom == atom_id ||
           a.id3.atom == atom_id;
  });
  vector_remove_if(cc.rt.torsions, [&](const Restraints::Torsion& t) {
    return t.id1.atom == atom_id || t.id2.atom == atom_id ||
           t.id3.atom == atom_id || t.id4.atom == atom_id;
  });
  vector_remove_if(cc.rt.chirs, [&](const Restraints::Chirality& c) {
    return c.id_ctr.atom == atom_id || c.id1.atom == atom_id ||
           c.id2.atom == atom_id || c.id3.atom == atom_id;
  });
  for (auto& plane : cc.rt.planes)
    vector_remove_if(plane.ids, [&](const Restraints::AtomId& a) {
      return a.atom == atom_id;
    });
  vector_remove_if(cc.rt.planes,
                   [](const Restraints::Plane& p) { return p.ids.size() < 3; });
}

// True if `idx`'s atom is hydrogen and has no heavy neighbour left.
bool is_orphan_hydrogen(const ChemComp& cc, const AceBondAdjacency& adj,
                        size_t idx) {
  if (!cc.atoms[idx].is_hydrogen())
    return false;
  for (const auto& nb : adj[idx])
    if (!cc.atoms[nb.idx].is_hydrogen())
      return false;
  return true;
}

// Heavy-atom bond-order sum at `idx` (ignores H neighbours).
float heavy_bond_order_sum(const ChemComp& cc, const AceBondAdjacency& adj,
                           size_t idx) {
  float sum = 0.0f;
  for (const auto& nb : adj[idx])
    if (!cc.atoms[nb.idx].is_hydrogen())
      sum += order_of(nb.type);
  return sum;
}

// Names of H atoms bonded to `idx`, sorted alphanumerically.
std::vector<std::string> bonded_hydrogen_ids(const ChemComp& cc,
                                             const AceBondAdjacency& adj,
                                             size_t idx) {
  std::vector<std::string> out;
  for (const auto& nb : adj[idx])
    if (cc.atoms[nb.idx].is_hydrogen())
      out.push_back(cc.atoms[nb.idx].id);
  std::sort(out.begin(), out.end());
  return out;
}

// Apply minimum-|charge| auto-H rule to one atom. Edits cc in place.
void apply_auto_h_one_atom(ChemComp& cc, const std::string& atom_id) {
  auto atom_index = cc.make_atom_index();
  auto it = atom_index.find(atom_id);
  if (it == atom_index.end())
    return;
  size_t idx = it->second;
  const ChemComp::Atom& a = cc.atoms[idx];

  if (a.is_hydrogen() || a.el.is_metal())
    return;
  int expected = expected_valence_for_nonmetal(a.el);
  if (expected == 0)
    return;

  AceBondAdjacency adj = build_bond_adjacency(cc, atom_index);
  int heavy = static_cast<int>(std::round(heavy_bond_order_sum(cc, adj, idx)));

  int target_charge;
  int target_H;
  if (heavy <= expected) {
    target_charge = 0;
    target_H = expected - heavy;
  } else {
    target_charge = heavy - expected;
    if (target_charge > 2)
      throw std::runtime_error(
          "auto-H: atom '" + atom_id + "' is over-valenced after link " +
          "insertion (heavy bond-order sum " + std::to_string(heavy) +
          " vs. expected valence " + std::to_string(expected) +
          "); a heavy-atom DELETE may be missing.");
    target_H = 0;
  }

  std::vector<std::string> hs = bonded_hydrogen_ids(cc, adj, idx);
  if (static_cast<int>(hs.size()) > target_H)
    for (size_t i = hs.size(); i > static_cast<size_t>(target_H); --i)
      purge_atom(cc, hs[i - 1]);

  atom_index = cc.make_atom_index();
  auto it2 = atom_index.find(atom_id);
  if (it2 != atom_index.end())
    cc.atoms[it2->second].charge = static_cast<float>(target_charge);
}

}  // anonymous namespace

// ─── make_joined_chemcomp ────────────────────────────────────────────────────

ChemComp make_joined_chemcomp(const ChemComp& cc1, const ChemComp& cc2,
                              const LinkSpec& spec) {
  if (cc1.find_atom(spec.atom1) == cc1.atoms.end())
    throw std::runtime_error(
        "Link atom '" + spec.atom1 + "' not found in component '" +
        cc1.name + "' (side 1).");
  if (cc2.find_atom(spec.atom2) == cc2.atoms.end())
    throw std::runtime_error(
        "Link atom '" + spec.atom2 + "' not found in component '" +
        cc2.name + "' (side 2).");

  ChemComp joined;
  joined.name = spec.id;
  joined.group = ChemComp::Group::NonPolymer;
  joined.has_coordinates = cc1.has_coordinates && cc2.has_coordinates;
  append_atoms_prefixed(joined, cc1, 1);
  append_atoms_prefixed(joined, cc2, 2);
  append_bonds_prefixed(joined, cc1, 1);
  append_bonds_prefixed(joined, cc2, 2);

  // Apply user heavy-atom DELETEs.
  for (const auto& d : spec.deletions) {
    if (d.first != 1 && d.first != 2)
      throw std::runtime_error(
          "DELETE side must be 1 or 2 (got " + std::to_string(d.first) + ")");
    std::string j = make_joined_atom_id(d.first, d.second);
    if (std::none_of(joined.atoms.begin(), joined.atoms.end(),
                     [&](const ChemComp::Atom& a) { return a.id == j; }))
      throw std::runtime_error(
          "DELETE: atom '" + d.second + "' not found on side " +
          std::to_string(d.first) + " (looked for joined id '" + j + "')");
    purge_atom(joined, j);
  }

  // Cascade-delete orphan H atoms (e.g. HO1 after DELETE O1).
  while (true) {
    auto atom_index = joined.make_atom_index();
    AceBondAdjacency adj = build_bond_adjacency(joined, atom_index);
    std::vector<std::string> to_drop;
    for (size_t i = 0; i < joined.atoms.size(); ++i)
      if (is_orphan_hydrogen(joined, adj, i))
        to_drop.push_back(joined.atoms[i].id);
    if (to_drop.empty()) break;
    for (const auto& id : to_drop)
      purge_atom(joined, id);
  }

  // Insert the link bond.
  std::string lj1 = make_joined_atom_id(1, spec.atom1);
  std::string lj2 = make_joined_atom_id(2, spec.atom2);
  if (std::none_of(joined.atoms.begin(), joined.atoms.end(),
                   [&](const ChemComp::Atom& a) { return a.id == lj1; }))
    throw std::runtime_error("Link atom side 1 was deleted: " + lj1);
  if (std::none_of(joined.atoms.begin(), joined.atoms.end(),
                   [&](const ChemComp::Atom& a) { return a.id == lj2; }))
    throw std::runtime_error("Link atom side 2 was deleted: " + lj2);

  Restraints::Bond link_bond;
  link_bond.id1 = {0, lj1};
  link_bond.id2 = {0, lj2};
  link_bond.type = spec.bond_type;
  link_bond.aromatic = is_aromatic_or_deloc(spec.bond_type);
  link_bond.value = NAN;
  link_bond.esd = NAN;
  link_bond.value_nucleus = NAN;
  link_bond.esd_nucleus = NAN;
  joined.rt.bonds.push_back(link_bond);

  // Auto-H + charge on the link atoms.
  apply_auto_h_one_atom(joined, lj1);
  apply_auto_h_one_atom(joined, lj2);

  // Validation: every heavy atom must still have at least one bond.
  auto atom_index = joined.make_atom_index();
  AceBondAdjacency adj = build_bond_adjacency(joined, atom_index);
  for (size_t i = 0; i < joined.atoms.size(); ++i)
    if (!joined.atoms[i].is_hydrogen() && adj[i].empty())
      throw std::runtime_error(
          "Joined graph: heavy atom '" + joined.atoms[i].id +
          "' has no bonds left after DELETEs; check the link spec.");

  return joined;
}

// ─── prepare_chemlink dispatch helpers (file-local) ─────────────────────────

namespace {

// Function codes used in ChemMod restraint rows: matches chem_mod_type() in
// monlib.cpp ('a' = add, 'c' = change, 'd' = delete).
constexpr int kMod_Add    = 'a';
constexpr int kMod_Change = 'c';
constexpr int kMod_Delete = 'd';

// Tolerance for considering two restraint values "the same" so we don't
// emit a `change` row for cosmetic noise. These are deliberately on the
// loose side: AceDRG emits 'change' only when the atom-type-tuple lookup
// key changes, so the resulting value almost always differs by more than
// noise. Tight tolerances here produce spurious changes for values that
// got re-derived with rounding drift.
constexpr double kBondValueTol  = 0.003;   // Angstrom
constexpr double kBondEsdTol    = 0.0015;  // sigma
constexpr double kAngleValueTol = 0.3;     // degree
constexpr double kAngleEsdTol   = 0.1;     // degree

// Re-tag every AtomId in `r` with `comp` (rebinds the meaning of the comp
// field). Use kMod_* for ChemMod rt rows, side numbers for ChemLink rt rows.
Restraints::Bond retag(const Restraints::Bond& src, int new_comp,
                       const std::string& a1, const std::string& a2) {
  Restraints::Bond out = src;
  out.id1 = {new_comp, a1};
  out.id2 = {new_comp, a2};
  return out;
}
Restraints::Angle retag(const Restraints::Angle& src, int new_comp,
                        const std::string& a1, const std::string& a2,
                        const std::string& a3) {
  Restraints::Angle out = src;
  out.id1 = {new_comp, a1};
  out.id2 = {new_comp, a2};
  out.id3 = {new_comp, a3};
  return out;
}
Restraints::Torsion retag(const Restraints::Torsion& src, int new_comp,
                          const std::string& a1, const std::string& a2,
                          const std::string& a3, const std::string& a4) {
  Restraints::Torsion out = src;
  out.id1 = {new_comp, a1};
  out.id2 = {new_comp, a2};
  out.id3 = {new_comp, a3};
  out.id4 = {new_comp, a4};
  return out;
}
Restraints::Chirality retag(const Restraints::Chirality& src, int new_comp,
                            const std::string& ac, const std::string& a1,
                            const std::string& a2, const std::string& a3) {
  Restraints::Chirality out = src;
  out.id_ctr = {new_comp, ac};
  out.id1 = {new_comp, a1};
  out.id2 = {new_comp, a2};
  out.id3 = {new_comp, a3};
  return out;
}

// True iff two bonds differ enough to warrant a 'change' mod row.
// We compare the electron-cloud value/sigma only — value_nucleus is a
// derived quantity that tends to drift cosmetically across pipeline runs.
bool bond_value_changed(const Restraints::Bond& a, const Restraints::Bond& b) {
  return std::abs(a.value - b.value) > kBondValueTol ||
         std::abs(a.esd   - b.esd)   > kBondEsdTol ||
         a.type != b.type;
}
bool angle_value_changed(const Restraints::Angle& a, const Restraints::Angle& b) {
  return std::abs(a.value - b.value) > kAngleValueTol ||
         std::abs(a.esd   - b.esd)   > kAngleEsdTol;
}
bool torsion_value_changed(const Restraints::Torsion& a, const Restraints::Torsion& b) {
  // Torsions are circular; treat 180 and -180 as the same. Allow esd/period
  // changes to trigger a mod though.
  double diff = std::abs(a.value - b.value);
  while (diff > 360.0) diff -= 360.0;
  if (diff > 180.0) diff = 360.0 - diff;
  return diff > 0.5 || std::abs(a.esd - b.esd) > 0.5 || a.period != b.period;
}

// Look-up bond / angle / torsion in a parent ChemComp's rt by atom names
// (order-insensitive for bonds, central-atom-fixed for angles, edge-fixed
// for torsions). Returns nullptr if absent.
const Restraints::Bond* find_bond(const Restraints& rt,
                                  const std::string& a, const std::string& b) {
  for (const auto& x : rt.bonds)
    if ((x.id1.atom == a && x.id2.atom == b) ||
        (x.id1.atom == b && x.id2.atom == a))
      return &x;
  return nullptr;
}
const Restraints::Angle* find_angle(const Restraints& rt,
                                    const std::string& a, const std::string& b,
                                    const std::string& c) {
  // Central atom is b; outers are unordered.
  for (const auto& x : rt.angles)
    if (x.id2.atom == b &&
        ((x.id1.atom == a && x.id3.atom == c) ||
         (x.id1.atom == c && x.id3.atom == a)))
      return &x;
  return nullptr;
}
const Restraints::Torsion* find_torsion(const Restraints& rt,
                                        const std::string& a, const std::string& b,
                                        const std::string& c, const std::string& d) {
  // Edge bc is fixed; (a,d) may be reversed direction.
  for (const auto& x : rt.torsions) {
    if (x.id2.atom == b && x.id3.atom == c &&
        x.id1.atom == a && x.id4.atom == d) return &x;
    if (x.id2.atom == c && x.id3.atom == b &&
        x.id1.atom == d && x.id4.atom == a) return &x;
  }
  return nullptr;
}
const Restraints::Chirality* find_chir(const Restraints& rt,
                                       const std::string& ctr,
                                       const std::string& a,
                                       const std::string& b,
                                       const std::string& c) {
  // Centre fixed; the three substituents may appear in any order.
  std::set<std::string> want{a, b, c};
  for (const auto& x : rt.chirs)
    if (x.id_ctr.atom == ctr &&
        std::set<std::string>{x.id1.atom, x.id2.atom, x.id3.atom} == want)
      return &x;
  return nullptr;
}

// Sides of every atom in a restraint, returned as a pair (set_of_sides, all_atom_names).
struct Membership {
  std::set<int> sides;        // 1, 2, or {1, 2}
  bool all_same() const { return sides.size() == 1; }
  int  single_side() const { return *sides.begin(); }
};
Membership classify_membership(std::initializer_list<std::string> joined_ids) {
  Membership m;
  for (const auto& jid : joined_ids) {
    auto s = split_joined_atom_id(jid);
    m.sides.insert(s.first);
  }
  return m;
}

// Emit one bond row in the right sink. Returns true if a row was added,
// so callers can track which monomer entries got 'change' vs 'unchanged'.
void dispatch_bond(const Restraints::Bond& jb,
                   const ChemComp& cc1, const ChemComp& cc2,
                   LinkGenerationResult& out) {
  auto s1 = split_joined_atom_id(jb.id1.atom);
  auto s2 = split_joined_atom_id(jb.id2.atom);
  if (s1.first == s2.first && s1.first != 0) {
    int side = s1.first;
    const ChemComp& parent = (side == 1) ? cc1 : cc2;
    ChemMod& mod = (side == 1) ? out.mod1 : out.mod2;
    const Restraints::Bond* orig = find_bond(parent.rt, s1.second, s2.second);
    if (!orig)
      mod.rt.bonds.push_back(retag(jb, kMod_Add, s1.second, s2.second));
    else if (bond_value_changed(*orig, jb))
      mod.rt.bonds.push_back(retag(jb, kMod_Change, s1.second, s2.second));
  } else {
    Restraints::Bond linkrow = jb;
    linkrow.id1 = {s1.first, s1.second};
    linkrow.id2 = {s2.first, s2.second};
    out.link.rt.bonds.push_back(linkrow);
  }
}

void dispatch_angle(const Restraints::Angle& ja,
                    const ChemComp& cc1, const ChemComp& cc2,
                    LinkGenerationResult& out) {
  auto s1 = split_joined_atom_id(ja.id1.atom);
  auto s2 = split_joined_atom_id(ja.id2.atom);
  auto s3 = split_joined_atom_id(ja.id3.atom);
  Membership m = classify_membership({ja.id1.atom, ja.id2.atom, ja.id3.atom});
  if (m.all_same()) {
    int side = m.single_side();
    const ChemComp& parent = (side == 1) ? cc1 : cc2;
    ChemMod& mod = (side == 1) ? out.mod1 : out.mod2;
    const Restraints::Angle* orig =
        find_angle(parent.rt, s1.second, s2.second, s3.second);
    if (!orig)
      mod.rt.angles.push_back(retag(ja, kMod_Add, s1.second, s2.second, s3.second));
    else if (angle_value_changed(*orig, ja))
      mod.rt.angles.push_back(retag(ja, kMod_Change, s1.second, s2.second, s3.second));
  } else {
    Restraints::Angle linkrow = ja;
    linkrow.id1 = {s1.first, s1.second};
    linkrow.id2 = {s2.first, s2.second};
    linkrow.id3 = {s3.first, s3.second};
    out.link.rt.angles.push_back(linkrow);
  }
}

void dispatch_torsion(const Restraints::Torsion& jt,
                      const ChemComp& cc1, const ChemComp& cc2,
                      LinkGenerationResult& out) {
  auto s1 = split_joined_atom_id(jt.id1.atom);
  auto s2 = split_joined_atom_id(jt.id2.atom);
  auto s3 = split_joined_atom_id(jt.id3.atom);
  auto s4 = split_joined_atom_id(jt.id4.atom);
  Membership m = classify_membership(
      {jt.id1.atom, jt.id2.atom, jt.id3.atom, jt.id4.atom});
  if (m.all_same()) {
    int side = m.single_side();
    const ChemComp& parent = (side == 1) ? cc1 : cc2;
    ChemMod& mod = (side == 1) ? out.mod1 : out.mod2;
    const Restraints::Torsion* orig = find_torsion(
        parent.rt, s1.second, s2.second, s3.second, s4.second);
    if (!orig)
      mod.rt.torsions.push_back(
          retag(jt, kMod_Add, s1.second, s2.second, s3.second, s4.second));
    else if (torsion_value_changed(*orig, jt))
      mod.rt.torsions.push_back(
          retag(jt, kMod_Change, s1.second, s2.second, s3.second, s4.second));
  } else {
    Restraints::Torsion linkrow = jt;
    linkrow.id1 = {s1.first, s1.second};
    linkrow.id2 = {s2.first, s2.second};
    linkrow.id3 = {s3.first, s3.second};
    linkrow.id4 = {s4.first, s4.second};
    out.link.rt.torsions.push_back(linkrow);
  }
}

void dispatch_chirality(const Restraints::Chirality& jc,
                        const ChemComp& cc1, const ChemComp& cc2,
                        LinkGenerationResult& out) {
  auto sc = split_joined_atom_id(jc.id_ctr.atom);
  auto s1 = split_joined_atom_id(jc.id1.atom);
  auto s2 = split_joined_atom_id(jc.id2.atom);
  auto s3 = split_joined_atom_id(jc.id3.atom);
  Membership m = classify_membership(
      {jc.id_ctr.atom, jc.id1.atom, jc.id2.atom, jc.id3.atom});
  if (m.all_same()) {
    int side = m.single_side();
    const ChemComp& parent = (side == 1) ? cc1 : cc2;
    ChemMod& mod = (side == 1) ? out.mod1 : out.mod2;
    const Restraints::Chirality* orig = find_chir(
        parent.rt, sc.second, s1.second, s2.second, s3.second);
    if (!orig || orig->sign != jc.sign)
      mod.rt.chirs.push_back(retag(jc, orig ? kMod_Change : kMod_Add,
                                   sc.second, s1.second, s2.second, s3.second));
  } else {
    Restraints::Chirality linkrow = jc;
    linkrow.id_ctr = {sc.first, sc.second};
    linkrow.id1 = {s1.first, s1.second};
    linkrow.id2 = {s2.first, s2.second};
    linkrow.id3 = {s3.first, s3.second};
    out.link.rt.chirs.push_back(linkrow);
  }
}

// Planes are special: each plane is a *set* of atoms, possibly spanning sides.
// We compare by atom set (not by label) — the joined-graph pipeline may
// label a plane differently from how the standalone monomer labelled the
// same chemistry (e.g. LYS's carboxylate plane comes back labelled "plan-2"
// rather than its original "plan-1").
void dispatch_plane(const Restraints::Plane& jp,
                    const ChemComp& cc1, const ChemComp& cc2,
                    LinkGenerationResult& out) {
  std::set<int> sides;
  for (const auto& a : jp.ids)
    sides.insert(split_joined_atom_id(a.atom).first);
  if (sides.size() == 1) {
    int side = *sides.begin();
    const ChemComp& parent = (side == 1) ? cc1 : cc2;
    ChemMod& mod = (side == 1) ? out.mod1 : out.mod2;
    // Collect this plane's atom set (unprefixed names).
    std::set<std::string> jp_atoms;
    for (const auto& a : jp.ids)
      jp_atoms.insert(split_joined_atom_id(a.atom).second);
    // Look for a parent plane whose atom set matches (set equality).
    const Restraints::Plane* orig = nullptr;
    for (const auto& p : parent.rt.planes) {
      std::set<std::string> p_atoms;
      for (const auto& a : p.ids) p_atoms.insert(a.atom);
      if (p_atoms == jp_atoms) { orig = &p; break; }
    }
    if (orig) {
      // Same plane already in the monomer. Only emit a 'change' if the
      // dist_esd differs noticeably. Tolerance matches kBondEsdTol.
      if (std::abs(orig->esd - jp.esd) > kBondEsdTol) {
        Restraints::Plane out_plane;
        out_plane.label = orig->label;
        out_plane.esd = jp.esd;
        for (const auto& a : jp.ids)
          out_plane.ids.push_back(
              {kMod_Change, split_joined_atom_id(a.atom).second});
        mod.rt.planes.push_back(out_plane);
      }
      return;
    }
    // New plane introduced by the link's re-typing.
    Restraints::Plane out_plane;
    out_plane.label = jp.label;
    out_plane.esd = jp.esd;
    for (const auto& a : jp.ids)
      out_plane.ids.push_back({kMod_Add, split_joined_atom_id(a.atom).second});
    mod.rt.planes.push_back(out_plane);
  } else {
    Restraints::Plane linkrow;
    linkrow.label = jp.label;
    linkrow.esd = jp.esd;
    for (const auto& a : jp.ids) {
      auto s = split_joined_atom_id(a.atom);
      linkrow.ids.push_back({s.first, s.second});
    }
    out.link.rt.planes.push_back(linkrow);
  }
}

// Emit 'delete' mod rows for every cc-side restraint that references an atom
// no longer in the joined CC (i.e. a deleted-by-user or auto-deleted atom).
void emit_delete_rows(const ChemComp& parent_cc, int side,
                      const std::set<std::string>& deleted_atoms,
                      ChemMod& mod) {
  auto touches = [&](const Restraints::AtomId& a) {
    return deleted_atoms.count(a.atom) != 0;
  };
  for (const auto& b : parent_cc.rt.bonds)
    if (touches(b.id1) || touches(b.id2))
      mod.rt.bonds.push_back(retag(b, kMod_Delete, b.id1.atom, b.id2.atom));
  for (const auto& a : parent_cc.rt.angles)
    if (touches(a.id1) || touches(a.id2) || touches(a.id3))
      mod.rt.angles.push_back(retag(a, kMod_Delete,
                                    a.id1.atom, a.id2.atom, a.id3.atom));
  for (const auto& t : parent_cc.rt.torsions)
    if (touches(t.id1) || touches(t.id2) || touches(t.id3) || touches(t.id4))
      mod.rt.torsions.push_back(retag(t, kMod_Delete,
                                      t.id1.atom, t.id2.atom,
                                      t.id3.atom, t.id4.atom));
  for (const auto& c : parent_cc.rt.chirs)
    if (touches(c.id_ctr) || touches(c.id1) || touches(c.id2) || touches(c.id3))
      mod.rt.chirs.push_back(retag(c, kMod_Delete,
                                   c.id_ctr.atom, c.id1.atom,
                                   c.id2.atom, c.id3.atom));
  // For planes: if a plane has at least one atom that got deleted, emit a
  // 'delete'-plane row with just the deleted atoms listed (AceDRG style).
  for (const auto& p : parent_cc.rt.planes) {
    Restraints::Plane delp;
    delp.label = p.label;
    delp.esd = p.esd;
    for (const auto& a : p.ids)
      if (touches(a))
        delp.ids.push_back({kMod_Delete, a.atom});
    if (!delp.ids.empty())
      mod.rt.planes.push_back(delp);
  }
  // Mark the deleted atoms themselves in atom_mods.
  for (const auto& aid : deleted_atoms) {
    auto src = parent_cc.find_atom(aid);
    ChemMod::AtomMod am{
        kMod_Delete, aid, "", Element(El::X), 0.0f, ""};
    if (src != parent_cc.atoms.end()) {
      am.el = src->el;
      am.charge = src->charge;
      am.chem_type = src->chem_type;
    }
    mod.atom_mods.push_back(am);
  }
  (void)side;
}

}  // anonymous namespace

// ─── prepare_chemlink ────────────────────────────────────────────────────────

LinkGenerationResult prepare_chemlink(const ChemComp& cc1, const ChemComp& cc2,
                                      const LinkSpec& spec,
                                      const AcedrgTables& tables) {
  // 1. Build joined skeleton (auto-H done, link bond inserted).
  ChemComp joined = make_joined_chemcomp(cc1, cc2, spec);

  // 2. Run the AceDRG monomer pipeline over the joined graph to fill in
  //    bond/angle values, derive torsions/chirs/planes, assign types.
  PrepareChemcompOptions opts;
  prepare_chemcomp(joined, tables, opts);

  // 3. Compute deleted-atom sets per side by name-set difference.
  std::set<std::string> deleted_side1, deleted_side2;
  for (const auto& a : cc1.atoms)
    if (joined.find_atom(make_joined_atom_id(1, a.id)) == joined.atoms.end())
      deleted_side1.insert(a.id);
  for (const auto& a : cc2.atoms)
    if (joined.find_atom(make_joined_atom_id(2, a.id)) == joined.atoms.end())
      deleted_side2.insert(a.id);

  // 4. Initialise the output containers.
  LinkGenerationResult out;
  out.link.id = spec.id;
  out.link.name = spec.id;
  out.link.side1.comp = cc1.name;
  out.link.side1.mod  = cc1.name + "m1";
  out.link.side1.group = cc1.group;
  out.link.side2.comp = cc2.name;
  out.link.side2.mod  = cc2.name + "m1";
  out.link.side2.group = cc2.group;
  out.mod1.id = cc1.name + "m1";
  out.mod1.name = cc1.name;
  out.mod1.comp_id = cc1.name;
  out.mod1.group_id = ChemComp::group_str(cc1.group);
  out.mod2.id = cc2.name + "m1";
  out.mod2.name = cc2.name;
  out.mod2.comp_id = cc2.name;
  out.mod2.group_id = ChemComp::group_str(cc2.group);

  // 5. Dispatch every derived restraint to one of the three sinks.
  for (const auto& b : joined.rt.bonds)   dispatch_bond(b, cc1, cc2, out);
  for (const auto& a : joined.rt.angles)  dispatch_angle(a, cc1, cc2, out);
  for (const auto& t : joined.rt.torsions) dispatch_torsion(t, cc1, cc2, out);
  for (const auto& c : joined.rt.chirs)   dispatch_chirality(c, cc1, cc2, out);
  for (const auto& p : joined.rt.planes)  dispatch_plane(p, cc1, cc2, out);

  // 6. Emit 'delete' rows for every parent restraint touching a deleted atom.
  emit_delete_rows(cc1, 1, deleted_side1, out.mod1);
  emit_delete_rows(cc2, 2, deleted_side2, out.mod2);

  return out;
}

}  // namespace gemmi
