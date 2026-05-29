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

#include "gemmi/ace_graph.hpp"      // expected_valence_for_nonmetal, build_bond_adjacency
#include "gemmi/acedrg_tables.hpp"  // AcedrgTables (used in prepare_chemlink stub)
#include "gemmi/util.hpp"           // vector_remove_if

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

// ─── prepare_chemlink — stub for now (Task 6) ────────────────────────────────

LinkGenerationResult prepare_chemlink(const ChemComp&, const ChemComp&,
                                      const LinkSpec&, const AcedrgTables&) {
  throw std::runtime_error("prepare_chemlink: not implemented yet");
}

}  // namespace gemmi
