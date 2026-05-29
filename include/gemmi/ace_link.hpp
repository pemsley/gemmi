// Copyright 2026 Global Phasing Ltd.
//
// Link-restraint generation for `gemmi drg --link`.
//
// AceDRG (Long et al., Acta Cryst. D73, 112-122, 2017) covers per-monomer
// restraint generation; this module extends the same atom-typing and
// table-lookup pipeline across a residue boundary so that a CCP4-style
// _chem_link block + two _chem_mod blocks can be derived for a covalent
// link between two ChemComps.
//
// The approach is:
//   1. Build a "joined ChemComp" containing atoms from both residues with
//      disambiguated IDs ("1/NZ", "2/C4A"), the link bond present and any
//      user-specified heavy-atom DELETEs applied. Implicit hydrogens and
//      formal charges are then normalised by valence.
//   2. Run the existing AceDRG classify/lookup machinery on the joined
//      object (it cares only about atom names and adjacency, not about
//      which "side" an atom belongs to).
//   3. Dispatch every derived restraint to one of three sinks based on
//      whether its atoms all belong to one side (-> ChemMod for that side)
//      or span both (-> ChemLink).

#ifndef GEMMI_ACE_LINK_HPP_
#define GEMMI_ACE_LINK_HPP_

#include <string>
#include <utility>
#include <vector>
#include "cifdoc.hpp"    // for cif::Document
#include "chemcomp.hpp"  // for ChemComp, BondType
#include "monlib.hpp"    // for ChemLink, ChemMod

namespace gemmi {

struct AcedrgTables;

/// User-supplied description of one chemical link to generate.
///
/// Modelled on AceDRG's one-liner instruction:
///   LINK: RES-NAME-1 LYS ATOM-NAME-1 NZ
///         RES-NAME-2 PLP ATOM-NAME-2 C4A
///         DELETE ATOM O4A 2  BOND-TYPE DOUBLE
struct LinkSpec {
  std::string id;            ///< Link identifier (e.g. "LYS-PLP")
  std::string comp1_name;    ///< Side 1 component three-letter code
  std::string comp2_name;    ///< Side 2 component three-letter code
  std::string atom1;         ///< Link atom on side 1 (bare name)
  std::string atom2;         ///< Link atom on side 2 (bare name)
  BondType bond_type = BondType::Single;  ///< Link bond order (authoritative)
  /// Heavy-atom deletions supplied by the user.
  /// Each pair is (side, atom_name). H atoms are *not* listed — they are
  /// cascade-deleted if orphaned, and surplus H's adjacent to the link are
  /// removed automatically by valence.
  std::vector<std::pair<int, std::string>> deletions;
};

/// Identifier prefix scheme for joined-graph atom IDs.
/// "1/NZ" means atom "NZ" of side 1; "2/C4A" means "C4A" of side 2.
GEMMI_DLL std::string make_joined_atom_id(int side, const std::string& bare);

/// Inverse of make_joined_atom_id. Returns {side, bare_name}.
/// `side` is 0 if `joined_id` has no recognised prefix.
GEMMI_DLL std::pair<int, std::string> split_joined_atom_id(
    const std::string& joined_id);

/// Build a single ChemComp representing both residues with the link bond
/// inserted and the user-supplied heavy-atom deletions applied (plus
/// orphaned-H cascade and minimum-|charge| auto-H normalisation).
///
/// Atom IDs in the returned ChemComp are prefixed (`"1/<name>"` or
/// `"2/<name>"`). Bonds copied from cc1/cc2 keep their original types and
/// the link bond is appended. ChemComp::name is set to `spec.id`.
///
/// On contract violation (over-valence, orphan atom, unresolvable charge)
/// this throws std::runtime_error with a diagnostic naming the bad atom.
GEMMI_DLL ChemComp make_joined_chemcomp(const ChemComp& cc1,
                                        const ChemComp& cc2,
                                        const LinkSpec& spec);

/// Output container for prepare_chemlink. The two ChemMod ids are
/// "<cc1.name>m1" and "<cc2.name>m1" per AceDRG convention.
struct LinkGenerationResult {
  ChemLink link;
  ChemMod  mod1;
  ChemMod  mod2;
};

/// Generate the ChemLink + per-side ChemMod blocks for a covalent link
/// between cc1 and cc2 using AceDRG-style atom typing on the joined graph.
///
/// Dispatch rule per derived restraint:
///   * all atoms on side 1 & type-tuple changed vs. standalone cc1 -> mod1
///   * all atoms on side 2 & type-tuple changed vs. standalone cc2 -> mod2
///   * all atoms on one side & type-tuple unchanged -> dropped (monomer
///     dictionary entry stands)
///   * atoms span both sides -> link
GEMMI_DLL LinkGenerationResult prepare_chemlink(const ChemComp& cc1,
                                                const ChemComp& cc2,
                                                const LinkSpec& spec,
                                                const AcedrgTables& tables);

/// Serialise the full link dictionary (CCP4 monomer-library format) into
/// `doc`: data_comp_list, data_mod_list, data_link_list summary blocks
/// followed by data_comp_<NAME> for each monomer, data_mod_<NAME>m1 for
/// each mod, and data_link_<id> for the link itself.
///
/// `instruction` is copied verbatim into a _CCP4_AceDRG_link_generation
/// .instruction text field on data_program_info (mirrors AceDRG's audit
/// trail). Pass an empty string to skip.
GEMMI_DLL void write_link_dictionary(const ChemComp& cc1,
                                     const ChemComp& cc2,
                                     const LinkGenerationResult& res,
                                     cif::Document& doc,
                                     const std::string& instruction = "");

} // namespace gemmi
#endif
