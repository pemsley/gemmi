// Copyright 2026 Global Phasing Ltd.
//
// SQLite backend for the heavy AceDRG bond and angle tables.
// One-time conversion via `gemmi drg --build-table-db` produces a
// single .sqlite file from the directory of ASCII tables; subsequent
// runs open that file and query per-molecule on demand instead of
// pre-loading ~1 GB into RAM.
//
// Falls back gracefully when gemmi was built without libsqlite3:
// build_acedrg_sqlite() throws std::runtime_error, callers can fall back
// to the in-memory or binary-cache loaders.

#include "gemmi/acedrg_tables.hpp"
#include "gemmi/fail.hpp"
#include "gemmi/util.hpp"

#include <cstdio>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string>

#ifdef GEMMI_HAS_SQLITE3
#include <sqlite3.h>
#endif

namespace gemmi {

#ifdef GEMMI_HAS_SQLITE3

namespace {

// Minimal RAII wrappers around sqlite3 and sqlite3_stmt.
struct SqliteDB {
  sqlite3* db = nullptr;
  ~SqliteDB() { if (db) sqlite3_close(db); }
  void open(const std::string& path, int flags) {
    if (sqlite3_open_v2(path.c_str(), &db, flags, nullptr) != SQLITE_OK) {
      std::string err = db ? sqlite3_errmsg(db) : "unknown";
      throw std::runtime_error("acedrg-db: cannot open " + path + ": " + err);
    }
  }
  void exec(const char* sql) {
    char* err = nullptr;
    if (sqlite3_exec(db, sql, nullptr, nullptr, &err) != SQLITE_OK) {
      std::string e = err ? err : "unknown";
      sqlite3_free(err);
      throw std::runtime_error("acedrg-db: SQL failed: " + e + " (" + sql + ")");
    }
  }
};

struct SqliteStmt {
  sqlite3_stmt* st = nullptr;
  ~SqliteStmt() { if (st) sqlite3_finalize(st); }
  void prepare(sqlite3* db, const char* sql) {
    if (sqlite3_prepare_v2(db, sql, -1, &st, nullptr) != SQLITE_OK)
      throw std::runtime_error(std::string("acedrg-db: prepare failed: ") +
                               sqlite3_errmsg(db));
  }
};

// Skip whitespace, return pointer to next non-blank.
inline const char* skip_blank_db(const char* p) {
  while (*p == ' ' || *p == '\t') ++p;
  return p;
}
inline const char* skip_word_db(const char* p) {
  while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') ++p;
  return p;
}
inline bool is_skip_line_db(const char* line) {
  const char* p = skip_blank_db(line);
  return *p == '\0' || *p == '#' || *p == '\n' || *p == '\r';
}

// Read every numbered N.table file in a directory; returns sorted ints.
std::vector<int> list_table_files(const std::string& dir) {
  std::vector<int> out;
  FILE* p = popen(("ls -1 " + dir + " 2>/dev/null").c_str(), "r");
  if (!p) return out;
  char line[256];
  while (std::fgets(line, sizeof(line), p)) {
    int n;
    if (std::sscanf(line, "%d.table", &n) == 1)
      out.push_back(n);
  }
  pclose(p);
  std::sort(out.begin(), out.end());
  return out;
}

// Load the coded -> full-type atom-type map (allAtomTypesFromMolsCoded.list).
std::unordered_map<std::string, std::string>
load_atom_codes(const std::string& path) {
  std::unordered_map<std::string, std::string> out;
  FILE* f = std::fopen(path.c_str(), "r");
  if (!f) return out;
  char line[2048];
  while (std::fgets(line, sizeof(line), f)) {
    if (is_skip_line_db(line)) continue;
    const char* p = line;
    const char* s = skip_blank_db(p); p = skip_word_db(s);
    std::string code(s, p - s);
    s = skip_blank_db(p); p = skip_word_db(s);
    std::string full(s, p - s);
    if (!code.empty() && !full.empty())
      out.emplace(std::move(code), std::move(full));
  }
  std::fclose(f);
  return out;
}

inline std::string prefix_before_db(const std::string& s, char c) {
  auto pos = s.find(c);
  return pos == std::string::npos ? s : s.substr(0, pos);
}

// --- bond table converter --------------------------------------------------

void convert_bond_tables(SqliteDB& db, const std::string& bond_dir,
                         const std::unordered_map<std::string, std::string>& codes) {
  db.exec(
    "DROP TABLE IF EXISTS bond_entries;"
    "CREATE TABLE bond_entries ("
    "  ha1 INTEGER, ha2 INTEGER,"
    "  hybr_comb TEXT, in_ring TEXT,"
    "  a1_nb2 TEXT, a2_nb2 TEXT, a1_nb TEXT, a2_nb TEXT,"
    "  a1_type_m TEXT, a2_type_m TEXT,"
    "  a1_type_f TEXT, a2_type_f TEXT,"
    "  value REAL, sigma REAL, count INTEGER,"
    "  value_1d REAL, sigma_1d REAL, count_1d INTEGER"
    ");"
  );

  SqliteStmt ins;
  ins.prepare(db.db,
    "INSERT INTO bond_entries (ha1, ha2, hybr_comb, in_ring,"
    "  a1_nb2, a2_nb2, a1_nb, a2_nb,"
    "  a1_type_m, a2_type_m, a1_type_f, a2_type_f,"
    "  value, sigma, count, value_1d, sigma_1d, count_1d)"
    " VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
  );

  db.exec("BEGIN");

  int n_files = 0, n_rows = 0;
  for (int file_num : list_table_files(bond_dir)) {
    std::string path = bond_dir + "/" + std::to_string(file_num) + ".table";
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) continue;
    ++n_files;

    char line[512];
    while (std::fgets(line, sizeof(line), f)) {
      if (is_skip_line_db(line)) continue;

      const char* p = line;
      const char* s;

      // ha1, ha2
      int ha1 = std::atoi(p);
      while (*p && *p != ' ' && *p != '\t') ++p;
      p = skip_blank_db(p);
      int ha2 = std::atoi(p);
      while (*p && *p != ' ' && *p != '\t') ++p;

      s = skip_blank_db(p); p = skip_word_db(s); std::string hybr_comb(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string in_ring(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a1_nb2(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a2_nb2(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a1_nb(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a2_nb(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string code1(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string code2(s, p - s);
      if (code2.empty()) continue;

      char* endp = nullptr;
      double value   = std::strtod(p, &endp); p = endp;
      double sigma   = std::strtod(p, &endp); p = endp;
      int    count   = std::strtol(p, &endp, 10); p = endp;
      double value2  = std::strtod(p, &endp); p = endp;
      double sigma2  = std::strtod(p, &endp); p = endp;
      int    count2  = std::strtol(p, &endp, 10); p = endp;

      auto it1 = codes.find(code1);
      auto it2 = codes.find(code2);
      std::string a1_type_f = it1 != codes.end() ? it1->second : std::string();
      std::string a2_type_f = it2 != codes.end() ? it2->second : std::string();
      std::string a1_type_m = prefix_before_db(a1_type_f, '{');
      std::string a2_type_m = prefix_before_db(a2_type_f, '{');

      sqlite3_reset(ins.st);
      sqlite3_bind_int   (ins.st, 1,  ha1);
      sqlite3_bind_int   (ins.st, 2,  ha2);
      sqlite3_bind_text  (ins.st, 3,  hybr_comb.c_str(), -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, 4,  in_ring.c_str(),   -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, 5,  a1_nb2.c_str(),    -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, 6,  a2_nb2.c_str(),    -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, 7,  a1_nb.c_str(),     -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, 8,  a2_nb.c_str(),     -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, 9,  a1_type_m.c_str(), -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, 10, a2_type_m.c_str(), -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, 11, a1_type_f.c_str(), -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, 12, a2_type_f.c_str(), -1, SQLITE_TRANSIENT);
      sqlite3_bind_double(ins.st, 13, value);
      sqlite3_bind_double(ins.st, 14, sigma);
      sqlite3_bind_int   (ins.st, 15, count);
      sqlite3_bind_double(ins.st, 16, value2);
      sqlite3_bind_double(ins.st, 17, sigma2);
      sqlite3_bind_int   (ins.st, 18, count2);
      if (sqlite3_step(ins.st) != SQLITE_DONE)
        throw std::runtime_error(std::string("acedrg-db: bond insert failed: ")
                                 + sqlite3_errmsg(db.db));
      ++n_rows;
    }
    std::fclose(f);
  }

  db.exec("COMMIT");
  std::fprintf(stderr, "    bond_entries: %d files, %d rows\n", n_files, n_rows);

  // Build indices last (faster than maintaining them during bulk insert).
  db.exec("CREATE INDEX idx_bond_hash ON bond_entries (ha1, ha2);");
  db.exec("CREATE INDEX idx_bond_full ON bond_entries"
          " (ha1, ha2, hybr_comb, in_ring,"
          "  a1_nb2, a2_nb2, a1_nb, a2_nb);");
}

// --- angle table converter -------------------------------------------------

void convert_angle_tables(SqliteDB& db, const std::string& angle_dir,
                          const std::unordered_map<std::string, std::string>& codes) {
  db.exec(
    "DROP TABLE IF EXISTS angle_entries;"
    "CREATE TABLE angle_entries ("
    "  ha1 INTEGER, ha2 INTEGER, ha3 INTEGER,"
    "  value_key TEXT,"
    "  a1_root TEXT, a2_root TEXT, a3_root TEXT,"
    "  a1_nb2 TEXT, a2_nb2 TEXT, a3_nb2 TEXT,"
    "  a1_nb  TEXT, a2_nb  TEXT, a3_nb  TEXT,"
    "  a1_type TEXT, a2_type TEXT, a3_type TEXT,"
    "  v1 REAL, s1 REAL, c1 INTEGER,"
    "  v2 REAL, s2 REAL, c2 INTEGER,"
    "  v3 REAL, s3 REAL, c3 INTEGER,"
    "  v4 REAL, s4 REAL, c4 INTEGER,"
    "  v5 REAL, s5 REAL, c5 INTEGER,"
    "  v6 REAL, s6 REAL, c6 INTEGER"
    ");"
  );

  SqliteStmt ins;
  ins.prepare(db.db,
    "INSERT INTO angle_entries VALUES ("
    "?,?,?,?,"             // ha1..3, value_key
    "?,?,?,?,?,?,?,?,?,"   // 3x (root, nb2, nb)
    "?,?,?,"               // types
    "?,?,?,?,?,?,?,?,?,"   // 3 levels of v/s/c
    "?,?,?,?,?,?,?,?,?)"   // 3 more levels
  );

  db.exec("BEGIN");

  int n_files = 0, n_rows = 0;
  for (int file_num : list_table_files(angle_dir)) {
    std::string path = angle_dir + "/" + std::to_string(file_num) + ".table";
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) continue;
    ++n_files;

    char line[1024];
    while (std::fgets(line, sizeof(line), f)) {
      if (is_skip_line_db(line)) continue;

      const char* p = line;
      const char* s;
      int ha1 = std::atoi(p); while (*p && *p != ' ' && *p != '\t') ++p; p = skip_blank_db(p);
      int ha2 = std::atoi(p); while (*p && *p != ' ' && *p != '\t') ++p; p = skip_blank_db(p);
      int ha3 = std::atoi(p); while (*p && *p != ' ' && *p != '\t') ++p;

      s = skip_blank_db(p); p = skip_word_db(s); std::string value_key(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a1_root(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a2_root(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a3_root(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a1_nb2(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a2_nb2(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a3_nb2(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a1_nb(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a2_nb(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string a3_nb(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string code1(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string code2(s, p - s);
      s = skip_blank_db(p); p = skip_word_db(s); std::string code3(s, p - s);
      if (code3.empty()) continue;

      double v[6]; double sg[6]; int c[6];
      char* endp = nullptr;
      for (int i = 0; i < 6; ++i) {
        v[i]  = std::strtod(p, &endp); p = endp;
        sg[i] = std::strtod(p, &endp); p = endp;
        c[i]  = std::strtol(p, &endp, 10); p = endp;
      }

      auto it1 = codes.find(code1);
      auto it2 = codes.find(code2);
      auto it3 = codes.find(code3);
      std::string a1_type = it1 != codes.end() ? prefix_before_db(it1->second, '{') : std::string();
      std::string a2_type = it2 != codes.end() ? prefix_before_db(it2->second, '{') : std::string();
      std::string a3_type = it3 != codes.end() ? prefix_before_db(it3->second, '{') : std::string();

      sqlite3_reset(ins.st);
      int b = 1;
      sqlite3_bind_int   (ins.st, b++, ha1);
      sqlite3_bind_int   (ins.st, b++, ha2);
      sqlite3_bind_int   (ins.st, b++, ha3);
      sqlite3_bind_text  (ins.st, b++, value_key.c_str(), -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a1_root.c_str(),   -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a2_root.c_str(),   -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a3_root.c_str(),   -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a1_nb2.c_str(),    -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a2_nb2.c_str(),    -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a3_nb2.c_str(),    -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a1_nb.c_str(),     -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a2_nb.c_str(),     -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a3_nb.c_str(),     -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a1_type.c_str(),   -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a2_type.c_str(),   -1, SQLITE_TRANSIENT);
      sqlite3_bind_text  (ins.st, b++, a3_type.c_str(),   -1, SQLITE_TRANSIENT);
      for (int i = 0; i < 6; ++i) {
        sqlite3_bind_double(ins.st, b++, v[i]);
        sqlite3_bind_double(ins.st, b++, sg[i]);
        sqlite3_bind_int   (ins.st, b++, c[i]);
      }
      if (sqlite3_step(ins.st) != SQLITE_DONE)
        throw std::runtime_error(std::string("acedrg-db: angle insert failed: ")
                                 + sqlite3_errmsg(db.db));
      ++n_rows;
    }
    std::fclose(f);
  }

  db.exec("COMMIT");
  std::fprintf(stderr, "    angle_entries: %d files, %d rows\n", n_files, n_rows);

  db.exec("CREATE INDEX idx_angle_hash ON angle_entries (ha1, ha2, ha3, value_key);");
}

}  // namespace

bool build_acedrg_sqlite(const std::string& tables_dir,
                         const std::string& sqlite_path) {
  // Erase any prior file so we don't append to old contents.
  std::remove(sqlite_path.c_str());

  SqliteDB db;
  db.open(sqlite_path,
          SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);

  // Pragmas to speed bulk insert.
  db.exec("PRAGMA journal_mode = OFF;");
  db.exec("PRAGMA synchronous  = OFF;");
  db.exec("PRAGMA temp_store   = MEMORY;");
  db.exec("PRAGMA cache_size   = -200000;");  // 200 MB page cache during build

  auto codes = load_atom_codes(tables_dir + "/allAtomTypesFromMolsCoded.list");
  if (codes.empty())
    throw std::runtime_error("acedrg-db: no atom-type codes found in " +
                             tables_dir + "/allAtomTypesFromMolsCoded.list");
  std::fprintf(stderr, "    atom-type codes: %zu entries\n", codes.size());

  convert_bond_tables(db, tables_dir + "/allOrgBondTables", codes);
  convert_angle_tables(db, tables_dir + "/allOrgAngleTables", codes);

  // Final analysis pass so the query planner has stats from the start.
  db.exec("ANALYZE;");

  return true;
}

#else  // !GEMMI_HAS_SQLITE3

bool build_acedrg_sqlite(const std::string&, const std::string&) {
  throw std::runtime_error(
      "gemmi was built without SQLite3 support; install libsqlite3-dev "
      "and reconfigure cmake.");
}

#endif

}  // namespace gemmi
