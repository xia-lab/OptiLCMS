#include <Rcpp.h>
#include "sqlite/sqlite3.h"
#include <string>
#include <vector>
#include <sstream>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix lipidomics_get_isotope_pattern_cpp(std::string db_path, int lipid_index) {
  if (db_path.empty() || lipid_index < 0) return NumericMatrix(0, 2);
  sqlite3* db = nullptr;
  int rc = sqlite3_open_v2(db_path.c_str(), &db, SQLITE_OPEN_READONLY, nullptr);
  if (rc != SQLITE_OK || db == nullptr) {
    if (db) sqlite3_close(db);
    return NumericMatrix(0, 2);
  }

  sqlite3_stmt* stmt = nullptr;
  std::vector<std::string> colnames;
  // get column names for isotope table via PRAGMA
  rc = sqlite3_prepare_v2(db, "PRAGMA table_info(isotope);", -1, &stmt, nullptr);
  if (rc == SQLITE_OK && stmt) {
    while (sqlite3_step(stmt) == SQLITE_ROW) {
      const unsigned char* name = sqlite3_column_text(stmt, 1); // name column
      if (name) colnames.emplace_back(reinterpret_cast<const char*>(name));
    }
    sqlite3_finalize(stmt);
  }

  if (colnames.size() < 2) {
    sqlite3_close(db);
    return NumericMatrix(0, 2);
  }

  std::string col_id = colnames[0];
  std::string col_val = colnames[1];
  std::string q = "SELECT \"" + col_val + "\" FROM isotope WHERE \"" + col_id + "\" = ? LIMIT 1;";

  rc = sqlite3_prepare_v2(db, q.c_str(), -1, &stmt, nullptr);
  if (rc != SQLITE_OK || stmt == nullptr) {
    if (stmt) sqlite3_finalize(stmt);
    sqlite3_close(db);
    return NumericMatrix(0, 2);
  }

  // bind as text to be safe
  sqlite3_bind_text(stmt, 1, std::to_string(lipid_index).c_str(), -1, SQLITE_TRANSIENT);
  NumericMatrix out(0, 2);
  if (sqlite3_step(stmt) == SQLITE_ROW) {
    const unsigned char* txt = sqlite3_column_text(stmt, 0);
    if (txt) {
      std::string s = reinterpret_cast<const char*>(txt);
      // parse lines of "mz intensity"
      std::istringstream iss(s);
      std::string line;
      std::vector<double> mzs;
      std::vector<double> ints;
      while (std::getline(iss, line)) {
        if (line.empty()) continue;
        std::istringstream ls(line);
        double a, b;
        if (ls >> a >> b) {
          if (std::isfinite(a) && std::isfinite(b)) {
            mzs.push_back(a);
            ints.push_back(b);
          }
        }
      }
      if (!mzs.empty()) {
        out = NumericMatrix(mzs.size(), 2);
        for (size_t i = 0; i < mzs.size(); ++i) {
          out(static_cast<int>(i), 0) = mzs[i];
          out(static_cast<int>(i), 1) = ints[i];
        }
      }
    }
  }
  sqlite3_finalize(stmt);
  sqlite3_close(db);
  return out;
}
