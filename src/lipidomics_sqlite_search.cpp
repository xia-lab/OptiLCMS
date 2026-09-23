#include <Rcpp.h>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include "sqlite/sqlite3.h"

using namespace Rcpp;

namespace {

std::string to_lower(const std::string& x) {
  std::string out = x;
  std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return out;
}

bool has_column(const std::unordered_map<std::string, int>& colmap, const std::string& key) {
  return colmap.find(to_lower(key)) != colmap.end();
}

std::string get_optional_text(sqlite3_stmt* stmt, const std::unordered_map<std::string, int>& colmap,
                              const std::string& key) {
  auto it = colmap.find(to_lower(key));
  if (it == colmap.end()) return "";
  const unsigned char* txt = sqlite3_column_text(stmt, it->second);
  return txt ? std::string(reinterpret_cast<const char*>(txt)) : std::string("");
}

double get_optional_double(sqlite3_stmt* stmt, const std::unordered_map<std::string, int>& colmap,
                           const std::string& key, double fallback = NA_REAL) {
  auto it = colmap.find(to_lower(key));
  if (it == colmap.end()) return fallback;
  if (sqlite3_column_type(stmt, it->second) == SQLITE_NULL) return fallback;
  return sqlite3_column_double(stmt, it->second);
}

int get_optional_int(sqlite3_stmt* stmt, const std::unordered_map<std::string, int>& colmap,
                     const std::string& key, int fallback = NA_INTEGER) {
  auto it = colmap.find(to_lower(key));
  if (it == colmap.end()) return fallback;
  if (sqlite3_column_type(stmt, it->second) == SQLITE_NULL) return fallback;
  return sqlite3_column_int(stmt, it->second);
}

NumericMatrix decode_peak_blob(sqlite3_stmt* stmt, int col_index) {
  if (sqlite3_column_type(stmt, col_index) == SQLITE_NULL) {
    return NumericMatrix(0, 2);
  }

  const void* blob = sqlite3_column_blob(stmt, col_index);
  const int bytes = sqlite3_column_bytes(stmt, col_index);
  if (blob == nullptr || bytes <= 0 || bytes % static_cast<int>(sizeof(double)) != 0) {
    return NumericMatrix(0, 2);
  }

  const std::size_t n_doubles = static_cast<std::size_t>(bytes / static_cast<int>(sizeof(double)));
  if (n_doubles == 0 || (n_doubles % 2) != 0) {
    return NumericMatrix(0, 2);
  }

  std::vector<double> values(n_doubles);
  std::memcpy(values.data(), blob, static_cast<std::size_t>(bytes));

  const std::size_t n_peaks = n_doubles / 2;
  std::vector<std::pair<double, double>> peaks(n_peaks);
  for (std::size_t i = 0; i < n_peaks; ++i) {
    peaks[i] = std::make_pair(values[i], values[i + n_peaks]);
  }

  std::sort(peaks.begin(), peaks.end(), [](const std::pair<double, double>& a,
                                           const std::pair<double, double>& b) {
    return a.first < b.first;
  });

  NumericMatrix spectrum(static_cast<int>(n_peaks), 2);
  for (std::size_t i = 0; i < n_peaks; ++i) {
    spectrum(static_cast<int>(i), 0) = peaks[i].first;
    spectrum(static_cast<int>(i), 1) = peaks[i].second;
  }

  return spectrum;
}

List search_lipidomics_precursor_blob_candidates(sqlite3* db,
                                                 const std::string& spectra_table,
                                                 double precursor_mz,
                                                 double ppm,
                                                 int max_candidates,
                                                 int max_per_table) {
  const int limit = std::min(max_candidates, max_per_table);
  const double mz_error = precursor_mz * ppm * 1e-6;
  const double min_mz = precursor_mz - mz_error;
  const double max_mz = precursor_mz + mz_error;

  std::string sql =
      "SELECT s.lipids_index, s.record_index, s.precursor_mz, s.adduct, s.DBID, s.DB_source, "
      "s.voltage, s.instrument, s.ref_spectra, s.spec_type, "
      "l.MolecularFormula, l.MolecularWeight, l.ExactMass, l.IUPACName, l.AbbreName, l.InChIKey, "
      "l.SMILES, l.InChI "
      "FROM \"" + spectra_table + "\" AS s "
      "LEFT JOIN lipids AS l ON s.lipids_index = l.lipids_index "
      "WHERE s.precursor_mz >= ?1 AND s.precursor_mz <= ?2 "
      "ORDER BY ABS(s.precursor_mz - ?3), s.record_index "
      "LIMIT ?4";

  sqlite3_stmt* row_stmt = nullptr;
  if (sqlite3_prepare_v2(db, sql.c_str(), -1, &row_stmt, nullptr) != SQLITE_OK) {
    std::string msg = sqlite3_errmsg(db);
    stop("Failed to prepare lipidomics blob query: " + msg);
  }

  sqlite3_bind_double(row_stmt, 1, min_mz);
  sqlite3_bind_double(row_stmt, 2, max_mz);
  sqlite3_bind_double(row_stmt, 3, precursor_mz);
  sqlite3_bind_int(row_stmt, 4, limit);

  const int ncols = sqlite3_column_count(row_stmt);
  std::unordered_map<std::string, int> row_colmap;
  for (int i = 0; i < ncols; ++i) {
    const char* cn = sqlite3_column_name(row_stmt, i);
    if (!cn) continue;
    row_colmap[to_lower(std::string(cn))] = i;
  }

  List candidates;
  while (sqlite3_step(row_stmt) == SQLITE_ROW) {
    NumericMatrix spectrum = decode_peak_blob(row_stmt, row_colmap[to_lower("ref_spectra")]);
    if (spectrum.nrow() == 0) {
      continue;
    }

    candidates.push_back(List::create(
      _["Table"] = spectra_table,
      _["LipidsIndex"] = get_optional_int(row_stmt, row_colmap, "lipids_index", NA_INTEGER),
      _["RecordIndex"] = get_optional_int(row_stmt, row_colmap, "record_index", NA_INTEGER),
      _["PrecursorMZ"] = get_optional_double(row_stmt, row_colmap, "precursor_mz", NA_REAL),
      _["Adduct"] = get_optional_text(row_stmt, row_colmap, "adduct"),
      _["DBID"] = get_optional_text(row_stmt, row_colmap, "dbid"),
      _["DBSource"] = get_optional_text(row_stmt, row_colmap, "db_source"),
      _["Voltage"] = get_optional_text(row_stmt, row_colmap, "voltage"),
      _["Instrument"] = get_optional_text(row_stmt, row_colmap, "instrument"),
      _["SpecType"] = get_optional_text(row_stmt, row_colmap, "spec_type"),
      _["MolecularFormula"] = get_optional_text(row_stmt, row_colmap, "molecularformula"),
      _["MolecularWeight"] = get_optional_double(row_stmt, row_colmap, "molecularweight", NA_REAL),
      _["ExactMass"] = get_optional_double(row_stmt, row_colmap, "exactmass", NA_REAL),
      _["IUPACName"] = get_optional_text(row_stmt, row_colmap, "iupacname"),
      _["AbbreName"] = get_optional_text(row_stmt, row_colmap, "abbrename"),
      _["InChIKey"] = get_optional_text(row_stmt, row_colmap, "inchikey"),
      _["SMILES"] = get_optional_text(row_stmt, row_colmap, "smiles"),
      _["InChI"] = get_optional_text(row_stmt, row_colmap, "inchi"),
      _["Spectrum"] = spectrum
    ));
  }

  sqlite3_finalize(row_stmt);
  return candidates;
}

} // namespace

// [[Rcpp::export]]
DataFrame lipidomics_db_search_precursor_cpp(
    std::string database_path,
    double precursor_mz,
    double ppm,
    int max_candidates = 300,
    int max_per_table = 100) {

  if (database_path.empty()) {
    stop("database_path is empty.");
  }
  if (precursor_mz <= 0) {
    stop("precursor_mz must be > 0.");
  }
  if (ppm <= 0) {
    stop("ppm must be > 0.");
  }
  if (max_candidates <= 0) {
    stop("max_candidates must be > 0.");
  }
  if (max_per_table <= 0) {
    stop("max_per_table must be > 0.");
  }

  sqlite3* db = nullptr;
  if (sqlite3_open(database_path.c_str(), &db) != SQLITE_OK) {
    std::string msg = sqlite3_errmsg(db);
    if (db) sqlite3_close(db);
    stop("Failed to open sqlite database: " + msg);
  }

  const double mz_error = precursor_mz * ppm * 1e-6;
  const double min_mz = precursor_mz - mz_error;
  const double max_mz = precursor_mz + mz_error;

  CharacterVector out_table;
  IntegerVector out_id;
  NumericVector out_precursor;
  CharacterVector out_name;
  CharacterVector out_inchikey;
  CharacterVector out_compound_class;
  CharacterVector out_comment;
  NumericVector out_rt;
  NumericVector out_ccs;
  CharacterVector out_ms2peaks;

  sqlite3_stmt* table_stmt = nullptr;
  const char* table_query =
      "SELECT name FROM sqlite_master "
      "WHERE type='table' AND name NOT LIKE 'sqlite_%' "
      "ORDER BY name";

  if (sqlite3_prepare_v2(db, table_query, -1, &table_stmt, nullptr) != SQLITE_OK) {
    std::string msg = sqlite3_errmsg(db);
    sqlite3_close(db);
    stop("Failed to list sqlite tables: " + msg);
  }

  int collected = 0;
  while (sqlite3_step(table_stmt) == SQLITE_ROW && collected < max_candidates) {
    const unsigned char* tname = sqlite3_column_text(table_stmt, 0);
    if (!tname) continue;
    std::string table_name(reinterpret_cast<const char*>(tname));

    std::string pragma_sql = "PRAGMA table_info(\"" + table_name + "\")";
    sqlite3_stmt* pragma_stmt = nullptr;
    if (sqlite3_prepare_v2(db, pragma_sql.c_str(), -1, &pragma_stmt, nullptr) != SQLITE_OK) {
      continue;
    }

    std::unordered_map<std::string, int> colmap;
    std::vector<std::string> colnames;
    while (sqlite3_step(pragma_stmt) == SQLITE_ROW) {
      const unsigned char* cn = sqlite3_column_text(pragma_stmt, 1);
      if (!cn) continue;
      std::string col(reinterpret_cast<const char*>(cn));
      colmap[to_lower(col)] = static_cast<int>(colnames.size());
      colnames.push_back(col);
    }
    sqlite3_finalize(pragma_stmt);

    std::string precursor_col;
    if (has_column(colmap, "PrecursorMZ")) precursor_col = "PrecursorMZ";
    else if (has_column(colmap, "PrecursorMz")) precursor_col = "PrecursorMz";
    else if (has_column(colmap, "precursor_mz")) precursor_col = "precursor_mz";

    std::string ms2_col;
    if (has_column(colmap, "MS2Peaks")) ms2_col = "MS2Peaks";
    else if (has_column(colmap, "ms2peaks")) ms2_col = "ms2peaks";
    else if (has_column(colmap, "MSMS")) ms2_col = "MSMS";

    if (precursor_col.empty() || ms2_col.empty()) {
      continue;
    }

    std::string sql = "SELECT * FROM \"" + table_name + "\" WHERE \"" + precursor_col +
      "\" >= ?1 AND \"" + precursor_col + "\" <= ?2 LIMIT " + std::to_string(max_per_table);

    sqlite3_stmt* row_stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &row_stmt, nullptr) != SQLITE_OK) {
      continue;
    }

    sqlite3_bind_double(row_stmt, 1, min_mz);
    sqlite3_bind_double(row_stmt, 2, max_mz);

    const int ncols = sqlite3_column_count(row_stmt);
    std::unordered_map<std::string, int> row_colmap;
    for (int i = 0; i < ncols; ++i) {
      const char* cn = sqlite3_column_name(row_stmt, i);
      if (!cn) continue;
      row_colmap[to_lower(std::string(cn))] = i;
    }

    while (sqlite3_step(row_stmt) == SQLITE_ROW && collected < max_candidates) {
      out_table.push_back(table_name);
      out_id.push_back(get_optional_int(row_stmt, row_colmap, "ID", collected + 1));
      out_precursor.push_back(get_optional_double(row_stmt, row_colmap, precursor_col, NA_REAL));

      std::string name = get_optional_text(row_stmt, row_colmap, "Name");
      if (name.empty()) name = get_optional_text(row_stmt, row_colmap, "CompoundName");
      if (name.empty()) name = get_optional_text(row_stmt, row_colmap, "Compounds");

      out_name.push_back(name);
      out_inchikey.push_back(get_optional_text(row_stmt, row_colmap, "InChIKey"));

      std::string cls = get_optional_text(row_stmt, row_colmap, "CompoundClass");
      if (cls.empty()) cls = get_optional_text(row_stmt, row_colmap, "Main_Class");
      out_compound_class.push_back(cls);

      out_comment.push_back(get_optional_text(row_stmt, row_colmap, "Comment"));

      double rt = get_optional_double(row_stmt, row_colmap, "RetentionTime", NA_REAL);
      if (R_IsNA(rt)) rt = get_optional_double(row_stmt, row_colmap, "RT", NA_REAL);
      out_rt.push_back(rt);

      double ccs = get_optional_double(row_stmt, row_colmap, "CollisionCrossSection", NA_REAL);
      if (R_IsNA(ccs)) ccs = get_optional_double(row_stmt, row_colmap, "CCS", NA_REAL);
      out_ccs.push_back(ccs);

      out_ms2peaks.push_back(get_optional_text(row_stmt, row_colmap, ms2_col));

      collected++;
    }

    sqlite3_finalize(row_stmt);
  }

  sqlite3_finalize(table_stmt);
  sqlite3_close(db);

  return DataFrame::create(
    _["Table"] = out_table,
    _["ID"] = out_id,
    _["PrecursorMZ"] = out_precursor,
    _["Name"] = out_name,
    _["InChIKey"] = out_inchikey,
    _["CompoundClass"] = out_compound_class,
    _["Comment"] = out_comment,
    _["RetentionTime"] = out_rt,
    _["CollisionCrossSection"] = out_ccs,
    _["MS2Peaks"] = out_ms2peaks,
    _["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
List lipidomics_db_search_precursors_blob_cpp(
    std::string database_path,
    std::string spectra_table,
    NumericVector precursor_mzs,
    double ppm,
    int max_candidates = 300,
    int max_per_table = 100) {

  if (database_path.empty()) {
    stop("database_path is empty.");
  }
  if (spectra_table != "spectra_pos" && spectra_table != "spectra_neg") {
    stop("spectra_table must be either 'spectra_pos' or 'spectra_neg'.");
  }
  if (ppm <= 0) {
    stop("ppm must be > 0.");
  }
  if (max_candidates <= 0) {
    stop("max_candidates must be > 0.");
  }
  if (max_per_table <= 0) {
    stop("max_per_table must be > 0.");
  }

  sqlite3* db = nullptr;
  if (sqlite3_open(database_path.c_str(), &db) != SQLITE_OK) {
    std::string msg = sqlite3_errmsg(db);
    if (db) sqlite3_close(db);
    stop("Failed to open sqlite database: " + msg);
  }

  List out(precursor_mzs.size());
  for (R_xlen_t i = 0; i < precursor_mzs.size(); ++i) {
    const double precursor_mz = precursor_mzs[i];
    if (!R_finite(precursor_mz) || precursor_mz <= 0) {
      out[i] = List::create();
      continue;
    }
    out[i] = search_lipidomics_precursor_blob_candidates(
      db, spectra_table, precursor_mz, ppm, max_candidates, max_per_table);
  }

  sqlite3_close(db);
  return out;
}
