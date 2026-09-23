#include "sqlite_lipid_utilities.h"
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <sstream>

bool detect_lipidomics_schema(sqlite3* db) {
  if (db == nullptr) return false;
  const char* q =
    "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name IN ('spectra_pos','spectra_neg','lipids');";
  sqlite3_stmt* local_stmt = nullptr;
  if (sqlite3_prepare_v2(db, q, -1, &local_stmt, nullptr) != SQLITE_OK) {
    if (local_stmt != nullptr) sqlite3_finalize(local_stmt);
    return false;
  }

  bool is_lipid = false;
  if (sqlite3_step(local_stmt) == SQLITE_ROW) {
    const int matched = sqlite3_column_int(local_stmt, 0);
    is_lipid = matched >= 2;
  }
  sqlite3_finalize(local_stmt);
  return is_lipid;
}

NumericMatrix decode_lipid_blob_to_matrix(const void* blob_ptr, int blob_size) {
  if (blob_ptr == nullptr || blob_size <= 0 || (blob_size % static_cast<int>(sizeof(double))) != 0) {
    return NumericMatrix(0, 2);
  }

  const int n_doubles = blob_size / static_cast<int>(sizeof(double));
  if (n_doubles < 2 || (n_doubles % 2) != 0) {
    return NumericMatrix(0, 2);
  }

  std::vector<double> values(n_doubles);
  std::memcpy(values.data(), blob_ptr, static_cast<size_t>(blob_size));

  const int n_pairs = n_doubles / 2;
  std::vector<std::pair<double, double>> peaks;
  peaks.reserve(static_cast<size_t>(n_pairs));

  for (int i = 0; i < n_pairs; ++i) {
    const double mz = values[i];
    const double intensity = values[i + n_pairs];
    if (std::isfinite(mz) && std::isfinite(intensity) && mz > 0.0 && intensity > 0.0) {
      peaks.emplace_back(mz, intensity);
    }
  }

  if (peaks.empty()) return NumericMatrix(0, 2);

  std::sort(peaks.begin(), peaks.end(), [](const std::pair<double,double>& a, const std::pair<double,double>& b){
    return a.first < b.first;
  });

  NumericMatrix spec(static_cast<int>(peaks.size()), 2);
  for (size_t i = 0; i < peaks.size(); ++i) {
    spec(static_cast<int>(i), 0) = peaks[i].first;
    spec(static_cast<int>(i), 1) = peaks[i].second;
  }
  return spec;
}
