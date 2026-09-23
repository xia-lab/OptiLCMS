#ifndef SQLITE_LIPID_UTILS_H
#define SQLITE_LIPID_UTILS_H

#include "sqlite/sqlite3.h"
#include <Rcpp.h>
#include <string>

using namespace Rcpp;

// Detect whether an opened sqlite DB follows the lipidomics schema
bool detect_lipidomics_schema(sqlite3* db);

// Decode a lipid BLOB (mz array then intensity array of doubles) to an Rcpp NumericMatrix (n x 2)
NumericMatrix decode_lipid_blob_to_matrix(const void* blob_ptr, int blob_size);

#endif
