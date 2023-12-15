#ifndef SQLITEU_H
#define SQLITEU_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <string.h>
#include "sqlite/sqlite3.h"

using namespace Rcpp;
using namespace std;

class SqliteDriver {
  
private:
  string database;
  string db_table;
  sqlite3* db;
  sqlite3_stmt* stmt;
  const unsigned char* MS2peaks;
  string MS2Peak_str;
  
  // These vectors are used to store the information of MS/MS records [they are paired]
  vector<int> IDs_vec;
  vector<string> FMs_vec;
  vector<string> MS2Peaks_vec;
  vector<string> cmpds_vec;
  vector<string> smiles_vec;
  vector<string> inchikeys_vec;
  vector<double> precmz_vec;
  vector<double> rt_vec;

  vector<string> all_DB;
  vector<string> all_expDB;

public:
  
  SqliteDriver(String dbase, string db_tb, int ion_mode);
  
  int setDatabase (string db_path);
  
  int setDB_table (string db_tb);
  
  int setEntireDatabase (int ion_mode);
  
  int create_connection(String database_Path);
  
  int db_Executor_core(const char* Query_statement);
  
  int disconnectDB();
  
  int extractMS2_with_ID(int ID);
  
  int extractMS2s_with_IDs(int IDs[]);
 
  vector<int> extractIDs_with_mzRange(double min_mz, double max_mz);
 
  int extractIDMS2_with_mzRange(double min_mz, double max_mz);
  
  int extractIDMS2_with_mzrtRange(double min_mz, double max_mz, double min_rt, double max_rt);
  
  vector<int> extractIDs_with_mzRange_entireDB(double min_mz, double max_mz);
  
  int extractIDMS2_with_mzRange_entireDB(double min_mz, double max_mz);
  
  int extractIDMS2_with_mzRange_expDB(double min_mz, double max_mz);
  
  int extractFMMS2_with_mzRange_entireDB(double min_mz, double max_mz);
  
  int extractALLMS2_with_mzRange(double min_mz, double max_mz);
  
  CharacterVector convertID2InChiKeys(IntegerVector IDs);
  
  CharacterVector convertID2CMPDNMs(IntegerVector IDs);
  
  CharacterVector convertID2Formulas(IntegerVector IDs);
  
  vector<CharacterVector> convertID2alls(IntegerVector IDs);

  CharacterVector convertID2MS2Peaks(IntegerVector IDs);
  
  vector<CharacterVector> extractClasses(IntegerVector IDs);
  
  bool clsTableExsiting();
  
  string getMS2Peaks();
  
  vector<int> getIDsVec();
  
  vector<string> getMS2PeaksVec();
  
  vector<string> getFMs();
  
  vector<double> getPrecMZVec();
  
  vector<double> getRTVec();
  
  vector<string> getCMPDsVec();
  
  vector<string> getSimlesVec();
  
  vector<string> getInchikeysVec();
  
};

#endif
