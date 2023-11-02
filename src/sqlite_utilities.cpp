#include "sqlite_utilities.h"

SqliteDriver::SqliteDriver(String dbase, string db_tb, int ion_mode = 0){
  database = dbase;
  db_table = db_tb;
  sqlite3* db;
  sqlite3_stmt* stmt;
  const unsigned char* MS2peaks;
  string MS2Peak_str;
  
  vector<int> IDs_vec;
  vector<string> FMs_vec;
  vector<string> MS2Peaks_vec;
  vector<string> cmpds_vec;
  vector<double> precmz_vec;
  vector<double> rt_vec;

  SqliteDriver::setEntireDatabase(ion_mode);

}

int SqliteDriver::setDatabase (string db_path){
  database = db_path;
  return 1;
}

int SqliteDriver::setDB_table (string db_tb){
  // used to set database table. e.g. HMDB_experimental_PosDB
  db_table = db_tb;
  return 1;
}

int SqliteDriver::setEntireDatabase (int ion_mode){
  if(ion_mode == 0){
    all_DB  = {"HMDB_experimental_NegDB","HMDB_predicted_NegDB",
               "GNPS_NegDB","MoNA_NegDB","MINEs_NegDB","MSDIAL_NegDB",
               "MassBank_NegDB", "RIKEN_NegDB", "ReSpect_NegDB", 
               "VaniyaNP_NegDB", "LipidBlast_NegDB"};
    all_expDB = {"HMDB_experimental_NegDB",
                 "GNPS_NegDB",
                 "MSDIAL_NegDB", 
                 "LipidBlast_NegDB"};
  } else if(ion_mode == 1){
    all_DB = {"HMDB_experimental_PosDB","HMDB_predicted_PosDB",
              "GNPS_PosDB","MoNA_PosDB","MINEs_PosDB","MSDIAL_PosDB",
              "MassBank_PosDB", "RIKEN_PosDB", "ReSpect_PosDB", 
              "VaniyaNP_PosDB", "LipidBlast_PosDB", "BMDMS_PosDB"};
    all_expDB = {"HMDB_experimental_PosDB",
                 "GNPS_PosDB", 
                 "MSDIAL_PosDB", 
                 "LipidBlast_PosDB"};
  }
  return 1;
};

int SqliteDriver::create_connection(String database_Path){
  if (sqlite3_open (database_Path.get_cstring(), &db) != SQLITE_OK) {
    fprintf(stderr, "Error opening database.\n");
    return 0;
  }

  return 1;
}

int SqliteDriver::disconnectDB(){
  sqlite3_close(db);
  return 1;
}


int SqliteDriver::db_Executor_core(const char* Query_statement){
  //cout << "Now, the Query_statement is : " << Query_statement << endl;
  sqlite3_prepare(db, Query_statement, -1, &stmt, NULL);
  int step_res  = sqlite3_step (stmt);
  return step_res;
}

int SqliteDriver::extractMS2_with_ID(int ID){
  // compose statement
  string q = "SELECT MS2Peaks FROM " + db_table + " WHERE ID = " + std::to_string(ID);
  // run the core 
  int res = db_Executor_core(q.c_str());
  // if the result returned
  if(res == SQLITE_ROW){
    // data found successfully
    MS2peaks  = sqlite3_column_text(stmt, 0);
    MS2Peak_str = (char*) MS2peaks;
    //cout << "Now MS2peaks --> " << MS2peaks << endl;
    //cout << "Now MS2Peak_str --> " << MS2Peak_str << endl;
  } else {
    // something else
    cout << "Now something wierd happening [xx0123opas] --> " << res << endl;
  }
  // finalize the statement (Destructor of the stmt)
  sqlite3_finalize(stmt);
  return 1;
}

int SqliteDriver::extractMS2s_with_IDs(int IDs[]){
  // TODO: enable this function later
  return 1;
}

// This function is used to extract IDs from a certain database [data table]
vector<int> SqliteDriver::extractIDs_with_mzRange(double min_mz, double max_mz){
  // Initiate result vector
  vector<int> res_ID;
  // compose statement
  string q = "SELECT ID FROM " + db_table + 
    " WHERE PrecursorMZ > " + std::to_string(min_mz) + 
    " AND PrecursorMZ < " + std::to_string(max_mz) ;
  // run the core 
  sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
  bool done = false;
  int res, ID;
  while(!done){
    // if the result returned
    res  = sqlite3_step (stmt);
    if(res == SQLITE_ROW){
      // data found successfully
      ID  = sqlite3_column_int(stmt, 0);
      res_ID.push_back(ID);
      //cout << "Now ID --> " << ID << endl;
      //cout << "Now MS2Peak_str --> " << MS2Peak_str << endl;
    } else if(res == SQLITE_DONE) {
      // all searching finished
      done = true;
      break;
    } else {
      // something else
      cout << "Now something wierd happening [xx0135opas] --> " << res << endl;
    }
  }
  // finalize the statement (Destructor of the stmt)
  sqlite3_finalize(stmt);
  return res_ID;
}

// This function is used to extract two columns table [ID + MS2Peak]
int SqliteDriver::extractIDMS2_with_mzRange(double min_mz, double max_mz){
  // Initiate result vector
  vector<int> res_ID;
  vector<string> res_MS2Peaks;
  string q = "";
  q = q + "SELECT ID, MS2Peaks FROM " + db_table +
    " WHERE PrecursorMZ > " + std::to_string(min_mz) +
    " AND PrecursorMZ < " + std::to_string(max_mz);

  // run the core to get all results with a while loop
  bool done = false;
  int res, ID;
  string ms2peak;
  sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
  while(!done){
    res  = sqlite3_step (stmt);
    // if the result returned
    if(res == SQLITE_ROW){
      // data found successfully
      ID = sqlite3_column_int(stmt, 0);
      ms2peak = (char*) sqlite3_column_text(stmt, 1);
      res_ID.push_back(ID);
      res_MS2Peaks.push_back(ms2peak);
    } else if(res == SQLITE_DONE) {
      // all searching finished
      done = true;
      break;
    } else {
      // something else
      cout << "Now something wierd happening [xx0152opas] --> " << res << endl;
      return 0;
    }
  }
  // finalize the statement (Destructor of the stmt)
  sqlite3_finalize(stmt);
  IDs_vec = res_ID;
  MS2Peaks_vec = res_MS2Peaks;
  return 1;
}

// This function is used to extract two columns table [ID + MS2Peak] based on mz + rt
int SqliteDriver::extractIDMS2_with_mzrtRange(double min_mz, double max_mz, double min_rt, double max_rt){
  // Initiate result vector
  vector<int> res_ID;
  vector<string> res_MS2Peaks;
  string q = "";
  q = q + "SELECT ID, MS2Peaks FROM " + db_table +
    " WHERE PrecursorMZ > " + std::to_string(min_mz) +
    " AND PrecursorMZ < " + std::to_string(max_mz) + 
    " AND RetentionTime > " + std::to_string(min_rt) +
    " AND RetentionTime < " + std::to_string(max_rt);
  
  // run the core to get all results with a while loop
  bool done = false;
  int res, ID;
  string ms2peak;
  sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
  while(!done){
    res  = sqlite3_step (stmt);
    // if the result returned
    if(res == SQLITE_ROW){
      // data found successfully
      ID = sqlite3_column_int(stmt, 0);
      ms2peak = (char*) sqlite3_column_text(stmt, 1);
      res_ID.push_back(ID);
      res_MS2Peaks.push_back(ms2peak);
    } else if(res == SQLITE_DONE) {
      // all searching finished
      done = true;
      break;
    } else {
      // something else
      cout << "Now something wierd happening [xx0195opas] --> " << res << endl;
      return 0;
    }
  }
  // finalize the statement (Destructor of the stmt)
  sqlite3_finalize(stmt);
  IDs_vec = res_ID;
  MS2Peaks_vec = res_MS2Peaks;
  return 1;
}

// This function is used to extract IDs from the entire database tables
vector<int> SqliteDriver::extractIDs_with_mzRange_entireDB(double min_mz, double max_mz){
  // Initiate result vector
  vector<int> res_ID;
  string q = "";
  for(string db_str : all_DB){
    q = q + "SELECT ID FROM " + db_str +
      " WHERE PrecursorMZ > " + std::to_string(min_mz) + 
      " AND PrecursorMZ < " + std::to_string(max_mz);
    if(db_str != all_DB[all_DB.size()-1]){
      q = q + " union ";
    }
  }
  
  // run the core to get all results with a while loop
  bool done = false;
  int res, ID;
  sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
  //int res = db_Executor_core(q.c_str());
  while(!done){
    res  = sqlite3_step (stmt);
    // if the result returned
    if(res == SQLITE_ROW){
      // data found successfully
      ID  = sqlite3_column_int(stmt, 0);
      res_ID.push_back(ID);
    } else if(res == SQLITE_DONE) {
      // all searching finished
      done = true;
      break;
    } else {
      // something else
      cout << "Now something wierd happening [xx0146opas] --> " << res << endl;
    }
  }
  // finalize the statement (Destructor of the stmt)
  sqlite3_finalize(stmt);
  return res_ID;
}

// This function is used to extract two columns table [ID + MS2Peak]
int SqliteDriver::extractIDMS2_with_mzRange_entireDB(double min_mz, double max_mz){
  // Initiate result vector
  vector<int> res_ID;
  vector<string> res_MS2Peaks;
  string q = "";
  
  for(string db_str : all_DB){
    q = q + "SELECT ID, MS2Peaks FROM " + db_str +
      " WHERE PrecursorMZ > " + std::to_string(min_mz) +
      " AND PrecursorMZ < " + std::to_string(max_mz);
    if(db_str != all_DB[all_DB.size()-1]){
      q = q + " union ";
    }
  }

  // run the core to get all results with a while loop
  bool done = false;
  int res, ID;
  string ms2peak;
  sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
  while(!done){
    res  = sqlite3_step (stmt);
    // if the result returned
    if(res == SQLITE_ROW){
      // data found successfully
      ID = sqlite3_column_int(stmt, 0);
      ms2peak = (char*) sqlite3_column_text(stmt, 1);
      res_ID.push_back(ID);
      res_MS2Peaks.push_back(ms2peak);
    } else if(res == SQLITE_DONE) {
      // all searching finished
      done = true;
      break;
    } else {
      // something else
      cout << "Now something wierd happening [xx0289opas] --> " << res << endl;
      // cout << "q.c_str() ==> " << q.c_str() << endl;
      // cout << "&stmt ===> " << stmt << endl;
      return 0;
    }
  }
  // finalize the statement (Destructor of the stmt)
  sqlite3_finalize(stmt);
  IDs_vec = res_ID;
  MS2Peaks_vec = res_MS2Peaks;
  return 1;
}

// This function is used to extract two columns table [ID + MS2Peak] from experimental db only
int SqliteDriver::extractIDMS2_with_mzRange_expDB(double min_mz, double max_mz){
  // Initiate result vector
  vector<int> res_ID;
  vector<string> res_MS2Peaks;
  string q = "";
  
  for(string db_str : all_expDB){
    q = q + "SELECT ID, MS2Peaks FROM " + db_str +
      " WHERE PrecursorMZ > " + std::to_string(min_mz) +
      " AND PrecursorMZ < " + std::to_string(max_mz);
    if(db_str != all_expDB[all_expDB.size()-1]){
      q = q + " union ";
    }
  }
  
  // run the core to get all results with a while loop
  bool done = false;
  int res, ID;
  int count = 0;
  string ms2peak;
  sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
  while(!done){
    res  = sqlite3_step (stmt);
    // if the result returned
    if(count > 50) {// stop to avoid too many results
      done = true;
      break;
    }
    if(res == SQLITE_ROW){
      //cout << "count -- > " << count << endl;
      // data found successfully
      ID = sqlite3_column_int(stmt, 0);
      ms2peak = (char*) sqlite3_column_text(stmt, 1);
      res_ID.push_back(ID);
      if(ms2peak.length() < 4000){
        //res_MS2Peaks.push_back("");
        res_MS2Peaks.push_back(ms2peak);
      } else{
        res_MS2Peaks.push_back("");
      }
      count++;
    } else if(res == SQLITE_DONE) {
      // all searching finished
      done = true;
      break;
    } else {
      // something else
      cout << "Now something wierd happening [xx0289opas] --> " << res << endl;
      return 0;
    }
  }
  // finalize the statement (Destructor of the stmt)
  sqlite3_finalize(stmt);
  IDs_vec = res_ID;
  MS2Peaks_vec = res_MS2Peaks;
  return 1;
}

// This function is used to extract two columns table [Formula + MS2Peak]
int SqliteDriver::extractFMMS2_with_mzRange_entireDB(double min_mz, double max_mz){
  // Initiate result vector
  vector<int> res_ID;
  vector<string> res_FM;
  vector<string> res_MS2Peaks;
  vector<double> res_precMZ;
  vector<double> res_rts;
  string q = "";
  for(string db_str : all_DB){
    q = q + "SELECT ID, PrecursorMZ, Formula, RetentionTime, MS2Peaks FROM " + db_str +
      " WHERE PrecursorMZ > " + std::to_string(min_mz) +
      " AND PrecursorMZ < " + std::to_string(max_mz);
    if(db_str != all_DB[all_DB.size()-1]){
      q = q + " union ";
    }
  }
 
  // run the core to get all results with a while loop
  bool done = false;
  int res, ID;
  string ms2peak, FM;
  double prec_mz, rt_value;
  sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
  while(!done){
    res  = sqlite3_step (stmt);
    // if the result returned
    if(res == SQLITE_ROW){
      // data found successfully
      ID = sqlite3_column_int(stmt, 0);
      prec_mz = sqlite3_column_double(stmt, 1);
      FM = (char*) sqlite3_column_text(stmt, 2);
      rt_value = sqlite3_column_double(stmt, 3);
      ms2peak = (char*) sqlite3_column_text(stmt, 4);
      
      res_ID.push_back(ID);
      res_rts.push_back(rt_value);
      res_precMZ.push_back(prec_mz);
      res_FM.push_back(FM);
      res_MS2Peaks.push_back(ms2peak);
    } else if(res == SQLITE_DONE) {
      // all searching finished
      done = true;
      break;
    } else {
      // something else
      cout << "Now something wierd happening [xx0346opas] --> " << res << endl;
      return 0;
    }
  }
  // finalize the statement (Destructor of the stmt)
  sqlite3_finalize(stmt);
  IDs_vec = res_ID;
  FMs_vec = res_FM;
  MS2Peaks_vec = res_MS2Peaks;
  precmz_vec = res_precMZ;
  rt_vec = res_rts;
  return 1;
}

vector<CharacterVector> SqliteDriver::convertID2alls(IntegerVector IDs){
  vector<CharacterVector> allRes;
  CharacterVector inchikey_res(IDs.size());
  CharacterVector compound_res(IDs.size());
  CharacterVector formula_res(IDs.size());
  CharacterVector dbrecord_res(IDs.size());
  CharacterVector ms2peaks_res(IDs.size()); // add the ms2peaks
  
  int thisid;
  string q, q0, q1;
  String db_name;
  for(int i = 0; i<IDs.size(); i++){
    thisid = IDs[i];
    
    // find specific table
    q1 = "SELECT DB_Tables FROM Index_table WHERE Min_ID <= " + 
      std::to_string(thisid) + " AND Max_ID >= " + std::to_string(thisid) + ";";
    
    bool done = false;
    int res;
    string db_table;
    sqlite3_prepare(db, q1.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      // if the result returned
      if(res == SQLITE_ROW){
        // data found successfully
        db_table = (char*) sqlite3_column_text(stmt, 0);
      } else if(res == SQLITE_DONE) {
        // all searching finished
        done = true;
        break;
      } else {
        // something else
        cout << "Now something wierd happening [xx0395opas] --> " << res << endl;
        return allRes;
      }
    }
    sqlite3_finalize(stmt);
    
    // find corresponding inchikeys from the table
    q = "SELECT CompoundName,InchiKey,Formula,MS2Peaks FROM " + db_table + " WHERE ID == " + std::to_string(thisid); // Step 2: Update SQL query to get MS2Peaks
    
    db_name = db_table.substr(0,db_table.size()-6);
    
    done = false;
    
    string inchikey_txt;
    string compound_txt;
    string formula_txt;
    string ms2peaks_txt; // for MS2Peaks
    sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      // if the result returned
      if(res == SQLITE_ROW){
        // data found successfully
        dbrecord_res[i] = db_name;
        compound_txt = (char*) sqlite3_column_text(stmt, 0);
        compound_res[i] = compound_txt;
        inchikey_txt = (char*) sqlite3_column_text(stmt, 1);
        inchikey_res[i] = inchikey_txt;
        formula_txt = (char*) sqlite3_column_text(stmt, 2);
        formula_res[i] = formula_txt;
        ms2peaks_txt = (char*) sqlite3_column_text(stmt, 3); // Step 3: Extract MS2Peaks data
        ms2peaks_res[i] = ms2peaks_txt;
      } else if(res == SQLITE_DONE) {
        // all searching finished
        done = true;
        break;
      } else {
        // something else
        cout << "Now something wierd happening [xx0419opas] --> " << res << endl;
        allRes.push_back(compound_res);
        allRes.push_back(inchikey_res);
        allRes.push_back(formula_res);
        allRes.push_back(dbrecord_res);
        allRes.push_back(ms2peaks_res); // Step 4: Add ms2peaks_res to allRes
        sqlite3_finalize(stmt);
        return allRes;
      }
    }
    // finalize the statement (Destructor of the stmt)
    sqlite3_finalize(stmt);
    
  }
  allRes.push_back(compound_res);
  allRes.push_back(inchikey_res);
  allRes.push_back(formula_res);
  allRes.push_back(dbrecord_res);
  allRes.push_back(ms2peaks_res); // Step 4: Add ms2peaks_res to allRes
  return allRes;
};

bool SqliteDriver::clsTableExsiting(){
  string qur = "SELECT EXISTS (SELECT name FROM sqlite_master WHERE name='Lipids_classification');";
  
  sqlite3_prepare(db, qur.c_str(), -1, &stmt, NULL);
  bool done = false;
  int res, ID=0;
  while(!done){
    // if the result returned
    res  = sqlite3_step (stmt);
    if(res == SQLITE_ROW){
      // data found successfully
      ID  = sqlite3_column_int(stmt, 0);
      if(ID == 1){
        return true;
      } else if (ID == 0){
        return false;
      }
    } else if(res == SQLITE_DONE) {
      // all searching finished
      done = true;
      break;
    } else {
      // something else
      cout << "Now something wierd happening [xx0135opas] --> " << res << endl;
    }
  }
  
  if(ID == 1){
    return true;
  } else if (ID == 0){
    return false;
  }
  return false;
}

vector<CharacterVector> SqliteDriver::extractClasses(IntegerVector IDs){
  vector<CharacterVector> allRes;
  CharacterVector supcls_res(IDs.size());
  CharacterVector maincls_res(IDs.size());
  CharacterVector subcls_res(IDs.size());
  
  int thisid;
  string q, q0, q1;
  for(int i = 0; i<IDs.size(); i++){
    thisid = IDs[i];
    bool done = false;
    int res;

    // find corresponding inchikeys from the table
    q = "SELECT superclass,mainclass,subclass FROM Lipids_classification WHERE ID == " + std::to_string(thisid);

    string supcls_txt = "Undefined";
    string maincls_txt = "Undefined";
    string subcls_txt = "Undefined";
    supcls_res[i] = supcls_txt;
    maincls_res[i] = maincls_txt;
    subcls_res[i] = subcls_txt;
    sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      // if the result returned
      if(res == SQLITE_ROW){
        // data found successfully
        supcls_txt = (char*) sqlite3_column_text(stmt, 0);
        supcls_res[i] = supcls_txt;
        maincls_txt = (char*) sqlite3_column_text(stmt, 1);
        maincls_res[i] = maincls_txt;
        subcls_txt = (char*) sqlite3_column_text(stmt, 2);
        subcls_res[i] = subcls_txt;
      } else if(res == SQLITE_DONE) {
        // all searching finished
        done = true;
        break;
      } else {
        // something else
        cout << "Now something wierd happening [xx0600008] --> " << res << endl;
        allRes.push_back(supcls_res);
        allRes.push_back(maincls_res);
        allRes.push_back(subcls_res);
        sqlite3_finalize(stmt);
        return allRes;
      }
    }
    // finalize the statement (Destructor of the stmt)
    sqlite3_finalize(stmt);
    
  }
  allRes.push_back(supcls_res);
  allRes.push_back(maincls_res);
  allRes.push_back(subcls_res);
  return allRes;
};

CharacterVector SqliteDriver::convertID2InChiKeys(IntegerVector IDs){
  CharacterVector inchikey_res(IDs.size());

  int thisid;
  string q, q0, q1;
  for(int i = 0; i<IDs.size(); i++){
    thisid = IDs[i];
    
    // find specific table
    q1 = "SELECT DB_Tables FROM Index_table WHERE Min_ID <= " + 
      std::to_string(thisid) + " AND Max_ID >= " + std::to_string(thisid) + ";";
    
    bool done = false;
    int res;
    string db_table;
    sqlite3_prepare(db, q1.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      // if the result returned
      if(res == SQLITE_ROW){
        // data found successfully
        db_table = (char*) sqlite3_column_text(stmt, 0);
      } else if(res == SQLITE_DONE) {
        // all searching finished
        done = true;
        break;
      } else {
        // something else
        cout << "Now something wierd happening [xx0395opas] --> " << res << endl;
        return inchikey_res;
      }
    }
    sqlite3_finalize(stmt);
    
    // find corresponding inchikeys from the table
    q = "SELECT InchiKey FROM " + db_table + " WHERE ID == " + std::to_string(thisid);
    
    done = false;

    string inchikey_txt;
    sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      // if the result returned
      if(res == SQLITE_ROW){
        // data found successfully
        inchikey_txt = (char*) sqlite3_column_text(stmt, 0);
        inchikey_res[i] = inchikey_txt;
      } else if(res == SQLITE_DONE) {
        // all searching finished
        done = true;
        break;
      } else {
        // something else
        cout << "Now something wierd happening [xx0419opas] --> " << res << endl;
        return inchikey_res;
      }
    }
    // finalize the statement (Destructor of the stmt)
    sqlite3_finalize(stmt);

  }
  return inchikey_res;
};

CharacterVector SqliteDriver::convertID2CMPDNMs(IntegerVector IDs){
  CharacterVector inchikey_res(IDs.size());
  
  int thisid;
  string q, q0, q1;
  for(int i = 0; i<IDs.size(); i++){
    thisid = IDs[i];
    //cout << "this id => " << thisid << endl;
    // find specific table
    q1 = "SELECT DB_Tables FROM Index_table WHERE Min_ID <= " + 
      std::to_string(thisid) + " AND Max_ID >= " + std::to_string(thisid) + ";";
    
    bool done = false;
    int res;
    string db_table;
    sqlite3_prepare(db, q1.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      // if the result returned
      if(res == SQLITE_ROW){
        // data found successfully
        db_table = (char*) sqlite3_column_text(stmt, 0);
      } else if(res == SQLITE_DONE) {
        // all searching finished
        done = true;
        break;
      } else {
        // something else
        cout << "Now something wierd happening [xx0395opas] --> " << res << endl;
        return inchikey_res;
      }
    }
    sqlite3_finalize(stmt);
    
    // find corresponding inchikeys from the table
    q = "SELECT CompoundName FROM " + db_table + " WHERE ID == " + std::to_string(thisid);
    
    done = false;
    
    string inchikey_txt;
    sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      // if the result returned
      if(res == SQLITE_ROW){
        // data found successfully
        inchikey_txt = (char*) sqlite3_column_text(stmt, 0);
        //cout << "inchikey_txt --> " << inchikey_txt << endl;
        inchikey_res[i] = inchikey_txt;
      } else if(res == SQLITE_DONE) {
        // all searching finished
        done = true;
        break;
      } else {
        // something else
        cout << "Now something wierd happening [xx0419opas] --> " << res << endl;
        return inchikey_res;
      }
    }
    // finalize the statement (Destructor of the stmt)
    sqlite3_finalize(stmt);
    
  }
  return inchikey_res;
};

CharacterVector SqliteDriver::convertID2Formulas(IntegerVector IDs){
  CharacterVector inchikey_res(IDs.size());
  
  int thisid;
  string q, q0, q1;
  for(int i = 0; i<IDs.size(); i++){
    thisid = IDs[i];
    //cout << "this id => " << thisid << endl;
    // find specific table
    q1 = "SELECT DB_Tables FROM Index_table WHERE Min_ID <= " + 
      std::to_string(thisid) + " AND Max_ID >= " + std::to_string(thisid) + ";";
    
    bool done = false;
    int res;
    string db_table;
    sqlite3_prepare(db, q1.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      // if the result returned
      if(res == SQLITE_ROW){
        // data found successfully
        db_table = (char*) sqlite3_column_text(stmt, 0);
      } else if(res == SQLITE_DONE) {
        // all searching finished
        done = true;
        break;
      } else {
        // something else
        cout << "Now something wierd happening [xx0395opas] --> " << res << endl;
        return inchikey_res;
      }
    }
    sqlite3_finalize(stmt);
    
    // find corresponding inchikeys from the table
    q = "SELECT Formula FROM " + db_table + " WHERE ID == " + std::to_string(thisid);
    
    done = false;
    
    string inchikey_txt;
    sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      // if the result returned
      if(res == SQLITE_ROW){
        // data found successfully
        inchikey_txt = (char*) sqlite3_column_text(stmt, 0);
        //cout << "inchikey_txt --> " << inchikey_txt << endl;
        inchikey_res[i] = inchikey_txt;
      } else if(res == SQLITE_DONE) {
        // all searching finished
        done = true;
        break;
      } else {
        // something else
        cout << "Now something wierd happening [xx0419opas] --> " << res << endl;
        return inchikey_res;
      }
    }
    // finalize the statement (Destructor of the stmt)
    sqlite3_finalize(stmt);
    
  }
  return inchikey_res;
};


string SqliteDriver::getMS2Peaks(){
  return MS2Peak_str;
}

vector<string> SqliteDriver::getFMs(){
  return FMs_vec;
}

vector<int> SqliteDriver::getIDsVec(){
  return IDs_vec;
}

vector<string> SqliteDriver::getMS2PeaksVec(){
  return MS2Peaks_vec;
}

vector<string> SqliteDriver::getCMPDsVec(){
  return cmpds_vec;
}

vector<double> SqliteDriver::getPrecMZVec(){
  return precmz_vec;
}

vector<double> SqliteDriver::getRTVec(){
  return rt_vec;
}

CharacterVector SqliteDriver::convertID2MS2Peaks(IntegerVector IDs){
  CharacterVector ms2peaks_res(IDs.size());
  
  int thisid;
  string q, q0, q1;
  for(int i = 0; i<IDs.size(); i++){
    thisid = IDs[i];
    
    // find specific table
    q1 = "SELECT DB_Tables FROM Index_table WHERE Min_ID <= " + 
      std::to_string(thisid) + " AND Max_ID >= " + std::to_string(thisid) + ";";
    
    bool done = false;
    int res;
    string db_table;
    sqlite3_prepare(db, q1.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      if(res == SQLITE_ROW){
        db_table = (char*) sqlite3_column_text(stmt, 0);
      } else if(res == SQLITE_DONE) {
        done = true;
        break;
      } else {
        cout << "Now something weird happening [xx0395opas] --> " << res << endl;
        return ms2peaks_res;
      }
    }
    sqlite3_finalize(stmt);
    
    // find corresponding MS2Peaks from the table
    q = "SELECT MS2Peaks FROM " + db_table + " WHERE ID == " + std::to_string(thisid);
    
    done = false;
    
    string ms2peaks_txt;
    sqlite3_prepare(db, q.c_str(), -1, &stmt, NULL);
    while(!done){
      res  = sqlite3_step (stmt);
      if(res == SQLITE_ROW){
        ms2peaks_txt = (char*) sqlite3_column_text(stmt, 0);
        ms2peaks_res[i] = ms2peaks_txt;
      } else if(res == SQLITE_DONE) {
        done = true;
        break;
      } else {
        cout << "Now something weird happening [xx0419opas] --> " << res << endl;
        return ms2peaks_res;
      }
    }
    sqlite3_finalize(stmt);
  }
  return ms2peaks_res;
};

// // [[Rcpp::export]]
// double test_SQLite_fun(int x) {
// 
// 
//   SqliteDriver sd;
//   sd.setDatabase("/data/COMPOUND_DBs/Curated_DB/MS2ID.sqlite");
//   sd.create_connection("/data/COMPOUND_DBs/Curated_DB/MS2ID.sqlite");
//   sd.setDB_table("HMDB_experimental_PosDB");
//   // // test-1
//   for(int i=26570; i< 26580;i++){
//     sd.extractMS2_with_ID(i);
//   }
//   string res = sd.getMS2Peaks();
// 
//   // test-2
//   for(int i=0; i <10;i++){
//     //vector<int> vec_int = sd.extractIDs_with_mzRange_entireDB(479+i, 480+i);
//     vector<int> vec_int = sd.extractIDs_with_mzRange(479+i, 480+i);
//     //cout << endl;
//     //for(int j : vec_int){
//     //  cout << "IDs --> " <<  j << endl;
//     //}
//   }
// 
//   // test-3
//   for(int i=0; i <10;i++){
//     int oo = sd.extractIDMS2_with_mzRange_entireDB(479+i*0.05,479.05+i*0.05);
//     //cout << endl;
//     vector<int> allIDs = sd.getIDsVec();
//     vector<string> allMSPeaks = sd.getMS2PeaksVec();
//     // cout << "size of allMSPeaks --> " << allMSPeaks.size() << endl;
//     // cout << "size of allIDs --> " << allIDs.size() << endl;
//     // for(int j : allIDs){
//     //   cout << "IDs --> " <<  j << endl;
//     // }
//     // for(string s : allMSPeaks){
//     //   cout << "allMSPeaks --> " <<  s << endl;
//     // }
//   }
//   
//   
//   
//   //cout << "Now the ms2 peak res is " << res << endl;
//   // cout << "q is --> " << q.get_cstring() << endl;
//   // cout << "sizeof q--> " << sizeof q << endl;
//   // sqlite3_prepare(db, q.get_cstring(), -1, &stmt, NULL);
//   // 
//   // cout << "stmt --> " << &stmt << endl;
//   // 
//   // bool done = false;
//   // while (!done) {
//   //   printf("In select while\n");
//   //   // cout << "SQLITE_ROW --> " << SQLITE_ROW << endl;
//   //   // cout << "SQLITE_DONE --> " << SQLITE_DONE << endl;
//   //   int a1  = sqlite3_step (stmt);
//   //   cout << "a1 --> " << a1 << endl;
//   //   
//   //   switch (a1) {
//   //   case SQLITE_ROW:
//   //     bytes = sqlite3_column_bytes(stmt, 0);
//   //     text  = sqlite3_column_text(stmt, 0);
//   //     printf ("count %d: %s (%d bytes)\n", row, text, bytes);
//   //     row++;
//   //     break;
//   //     
//   //   case SQLITE_DONE:
//   //     done = true;
//   //     break;
//   //     
//   //   default:
//   //     fprintf(stderr, "Failed.\n");
//   //   return 1;
//   //   }
//   // }
//   // 
//   // sqlite3_finalize(stmt);
//   // 
//   sd.disconnectDB();
//   return 0;
// }
// 
// 
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically 
// // run after the compilation.
// //
// 
// /*** R
// system.time(
// a1 <- test_SQLite_fun(42)
// )
// */
