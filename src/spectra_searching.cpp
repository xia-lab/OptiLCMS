#include <RcppArmadillo.h>
#include "dda_utilities.h"
#include "sqlite_utilities.h"
#include "utilities.h"
#include "rules.h"

using namespace Rcpp;
using namespace std;

List neutralLoss_matching_Live(SqliteDriver &SQLiteObj, 
                          vector<string> rules, vector<int> rule_dirs, vector<double> rule_ms,
                          double ppm2, double prec_exp,
                          NumericVector exp_mzs, NumericVector exp_ints){
  List res;
  
  // predict neutral loss of experimental spectrum first
  NumericVector nl_exp_mzs;
  NumericVector nl_exp_ints;
  for(int i=0; i<exp_mzs.size(); i++){
    if((prec_exp - exp_mzs[i]) > 2){ // at least H2 is needed to be considered as a "neutral loss"
      nl_exp_mzs.push_back(prec_exp - exp_mzs[i]);
      nl_exp_ints.push_back(exp_ints[i]);
    }
  }
  
  // cout << "nl_exp_mzsv -> " << nl_exp_mzs << endl;
  vector<int> IDs;
  vector<double> scores;
  vector<string> rulesVec;
  
  if(nl_exp_mzs.size() == 0){
    res = List::create(Named("IDs") = IDs, 
                       Named("Scores") = scores,
                       Named("rules") = rulesVec);
    return res;
  }
  
  // extract and predict reference from database
  // propagate with rules
  int tmp_dir, x;
  double diff_ms, tmp_mz, mass_error;
  for(int r=0; r<rule_dirs.size(); r++){
    tmp_dir = rule_dirs[r];
    diff_ms = rule_ms[r];
    if(tmp_dir == 0){
      // + 1
      tmp_mz = prec_exp + diff_ms;
      mass_error = tmp_mz*ppm2*1e-6;
      x = SQLiteObj.extractIDMS2_with_mzRange_entireDB(tmp_mz - mass_error, tmp_mz + mass_error);
      if(x == 0){
        warning("Database searching failed for m/z = " + std::to_string(tmp_mz) + " !\n");
      }
      vector<string> allMS2refs = SQLiteObj.getMS2PeaksVec();
      vector<int> allIDs = SQLiteObj.getIDsVec();
      if(allMS2refs.size() > 0){
        double best_score = 0.2, tmp_score;
        int best_id = -1;
        //for(string ms2ref : allMS2refs){
        for(int s=0; s<allMS2refs.size(); s++){
          string ms2ref = allMS2refs[s];
          NumericMatrix ref_spec_msms = ms2peak_parse(ms2ref);
          NumericVector ref_mzs = ref_spec_msms(_,0);
          NumericVector ref_ints = ref_spec_msms(_,1);
          // calculate neutral loss
          NumericVector nl_ref_mzs;
          NumericVector nl_ref_ints;
          for(int i=0; i<ref_mzs.size(); i++){
            if((tmp_mz - ref_mzs[i]) > 2){ // at least H2 is needed to be considered as a "neutral loss"
              nl_ref_mzs.push_back(tmp_mz - ref_mzs[i]);
              nl_ref_ints.push_back(ref_ints[i]);
            }
          }
          if(nl_ref_mzs.size() > 1){
            tmp_score = neutral_loss_similarity(nl_exp_mzs, nl_exp_ints, nl_ref_mzs, nl_ref_ints, ppm2);
            if(tmp_score > best_score){
              best_score = tmp_score;
              best_id = allIDs[s];
            }
          }
        }
        if(best_id != -1){
          // found a useful candidate to return
          rulesVec.push_back(rules[r]);
          scores.push_back(best_score);
          IDs.push_back(best_id);
        }
      }

      // -1
      tmp_mz = prec_exp - diff_ms;
      mass_error = tmp_mz*ppm2*1e-6;
      x = SQLiteObj.extractIDMS2_with_mzRange_entireDB(tmp_mz - mass_error, tmp_mz + mass_error);
      if(x == 0){
        warning("Database searching failed for m/z = " + std::to_string(tmp_mz) + " !\n");
      }
      allMS2refs = SQLiteObj.getMS2PeaksVec();
      allIDs = SQLiteObj.getIDsVec();
      if(allMS2refs.size() > 0){
        double best_score = 0.2, tmp_score; // minimum score is 0.2
        int best_id = -1;
        //for(string ms2ref : allMS2refs){
        for(int s=0; s<allMS2refs.size(); s++){
          string ms2ref = allMS2refs[s];
          NumericMatrix ref_spec_msms = ms2peak_parse(ms2ref);
          NumericVector ref_mzs = ref_spec_msms(_,0);
          NumericVector ref_ints = ref_spec_msms(_,1);
          // calculate neutral loss
          NumericVector nl_ref_mzs;
          NumericVector nl_ref_ints;
          for(int i=0; i<ref_mzs.size(); i++){
            if((tmp_mz - ref_mzs[i]) > 2){ // at least H2 is needed to be considered as a "neutral loss"
              nl_ref_mzs.push_back(tmp_mz - ref_mzs[i]);
              nl_ref_ints.push_back(ref_ints[i]);
            }
          }
          if(nl_ref_mzs.size() > 1){
            tmp_score = neutral_loss_similarity(nl_exp_mzs, nl_exp_ints, nl_ref_mzs, nl_ref_ints, ppm2);
            if(tmp_score > best_score){
              best_score = tmp_score;
              best_id = allIDs[s];
            }
          }
        }
        if(best_id != -1){
          // found a useful candidate to return
          rulesVec.push_back(rules[r]);
          scores.push_back(best_score);
          IDs.push_back(best_id);
        }
      }
    
    } else {
      
      tmp_mz = prec_exp + diff_ms*tmp_dir;
      mass_error = tmp_mz*ppm2*1e-6;
      x = SQLiteObj.extractIDMS2_with_mzRange_entireDB(tmp_mz - mass_error, tmp_mz + mass_error);
      if(x == 0){
        warning("Database searching failed for m/z = " + std::to_string(tmp_mz) + " !\n");
      }
      vector<string> allMS2refs = SQLiteObj.getMS2PeaksVec();
      vector<int> allIDs = SQLiteObj.getIDsVec();
      if(allMS2refs.size() > 0){
        double best_score = 0.2, tmp_score;
        int best_id = -1;
        //for(string ms2ref : allMS2refs){
        for(int s=0; s<allMS2refs.size(); s++){
          string ms2ref = allMS2refs[s];
          NumericMatrix ref_spec_msms = ms2peak_parse(ms2ref);
          NumericVector ref_mzs = ref_spec_msms(_,0);
          NumericVector ref_ints = ref_spec_msms(_,1);
          // calculate neutral loss
          NumericVector nl_ref_mzs;
          NumericVector nl_ref_ints;
          for(int i=0; i<ref_mzs.size(); i++){
            if((tmp_mz - ref_mzs[i]) > 2){ // at least H2 is needed to be considered as a "neutral loss"
              nl_ref_mzs.push_back(tmp_mz - ref_mzs[i]);
              nl_ref_ints.push_back(ref_ints[i]);
            }
          }
          if(nl_ref_mzs.size() > 1){
            tmp_score = neutral_loss_similarity(nl_exp_mzs, nl_exp_ints, nl_ref_mzs, nl_ref_ints, ppm2);
            if(tmp_score > best_score){
              best_score = tmp_score;
              best_id = allIDs[s];
            }
          }
        }
        if(best_id != -1){
          // found a useful candidate to return
          rulesVec.push_back(rules[r]);
          scores.push_back(best_score);
          IDs.push_back(best_id);
        }
      }
      
    }
    
  }
  
  
  // Organize results a [List] with three vectors
  // <1>. vector of ID <- database int | vector<int>;
  // <2>. vector of score <- matching score (neutral loss matching : dot product) | vector<double>;
  // <3>. vector of rules <- rules name | vector<string>;
  res = List::create(Named("IDs") = IDs, 
                     Named("Scores") = scores,
                     Named("rules") = rulesVec);
  return res;
}



List neutralLoss_matching(SqliteDriver &SQLiteObjNL, 
                          vector<string> rules, vector<int> rule_dirs, vector<double> rule_ms,
                          double ppm2, double prec_exp,
                          NumericVector exp_mzs, NumericVector exp_ints){
  List res;
  
  // predict neutral loss of experimental spectrum first
  NumericVector nl_exp_mzs;
  NumericVector nl_exp_ints;
  for(int i=0; i<exp_mzs.size(); i++){
    if((prec_exp - exp_mzs[i]) > 2){ // at least H2 is needed to be considered as a "neutral loss"
      nl_exp_mzs.push_back(prec_exp - exp_mzs[i]);
      nl_exp_ints.push_back(exp_ints[i]);
    }
  }
  
  // cout << "nl_exp_mzsv -> " << nl_exp_mzs << endl;
  vector<int> IDs;
  vector<double> scores;
  vector<string> rulesVec;
  
  if(nl_exp_mzs.size() == 0){
    res = List::create(Named("IDs") = IDs, 
                       Named("Scores") = scores,
                       Named("rules") = rulesVec);
    return res;
  }
  
  // extract and predict reference from database
  // propagate with rules
  int tmp_dir, x;
  double diff_ms, tmp_mz, mass_error;
  for(int r=0; r<rule_dirs.size(); r++){
    tmp_dir = rule_dirs[r];
    diff_ms = rule_ms[r];
    if(tmp_dir == 0){
      // + 1
      tmp_mz = prec_exp + diff_ms;
      mass_error = tmp_mz*ppm2*1e-6;
      x = SQLiteObjNL.extractIDMS2_with_mzRange_entireDB(tmp_mz - mass_error, tmp_mz + mass_error);
      if(x == 0){
        warning("Database searching failed for m/z = " + std::to_string(tmp_mz) + " !\n");
      }
      vector<string> allMS2refs = SQLiteObjNL.getMS2PeaksVec();
      vector<int> allIDs = SQLiteObjNL.getIDsVec();
      if(allMS2refs.size() > 0){
        double best_score = 0.2, tmp_score;
        int best_id = -1;
        for(int s=0; s<allMS2refs.size(); s++){
          string ms2ref = allMS2refs[s];
          NumericMatrix ref_spec_msms = ms2peak_parse(ms2ref);
          NumericVector nl_ref_mzs = ref_spec_msms(_,0);
          NumericVector nl_ref_ints = ref_spec_msms(_,1);
          // calculate neutral loss similarity
          if(nl_ref_mzs.size() > 1){
            tmp_score = neutral_loss_similarity(nl_exp_mzs, nl_exp_ints, nl_ref_mzs, nl_ref_ints, ppm2);
            if(tmp_score > best_score){
              best_score = tmp_score;
              best_id = allIDs[s];
            }
          }
        }
        if(best_id != -1){
          // found a useful candidate to return
          rulesVec.push_back(rules[r]);
          scores.push_back(best_score);
          IDs.push_back(best_id);
        }
      }
      
      // -1
      tmp_mz = prec_exp - diff_ms;
      mass_error = tmp_mz*ppm2*1e-6;
      x = SQLiteObjNL.extractIDMS2_with_mzRange_entireDB(tmp_mz - mass_error, tmp_mz + mass_error);
      if(x == 0){
        warning("Database searching failed for m/z = " + std::to_string(tmp_mz) + " !\n");
      }
      allMS2refs = SQLiteObjNL.getMS2PeaksVec();
      allIDs = SQLiteObjNL.getIDsVec();
      if(allMS2refs.size() > 0){
        double best_score = 0.2, tmp_score; // minimum score is 0.2
        int best_id = -1;
        for(int s=0; s<allMS2refs.size(); s++){
          string ms2ref = allMS2refs[s];
          NumericMatrix ref_spec_msms = ms2peak_parse(ms2ref);
          NumericVector nl_ref_mzs = ref_spec_msms(_,0);
          NumericVector nl_ref_ints = ref_spec_msms(_,1);
          // calculate neutral loss similarity
          if(nl_ref_mzs.size() > 1){
            tmp_score = neutral_loss_similarity(nl_exp_mzs, nl_exp_ints, nl_ref_mzs, nl_ref_ints, ppm2);
            if(tmp_score > best_score){
              best_score = tmp_score;
              best_id = allIDs[s];
            }
          }
        }
        if(best_id != -1){
          // found a useful candidate to return
          rulesVec.push_back(rules[r]);
          scores.push_back(best_score);
          IDs.push_back(best_id);
        }
      }
      
    } else {
      
      tmp_mz = prec_exp + diff_ms*tmp_dir;
      mass_error = tmp_mz*ppm2*1e-6;
      x = SQLiteObjNL.extractIDMS2_with_mzRange_entireDB(tmp_mz - mass_error, tmp_mz + mass_error);
      if(x == 0){
        warning("Database searching failed for m/z = " + std::to_string(tmp_mz) + " !\n");
      }
      vector<string> allMS2refs = SQLiteObjNL.getMS2PeaksVec();
      vector<int> allIDs = SQLiteObjNL.getIDsVec();
      if(allMS2refs.size() > 0){
        double best_score = 0.2, tmp_score;
        int best_id = -1;
        for(int s=0; s<allMS2refs.size(); s++){
          string ms2ref = allMS2refs[s];
          NumericMatrix ref_spec_msms = ms2peak_parse(ms2ref);
          NumericVector nl_ref_mzs = ref_spec_msms(_,0);
          NumericVector nl_ref_ints = ref_spec_msms(_,1);
          // calculate neutral loss similarity
          if(nl_ref_mzs.size() > 1){
            tmp_score = neutral_loss_similarity(nl_exp_mzs, nl_exp_ints, nl_ref_mzs, nl_ref_ints, ppm2);
            if(tmp_score > best_score){
              best_score = tmp_score;
              best_id = allIDs[s];
            }
          }
        }
        if(best_id != -1){
          // found a useful candidate to return
          rulesVec.push_back(rules[r]);
          scores.push_back(best_score);
          IDs.push_back(best_id);
        }
      }
    }
  }
  
  
  // Organize results a [List] with three vectors
  // <1>. vector of ID <- database int | vector<int>;
  // <2>. vector of score <- matching score (neutral loss matching : dot product) | vector<double>;
  // <3>. vector of rules <- rules name | vector<string>;
  res = List::create(Named("IDs") = IDs, 
                     Named("Scores") = scores,
                     Named("rules") = rulesVec);
  return res;
}




// [[Rcpp::export]]
List SpectraSearching(List ConsensusRes, 
                      IntegerVector idxs,
                      NumericMatrix peak_matrix,
                      double ppm_ms1,
                      double ppm_ms2,
                      double rt_tol,
                      List rt_ms1,
                      List scan_ms1,
                      int ion_mode,
                      std::string database_path = "",
                      bool use_rt = false,
                      bool enableNL = false,
                      std::string NLdatabase_path = "",
                      bool useEntropy = false) {
  // idxs are the indexs number of all consensus results (0,1,2,...)
  // ion mode, int: 0 is negative; 1 is positive.
  
  List SearchingRes(idxs.size());
  
  SqliteDriver SQLiteObj(database_path, "HMDB_experimental_PosDB", ion_mode);
  SQLiteObj.create_connection(database_path);
  
  SqliteDriver SQLiteObjNL(NLdatabase_path, "HMDB_experimental_PosDB", ion_mode);
  
  if(enableNL) {
    SQLiteObjNL.create_connection(NLdatabase_path);
  }
 
  IntegerVector ft_idxs = ConsensusRes[0];
  List all_Consensus_spec = ConsensusRes[1];
  
  int searching_size = idxs.size();
  
  double mz_max, mz_min, rt_max, rt_min, mz_med, rt_med;
  int this_idx, spec_idx;
  for(int i=0; i<searching_size; i++) {
    this_idx = ft_idxs[idxs[i]];
    spec_idx = idxs[i];
    //cout << i <<" this_idx -> " << this_idx << endl;
    cout << ".";
    mz_min = peak_matrix(this_idx, 0);
    mz_min = mz_min - mz_min*ppm_ms1*1e-6;
    mz_max = peak_matrix(this_idx, 1);
    mz_max = mz_max + mz_max*ppm_ms1*1e-6;
    mz_med = (mz_min + mz_max)/2.0;
    rt_min = peak_matrix(this_idx, 2);
    rt_max = peak_matrix(this_idx, 3);
    rt_med = (rt_min + rt_max)/2.0;
    int x;
    if(use_rt){
      x = SQLiteObj.extractIDMS2_with_mzrtRange(mz_min, mz_max, rt_min, rt_max);
    } else {
      x = SQLiteObj.extractFMMS2_with_mzRange_entireDB(mz_min, mz_max);
    }
    
    if(x == 0){
      warning("Database searching failed for m/z = " + std::to_string(mz_min) + " !\n");
    }
    vector<string> allMS2refs = SQLiteObj.getMS2PeaksVec();
    vector<string> allFormula = SQLiteObj.getFMs();
    vector<double> allPrecMZs = SQLiteObj.getPrecMZVec();
    vector<double> allRTs = SQLiteObj.getRTVec();
    vector<int> allIDs = SQLiteObj.getIDsVec();
    //cout << mz_min << " <- mz_min -- | Found ---> " << allMS2refs.size() << " | allFormula -> " << allFormula.size() << endl;
    //cout << "mz_med --> " << mz_med << endl;
    
    // Initialize vairiable
    vector<int> formula;
    double mz_score, rt_score, msms_score, iso_score, all_score; // these four values are the similarity score of the four aspects (initialization).
    double mz_exponent, rt_exponent;
    
    // Prepare isotope identification data/variables
    // <1> define some const ratio
    const double C13_12_ratio = 0.01112;
    const double H2_H1_ratio = 0.0001560243;
    const double O17_O16_ratio = 0.0003809256;
    const double O18_O16_ratio = 0.002054994;
    const double N15_N14_ratio = 0.008530711;
    const double S33_S32_ratio = 0.007893075;
    const double S34_S32_ratio = 0.04430646;
    // <2> calculate isotope ratio in original data [M+2]/[M] and [M+1]/[M]
    NumericVector exp_ratio_mp1, exp_ratio_mp2;
    double exp_ratio_mp1_mean, exp_ratio_mp2_mean, tmp_val0 = 0, tmp_val1 = 0, tmp_val2 = 0;
    for(int f=0; f<rt_ms1.size(); f++){ // file level
      NumericVector all_rts = rt_ms1[f];
      List all_scans = scan_ms1[f];
      IntegerVector idx_rts = whichTrue((all_rts > rt_min) & (all_rts < rt_max));
      for(int ix : idx_rts){ // scan level
        NumericMatrix ms1_mtx = all_scans[ix];
        IntegerVector idx_nr = whichTrue(abs(ms1_mtx(_,0) - mz_med)/(mz_med) < ppm_ms1*1e-6); // centroid level
        IntegerVector idx_nr_iso = whichTrue(abs(ms1_mtx(_,0) - mz_med - 1.003)/(mz_med + 1.003) < ppm_ms1*1e-6);
        IntegerVector idx_nr_iso2 = whichTrue(abs(ms1_mtx(_,0) - mz_med - 1.003*2)/(mz_med + 1.003*2) < ppm_ms1*1e-6);
        if(idx_nr.size() > 0){
          NumericVector ints_vec = ms1_mtx(_,1);
          NumericVector ints_vec_sub = ints_vec[idx_nr];
          NumericVector ints_vec_iso_sub = ints_vec[idx_nr_iso];
          NumericVector ints_vec_iso2_sub = ints_vec[idx_nr_iso2];
          tmp_val0 = sum(ints_vec_sub)/idx_nr.size();
          if(idx_nr_iso.size()>0){
            tmp_val1 = sum(ints_vec_iso_sub)/idx_nr_iso.size();
          }
          if(idx_nr_iso2.size()>0){
            tmp_val2 = sum(ints_vec_iso2_sub)/idx_nr_iso2.size();
          }
          // cout << "tmp 0 - > " << tmp_val0 << endl;
          // cout << "tmp 1 - > " << tmp_val1 << endl;
          // cout << "tmp 2 - > " << tmp_val2 << endl << endl;
          exp_ratio_mp1.push_back(tmp_val1/tmp_val0);
          exp_ratio_mp2.push_back(tmp_val2/tmp_val0);
        }
      }
    }
    if(exp_ratio_mp1.size() > 0){
      exp_ratio_mp1_mean = mean(exp_ratio_mp1);
    }
    if(exp_ratio_mp2.size() > 0){
      exp_ratio_mp2_mean = mean(exp_ratio_mp2);
    }
    
    List current_conss_spec = all_Consensus_spec[spec_idx];
    //cout << "current_conss_spec size ->" << current_conss_spec.size() << endl;
    //current_conss_spec size is the number of isomerics size
    List score_list(current_conss_spec.size());
    List dot_list(current_conss_spec.size());
    List neutral_loss_list(current_conss_spec.size());
    for(int k=0; k<current_conss_spec.size(); k++){
      NumericVector score_vec(allMS2refs.size());
      NumericVector dot_vec(allMS2refs.size());
      score_list[k] = score_vec;
      dot_list[k] = dot_vec;
      neutral_loss_list[k] = List();
    }
    
    // cout << "allMS2refs.size() --> " << allMS2refs.size() << endl;
    // vector<double> score_vec(allMS2refs.size()); // this vector is used to store the similarity values of all matching results from DB
    for(int s=0; s<allMS2refs.size(); s++){
      // calculate scores (scoring) within this for-loop
      
      // 1). mz_score
      mz_exponent = -0.5*(pow((allPrecMZs[s] - mz_med)/(mz_med*ppm_ms1*1e-6), 2));
      mz_score = std::exp(mz_exponent);
      
      // 2). rt_socre (optional)
      if(use_rt){
        rt_exponent = -0.5*(pow((allRTs[s] - rt_med)/(rt_tol), 2));
        rt_score = std::exp(rt_exponent);
      }
      
      // 3). isotope score
      formula = parse_formula(allFormula[s]); // formula vector
      // isotope ratio of [M+1]/[M]
      double ratio_mp1 = C13_12_ratio*formula[0] + 
                         H2_H1_ratio*formula[1] + 
                         O17_O16_ratio*formula[2] + 
                         N15_N14_ratio*formula[3] + 
                         S33_S32_ratio*formula[5];
      // isotope ratio of [M+2]/[M]
      double ratio_mp2 = C13_12_ratio*formula[0]*H2_H1_ratio*formula[1] + // 1 C13 + 1 H2
                         C13_12_ratio*formula[0]*O17_O16_ratio*formula[2] + // 1 C13 + 1 O17
                         C13_12_ratio*formula[0]*N15_N14_ratio*formula[3] + // 1 C13 + 1 N15
                         C13_12_ratio*formula[0]*S33_S32_ratio*formula[5] + // 1 C13 + 1 S33
                         H2_H1_ratio*formula[1]*O17_O16_ratio*formula[2] + // 1 H2 + 1 O17
                         H2_H1_ratio*formula[1]*N15_N14_ratio*formula[3] + // 1 H2 + 1 N15
                         H2_H1_ratio*formula[1]*S33_S32_ratio*formula[5] + // 1 H2 + 1 S33
                         O17_O16_ratio*formula[2]*N15_N14_ratio*formula[3] + // 1 O17 + 1 N15
                         O17_O16_ratio*formula[2]*S33_S32_ratio*formula[5] + // 1 O17 + 1 N15
                         N15_N14_ratio*formula[3]*S33_S32_ratio*formula[5] + // 1 S33 + 1 N15
                         O18_O16_ratio*formula[2] + // 1 O18
                         S34_S32_ratio*formula[5]; // 1 S34
      
      iso_score = 1 - abs(exp_ratio_mp1_mean - ratio_mp1) - abs(exp_ratio_mp2_mean - ratio_mp2);
      if(iso_score>1 | iso_score<0){
        iso_score = 0;
      }
      // 
      // cout << "exp_ratio_mp1_mean --> " << exp_ratio_mp1_mean << " ratio_mp1 -> " << ratio_mp1 << endl;
      // cout << "exp_ratio_mp2_mean --> " << exp_ratio_mp2_mean << " ratio_mp2 -> " << ratio_mp2 << endl;
      // cout << "iso_score -> " << iso_score << endl;
      
      // 4). ms/ms score -> allMS2refs
      // prepare consensus spectrum of current "spec_idx"
      string tmp_msms = allMS2refs[s];
      NumericMatrix ref_msms_spec = ms2peak_parse(tmp_msms);
      NumericVector reff_mzs= ref_msms_spec(_,0);
      NumericVector reff_intss= ref_msms_spec(_,1);
      double sim_dot = 0, matched_ratio = 0;
      
      for(int l=0; l<current_conss_spec.size(); l++) { // different isomerics
        NumericMatrix conss_mtx_spec = current_conss_spec[l];
        NumericVector conss_mzs= conss_mtx_spec(_,0);
        NumericVector conss_ints= conss_mtx_spec(_,1);
        NumericVector all_mzs4dot = clone(conss_mzs);
        for(double m : reff_mzs){
          if(is_true(all(abs(m-conss_mzs)/conss_mzs > ppm_ms2*1e-6))){
            all_mzs4dot.push_back(m);
          }
        }
        all_mzs4dot.sort();
        NumericVector ints_exp(all_mzs4dot.size());
        NumericVector ints_ref(all_mzs4dot.size());
        int tmp_ix1, tmp_ix2, matched_count = 0;
        
        for(int n=0; n<all_mzs4dot.size();n++){
          tmp_ix1 = whichTrue1(abs(all_mzs4dot[n]-conss_mzs)/conss_mzs < ppm_ms2*1e-6);
          tmp_ix2 = whichTrue1(abs(all_mzs4dot[n]-reff_mzs)/reff_mzs < ppm_ms2*1e-6);
          if(tmp_ix1 != -1){
            ints_exp[n] = conss_ints[tmp_ix1];
          }
          if(tmp_ix2 != -1){
            ints_ref[n] = reff_intss[tmp_ix2];
          }
          if((ints_exp[n] != 0) & (ints_ref[n] != 0)) {
            matched_count++;
          }
        }
        ints_ref = ints_ref/max(ints_ref);
        // cout << "conss_mzs   -> " << conss_mzs << endl;
        // cout << "conss_ints  -> " << conss_ints << endl;
        // cout << "reff_mzs    -> " << reff_mzs << endl;
        // cout << "reff_intss  -> " << reff_intss << endl << endl;
        // cout << "all_mzs4dot -> " << all_mzs4dot << endl;
        // cout << "ints_exp    -> " << ints_exp << endl;
        // cout << "ints_ref    -> " << ints_ref << endl << endl;
        // cout << " ------------ " << endl;
        if(useEntropy) {
          sim_dot = entropy(conss_mtx_spec, ref_msms_spec);
        } else {
          sim_dot = dot(ints_exp, ints_ref, true);
        }
        
        matched_ratio = (double)matched_count/(double)conss_mzs.size();
        //cout << "sim dot is --> " << sim_dot << " | matched_ratio--> " << matched_ratio << endl;
        
        // similarity score
        msms_score = (sim_dot + matched_ratio)/2.0;
        if(use_rt){
          all_score = (msms_score + mz_score + rt_score + iso_score*0.5)/3.5*100;
        } else {
          all_score = (msms_score + mz_score + iso_score*0.5)/2.5*100;
        }
        
        NumericVector tmp_vec = score_list[l];
        tmp_vec[s] = all_score;
        score_list[l] = tmp_vec;
        
        NumericVector tmp_vex = dot_list[l];
        tmp_vex[s] = sim_dot;
        dot_list[l] = tmp_vex;
        
      }
    }
    
    // process to see if neutral loss matching is needed
    if(enableNL){
      //cout << this_idx << " <- this_idx | i -> " << i << endl;
      for(int k=0; k<neutral_loss_list.size();k++){
        NumericVector tmp_vec = score_list[k];
        if((max(tmp_vec)<=10) | (tmp_vec.size() == 0)) {
          // run neutral loss matching and return results into "neutral_loss_list[l];" (l == k)
          // cout << "RUN NL --> Now max score is --> " << max(tmp_vec) << endl;
          // cout << "tmp_vec --> " << tmp_vec << " | k--> " << k << endl;
          // here we only use all rules by default
          vector<string> rules = getAll_rules_name();
          vector<int> rules_dir = getAll_dir_change();
          vector<double> rules_ms = getAll_ms_change();
          NumericMatrix conss_mtx_spec = current_conss_spec[k];
          NumericVector conss_mzs= conss_mtx_spec(_,0);
          NumericVector conss_ints= conss_mtx_spec(_,1);
          List nl_list = neutralLoss_matching(SQLiteObjNL, 
                                              rules, rules_dir, rules_ms,
                                              ppm_ms2, mz_med, 
                                              conss_mzs, conss_ints);
          neutral_loss_list[k] = nl_list;
        }
      }
    }

    // 
    // if(enableNL & (max(tmp_vec)<=10)) {
    //   // run neutral loss matching and return results into "neutral_loss_list[l];"
    //   cout << "RUN NL --> Now all_score is --> " << all_score << endl;
    //   vector<string> rules = getAll_rules_name();
    //   vector<int> rules_dir = getAll_dir_change();
    //   vector<double> rules_ms = getAll_ms_change();
    //   List nl_list = neutralLoss_matching(SQLiteObj, rules, rules_dir, rules_ms, 
    //                                       ppm_ms2, mz_med, conss_mzs, conss_ints);
    // }

    
    List search_res_sub = List::create(Named("IDs") = allIDs, 
                                       Named("Scores") = score_list,
                                       Named("dot_product") = dot_list,
                                       Named("Neutral_loss") = neutral_loss_list);
    
    SearchingRes[i] = search_res_sub;
  }
  
  SQLiteObj.disconnectDB();
  if(enableNL){
    SQLiteObjNL.disconnectDB();
  }
  
  return SearchingRes;
}



// [[Rcpp::export]]
List annotation_export(List searching_res, 
                       int type = 0, 
                       int topN = 10, 
                       int ion_mode = 0, 
                       std::string database_path = "", 
                       bool lipidsClass = false){
  /*
   *  searching_res, is a list returned by function, "SpectraSearching"
   *  type, integer of returned matching results; 
   *    0, return all; [compound name, inchikey, formula]
   *    1, compound name; 
   *    2, compound formula; 
   *    3, InChiKeys;
   *  topN, integer of the number to return (highest score);
   *  
   *  returned results is a List1, which contains a vector of List0. 
   *  The List0 has the same length of isomeric;
   *  The list1 has the same length as input:
   *    0: ID; 
   *    1: matching res, as needed; 
   *    2: matching score; 
   *    3: is neutral loss matching (1) or not (0);  
   *        
   */
  List annots_list(searching_res.size());
  
  SqliteDriver SQLiteObj(database_path, "HMDB_experimental_PosDB", ion_mode);
  SQLiteObj.create_connection(database_path);
  
  if(lipidsClass) {
    if(!SQLiteObj.clsTableExsiting()){
      lipidsClass = false; // forcedly change to false
      warning("No a valid lipidomics database used! Will skip lipids classfication information exporting!");
    }
  }
  
  for(int i=0; i<searching_res.size(); i++){//searching_res.size()
    //cout << "Now converting --> " << i << endl;
    cout << ".";
    List this_item = searching_res[i];
    IntegerVector allIDs = this_item[0];
    List allScores = this_item[1];
    List alldots = this_item[2];
    vector<List> this_res(allScores.size());
    
    for(int s=0; s<allScores.size(); s++){
      // this score vector
      NumericVector this_all_scores = allScores[s];
      NumericVector all_scores_sorted = clone(this_all_scores);//this_all_scores.sort(true);
      all_scores_sorted.sort(true);
      int topNx;
      NumericVector all_scores_topN0, all_scores_topN;
      if(all_scores_sorted.size() > topN){
        topNx = topN;
        all_scores_topN0 = all_scores_sorted[seq(0, topN-1)];
        for(double sc : all_scores_topN0){
          if(sc > 0){
            all_scores_topN.push_back(sc);
          }
        }
      } else {
        topNx = all_scores_sorted.size();
        all_scores_topN0 = all_scores_sorted;
        for(double sc : all_scores_topN0){
          if(sc > 0){
            all_scores_topN.push_back(sc);
          }
        }
      }

      IntegerVector allIDs_sub;
      NumericVector allScore_sub;
      NumericVector allDots_sub;
      
      NumericVector this_all_dots = alldots[s];
      int kc=0;
      IntegerVector includedVec;
      bool IDadded = false;
      for(int u=0; u<topNx; u++){
        for(int k=0; k<this_all_scores.size(); k++){
          IDadded = false;
          if(this_all_scores[k] == all_scores_topN[u]){
            // to see if this ID has already been added
            if(includedVec.size()>0){
              for(int nnn : includedVec){
                if(nnn==k){
                  IDadded = true;
                  break;
                }
              }
              if(IDadded){
                continue;
              }
            }
            includedVec.push_back(k);
            allIDs_sub.push_back(allIDs[k]);
            allScore_sub.push_back(all_scores_topN[u]);
            allDots_sub.push_back(this_all_dots[k]);
            kc++;
            break;
            //if(kc>=topNx){break;}
          }
        }
        if(kc>=topNx){break;}
      }

      if(type == 0){
        vector<CharacterVector> allResNMs = SQLiteObj.convertID2alls(allIDs_sub);
        
        List res_tmp;
        if(lipidsClass){
          // extract lipid class information of all IDs
          vector<CharacterVector> lipidList = SQLiteObj.extractClasses(allIDs_sub);
          res_tmp = List::create(_["IDs"] = allIDs_sub,
                                 _["Compounds"] = allResNMs[0],
                                 _["InchiKeys"] = allResNMs[1],
                                 _["Formula"] = allResNMs[2],
                                 _["Database"] = allResNMs[3],
                                 _["Scores"] = allScore_sub,    
                                 _["Dot_Similarity"] = allDots_sub,
                                 _["Super_Class"] = lipidList[0],
                                 _["Main_Class"] = lipidList[1],
                                 _["Sub_Class"] = lipidList[2]);
        } else {
          res_tmp = List::create(_["IDs"] = allIDs_sub,
                                 _["Compounds"] = allResNMs[0],
                                 _["InchiKeys"] = allResNMs[1],
                                 _["Formula"] = allResNMs[2],
                                 _["Database"] = allResNMs[3],
                                 _["Scores"] = allScore_sub,
                                 _["Dot_Similarity"] = allDots_sub);
        }
        this_res[s] = res_tmp;
      } else if(type == 1){
        //cout << "allIDs_sub -> " << allIDs_sub << endl;
        CharacterVector allCMPDNMs = SQLiteObj.convertID2CMPDNMs(allIDs_sub);
        List res_tmp;
        if(lipidsClass){
          // extract lipid class information of all IDs
          vector<CharacterVector> lipidList = SQLiteObj.extractClasses(allIDs_sub);
          res_tmp = List::create(_["IDs"] = allIDs_sub,
                                 _["Compounds"] = allCMPDNMs,
                                 _["Scores"] = allScore_sub,   
                                 _["Dot_Similarity"] = allDots_sub,
                                 _["Super_Class"] = lipidList[0],
                                 _["Main_Class"] = lipidList[1],
                                 _["Sub_Class"] = lipidList[2]);
        } else {
          res_tmp = List::create(_["IDs"] = allIDs_sub,
                                 _["Compounds"] = allCMPDNMs,
                                 _["Scores"] = allScore_sub, 
                                 _["Dot_Similarity"] = allDots_sub);
        }
        this_res[s] = res_tmp;
        
      } else if(type == 2){
        CharacterVector allFormulas = SQLiteObj.convertID2Formulas(allIDs_sub);
        List res_tmp;
        if(lipidsClass){
          // extract lipid class information of all IDs
          vector<CharacterVector> lipidList = SQLiteObj.extractClasses(allIDs_sub);
          res_tmp = List::create(_["IDs"] = allIDs_sub,
                                 _["Formula"] = allFormulas,
                                 _["Scores"] = allScore_sub,   
                                 _["Dot_Similarity"] = allDots_sub,
                                 _["Super_Class"] = lipidList[0],
                                 _["Main_Class"] = lipidList[1],
                                 _["Sub_Class"] = lipidList[2]);
        } else {
          res_tmp = List::create(_["IDs"] = allIDs_sub,
                                 _["Formula"] = allFormulas,
                                 _["Scores"] = allScore_sub, 
                                 _["Dot_Similarity"] = allDots_sub);
        }
        this_res[s] = res_tmp;
        
      } else if(type == 3){
        CharacterVector allInCKys = SQLiteObj.convertID2InChiKeys(allIDs_sub);
        List res_tmp;
        if(lipidsClass){
          // extract lipid class information of all IDs
          vector<CharacterVector> lipidList = SQLiteObj.extractClasses(allIDs_sub);
          res_tmp = List::create(_["IDs"] = allIDs_sub,
                                 _["InchiKeys"] = allInCKys,
                                 _["Scores"] = allScore_sub,   
                                 _["Dot_Similarity"] = allDots_sub,
                                 _["Super_Class"] = lipidList[0],
                                 _["Main_Class"] = lipidList[1],
                                 _["Sub_Class"] = lipidList[2]);
        } else {
          res_tmp = List::create(_["IDs"] = allIDs_sub,
                                 _["InchiKeys"] = allInCKys,
                                 _["Scores"] = allScore_sub, 
                                 _["Dot_Similarity"] = allDots_sub);
        }
        this_res[s] = res_tmp;
      }
      
    }

    annots_list[i] = this_res;
    
  }
  
  SQLiteObj.disconnectDB();
  return annots_list;
}




/*** R
# system.time(res <- SpectraSearching(a1, 238, peak_matrix, 10, 25, 5, rt_ms1, scan_ms1, ion_mode = 1,
#                                     database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite", enableNL = T))
# load("/data/ms2_benchmark/MTBLS2207/IROA/optilcms_neg_dda/annotated_res_fullDB.rda")
# system.time(annote_res <- annotation_export(res, 3, ion_mode = 1,
#                                             database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite"))
*/
