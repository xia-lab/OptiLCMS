#include <RcppArmadillo.h>
#include "utilities.h"
#include "dda_utilities.h"
#include "sqlite_utilities.h"

using namespace Rcpp;
using namespace std;


// This script is going to be used by DDA and DIA to consensus the spectrum data from replicates

// [[Rcpp::export]]
List SpectrumConsensus(List DecResList, 
                       NumericMatrix peak_mtx, 
                       double ppm_ms2,
                       double concensus_frac = 0.5, 
                       bool db_correction = false,
                       std::string database_path = "",
                       bool use_rt = false,
                       int ion_mode = 0,
                       bool useEntropy = false) {
  // DecResList: is a list containing the deconvolution results <- DDA
  List res;
  // Detect all features
  IntegerVector this_ft_idx;
  List thisRes;
  int n = 0, m = 0;
  for(int i=0;i<DecResList.size();i++){
    thisRes = DecResList[i];
    this_ft_idx = thisRes[2];
    n = n + this_ft_idx.size();
  }
  IntegerVector all_ft_idx(n, -1);
  
  for(int i=0;i<DecResList.size();i++){
    thisRes = DecResList[i];
    this_ft_idx = thisRes[2];
    for(int j=0;j<this_ft_idx.size();j++){
      all_ft_idx[m] = this_ft_idx[j];
      m++;
    }
  }
  
  all_ft_idx = unique(all_ft_idx);
  all_ft_idx.sort();
  
  // to generate consensus spectrum
  List thisRes_spec, thisRes_indic, tmpL;
  List concensus_res(all_ft_idx.size());
  
  bool no_isomeric;
  
  for(int i=0; i<all_ft_idx.size(); i++){ // for-loop for all features // all_ft_idx.size()
    //cout << "i -> " << i << endl;
    int thisFt_idx, thisFt_idx0;
    thisFt_idx = all_ft_idx[i];
    vector<int> idx_vec(DecResList.size());
    vector<List> spec_vec(DecResList.size());
    //vector<IntegerVector> idc_vec(DecResList.size());
    IntegerVector idc_vec_int;

    no_isomeric = true;
    for(int j=0; j<DecResList.size(); j++){ // for-loop for all files
      thisRes = DecResList[j];
      this_ft_idx = thisRes[2];
      thisRes_spec = thisRes[0];
      thisRes_indic = thisRes[1];
      thisFt_idx0 = whichTrue1(this_ft_idx == thisFt_idx);
      //cout << "thisFt_idx0 -> "<< thisFt_idx0 << endl;
      if(thisFt_idx0 != -1){ // if exists
        idx_vec[j] = thisFt_idx0;
        spec_vec[j] = thisRes_spec[thisFt_idx0];
        idc_vec_int = thisRes_indic[thisFt_idx0];
        //idc_vec[j] = idc_vec_int;
        if(idc_vec_int.size()>1){
          no_isomeric = false;
        }
      }
    }
    
    // starting concensus..
    vector<NumericMatrix> concensus_res_vec;
    if(no_isomeric){
      // a simple case: there is no isomeric found in all files for this feature
      // spec_vec, idc_vec, idx_vec;
      // NOT consider indicator at current stage
      NumericVector concensus_mzs;
      LogicalVector logVec;
      for(int d=0; d<DecResList.size(); d++){
        tmpL = spec_vec[d];
        if(tmpL.size() == 0){
          continue;
        }

        NumericMatrix num_mtx = tmpL[0];

        if(concensus_mzs.size() == 0){
          concensus_mzs = num_mtx(_,0);
        } else {
          for(int n=0; n<num_mtx.nrow(); n++) {
            logVec = abs(concensus_mzs - num_mtx(n,0))/num_mtx(n,0) > ppm_ms2*1e-6;
            if(is_true(all(logVec))){
              double new_mz = num_mtx(n,0);
              concensus_mzs.push_back(new_mz);
            }
          }
        }
      }
      concensus_mzs.sort();
      
      NumericMatrix concensus_dt(concensus_mzs.size(),0);
      concensus_dt = cbind(concensus_dt, concensus_mzs);
      for(int d=0; d<DecResList.size(); d++) {
        tmpL = spec_vec[d];
        if(tmpL.size() == 0){
          continue;
        }
        NumericVector concensus_ints(concensus_mzs.size());
        NumericMatrix num_mtx = tmpL[0];
        for(int m=0; m<concensus_mzs.size(); m++){
          for(int n=0; n<num_mtx.nrow(); n++){
            if(abs(concensus_mzs[m] - num_mtx(n,0))/num_mtx(n,0) < ppm_ms2*1e-6) {
              concensus_ints[m] = num_mtx(n,1);
              break;
            }
          }
        }
        concensus_dt = cbind(concensus_dt, concensus_ints);
      }

      // The first column [column 0] of concensus_dt is mzs vector
      concensus_res_vec.push_back(concensus_dt);
      concensus_res[i] = concensus_res_vec;
    } else {
      // a complicated case: there are isomerics found in at least one file for this feature
      // cout << "Found isomeric case !\n";
      // If there are isomerics found, will calculate the similarity and merge the paired ones with highest similarity
      // at first, determine how many pairs <- spectrumSimilarity
      // spec_vec is a vector of List including the matrixs' List (one or more member) from multiple files
      int max_num = 0, max_idx = 0; //max_num is the maximum of isomerics
      for(int u=0; u<spec_vec.size(); u++){
        List dl1 = spec_vec[u];
        if(dl1.size()>max_num){
          max_num = dl1.size();
          max_idx = u;
        }
      }
      //cout << "now the max_num is -> " << max_num << endl;
      List dl1 = spec_vec[max_idx];
      vector<List> newSpecList(max_num); // vector level is the isomerics; members in List are the matrixs from multiple files
      for(int p=0; p<max_num; p++){
        List tmpList;
        NumericMatrix nmt_dl1 = dl1[p];
        tmpList.push_back(nmt_dl1);
        newSpecList[p] = tmpList;
      }
      
      for(int u=0; u<spec_vec.size(); u++){ // multiple file
        if(u != max_idx){
          List dl2 = spec_vec[u];
          for(int q=0; q<dl2.size(); q++){ // multiple isomeric in file 2
            NumericVector sim_vec;
            NumericMatrix mtx_dl2 = dl2[q];
            for(int p =0; p<max_num; p++){ // multiple isomeric in file 0
              NumericMatrix mtx_dl1 = dl1[p];
              double simval = 0.0;
              if(useEntropy){
                // entropySimilarity
                simval = entropySimilarity(mtx_dl2, mtx_dl1, ppm_ms2);
              } else {
                simval = spectrumSimilarity(mtx_dl2, mtx_dl1, ppm_ms2);
              }
              sim_vec.push_back(simval);
            }
            int iso_idx = which_max(sim_vec);
            List tmpListx = newSpecList[iso_idx];
            tmpListx.push_back(mtx_dl2);
            newSpecList[iso_idx] = tmpListx;
          }
        }
      }

      // then, run spectrum concensus
      for(int w=0; w < newSpecList.size(); w++){
        List lt1 = newSpecList[w];
        
        // running...
        NumericVector concensus_mzs;
        LogicalVector logVec;
        for(int d=0; d<lt1.size(); d++){
          NumericMatrix num_mtx = lt1[d];
          if(concensus_mzs.size() == 0){
            concensus_mzs = num_mtx(_,0);
          } else {
            for(int n=0; n<num_mtx.nrow(); n++) {
              logVec = abs(concensus_mzs - num_mtx(n,0))/num_mtx(n,0) > ppm_ms2*1e-6;
              if(is_true(all(logVec))){
                double new_mz = num_mtx(n,0);
                concensus_mzs.push_back(new_mz);
              }
            }
          }
        }
        concensus_mzs.sort();
        
        //cout << "concensus_mzs -> " << concensus_mzs << endl;
        NumericMatrix concensus_dt(concensus_mzs.size(),0);
        concensus_dt = cbind(concensus_dt, concensus_mzs);
        for(int d=0; d<lt1.size(); d++) {
          NumericVector concensus_ints(concensus_mzs.size());
          NumericMatrix num_mtx = lt1[d];
          for(int m=0; m<concensus_mzs.size(); m++){
            for(int n=0; n<num_mtx.nrow(); n++){
              if(abs(concensus_mzs[m] - num_mtx(n,0))/num_mtx(n,0) < ppm_ms2*1e-6) {
                concensus_ints[m] = num_mtx(n,1);
                break;
              }
            }
          }
          concensus_dt = cbind(concensus_dt, concensus_ints);
        }
        // The first column [column 0] of concensus_dt is mzs vector
        concensus_res_vec.push_back(concensus_dt);
      }
      concensus_res[i] = concensus_res_vec;
    }
  }
  
  const int min_frac = floor(DecResList.size()*concensus_frac);
  
  SqliteDriver SQLiteObj(database_path, "HMDB_experimental_PosDB", ion_mode);
  if(db_correction){
    SQLiteObj.create_connection(database_path);
  }
  //cout << "Now the min_frac is --> " << min_frac << endl;
  
  //concensus_res size == all_ft_idx size. They are corresponding to each other.
  int tmp_ft_idx;
  for(int k=0; k<concensus_res.size(); k++){
    vector<NumericMatrix> this_vec_list = concensus_res[k];
    for(int l=0; l<this_vec_list.size(); l++){
      NumericMatrix this_nmtx = this_vec_list[l];
      if(this_nmtx.ncol()<=2){ 
        // it there is only one file contained this feature, skip concensus. keep it!
        continue;
      }
      NumericVector mz_col = this_nmtx(_,0);
      NumericVector new_consus_int;
      NumericVector new_consus_mz;
      for(int m=0; m<this_nmtx.nrow(); m++){
        NumericVector numv_row = this_nmtx(m,_);
        int cout_int = 0;
        double ints_val = 0;
        for(int n=1; n<numv_row.size(); n++){
          if(numv_row[n]>0) {
            ints_val = ints_val + numv_row[n];
            cout_int++;
          }
        }
        if(cout_int>min_frac){
          new_consus_mz.push_back(numv_row[0]);
          new_consus_int.push_back(ints_val/(double)cout_int);
        } else if(cout_int>1) {
          // ar least two fragments needed to enable a database assistance
          if(db_correction) {
            string ms2_ref_str;
            tmp_ft_idx = all_ft_idx[k];
            double min_mz = peak_mtx(tmp_ft_idx,0);
            double max_mz = peak_mtx(tmp_ft_idx,1);
            double min_rt = peak_mtx(tmp_ft_idx,2);
            double max_rt = peak_mtx(tmp_ft_idx,3);
            int x;
            if(use_rt){
              x = SQLiteObj.extractIDMS2_with_mzrtRange(min_mz, max_mz, min_rt, max_rt);
            } else {
              x = SQLiteObj.extractIDMS2_with_mzRange(min_mz, max_mz);
            }

            if(x == 0){
              warning("Database searching failed for m/z = " + std::to_string(min_mz) + " !\n");
            }
            vector<string> allMS2refs = SQLiteObj.getMS2PeaksVec();
            bool controller0=false;
            for(int v= 0; v<allMS2refs.size(); v++){
              NumericMatrix db_ref_msms = ms2peak_parse(allMS2refs[v]);
              NumericVector this_db_ref_msms_mzs = db_ref_msms(_,0);
              if(is_true(any(abs(numv_row[0] - this_db_ref_msms_mzs)/numv_row[0] < ppm_ms2*1e-6))){
                // this means there is a fragment of specific m/s found in db -> keep this fragment
                //cout << "Found a stuff from db\n";
                controller0 = true;
                break;
              }
            }
            if(controller0){
              new_consus_mz.push_back(numv_row[0]);
              new_consus_int.push_back(ints_val/(double)cout_int);
            }
          }
        }
      }
      new_consus_int = new_consus_int/(max(new_consus_int));
      this_nmtx = cbind(new_consus_mz, new_consus_int);
      this_vec_list[l] = this_nmtx;
    }
    concensus_res[k] = this_vec_list;
  }
  
  if(db_correction){
    SQLiteObj.disconnectDB();
  }
  res.push_back(all_ft_idx);
  res.push_back(concensus_res);
  return res;
}




/*** R
system.time(a1 <- SpectrumConsensus(DecRes, peak_matrix, 25, 0.5, FALSE, ion_mode = 0L))
*/
