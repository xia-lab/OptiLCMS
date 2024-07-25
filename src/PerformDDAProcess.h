#ifndef DDA_H
#define DDA_H

#include <RcppArmadillo.h>
#include "dda_utilities.h"
#include "hclust_ultrafast.h"
#include "utilities.h"
#include "sqlite_utilities.h"
#include "rules.h"
#include "linear_regression.h"

using namespace Rcpp;
using namespace std;

class ddaprocess {
  
private:
  NumericMatrix peak_matrix;
  NumericVector scanrt_ms1;
  NumericVector scanrt_ms2;
  List scan_ms1;
  List scan_ms2;
  NumericMatrix precursors_mzs;
  double isolation_window_size;
  double inclusion_inten_thre;
  double ppm_ms1;
  double ppm_ms2;
  double sn_thres;
  double filter_thres;
  double rt_size;
  
  bool useRT = false;
  bool enable_propagation = true;
  bool useEntropy = false;
  bool showOutput;
  
  double C13_12_ratio = 0.01112;
  double C13_diff = 1.003;
  
  int propagation_opt = 0; //0, all; 1, bio; 2, adduct; 3, fragment.
  int ion_mode; //0, negative; 1, positive.
  int thread_int;
  string database_path;
  string file_name;
  
  List results;
  // results list includes two sub list, data structure designed as below:
  // List 1: all spectrum results (including clean/deconvoluted/contaminated spectrum)
  //         Each item in the list is a vector of NumericMatrix (single member vector means no isomeric; multiple member vector means isomeric included);
  // List 2: this list include the information related to the purity of the spectrum
  //         Each item in the list is a vector of int (single integer means no isomerics; multiple int member means isomerics included)
  //         int mean: 0, clean; 1, deconvoluted; 2, partly deconvoluted; 3, contaminated (not convoluted).
  List GroupList, convolved_prec_List, decoResList; 
  // GroupList stores the indexes of precursors grouped as the same group
  // convolved_prec_List stores potential convolved precursor's info [m/z], which is corresponding to the list of ContmGroupIdxVec
  // decoResList stores the deconvolution results [SpecList + SpecIdex] <--> ContmGroupIdxVec
  IntegerVector featureIdxVec; // this vector to store the index integer of MS1 features <--matched--> precurser groups
  
  IntegerVector CleanGroupIdxVec; // this vector is used to store the indexes of clean MS/MS spectrum groups
  IntegerVector ContmGroupIdxVec; // this vector is used to store the indexes of (potentially) contaminated MS/MS spectrum groups
  
  int formatPeakMatrix() {
    // this function is used to format peak_matrix as a matrix with four columns
    if((peak_matrix.nrow() == 1) & (peak_matrix.ncol() == 1)){
      return 0;
    }
    NumericVector RT_min, RT_max;
    
    if(peak_matrix.ncol() < 2){
      stop("At least two columns (m/z and RT) are required for MS1 peak matrix.");
    } else if(peak_matrix.ncol() == 2){
      // only two columns: m/z and RT
      if(rt_size <= 0){
        stop("Error: RT size must be greater than 0 for the MS1 peak matrix with only two columns (m/z and RT)");
      } else {
        string rs = std::to_string(rt_size);
        warning("RT ranges are highly recommended to be provided, otherwise, a hard value " + rs + " will be used.");
        RT_min = peak_matrix(_,1) - rt_size;
        RT_max = peak_matrix(_,1) + rt_size;
      }
      
      NumericVector mz_min = peak_matrix(_,0) - peak_matrix(_,0)*ppm_ms1*2e-6;
      //2e-6: using 1 + 1, first 1 is used to generate range, second 1 is used to generate variation
      //NumericVector mzmint = mz_min[Range(0, 10)];
      //cout << ppm_ms1 << " <--- ppm_ms1 | mzmint --> " << mzmint << endl;
      NumericVector mz_max = peak_matrix(_,0) + peak_matrix(_,0)*ppm_ms1*2e-6;
      peak_matrix = cbind(mz_min, mz_max, RT_min, RT_max);
    } else if(peak_matrix.ncol() == 3) {
      NumericVector mz_min = peak_matrix(_,0) - peak_matrix(_,0)*ppm_ms1*2e-6;
      NumericVector mz_max = peak_matrix(_,0) + peak_matrix(_,0)*ppm_ms1*2e-6;
      RT_min = peak_matrix(_,1);
      RT_max = peak_matrix(_,2);
      peak_matrix = cbind(mz_min, mz_max, RT_min, RT_max);
    } else if(peak_matrix.ncol() == 4){
      NumericVector mz_min = peak_matrix(_,0) - peak_matrix(_,0)*ppm_ms1*1e-6;
      NumericVector mz_max = peak_matrix(_,1) + peak_matrix(_,1)*ppm_ms1*1e-6;
      RT_min = peak_matrix(_,2);
      RT_max = peak_matrix(_,3);
      peak_matrix = cbind(mz_min, mz_max, RT_min, RT_max);
    } else {
      stop("More than four columns for MS1 peak matrix is not accepted! Please correct!");
    }
    if(showOutput){
      cout << "A total of " << peak_matrix.nrow() << " MS target peaks have been formatted." << endl;
    }
    return 1;
  }
  
  int precursorsGrouping(){
    if(showOutput){
      cout << scan_ms1.size()  << " MS1 scans and " << scan_ms2.size() << " MS2 scans have been included."<< endl;
    }
    if((peak_matrix.ncol() > 1) & (peak_matrix.nrow() >= 1)){
      // there are some targets for detection
      double tf_rt_min, tf_rt_max, tf_mz_min, tf_mz_max;
      vector<int> prec_idxs;
      for(int i = 0; i < peak_matrix.nrow(); i++){
        tf_mz_min = peak_matrix(i,0);
        tf_mz_max = peak_matrix(i,1);
        tf_rt_min = peak_matrix(i,2);
        tf_rt_max = peak_matrix(i,3);
        prec_idxs.clear();
        
        NumericVector allMzs_vals = precursors_mzs(_,0);
        for(int j = 0; j < allMzs_vals.size(); j++){
          if((allMzs_vals[j] < tf_mz_min) | (allMzs_vals[j] > tf_mz_max)){
            continue;
          }
          if((scanrt_ms2[j] > tf_rt_min) & (scanrt_ms2[j] < tf_rt_max)){
            prec_idxs.push_back(j);
          }
        }
        if(prec_idxs.size() > 0){
          GroupList.push_back(prec_idxs);
          featureIdxVec.push_back(i);
        }
      }
      
    } else if(scan_ms1.size() > 0){
      //} else if(false){
      // there is no targets but MS1 has been detected, so do ms1 data processing first
      if(showOutput){
        cout << "No target feature matrix provided. Feature detection is starting.\nSeveral minutes may be needed.";
      }
      
      // Let's do a data filtration to filter low intensity signal
      
      for(int i=0; i < scan_ms1.size(); i++){
        NumericMatrix tmpMtx = scan_ms1[i];
        IntegerVector rows2rm;
        for(int j = 0; j<tmpMtx.nrow(); j++){
          if(tmpMtx(j,1) < filter_thres){
            rows2rm.push_back(j);
          }
        }
        scan_ms1[i] = row_erase(tmpMtx, rows2rm);
      }
      
      int sumMZ = 0;
      for(int i=0; i < scan_ms1.size(); i++){
        NumericMatrix tmpM = scan_ms1[i];
        sumMZ = sumMZ + tmpM.nrow();
      }
      NumericVector allMzs(sumMZ);
      int k = 0;
      for(int i = 0; i <  scan_ms1.size(); i++){
        NumericMatrix tmpM = scan_ms1[i];
        for(int j =0; j < tmpM.nrow(); j++){
          allMzs[k] = tmpM(j,0);
          k++;
        }
      }
      if(showOutput){
        cout << ".";
      }
      allMzs.sort();
      allMzs = allMzs[allMzs > 1];//lowest m/z is 1
      List allMzsx = MSCentroidsGrouping(allMzs);
      if(showOutput){
        cout << ".";
      }
      //cout << "allMzsx size --> " << allMzsx.size() << endl;
      List allRtsx;
      NumericVector mzVals, tmpV, rtVals;
      double min_mz, max_mz, min_rt, max_rt;
      vector<double> rts;
      for(int i = 0; i < allMzsx.size(); i++){ //allMzsx.size()
        //cout << "i --> " << i << endl;
        rts.clear();
        mzVals = allMzsx[i];
        min_mz = min(mzVals);
        max_mz = max(mzVals);
        
        for(int j = 0; j < scan_ms1.size(); j++){
          NumericMatrix tmpM = scan_ms1[j];
          tmpV = tmpM(_,0);
          for(int t = 0; t < tmpV.size(); t++){
            if ((tmpV[t] > min_mz) & (tmpV[t] < max_mz)){
              rts.push_back(scanrt_ms1[j]);
              break;
            }
          }
        }
        allRtsx.push_back(rts);
      }
      if(showOutput){
        cout << ".";
      }
      List RTClusters(allRtsx.size());
      for(int i=0; i < allRtsx.size(); i++){
        //cout << "Processing cluster --> " << i << endl;
        NumericVector thisRts = allRtsx[i];
        if(thisRts.size() > 2){
          NumericVector thisClsRes = auto_hclust_median(thisRts);
          RTClusters[i] = thisClsRes;
        }
      }
      // GroupList.push_back(allRtsx);
      // GroupList.push_back(allMzsx);
      // GroupList.push_back(RTClusters);
      
      //peak_matrix generation
      NumericVector mz_minVecs;
      NumericVector mz_maxVecs;
      NumericVector rt_minVecs;
      NumericVector rt_maxVecs;
      if(showOutput){
        cout << "." << endl;
      }
      for(int f=0; f<allRtsx.size();f++){
        IntegerVector thisCls;
        if(RTClusters[f] != R_NilValue){
          thisCls = RTClusters[f];
        } else {
          continue;
        }
        if(thisCls.size() == 0){
          continue;
        }
        IntegerVector UniCls = unique(thisCls);
        int clsNum = UniCls.size();
        mzVals = allMzsx[f];
        min_mz = min(mzVals) - min(rtVals)*ppm_ms2*1e-6;
        max_mz = max(mzVals) + max(rtVals)*ppm_ms2*1e-6;
        
        if(clsNum == 1){
          rtVals = allRtsx[f];
          min_rt = min(rtVals);
          max_rt = max(rtVals);
          mz_minVecs.push_back(min_mz);
          mz_maxVecs.push_back(max_mz);
          rt_minVecs.push_back(min_rt-1);
          rt_maxVecs.push_back(max_rt+1);
        } else {
          // Case for Multiple Clusters
          rtVals = allRtsx[f];
          LogicalVector logVec;
          IntegerVector idxVec;
          NumericVector rtValsSub;
          for(int m : UniCls){
            NumericVector tmpRTvls;
            //for(int n )
            logVec = m == thisCls;
            idxVec = whichTrue(logVec);
            rtValsSub = rtVals[idxVec];
            
            min_rt = min(rtValsSub);
            max_rt = max(rtValsSub);
            mz_minVecs.push_back(min_mz);
            mz_maxVecs.push_back(max_mz);
            rt_minVecs.push_back(min_rt-1);
            rt_maxVecs.push_back(max_rt+1);
          }
        }
      }
      
      peak_matrix = cbind(mz_minVecs, mz_maxVecs, rt_minVecs, rt_maxVecs);
      if(showOutput){
        cout << peak_matrix.nrow() << " peaks have been detected." << endl;
      }
      double tf_rt_min, tf_rt_max, tf_mz_min, tf_mz_max;
      vector<int> prec_idxs;
      for(int i = 0; i < peak_matrix.nrow(); i++){
        tf_mz_min = peak_matrix(i,0);
        tf_mz_max = peak_matrix(i,1);
        tf_rt_min = peak_matrix(i,2);
        tf_rt_max = peak_matrix(i,3);
        prec_idxs.clear();
        NumericVector allMzs_vals = precursors_mzs(_,0);
        for(int j = 0; j < allMzs_vals.size(); j++){
          if((allMzs_vals[j] < tf_mz_min) | (allMzs_vals[j] > tf_mz_max)){
            continue;
          }
          if((scanrt_ms2[j] > tf_rt_min) & (scanrt_ms2[j] < tf_rt_max)){
            prec_idxs.push_back(j);
          }
        }
        if(prec_idxs.size() > 0){
          GroupList.push_back(prec_idxs);
          featureIdxVec.push_back(i);
        }
      }
      
    } else {
      // here, no targets, no ms1 full scans, process precursors grouping directly
      if(showOutput){
        cout << "Neither target feature matrix nor MS1 full scan detected. All precursors will be included.";
      }
      NumericVector allMzs = precursors_mzs(_,0);
      NumericVector allMzs0 = clone(allMzs);
      allMzs.sort();
      
      List allMzsx;
      double mzdiff, thisppm;
      vector<double> mz4samePeak;
      for(int i = 1; i < allMzs.size(); i++){
        mzdiff = allMzs[i] - allMzs[i-1];
        thisppm = mzdiff/allMzs[i]*1e6;
        if(thisppm <= ppm_ms1){
          mz4samePeak.push_back(allMzs[i-1]);
        } else {
          mz4samePeak.push_back(allMzs[i-1]);
          allMzsx.push_back(mz4samePeak);
          mz4samePeak.clear();
        }
        if(i == (allMzs.size()-1)){
          if(thisppm <= ppm_ms1){
            mz4samePeak.push_back(allMzs[i]);
            allMzsx.push_back(mz4samePeak);
          } else {
            allMzsx.push_back(mz4samePeak);
            mz4samePeak.clear();
            mz4samePeak.push_back(allMzs[i]);
            allMzsx.push_back(mz4samePeak);
          }
        }
      }
      
      List allRtsx;
      NumericVector mzVals, rtVals;
      double min_mz, max_mz, min_rt, max_rt;
      for(int i = 0; i < allMzsx.size(); i++){
        NumericVector mzs4SmPeak = allMzsx[i];
        min_mz = min(mzs4SmPeak);
        max_mz = max(mzs4SmPeak);
        vector<double> rts4thisPeak;
        for(int j = 0; j < allMzs0.size(); j++){
          if((allMzs0[j] >= min_mz) & (allMzs0[j] <= max_mz)){
            rts4thisPeak.push_back(scanrt_ms2[j]);
          }
        }
        allRtsx.push_back(rts4thisPeak);
        rts4thisPeak.clear();
      }
      
      List RTClusters(allRtsx.size());
      for(int i=0; i < allRtsx.size(); i++){
        //cout << "Processing cluster --> " << i << endl;
        NumericVector thisRts = allRtsx[i];
        if(thisRts.size() > 2){
          NumericVector thisClsRes = auto_hclust_median(thisRts);
          RTClusters[i] = thisClsRes;
        } else {
          NumericVector thisClsRes(thisRts.size());
          RTClusters[i] = thisClsRes;
        }
      }
      // GroupList.push_back(allRtsx);
      // GroupList.push_back(allMzsx);
      // GroupList.push_back(RTClusters);
      
      
      //peak_matrix generation
      NumericVector mz_minVecs;
      NumericVector mz_maxVecs;
      NumericVector rt_minVecs;
      NumericVector rt_maxVecs;
      
      for(int f=0; f<allRtsx.size();f++){
        IntegerVector thisCls;
        if(RTClusters[f] != R_NilValue){
          thisCls = RTClusters[f];
        } else {
          continue;
        }
        if(thisCls.size() == 0){
          continue;
        }
        //cout << "Running into line 342 --> " << f << endl;
        IntegerVector UniCls = unique(thisCls);
        int clsNum = UniCls.size();
        mzVals = allMzsx[f];
        min_mz = min(mzVals);
        max_mz = max(mzVals);
        
        if(clsNum == 1){
          rtVals = allRtsx[f];
          min_rt = min(rtVals) - min(rtVals)*ppm_ms2*1e-6;
          max_rt = max(rtVals) + max(rtVals)*ppm_ms2*1e-6;
          mz_minVecs.push_back(min_mz);
          mz_maxVecs.push_back(max_mz);
          rt_minVecs.push_back(min_rt-1);
          rt_maxVecs.push_back(max_rt+1);
        } else {
          // Case for Multiple Clusters
          rtVals = allRtsx[f];
          LogicalVector logVec;
          IntegerVector idxVec;
          NumericVector rtValsSub;
          for(int m : UniCls){
            NumericVector tmpRTvls;
            logVec = m == thisCls;
            idxVec = whichTrue(logVec);
            rtValsSub = rtVals[idxVec];
            
            min_rt = min(rtValsSub) - min(rtVals)*ppm_ms2*1e-6;
            max_rt = max(rtValsSub) + max(rtVals)*ppm_ms2*1e-6;
            mz_minVecs.push_back(min_mz);
            mz_maxVecs.push_back(max_mz);
            rt_minVecs.push_back(min_rt-1);
            rt_maxVecs.push_back(max_rt+1);
          }
        }
      }
      
      peak_matrix = cbind(mz_minVecs, mz_maxVecs, rt_minVecs, rt_maxVecs);
      if(showOutput){
        cout << peak_matrix.nrow() << " peaks from precursors have been detected." << endl;
      }
      double tf_rt_min, tf_rt_max, tf_mz_min, tf_mz_max;
      vector<int> prec_idxs;
      for(int i = 0; i < peak_matrix.nrow(); i++){
        tf_mz_min = peak_matrix(i,0);
        tf_mz_max = peak_matrix(i,1);
        tf_rt_min = peak_matrix(i,2);
        tf_rt_max = peak_matrix(i,3);
        prec_idxs.clear();
        NumericVector allMzs_vals = precursors_mzs(_,0);
        for(int j = 0; j < allMzs_vals.size(); j++){
          if((allMzs_vals[j] < tf_mz_min) | (allMzs_vals[j] > tf_mz_max)){
            continue;
          }
          if((scanrt_ms2[j] > tf_rt_min) & (scanrt_ms2[j] < tf_rt_max)){
            prec_idxs.push_back(j);
          }
        }
        if(prec_idxs.size() > 0){
          GroupList.push_back(prec_idxs);
          featureIdxVec.push_back(i);
        }
      }
      
    }
    if(showOutput){
      //cout << GroupList.size() << " MS/MS features groups have been detected." << endl;
    }
    return 1;
  }
  
  int chimericSpectraDetection(){
    // This function is used to detect from the GroupList to find out which Feature Targets are contaminated.
    if(GroupList.size() == 0){
      warning("No grouped peaks are detected. Chimeric detection skipped!");
      return 0;
    }
    if(showOutput){
      cout << "The isolation window of DDA acquisition is " << isolation_window_size << endl;
    }
    if(scan_ms1.size() > 0){
      // if MS1 full scan existing
      NumericVector prec_mzs_vals = precursors_mzs(_,0);
      for(int i=0; i < GroupList.size(); i++){//GroupList.size()
        vector<int> thisGrpIdxs = GroupList[i];
        int fidx = featureIdxVec[i];
        double thisRt = 0, thisMz = 0;
        int k = 0;
        for(int thisIdx : thisGrpIdxs){
          thisRt = thisRt + scanrt_ms2[thisIdx];
          thisMz = thisMz + prec_mzs_vals[thisIdx];
          k++;
        }
        thisRt = thisRt/((double)k); // mean value of retention times of this group
        thisMz = thisMz/((double)k); // mean value of m/z values of this group 
        //cout << thisMz << " <-- thisMz | thisRt --> " << thisRt << endl;
        
        NumericVector resRt = abs(scanrt_ms1 - thisRt);
        
        double mz_range_min = peak_matrix(fidx, 0);
        double mz_range_max = peak_matrix(fidx, 1);
        double rt_range_min = peak_matrix(fidx, 2);
        double rt_range_max = peak_matrix(fidx, 3);
        
        int idx_near = which_min(resRt);
        // idx_near: is the index of the nearest MS1 scan
        // cout << idx_near << " <-- idx_near" << endl;
        if((scanrt_ms1[idx_near] < rt_range_min) | (scanrt_ms1[idx_near] > rt_range_max)){
          // cout << "reaching here --> rt_range_min: " << rt_range_min << " | rt_range_max: " << rt_range_max << "|| -> " << scanrt_ms1[idx_near] << endl;
          CleanGroupIdxVec.push_back(i);
          continue;
        }
        NumericMatrix thisMS1Scan = scan_ms1[idx_near];
        bool isClean = true;
        for(int j = 0; j < thisMS1Scan.nrow(); j++){
          if((abs(thisMS1Scan(j,0) - thisMz) < isolation_window_size/2.0) & 
             (thisMS1Scan(j,0) > mz_range_max | thisMS1Scan(j,0) < mz_range_min) &
             (thisMS1Scan(j,1) > inclusion_inten_thre)){
            isClean = false;
            break;
          }
        }
        if(idx_near>0){
          int idx_near1 = idx_near - 1;
          if((scanrt_ms1[idx_near1] > rt_range_min) & (scanrt_ms1[idx_near1] < rt_range_max)){
            NumericMatrix thisMS1Scan = scan_ms1[idx_near1];
            for(int j = 0; j < thisMS1Scan.nrow(); j++){
              if((abs(thisMS1Scan(j,0) - thisMz) < isolation_window_size/2.0) & 
                 (thisMS1Scan(j,0) > mz_range_max | thisMS1Scan(j,0) < mz_range_min) &
                 (thisMS1Scan(j,1) > inclusion_inten_thre)){
                isClean = false;
                break;
              }
            }
          }
        }
        if(idx_near < scan_ms1.size()-1){
          int idx_near2 = idx_near + 1;
          if((scanrt_ms1[idx_near2] > rt_range_min) & (scanrt_ms1[idx_near2] < rt_range_max)){
            NumericMatrix thisMS1Scan = scan_ms1[idx_near2];
            for(int j = 0; j < thisMS1Scan.nrow(); j++){
              if((abs(thisMS1Scan(j,0) - thisMz) < isolation_window_size/2.0) & 
                 (thisMS1Scan(j,0) > mz_range_max | thisMS1Scan(j,0) < mz_range_min) &
                 (thisMS1Scan(j,1) > inclusion_inten_thre)){
                isClean = false;
                break;
              }
            }
          }
        }
        
        if(isClean){
          CleanGroupIdxVec.push_back(i);
        } else {
          ContmGroupIdxVec.push_back(i);
        }
      }
    } else {
      // if no MS1 full scan
      // To see if the nearest precursors fall into the isolation window
      NumericVector prec_mzs_vals = precursors_mzs(_,0);
      IntegerVector thisGrp;
      for(int i=0; i < GroupList.size(); i++){
        thisGrp = GroupList[i];
        int ths_scan_idx = thisGrp[0];
        int pre_scan_idx = min(thisGrp) - 1;
        int nxt_scan_idx = max(thisGrp) + 1;
        bool isClean = true;
        if((pre_scan_idx > -1) & (nxt_scan_idx < prec_mzs_vals.size())){
          if(abs(prec_mzs_vals[pre_scan_idx] - prec_mzs_vals[ths_scan_idx]) < isolation_window_size/2.0){
            isClean = false;
          }
          if(abs(prec_mzs_vals[nxt_scan_idx] - prec_mzs_vals[ths_scan_idx]) < isolation_window_size/2.0){
            isClean = false;
          }
        } else if(pre_scan_idx == -1) {
          if(abs(prec_mzs_vals[nxt_scan_idx] - prec_mzs_vals[ths_scan_idx]) < isolation_window_size/2.0){
            isClean = false;
          }
        } else if(nxt_scan_idx == prec_mzs_vals.size()){
          if(abs(prec_mzs_vals[pre_scan_idx] - prec_mzs_vals[ths_scan_idx]) < isolation_window_size/2.0){
            isClean = false;
          }
        }
        if(isClean){
          CleanGroupIdxVec.push_back(i);
        } else {
          ContmGroupIdxVec.push_back(i);
        }
      }
    }
    if(showOutput){
      cout << "Inclusion intensity threshold: " << inclusion_inten_thre << endl;
    }
    cout << ".........." << endl;
    cout << CleanGroupIdxVec.size() << "    clean feature groups found in thread - " << thread_int << " from file: " << file_name << endl;
    cout << ContmGroupIdxVec.size() << " chimeric feature groups found in thread - " << thread_int << " from file: " << file_name << endl;
    
    return 1;
  }
  
  int findContaminationList(){
    // this function is designed to find all precursors [convolved_prec_List], which convolve with the MS/MS target precurosrs
    // this function is using the "ContmGroupIdxVec" to find out all precursors and to see if they are orphanisotopes from another precursor in the same MS1 scan
    // will not detect orphan isotopes if there is not MS1 scan included.
    double C13_diff_min = C13_diff-ppm_ms1*1e-6;
    double C13_diff_max = C13_diff+ppm_ms1*1e-6;
    int kcount = 0;
    Rcpp::List resList(ContmGroupIdxVec.size());
    
    if(scan_ms1.size()>0){
      // if there is MS1 scan
      NumericVector prec_mzs_vals = precursors_mzs(_,0);
      for(int idx : ContmGroupIdxVec){
        IntegerVector thisGrp = GroupList[idx];
        int fidx = featureIdxVec[idx];
        double thisRt = 0, thisMz = 0;
        int k = 0;
        for(int idxg : thisGrp){
          // cout << "scanrt_ms2[idxg]-> " << scanrt_ms2[idxg] << endl;
          // cout << "prec_mzs_vals[idxg]-> " << prec_mzs_vals[idxg] << endl;
          thisRt = thisRt + scanrt_ms2[idxg];
          thisMz = thisMz + prec_mzs_vals[idxg];
          k++;
        }
        
        thisRt = thisRt/((double)k); // mean value of retention times of this group
        thisMz = thisMz/((double)k); // mean value of m/z values of this group 
        double mz_range_min = peak_matrix(fidx, 0);
        double mz_range_max = peak_matrix(fidx, 1);
        double rt_range_min = peak_matrix(fidx, 2);
        double rt_range_max = peak_matrix(fidx, 3);
        NumericVector resRt = abs(scanrt_ms1 - thisRt);
        int idx_near = which_min(resRt);
        
        NumericMatrix allConvolutedPrcs(0,4); // this matrix include 4 columns [m/z, rt, intensity, orphanisotope (yes=1, no=0)]
        NumericVector mzs, rts, ints, orphs;
        
        NumericMatrix thisMS1Scan = scan_ms1[idx_near];
        for(int j = 0; j < thisMS1Scan.nrow(); j++){
          if((abs(thisMS1Scan(j,0) - thisMz) < isolation_window_size/2.0) & 
             (thisMS1Scan(j,0) > mz_range_max | thisMS1Scan(j,0) < mz_range_min) &
             (thisMS1Scan(j,1) > inclusion_inten_thre)){
            mzs.push_back(thisMS1Scan(j,0));
            rts.push_back(scanrt_ms1[idx_near]);
            ints.push_back(thisMS1Scan(j,1));
            //cout << "REaching here -->1\n";
            //if there is M+0 [mzdiff + ratio] exits
            LogicalVector bxv1 = ((thisMS1Scan(j,0)-thisMS1Scan(_,0) > C13_diff_min));
            LogicalVector bxv2 = ((thisMS1Scan(j,0)-thisMS1Scan(_,0) < C13_diff_max));
            bool bl1 = is_true(any(bxv1 & bxv2));
            bool bl_final = false;
            //cout << "bxv --> " << bxv1 << endl;
            // numC <- abs(round(theo.mass / 12)); #max. number of C in molecule
            // inten.max <- int.c12 * numC * 0.011; #highest possible intensity
            if(bl1){
              // cout << "REaching here found orphan isotopes --> FOUND [M+0] <-- \n";
              IntegerVector idxs_iso = whichTrue(bxv1 & bxv2);
              for(int ii : idxs_iso){
                double tmpInts = thisMS1Scan(ii,1);
                int numC = floor(thisMS1Scan(ii,0)/12); // maximum possible C number
                if((thisMS1Scan(j,1)/tmpInts > C13_12_ratio*0.95) & (thisMS1Scan(j,1)/tmpInts < C13_12_ratio*numC*1.05)){ //allow a variance of 5% at most
                  //cout <<kcount << " | " << idx << "<- group idx || REaching here found orphan isotopes --> truly FOUND [M+0] <-- \n";
                  bl_final = true;
                }
              }
            }
            
            if(bl_final){
              orphs.push_back(1);
            } else {
              orphs.push_back(0);
            }
          }
        }
        if(idx_near>0){
          int idx_near1 = idx_near - 1;
          if((scanrt_ms1[idx_near1] > rt_range_min) & (scanrt_ms1[idx_near1] < rt_range_max)){
            NumericMatrix thisMS1Scan = scan_ms1[idx_near1];
            for(int j = 0; j < thisMS1Scan.nrow(); j++){
              if((abs(thisMS1Scan(j,0) - thisMz) < isolation_window_size/2.0) & 
                 (thisMS1Scan(j,0) > mz_range_max | thisMS1Scan(j,0) < mz_range_min) &
                 (thisMS1Scan(j,1) > inclusion_inten_thre)){
                
                //cout << "REaching here -->2\n";
                mzs.push_back(thisMS1Scan(j,0));
                rts.push_back(scanrt_ms1[idx_near1]);
                ints.push_back(thisMS1Scan(j,1));
                //cout << "REaching here -->1\n";
                //if there is M+0 [mzdiff + ratio] exits
                LogicalVector bxv1 = ((thisMS1Scan(j,0)-thisMS1Scan(_,0) > C13_diff_min));
                LogicalVector bxv2 = ((thisMS1Scan(j,0)-thisMS1Scan(_,0) < C13_diff_max));
                bool bl1 = is_true(any(bxv1 & bxv2));
                bool bl_final = false;
                //cout << "bxv --> " << bxv1 << endl;
                if(bl1){
                  //cout << "REaching here found orphan isotopes --> FOUND [M+0] <-- \n";
                  IntegerVector idxs_iso = whichTrue(bxv1 & bxv2);
                  for(int ii : idxs_iso){
                    double tmpInts = thisMS1Scan(ii,1);
                    int numC = floor(thisMS1Scan(ii,0)/12); // maximum possible C number
                    if((thisMS1Scan(j,1)/tmpInts > C13_12_ratio*0.95) & (thisMS1Scan(j,1)/tmpInts < C13_12_ratio*numC*1.05)){ //allow a variance of 5% at most
                      // cout <<kcount << " | " << idx << "<- idx || REaching here found orphan isotopes --> truly FOUND [M+0] <-- \n";
                      bl_final = true;
                    }
                  }
                }
                
                if(bl_final){
                  orphs.push_back(1);
                } else {
                  orphs.push_back(0);
                }
                
              }
            }
          }
        }
        if(idx_near < scan_ms1.size()-1){
          int idx_near2 = idx_near + 1;
          if((scanrt_ms1[idx_near2] > rt_range_min) & (scanrt_ms1[idx_near2] < rt_range_max)){
            NumericMatrix thisMS1Scan = scan_ms1[idx_near2];
            for(int j = 0; j < thisMS1Scan.nrow(); j++){
              if((abs(thisMS1Scan(j,0) - thisMz) < isolation_window_size/2.0) & 
                 (thisMS1Scan(j,0) > mz_range_max | thisMS1Scan(j,0) < mz_range_min) &
                 (thisMS1Scan(j,1) > inclusion_inten_thre)){
                
                //cout << "REaching here -->3\n";
                mzs.push_back(thisMS1Scan(j,0));
                rts.push_back(scanrt_ms1[idx_near2]);
                ints.push_back(thisMS1Scan(j,1));
                //cout << "REaching here -->1\n";
                //if there is M+0 [mzdiff + ratio] exits
                LogicalVector bxv1 = ((thisMS1Scan(j,0)-thisMS1Scan(_,0) > C13_diff_min));
                LogicalVector bxv2 = ((thisMS1Scan(j,0)-thisMS1Scan(_,0) < C13_diff_max));
                bool bl1 = is_true(any(bxv1 & bxv2));
                bool bl_final = false;
                //cout << "bxv --> " << bxv1 << endl;
                if(bl1){
                  //cout << "REaching here found orphan isotopes --> FOUND [M+0] <-- \n";
                  IntegerVector idxs_iso = whichTrue(bxv1 & bxv2);
                  for(int ii : idxs_iso){
                    double tmpInts = thisMS1Scan(ii,1);
                    int numC = floor(thisMS1Scan(ii,0)/12); // maximum possible C number
                    if((thisMS1Scan(j,1)/tmpInts > C13_12_ratio*0.95) & (thisMS1Scan(j,1)/tmpInts < C13_12_ratio*numC*1.05)){ //allow a variance of 5% at most
                      // cout << kcount << " | " << idx << "<- idx || REaching here found orphan isotopes --> truly FOUND [M+0] <-- \n";
                      bl_final = true;
                    }
                  }
                }
                
                if(bl_final){
                  orphs.push_back(1);
                } else {
                  orphs.push_back(0);
                }
                
              }
            }
          }
        }
        allConvolutedPrcs = cbind(mzs, rts, ints, orphs);
        resList[kcount] = allConvolutedPrcs;
        kcount++;
        
        //cout << endl;
      }
    } else {
      // If there is no MS1 scan
      NumericVector prec_mzs_vals = precursors_mzs(_,0);
      NumericVector prec_ints_vals = precursors_mzs(_,1);
      
      for(int idx : ContmGroupIdxVec){
        IntegerVector thisGrp = GroupList[idx]; // thisGrp is a vector of indexs, which include scan indexs, which are contaminated by other ions
        
        int ths_scan_idx = thisGrp[0];
        int pre_scan_idx = min(thisGrp) - 1;
        int nxt_scan_idx = max(thisGrp) + 1;
        
        NumericMatrix allConvolutedPrcs(0,4); // this matrix include 4 columns [m/z, rt, intensity, orphanisotope (yes=1, no=0)]
        NumericVector mzs, rts, ints, orphs;
        
        if((pre_scan_idx > -1) & (nxt_scan_idx < prec_mzs_vals.size())){
          if(abs(prec_mzs_vals[pre_scan_idx] - prec_mzs_vals[ths_scan_idx]) < isolation_window_size/2.0){
            mzs.push_back(prec_mzs_vals[pre_scan_idx]);
            rts.push_back(scanrt_ms2[pre_scan_idx]);
            ints.push_back(prec_ints_vals[pre_scan_idx]);
            orphs.push_back(0);
          }
          if(abs(prec_mzs_vals[nxt_scan_idx] - prec_mzs_vals[ths_scan_idx]) < isolation_window_size/2.0){
            mzs.push_back(prec_mzs_vals[nxt_scan_idx]);
            rts.push_back(scanrt_ms2[nxt_scan_idx]);
            ints.push_back(prec_ints_vals[nxt_scan_idx]);
            orphs.push_back(0);
          }
        } else if(pre_scan_idx == -1) {
          if(abs(prec_mzs_vals[nxt_scan_idx] - prec_mzs_vals[ths_scan_idx]) < isolation_window_size/2.0){
            mzs.push_back(prec_mzs_vals[nxt_scan_idx]);
            rts.push_back(scanrt_ms2[nxt_scan_idx]);
            ints.push_back(prec_ints_vals[nxt_scan_idx]);
            orphs.push_back(0);
          }
        } else if(nxt_scan_idx == prec_mzs_vals.size()){
          if(abs(prec_mzs_vals[pre_scan_idx] - prec_mzs_vals[ths_scan_idx]) < isolation_window_size/2.0){
            mzs.push_back(prec_mzs_vals[pre_scan_idx]);
            rts.push_back(scanrt_ms2[pre_scan_idx]);
            ints.push_back(prec_ints_vals[pre_scan_idx]);
            orphs.push_back(0);
          }
        }
        
        allConvolutedPrcs = cbind(mzs, rts, ints, orphs);
        resList[kcount] = allConvolutedPrcs;
        kcount++;
      }
      
    }
    
    convolved_prec_List = resList;
    return 1;
  }
  
  /*
   * 
   this function is a general function to deconvolve all convolved (chimeric) spectrum
   convolved_prec_List <---> ContmGroupIdxVec
   
   General idea for this function: 
   I: prepare data
   1). extract the ms/ms spectrum for deconvolution from raw data list;                              --> spectrum2Deco
   2). check isomerics within the same precursor group (over 5 mz centroids) based on similarity (cosine) of ms/ms peaks;
   3). extract the ms/ms reference spectrum as main candidate from database (with m/z);              --> spectrum_main
   4). extract the ms/ms reference spectrum of contamination (chimeras) as candidates from database (with m/z); --> spectrum_candis
   II: candidates selection
   5). find the most similar spectrum (to spectrum2Deco) from main candidate list;
   6). find the most similar spectrum (to spectrum2Deco) from contamination candidate lists (which are not an orphan isotope);
   7). predict the spectrum of orphan isotopes based on the spectrum from their corresponding [M+0] ions, and select the most similar one as the candidate;
   III: perform deconvolution
   8). summary main candidate and other candidates [non-prphan isotopes + orphan isotopes, knowns + unknowns] as model spectrum;
   9). run Elastic Net model to simulate the "spectrum2Deco";
   10). pseudo spectrum generation for main spectrum;
   11). if main spectrum is also contaminating others, use the pseudo spectrum to deconvolve others ("unknown" without matching to db, or "knowns" with match to db);
   IV: redundancy check and scoring
   12). check redundancy towards the database
   13). score the results based on the similarity and retention time (if applicable based on users' input)
   14). return the deconvolution results for consensus analysis
   */
  List tmpmtx;
  List performDeco(){
    
    // 1) extract raw ms/ms spectrum + 2) isomeric distinguishing
    List rawMS2Deco(ContmGroupIdxVec.size());
    for(int i=0; i<ContmGroupIdxVec.size();i++){
      //variable 2 use: scan_ms2, GroupList, convolved_prec_List
      int &idx = ContmGroupIdxVec[i];
      // cout << "&idx is --> " << &idx << endl;
      //cout << "idx is --> " << idx << endl;
      vector<int> idxs = GroupList[idx];
      //cout << "idx is --> " << idxs << endl;
      //rawMS2Deco[i] = scan_ms2;
      NumericMatrix res = spectrumPreparation(idxs);
      //cout << "spectrumPreparation DONE is --> " << idx << endl;
      rawMS2Deco[i] = res;
      
      // if(idxs.size()>3){
      //   NumericVector v1 = res(_,0);
      //   NumericVector v2 = res(_,1);
      //   NumericVector v3 = res(_,2);
      //   cout << "this res(_,0) --> " << v1 << endl;
      //   cout << "this res(_,1) --> " << v2 << endl;
      //   cout << "this res(_,2) --> " << v3 << endl << endl;
      // }
      //cout << i << " <- i | res.nrow() --> " << res.nrow() << endl;
      idxs.clear();
    }
    if(showOutput){
      //cout << "rawMS2Deco we have got is --> " << rawMS2Deco.size() << endl;
    }
    // 3) extract ms/ms of main candidate + 4) contamination candidate from DB 
    // [on-the-fly with deconvolution]
    NumericVector prec_mzs_vals = precursors_mzs(_,0);
    SqliteDriver SQLiteObj(database_path, "HMDB_experimental_PosDB", ion_mode);
    SQLiteObj.create_connection(database_path);
    
    int nn =0, mm=0;
    
    List allSpectraList, allSpectraList_res, allSpectraList_idx;
    // allSpectraList is a part for List results (for contaminated and deconvoluted spectrum);
    // allSpectraList_res is the 1st List; allSpectraList_idx is the 2nd List;
    
    for(int i=0; i<ContmGroupIdxVec.size(); i++) { //ContmGroupIdxVec.size()
      //cout << "Now running deco dda i -> " << i << endl;
      cout << ".";
      int ii = ContmGroupIdxVec[i];
      vector<int> thisGrpIdxs = GroupList[ii];
      
      double thisRt = 0, thisMz = 0;
      int k = 0;
      for(int thisIdx : thisGrpIdxs){
        thisRt = thisRt + scanrt_ms2[thisIdx];
        thisMz = thisMz + prec_mzs_vals[thisIdx];
        k++;
      }
      thisRt = thisRt/((double)k); // mean value of retention times of this group
      thisMz = thisMz/((double)k); // mean value of m/z values of this group 
      // cout << "Now the mz of main candidate is --> " << thisMz << endl;
      // cout << "Now the rt of main candidate is --> " << thisRt << endl;
      // 
      // Now let's start extration --> thisMz <--
      double &prec_rt = thisRt;
      if(!useRT){
        prec_rt = -1.0;
      }
      
      NumericMatrix thisSpectra = rawMS2Deco[i];
      
      vector<NumericMatrix> main_candidate_msms = extractMainCandidates(SQLiteObj, thisSpectra, thisMz, 
                                                                        ppm_ms1, ppm_ms2, rt_size, prec_rt,
                                                                        useEntropy);
      List chimeras_candidate_msms, chimeras_candidate_msms_prec;
      if (main_candidate_msms.size()>0) {
        // (1). will run the contamination candidate extraction ONLY if the main candidate exists
        NumericMatrix chimericPrecrsrs = convolved_prec_List[i];
        chimeras_candidate_msms = extractContamCandidates(SQLiteObj, thisSpectra, chimericPrecrsrs,
                                                          ppm_ms1, ppm_ms2, rt_size, useRT, useEntropy);
        // (2). next step, using clean (pure) DDA spectrum to see if any of them can be used to help deconvolve the
        // chimera with main candidate but without contam candidate (List empty)
        // Check if the clean spectrum can be helpful for deconvolution (as candidate)
        double min_rtc, max_rtc, min_mzc, max_mzc, ch_mz, ch_rt, mass_error;
        bool cond1, cond2;
        NumericVector prec_ints_vals = precursors_mzs(_,1);
        vector<NumericMatrix> cont_ms_vec;
        int idxxx;
        for(int cl=0; cl<CleanGroupIdxVec.size();cl++){
          cond1 = false;
          cond2 = false;
          cont_ms_vec.clear();
          IntegerVector scan_idxs = GroupList[cl];
          NumericVector rt_vec = scanrt_ms2[scan_idxs];
          NumericVector mz_vec = prec_mzs_vals[scan_idxs];
          NumericVector ints_vec = prec_ints_vals[scan_idxs];
          min_rtc = min(rt_vec);
          max_rtc = max(rt_vec);
          min_mzc = min(mz_vec);
          max_mzc = max(mz_vec);
          mass_error = max_rtc*ppm_ms1*1e-6;
          for(int ch=0; ch<chimericPrecrsrs.nrow();ch++){
            ch_mz = chimericPrecrsrs(ch,0);
            ch_rt = chimericPrecrsrs(ch,1);
            if((ch_mz > min_mzc-mass_error) & (ch_mz < max_mzc+mass_error)){cond1 = true;}
            if((ch_rt > min_rtc) & (ch_rt < max_rtc)){cond2 = true;}
          }
          if(cond1 & cond2){
            idxxx = which_max(ints_vec);
            idxxx = scan_idxs[idxxx];
            for(int s=0; s < main_candidate_msms.size(); s++){
              cont_ms_vec.push_back(scan_ms2[idxxx]); // push multiple times because the spectrum is clean
            }
            chimeras_candidate_msms.push_back(cont_ms_vec);
          }
        }
        
        // (3). if chimeras_candidate_msms.size() <= 1 (still) -> have to find another approach (leverage concept from NetID)
        // Searching homologues/or similar structure from database to generate pseudo spectrum for deconvolution 
        // [NedID approach]:  May also search bio-transform/abiotic transformation results.
        if((chimeras_candidate_msms.size() <= 1) & enable_propagation){
          chimeras_candidate_msms_prec = predictRef_Spectrum(SQLiteObj, 
                                                             thisSpectra, 
                                                             chimericPrecrsrs,
                                                             ppm_ms1, 
                                                             ppm_ms2, 
                                                             rt_size, 
                                                             useRT, 
                                                             propagation_opt,
                                                             useEntropy);
          // if(chimeras_candidate_msms_prec.size()>0){
          //   //cout << "find some predicted contamination spectra\n";
          // }
        }
      }
      
      // Now let's perform deconvolution with Elastic Net linear regression model (if main_candidate_msms.size > 0)
      // If main_candidate_msms size is 0, the spectra will be marked as Non-deconvoluted directly.
      // The linear regression model is ready in the linear_regression.h/cpp script
      if (main_candidate_msms.size()>0) {
        // The main candiate exits, no matter if there is any chimeric candidate (predicted or not);
        // For linear regression, we need two or three arguments, (components matrix, response vector, and [optional] loading penalty)
        // <1>. components matrix consists of main candiates [main_candidate_msms], chimeric candidate [chimeras_candidate_msms] and predicted chimeric candidates [chimeras_candidate_msms_prec]
        // <2>. response vector consists of the (convolved) experimental spectrum --> thisSpectra;
        // <3>. (Optional) assign penalty to the predicted chimeric candidate based on the similarity of towards to raw experimental spectrum
        
        // It is noted that: [main_candidate_msms] is a vector, which includes one or more member(s) from different isomerics (a1);
        // [chimeras_candidate_msms] is a list, which include one or more vector<NumericMatrix>, the number of "vector<NumericMatrix>" is consistent to [main_candidate_msms] (a1)
        // one or more members (a2, a2 may != a1) are included in each vector;
        // [chimeras_candidate_msms_prec] includes ONLY one most similar predicted spectrum, which is not corresponding to the number (a1);
        
        NumericMatrix prec_msms_candi0, main_msms_candi, exp_msms_data;
        
        if(chimeras_candidate_msms_prec.size() > 0){
          NumericMatrix prec_msms_candi = chimeras_candidate_msms_prec[0]; // there is only one predicted spectrum to deconvolute all
          prec_msms_candi0 = prec_msms_candi;
        }
        // thisSpectra includes three columns: 0. mz; 1, int; 2, class;
        vector<NumericMatrix> allChimCandiMtxs;
        vector<NumericMatrix> deco_exp_spectrum(main_candidate_msms.size());
        vector<int> deco_idxs(main_candidate_msms.size());
        
        for(int i=0; i<main_candidate_msms.size();i++){
          // the idx of main_candidate_msms is consistent with the class label in "this spectra"
          exp_msms_data = extractSingleSpectrum(thisSpectra, i); // to deconvolve
          main_msms_candi = main_candidate_msms[i]; 
          allChimCandiMtxs.clear();
          if (chimeras_candidate_msms.size() != 0) {
            for(int j=0; j<chimeras_candidate_msms.size();j++){
              vector<NumericMatrix> ChimCandiVec = chimeras_candidate_msms[j];
              NumericMatrix thisChimCandiMtx = ChimCandiVec[i];
              allChimCandiMtxs.push_back(thisChimCandiMtx);
            }
          }
          // summarize all candidates
          NumericVector allMzs;
          // I. process experimental data:
          double thismz;
          for(int k=0; k<exp_msms_data.nrow();k++){
            thismz = exp_msms_data(k,0);
            allMzs.push_back(thismz);
          }
          
          // II. process main candidate data
          bool blx1;
          for(int kk=0; kk<main_msms_candi.nrow();kk++){
            blx1 = false;
            thismz = main_msms_candi(kk,0);
            for(int kl=0; kl<allMzs.size();kl++){
              if(abs(allMzs[kl] - thismz)/thismz < 1e-6*ppm_ms2){
                blx1 = true;
              }
            }
            if(!blx1){
              // Found a fragment not in the experimental pattern
              allMzs.push_back(thismz);
            }
          }
          
          // III. process chimeric candidate (if exits)
          if((allChimCandiMtxs.size()>0) & (chimeras_candidate_msms.size() != 0)){
            bool blx2;
            for(int km=0; km<allChimCandiMtxs.size();km++){
              NumericMatrix thisChimCandiMtx = allChimCandiMtxs[km];
              for(int kk=0; kk<thisChimCandiMtx.nrow();kk++){
                blx2 = false;
                thismz = thisChimCandiMtx(kk,0);
                for(int kl=0; kl<allMzs.size();kl++){
                  if(abs(allMzs[kl] - thismz)/thismz < 1e-6*ppm_ms2){
                    blx2 = true;
                  }
                }
                if(!blx2){
                  // Found a fragment not in current pattern
                  allMzs.push_back(thismz);
                }
              }
            }
          }
          
          // IV. process predicted candidate (if exits)
          if(prec_msms_candi0.nrow()>0){
            bool blx3;
            for(int kk=0; kk<prec_msms_candi0.nrow();kk++){
              blx3 = false;
              thismz = prec_msms_candi0(kk,0);
              for(int kl=0; kl<allMzs.size();kl++){
                if(abs(allMzs[kl] - thismz)/thismz < 1e-6*ppm_ms2){
                  blx3 = true;
                }
              }
              if(!blx3){
                // Found a fragment not in all current pattern
                allMzs.push_back(thismz);
              }
            }
          }
          allMzs.sort();
          // V. process the intensity normalization and NumericMatrix generation
          // Now all mzs are ready, start process int for all candidates (as a NumericMatrix)
          
          NumericVector response_vec(allMzs.size()), loading_penalties; 
          // V_I. let's generate "response_vec"
          double thisInt;
          double max_int = max(exp_msms_data(_,1));
          NumericVector nmv1 = exp_msms_data(_,0);
          
          for(int l=0; l<allMzs.size();l++){
            for(int k=0; k<exp_msms_data.nrow();k++){
              thismz = exp_msms_data(k,0);
              thisInt = exp_msms_data(k,1);
              if (abs(thismz - allMzs[l])/thismz < 1e-6*ppm_ms2) {
                response_vec[l] = thisInt/max_int;
              }
            }
          }
          
          // V_II. let's generate "all_candidates_mtx";
          NumericMatrix all_candidates_mtx(allMzs.size(),0);
          NumericVector cand_vec(allMzs.size());
          // main
          max_int = max(main_msms_candi(_,1));
          for(int l=0; l<allMzs.size();l++){
            for(int kk=0; kk<main_msms_candi.nrow();kk++){
              thismz = main_msms_candi(kk,0);
              thisInt = main_msms_candi(kk,1);
              if(abs(allMzs[l] - thismz)/thismz < 1e-6*ppm_ms2){
                cand_vec[l] = thisInt/max_int;
              }
            }
          }
          all_candidates_mtx = cbind(all_candidates_mtx, cand_vec);
          loading_penalties.push_back(1);
          // chimeric
          if((allChimCandiMtxs.size()>0) & (chimeras_candidate_msms.size() != 0)){
            for(int km=0; km<allChimCandiMtxs.size();km++){
              NumericVector cand_vec(allMzs.size());
              NumericMatrix thisChimCandiMtx = allChimCandiMtxs[km];
              max_int = max(thisChimCandiMtx(_,1));
              for(int l=0; l<allMzs.size();l++){
                for(int kk=0; kk<thisChimCandiMtx.nrow();kk++){
                  thismz = thisChimCandiMtx(kk,0);
                  thisInt = thisChimCandiMtx(kk,1);
                  if(abs(allMzs[l] - thismz)/thismz < 1e-6*ppm_ms2){
                    cand_vec[l] = thisInt/max_int;
                  }
                }
              }
              all_candidates_mtx = cbind(all_candidates_mtx, cand_vec);
              loading_penalties.push_back(1);
            }
          }
          // predicted_chimeric
          if(prec_msms_candi0.nrow()>0){
            max_int = max(prec_msms_candi0(_,1));
            NumericVector cand_vec(allMzs.size());
            for(int l=0; l<allMzs.size();l++){
              for(int kk=0; kk<prec_msms_candi0.nrow();kk++){
                thismz = prec_msms_candi0(kk,0);
                thisInt = prec_msms_candi0(kk,1);
                if(abs(allMzs[l] - thismz)/thismz < 1e-6*ppm_ms2){
                  cand_vec[l] = thisInt/max_int;
                }
              }
            }
            all_candidates_mtx = cbind(all_candidates_mtx, cand_vec);
            double sim_val = 0.0;
            // calculate similarity again Y_Y <spectrumSimilarity>
            if(useEntropy){
              sim_val = entropySimilarity(prec_msms_candi0, exp_msms_data, ppm_ms2);
            } else {
              sim_val = spectrumSimilarity(prec_msms_candi0, exp_msms_data, ppm_ms2);
            }
            loading_penalties.push_back(10 - sim_val*9);
            //cout << "sim_val -> " << sim_val << endl;
            //cout << "Now the loading_penalties -> " << loading_penalties << endl;
          }
          
          // VI. Run Linear Regression [Penalized Elastic Net Analysis]
          // response_vec, all_candidates_mtx, loading_penalties;
          NumericVector betaVec;
          if(loading_penalties.size() > 1){
            betaVec = PerformLinearRegress(all_candidates_mtx, response_vec, loading_penalties);
          } else {
            NumericVector dbl_empty;
            betaVec = PerformLinearRegress(all_candidates_mtx, response_vec, dbl_empty);
          }
          
          // NumericVector tmvv = all_candidates_mtx(_,0);
          //cout << "db spect_vec -> " << tmvv << endl;
          //cout << "response_vec -> " <<response_vec << endl;
          //cout << all_candidates_mtx.ncol() << " <- all_candidates_mtx | betaVec --> " << betaVec << endl << endl;
          
          // VII. Summarize deconvolution results
          vector<int> keepIdxes;
          //cout << " betaVec --> " << betaVec << endl << endl;
          for(int b=0; b<betaVec.size();b++){
            if(betaVec[b] > 0){
              keepIdxes.push_back(b);
            }
          }
          
          if(keepIdxes.size() == 0){
            // all beta equals or less than 0 --> CANNOT be deconvoluted at all ===>> return orignal spectrum directly [3]
            // allSpectraList_res | allSpectraList_idx <-- two sub lists
            NumericVector all_exp_ints = exp_msms_data(_,1);
            NumericVector all_exp_mzs = exp_msms_data(_,0);
            double max_exp_ints = max(all_exp_ints);
            all_exp_ints = all_exp_ints / max_exp_ints;
            NumericMatrix exp_msms_data_normalized = cbind(all_exp_mzs, all_exp_ints);
            deco_exp_spectrum[i] = exp_msms_data_normalized; // normalized
            deco_idxs[i] = 3;
          } else if(betaVec[0] == 0) {
            // the main candidate does not contribute, other contaminating ions contribute --> partly deconvoluted [2]
            // exclude convolved spectrum ===>> remove and return the remaining spectrum as results
            deco_exp_spectrum[i] = reconstructSpectrum(response_vec, allMzs, exp_msms_data, all_candidates_mtx, betaVec);
            deco_idxs[i] = 2;
          } else {
            // successfully deconvoluted ===>> return the deconvoluted spectrum [1]
            deco_exp_spectrum[i] = reconstructSpectrum(response_vec, allMzs, exp_msms_data, all_candidates_mtx, betaVec);
            deco_idxs[i] = 1;
          }
          
        }
        
        allSpectraList_res.push_back(deco_exp_spectrum);
        allSpectraList_idx.push_back(deco_idxs);
        
      } else {
        // No main candidate exists, will not deconvolve this spectrum
        // cout << "No main candidate found, return the raw spectrum directly" << endl;
        
        // allSpectraList_res | allSpectraList_idx <-- two sub lists
        // prepare spectrum/spectra data
        NumericVector thisclss = thisSpectra(_,2);
        NumericVector uni_thisClss = unique(thisclss);
        //cout << "now this uni_thisClss is --> " << uni_thisClss << endl;
        vector<NumericMatrix> undeco_exp_spectrum(uni_thisClss.size());
        vector<int> undeco_idxs(uni_thisClss.size());
        if(uni_thisClss.size() == 1){
          NumericVector specmz = thisSpectra(_,0);
          NumericVector specit = thisSpectra(_,1);
          specit = specit / max(specit);
          undeco_exp_spectrum[0] = cbind(specmz, specit);
          undeco_idxs[0] = 3;
        } else {
          int min_bd, max_bd, w=0;
          for(double u : uni_thisClss){
            IntegerVector idxVec = whichTrue(u == thisclss);
            min_bd = min(idxVec); // lower limit
            max_bd = max(idxVec); // upper limit
            NumericMatrix m2 = thisSpectra(Range(min_bd,max_bd) ,_);
            NumericVector specmz2 = m2(_,0);
            NumericVector specit2 = m2(_,1);
            specit2 = specit2 / max(specit2);
            undeco_exp_spectrum[w] = cbind(specmz2, specit2);;
            undeco_idxs[w] = 3; // Non-deconvoluted at all
            w++;
          }
        }
        
        allSpectraList_res.push_back(undeco_exp_spectrum);
        allSpectraList_idx.push_back(undeco_idxs);
      }
      
      
    }
    
    SQLiteObj.disconnectDB();
    allSpectraList = List::create(_["SpecList"] = allSpectraList_res, 
                                  _["SpecIdex"] = allSpectraList_idx);
    decoResList = allSpectraList;
    return allSpectraList;
  }
  
  NumericMatrix spectrumPreparation(vector<int> &idxs){
    NumericMatrix res(0,3); // 3 columns as results [1st, mz; 2nd, int; 3rd, cluster(isomerics)]
    if(idxs.size() == 1){
      // 1, simply return the spetrum
      NumericMatrix thisNMtx = scan_ms2[idxs[0]];
      NumericVector thisVecs (thisNMtx.nrow());
      thisNMtx = cbind(thisNMtx, thisVecs);
      return thisNMtx;
    } else if(idxs.size() >= 4) {
      // 4 or more, need to cluster based on cosine similarity of MS/MS spectrum to distinguish isomerics
      // generating matrix of ms/ms spectra
      NumericVector mzs;
      //int cc = 0;
      for(int id : idxs){
        NumericMatrix df = scan_ms2[id];
        if(id == idxs[0]){
          mzs = df(_,0);
          //  cc = df.nrow();
        } else {
          NumericVector thismzs = df(_,0);
          //cc = cc + df.nrow();
          for(double t : thismzs){
            bool included = false;
            for(double m : mzs){
              if(abs(t-m)/m < ppm_ms2*1e-6){
                included = true;
                break;
              }
            }
            if(!included){
              mzs.push_back(t);
            }
          }
        }
      }
      // calculate weights of different spectrum
      NumericVector these_prec_ints, weights_vec;
      for(int id : idxs){
        double this_prc_int = precursors_mzs(id,1);
        these_prec_ints.push_back(this_prc_int);
      }
      if(max(these_prec_ints) == 0){
        weights_vec = rep(1, these_prec_ints.size());
        //cout << "weights_vec --> " << weights_vec << endl;
      } else {
        weights_vec = these_prec_ints/max(these_prec_ints);
        //
      }
      //cout << "weights_vec --> " << weights_vec << endl;
      //cout << idxs.size() << " <-- idxs.size() | now the mzs number is --> " << mzs.size() << " | cc --> " << cc << endl;
      // now all mzs are ready -> re-run the for loop to generate a matrix for cosine similarity
      mzs.sort();
      NumericMatrix dataMatrix(mzs.size(), idxs.size());
      
      int w = 0;
      for(int id : idxs){
        NumericMatrix df = scan_ms2[id];
        NumericVector thismzs = df(_,0);
        NumericVector thisInts = df(_,1);
        NumericVector thisCol(mzs.size());
        for(int j=0; j < mzs.size(); j++){
          double &m = mzs[j];
          for(int i=0; i < thismzs.size(); i++){
            double &t = thismzs[i];
            if(abs(t-m)/m < ppm_ms2*1e-6){
              thisCol[j] = thisInts[i];
              break;
            }
          }
        }
        dataMatrix(_,w) = thisCol;
        w++;
      }
      NumericMatrix mt = cosineSimilarity(dataMatrix);
      NumericVector cls = matrix_hclust(mt);
      // now we have clusters infomation
      // cout << "Now cl is -> " << cls << endl;
      NumericVector cls_unique = unique(cls);
      cls_unique.sort();
      NumericVector allMzs, allInts, allCls;
      for(int k=0; k<cls_unique.size(); k++){
        NumericMatrix res(0,3);
        double this_cls = cls_unique[k];
        for(int l=0; l<cls.size();l++){
          double cl = cls[l];
          if(cl != this_cls){
            continue;
          }
          int id = idxs[l];
          NumericMatrix df = scan_ms2[id];
          //cout << id << " <- id | now DF nrow is " << df.nrow() << endl;
          //cout << "weights_vec[l] -> " << weights_vec[l] << endl;
          if(res.nrow() == 0){
            NumericVector specInts = df(_,1);
            NumericVector specMzs = df(_,0);
            specInts = specInts*weights_vec[l];
            res = cbind(specMzs, specInts);
            //res = df;
          } else {
            // merge iteratively
            NumericVector newMzs = res(_,0), newInts = res(_,1);
            for(int d=0; d < df.nrow(); d++){
              bool includ = false; int rr = 0;
              for(int r=0; r < res.nrow(); r++){
                if(abs(res(r,0)-df(d,0))/res(r,0) < ppm_ms2*1e-6){
                  includ = true;
                  rr = r;
                  break;
                }
              }
              if(includ){
                newMzs[rr] = (res(rr,0)+df(d,0))/2.0;
                newInts[rr] = res(rr,1)+df(d,1)*weights_vec[l];
              } else {
                newMzs.push_back(df(d,0));
                newInts.push_back(df(d,1)*weights_vec[l]);
              }
            }
            res = cbind(newMzs, newInts);
          }
        }
        
        for(int rs=0; rs < res.nrow(); rs++){
          allMzs.push_back(res(rs,0));
          allInts.push_back(res(rs,1));
          allCls.push_back(this_cls);
        }
        
      }
      res = cbind(allMzs, allInts, allCls);
      return res;
      
    } else {
      // 2-3, merge all of them for further deconvolution
      List these_scan_ms2;
      NumericVector these_prec_ints;
      for(int id : idxs){
        NumericMatrix dfxxx = scan_ms2[id];
        double this_prc_int = precursors_mzs(id,1);
        these_scan_ms2.push_back(dfxxx);
        these_prec_ints.push_back(this_prc_int);
      }
      res = spectrumMerging(these_prec_ints, these_scan_ms2, ppm_ms2);
      
      // below is the old code chunk <<-- replaced with a static function (spectrumMerging)
      /*
       * 
       for(int id : idxs){
       NumericMatrix df = scan_ms2[id];
       //cout << id << " <- id | now DF nrow is " << df.nrow() << endl;
       if(res.nrow() == 0){
       //cout << "now res nrow is 0\n";
       res = df;
       } else {
       // merge iteratively
       NumericVector newMzs = res(_,0), newInts = res(_,1);
       for(int d=0; d < df.nrow(); d++){
       bool includ = false; int rr = 0;
       for(int r=0; r < res.nrow(); r++){
       if(abs(res(r,0)-df(d,0))/res(r,0) < ppm_ms2*1e-6){
       includ = true;
       rr = r;
       break;
       }
       }
       if(includ){
       newMzs[rr] = (res(rr,0)+df(d,0))/2.0;
       newInts[rr] = res(rr,1)+df(d,1);
       } else {
       newMzs.push_back(df(d,0));
       newInts.push_back(df(d,1));
       }
       }
       
       res = cbind(newMzs, newInts);
       }
       }
       */
      
      NumericVector thisVecs (res.nrow());
      res = cbind(res, thisVecs);
    }
    return res;
  }
  
  int MS2listSummarize(){
    /// this function we summarize clean and deconvoluted (as well as undeconvoluted) spectrum
    // generate class variable -> ["results"]
    List Spectra_res(GroupList.size());
    List Indicator_res(GroupList.size());
    
    // Process Clean spectra first - merging
    int cln_idx;
    NumericVector prec_mzs = precursors_mzs(_,0);
    NumericVector prec_int = precursors_mzs(_,1);
    
    IntegerVector grpIdxs;
    int gpidx;
    IntegerVector intv0 = {0};
    for(int cc=0; cc<CleanGroupIdxVec.size();cc++) {
      cln_idx = CleanGroupIdxVec[cc];
      grpIdxs = GroupList[cln_idx];
      
      NumericVector prec_int_vec(grpIdxs.size());
      List spectList(grpIdxs.size());
      for(int g=0; g < grpIdxs.size(); g++) {
        gpidx = grpIdxs[g];
        prec_int_vec[g] = prec_int[gpidx];
        spectList[g] = scan_ms2[gpidx];
      }
      NumericMatrix nmtx = spectrumMerging(prec_int_vec, spectList, ppm_ms2);
      Spectra_res[cln_idx] = List::create(nmtx);
      Indicator_res[cln_idx] = intv0; // a fixed value, include only a "0";
    }
    
    // Adding deconvoluted spectra into the list from decoResList
    int dec_idx, tmp_idx = 0;
    List spec_decoed = decoResList[0];
    List indc_decoed = decoResList[1];
    for(int dd=0; dd<ContmGroupIdxVec.size(); dd++){
      dec_idx = ContmGroupIdxVec[dd];
      Spectra_res[dec_idx] = spec_decoed[dd];
      Indicator_res[dec_idx] = indc_decoed[dd];
    }
    
    // push back all results into [results]
    results = List::create(_["Spectra"] = Spectra_res , 
                           _["Indicator"] = Indicator_res);
    
    // cout << "CleanGroupIdxVec -> " <<  CleanGroupIdxVec.size() << endl;
    // cout << "ContmGroupIdxVec -> " <<  ContmGroupIdxVec.size() << endl;
    return 1;
  }
  
  int MS2listSummarize_noDeco(){
    /// this function we summarize all undeconvoluted spectrum
    // generate class variable -> ["results"]
    List Spectra_res(GroupList.size());
    List Indicator_res(GroupList.size());
    
    NumericVector prec_int = precursors_mzs(_,1);
    
    IntegerVector grpIdxs;
    int gpidx;
    IntegerVector intv0 = {4}; // 4 here means no neither chimeric checking nor deconvolution
    for(int cc=0; cc<GroupList.size();cc++) {
      grpIdxs = GroupList[cc];
      NumericVector prec_int_vec(grpIdxs.size());
      List spectList(grpIdxs.size());
      for(int g=0; g < grpIdxs.size(); g++) {
        gpidx = grpIdxs[g];
        prec_int_vec[g] = prec_int[gpidx];
        spectList[g] = scan_ms2[gpidx];
      }
      
      NumericMatrix nmtx = spectrumMerging(prec_int_vec, spectList, ppm_ms2);
      Spectra_res[cc] = List::create(nmtx);
      Indicator_res[cc] = intv0; // a fixed value, include only a "4" here;
    }
    
    // push back all results into [results]
    results = List::create(_["Spectra"] = Spectra_res , 
                           _["Indicator"] = Indicator_res);
    
    return 1;
  }
  
  int MS2listSummarize_noMerge(){
    /// this function we summarize all undeconvoluted spectrum
    // generate class variable -> ["results"]
    List Spectra_res(GroupList.size());
    List Indicator_res(GroupList.size());
    
    NumericVector prec_int = precursors_mzs(_,1);
    
    IntegerVector grpIdxs;
    int gpidx, maxidx = 0;
    IntegerVector intv0 = {4}; // 4 here means no neither chimeric checking nor deconvolution
    for(int cc=0; cc<GroupList.size();cc++) {
      grpIdxs = GroupList[cc];
      NumericVector prec_int_vec(grpIdxs.size());
      List spectList(grpIdxs.size());
      for(int g=0; g < grpIdxs.size(); g++) {
        gpidx = grpIdxs[g];
        prec_int_vec[g] = prec_int[gpidx];
        spectList[g] = scan_ms2[gpidx];
      }
      maxidx = which_max(prec_int_vec);
      NumericMatrix nmtx = spectList[maxidx]; //spectrumMerging(prec_int_vec, spectList, ppm_ms2);
      Spectra_res[cc] = List::create(nmtx);
      Indicator_res[cc] = intv0; // a fixed value, include only a "4" here;
    }
    
    // push back all results into [results]
    results = List::create(_["Spectra"] = Spectra_res , 
                           _["Indicator"] = Indicator_res);
    
    return 1;
  }
  
  NumericMatrix static spectrumMerging(NumericVector prec_ints, List SpectraLists, double ppm_val){
    NumericMatrix mergedSpectrum, res;
    // cout << "prec_ints  -->" << prec_ints.size() << endl;
    // cout << "SpectraLists->" << SpectraLists.size() << endl << endl;
    
    NumericVector weights;// = prec_ints/max(prec_ints);
    if(max(prec_ints) == 0){
      weights = rep(1, prec_ints.size());
    } else {
      weights = prec_ints/max(prec_ints);
    }
    
    for(int i=0; i<prec_ints.size(); i++){
      NumericMatrix df = SpectraLists[i];
      if(res.nrow() == 0){
        NumericVector specInts = df(_,1);
        NumericVector specMzs = df(_,0);
        specInts = specInts*weights[i];
        res = cbind(specMzs, specInts);
      } else {
        // merge iteratively
        NumericVector newMzs = res(_,0), newInts = res(_,1);
        for(int d=0; d < df.nrow(); d++){
          bool includ = false; int rr = 0;
          for(int r=0; r < res.nrow(); r++){
            if(abs(res(r,0)-df(d,0))/res(r,0) < ppm_val*1e-6){
              includ = true;
              rr = r;
              break;
            }
          }
          if(includ){
            newMzs[rr] = (res(rr,0)+df(d,0))/2.0;
            newInts[rr] = res(rr,1)+df(d,1)*weights[i];
          } else {
            newMzs.push_back(df(d,0));
            newInts.push_back(df(d,1)*weights[i]);
          }
        }
        res = cbind(newMzs, newInts);
      }
    }
    //NumericVector thisVecs (res.nrow());
    //mergedSpectrum = cbind(res, thisVecs);
    mergedSpectrum = res;
    return mergedSpectrum;
  }
  
  static vector<NumericMatrix> extractMainCandidates(SqliteDriver SQLiteObj, 
                                                     NumericMatrix thisSpectrum, 
                                                     double precrsr_mz, 
                                                     double ppm_ms1,
                                                     double ppm_ms2, 
                                                     double rt_tol,
                                                     double precrsr_rt = -1.0,
                                                     bool useEntropy = false){
    vector<NumericMatrix> mainCandidates;
    bool use_rt = false;
    if(precrsr_rt >= 0.0){
      use_rt = true;
    }
    
    // prepare spectrum/spectra data
    NumericVector thisclss = thisSpectrum(_,2);
    NumericVector uni_thisClss = unique(thisclss);
    //cout << "now this uni_thisClss is --> " << uni_thisClss << endl;
    vector<NumericMatrix> main_exp_spectrum(uni_thisClss.size());
    if(uni_thisClss.size() == 1){
      main_exp_spectrum[0] = thisSpectrum;
    } else {
      int n=0, min_bd, max_bd;
      for(double u : uni_thisClss){
        //cout << "n --> " << n << endl;
        IntegerVector idxVec = whichTrue(u == thisclss);
        min_bd = min(idxVec); // lower limit
        max_bd = max(idxVec); // upper limit
        NumericMatrix m2 = thisSpectrum(Range(min_bd,max_bd) ,_);
        main_exp_spectrum[n] = m2;
        n++;
      }
    }
    
    // extract from database tables
    string ms2_ref_str;
    NumericMatrix ms2_ref_mtx;
    double ms_error = precrsr_mz*ppm_ms1*1e-6;
    int x;
    if(use_rt){
      x = SQLiteObj.extractIDMS2_with_mzrtRange(precrsr_mz - ms_error, precrsr_mz + ms_error, 
                                                precrsr_rt - rt_tol, precrsr_rt + rt_tol);
    } else {
      x = SQLiteObj.extractIDMS2_with_mzRange_entireDB(precrsr_mz - ms_error, precrsr_mz + ms_error);
    }
    
    if(x == 0){
      warning("Database searching failed for m/z = " + std::to_string(precrsr_mz) + " !\n");
    }
    vector<int> allIDs = SQLiteObj.getIDsVec();
    vector<string> allMS2refs = SQLiteObj.getMS2PeaksVec();
    // cout << precrsr_mz << " <- precrsr_mz | allMS2refs size --> " << allMS2refs.size() << endl;
    // if(allIDs.size() > 2){
    //   cout << "allIDs[0] -> " << allIDs[0] << endl;
    //   cout << "allIDs[1] -> " << allIDs[1] << endl;
    //   cout << "allIDs[2] -> " << allIDs[2] << endl;
    // }
    
    // Calculate matrix similarity by using function spectrumSimilarity or entropySimilarity
    if (allIDs.size() != 0) {
      //cout << "main_exp_spectrum.size() --> " << main_exp_spectrum.size() << endl;
      for(int p=0; p<main_exp_spectrum.size();p++){
        NumericMatrix best_mtx;
        double curr_simi = 0.0, last_simi = -0.1;
        for(string ms2_ref_str : allMS2refs){
          //ms2_ref_str = allMS2refs[0];
          ms2_ref_mtx = ms2peak_parse(ms2_ref_str);
          // cout << "extracting " << precrsr_mz << ", and got num of ID is -> " << allIDs.size() << endl;
          // NumericVector t1 = ms2_ref_mtx(_,0);
          // NumericVector t2 = ms2_ref_mtx(_,1);
          // cout << "mzs:-->" << t1 << endl;
          // cout << "ints:-->" << t2 << endl;
          if(useEntropy){
            curr_simi = entropySimilarity(main_exp_spectrum[p], ms2_ref_mtx, ppm_ms2);
          } else {
            curr_simi = spectrumSimilarity(main_exp_spectrum[p], ms2_ref_mtx, ppm_ms2);
          }
          if(curr_simi > last_simi){
            last_simi = curr_simi;
            best_mtx = ms2_ref_mtx;
          }
          //cout << "ok, now the sim num is --> " << curr_simi << endl;
        }
        mainCandidates.push_back(best_mtx);
      }
    } else {
      // No reference spectrum found from database, do nothing.
    }
    
    return mainCandidates;
  }
  
  vector<int> static assignSubFormula_CarbonLimited(double molecula_weight, double ppm, int C_max){
    vector<int> formula{0,0,0,0,0,0};
    // this function is used to assign formulas to one ms/ms fragment
    // Possible chemical elements: 
    // C (12.00), H (1.007825), O (15.994915), N (14.003074), P (30.973763), S (31.972072)
    
    int O_max = floor(molecula_weight/15.995);
    int N_max = floor(molecula_weight/14.003);
    int P_max = floor(molecula_weight/30.974);
    int S_max = floor(molecula_weight/31.972);
    
    double tmp_mw;
    for(int s=0; s<= S_max;s++) {
      for(int p=0; p<= P_max;p++) {
        for(int n=0; n<= N_max;n++) {
          for(int o=0; o<= O_max;o++) {
            for(int c=0; c<= C_max;c++) {
              for(int h=0; h<=c*4+o*2+n*3+p*3+s*2; h++) {
                tmp_mw = c*12.00+h*1.007825+o*15.994915+n*14.003074+p*30.973763+s*31.972072;
                if((abs(tmp_mw - molecula_weight)/molecula_weight) < ppm*1e-6){
                  formula[0] = c;
                  formula[1] = h;
                  formula[2] = o;
                  formula[3] = n;
                  formula[4] = p;
                  formula[5] = s;
                  return formula;
                }
              } 
            } 
          }
        }
      }
    }
    
    return formula;
  }
  
  vector<int> static assignFullFormula(double molecula_weight, double ppm){
    vector<int> formula{0,0,0,0,0,0};
    // this function is used to assign formulas to one ms/ms fragment
    // Possible chemical elements: 
    // C (12.00), H (1.007825), O (15.994915), N (14.003074), P (30.973763), S (31.972072)
    int C_max = floor(molecula_weight/12.00);
    //int H_max = C_max*4;
    int O_max = floor(molecula_weight/15.995);
    int N_max = floor(molecula_weight/14.003);
    int P_max = floor(molecula_weight/30.974);
    int S_max = floor(molecula_weight/31.972);
    
    double tmp_mw;
    for(int s=0; s<= S_max;s++) {
      for(int p=0; p<= P_max;p++) {
        for(int n=0; n<= N_max;n++) {
          for(int o=0; o<= O_max;o++) {
            for(int c=0; c<= C_max;c++) {
              for(int h=floor((c+o+n+p+s)*0.5); h<=c*4+o*2+n*3+p*3+s*2; h++) {
                tmp_mw = c*12.00+h*1.007825+o*15.994915+n*14.003074+p*30.973763+s*31.972072;
                if((abs(tmp_mw - molecula_weight)/molecula_weight) < ppm*1e-6){
                  formula[0] = c;
                  formula[1] = h;
                  formula[2] = o;
                  formula[3] = n;
                  formula[4] = p;
                  formula[5] = s;
                  return formula;
                }
              } 
            } 
          }
        }
      }
    }
    
    return formula;
  }
  
  vector<int> static assignFullFormula_accurate(double molecula_weight, double ppm){
    vector<int> formula{0,0,0,0,0,0};
    // this function is used to assign formulas to one ms/ms fragment
    // Possible chemical elements: 
    // C (12.00), H (1.007825), O (15.994915), N (14.003074), P (30.973763), S (31.972072)
    int C_max = floor(molecula_weight/12.00);
    //int H_max = C_max*4;
    int O_max = floor(molecula_weight/15.995);
    int N_max = floor(molecula_weight/14.003);
    int P_max = floor(molecula_weight/30.974);
    int S_max = floor(molecula_weight/31.972);
    
    double tmp_mw, tmp_val, best_ppm = ppm;
    for(int s=0; s<= S_max;s++) {
      for(int p=0; p<= P_max;p++) {
        for(int n=0; n<= N_max;n++) {
          for(int o=0; o<= O_max;o++) {
            for(int c=0; c<= C_max;c++) {
              for(int h=floor((c+o+n+p+s)*0.5); h<=c*4+o*2+n*3+p*3+s*2; h++) {
                tmp_mw = c*12.00+h*1.007825+o*15.994915+n*14.003074+p*30.973763+s*31.972072;
                tmp_val = (abs(tmp_mw - molecula_weight)/molecula_weight)*1e6;
                if(tmp_val < best_ppm){
                  formula[0] = c;
                  formula[1] = h;
                  formula[2] = o;
                  formula[3] = n;
                  formula[4] = p;
                  formula[5] = s;
                  best_ppm = tmp_val;
                }
              } 
            } 
          }
        }
      }
    }
    //cout << "best_ppm -> " << best_ppm << endl;
    return formula;
  }
  
  List static assignSubFormula_all(double molecula_weight, double ppm){
    List formulaList;
    // this function is used to all possible assign formulas to one ms/ms fragment within the ppm range
    // Possible chemical elements: 
    // C (12.00), H (1.007825), O (15.994915), N (14.003074), P (30.973763), S (31.972072)
    int C_max = floor(molecula_weight/12.00);
    int O_max = floor(molecula_weight/15.995);
    int N_max = floor(molecula_weight/14.003);
    int P_max = floor(molecula_weight/30.974);
    int S_max = floor(molecula_weight/31.972);
    vector<int> formula;
    double tmp_mw, tmp_val;
    for(int s=0; s<= S_max;s++) {
      for(int p=0; p<= P_max;p++) {
        for(int n=0; n<= N_max;n++) {
          for(int o=0; o<= O_max;o++) {
            for(int c=0; c<= C_max;c++) {
              for(int h=0; h<=c*4+o*2+n*3+p*3+s*2; h++) {
                tmp_mw = c*12.00+h*1.007825+o*15.994915+n*14.003074+p*30.973763+s*31.972072;
                tmp_val = (abs(tmp_mw - molecula_weight)/molecula_weight)*1e6;
                if(tmp_val < ppm){
                  formula = {0,0,0,0,0,0};
                  formula[0] = c;
                  formula[1] = h;
                  formula[2] = o;
                  formula[3] = n;
                  formula[4] = p;
                  formula[5] = s;
                  formulaList.push_back(formula);
                  if(formulaList.size() > 3000) {
                    return formulaList;
                  }
                  //best_ppm = tmp_val;
                }
              } 
            } 
          }
        }
      }
    }
    //cout << "best_ppm -> " << best_ppm << endl;
    return formulaList;
  }
  
  NumericMatrix static orphanIsotopePredict(NumericMatrix ms2_mtx, 
                                            double prec_mz, 
                                            double ppm1, 
                                            double ppm2){
    // This function is used to predict fragments of orphan isotopes based on database reference spectrum
    
    // orphan isotopes contamination is from the M+1, whose M+0 is not within the isolation window but the M+1 falls into the 
    // window and the intensity of M+1 is over the threshold.
    
    // Possible chemical elements: C (12.00), H (1.007825), O (15.994915), N (14.003074), P (30.973763), S (31.972072)
    
    NumericVector mzs = ms2_mtx(_,0);
    NumericVector ints = ms2_mtx(_,1);
    
    // cout << "ms2_mtx 1 -> " << mzs <<endl;
    // cout << "ms2_mtx 2 -> " << ints <<endl;
    //
    vector<int> Complete_formula = assignFullFormula(prec_mz, ppm1);
    int C_num_total = Complete_formula[0];
    
    double mz_iso_val, mz_val, int_val, int_iso_val, int_noniso_val;
    vector<int> sub_formula;
    int C_num_sub;
    double nonIso_ratio, iso_ratio;
    
    NumericVector new_mzs, new_ints;
    
    for(int i=0; i<mzs.size();i++){
      sub_formula.clear();
      mz_val = mzs[i];
      mz_iso_val = mz_val + 1.003;
      int_val = ints[i];
      sub_formula = assignSubFormula_CarbonLimited(mz_val, ppm2, C_num_total);
      C_num_sub = sub_formula[0];
      
      iso_ratio = (double) C_num_sub/C_num_total;
      nonIso_ratio = 1 - iso_ratio;
      
      int_iso_val = int_val*iso_ratio;
      int_noniso_val = int_val*nonIso_ratio;
      
      if(int_iso_val > 0.0){
        new_mzs.push_back(mz_iso_val);
        new_ints.push_back(int_iso_val);
      }
      if(int_noniso_val > 0.0){
        new_mzs.push_back(mz_val);
        new_ints.push_back(int_noniso_val);
      }
      // cout << iso_ratio << " <--- iso_ratio | nonIso_ratio --> " << nonIso_ratio << endl;
    }
    NumericMatrix predictedPattern = cbind(new_mzs, new_ints);
    
    return predictedPattern;
  }
  
  static List extractContamCandidates(SqliteDriver SQLiteObj, 
                                      NumericMatrix thisSpectrum, 
                                      NumericMatrix precrsr_mtx, 
                                      double ppm_ms1,
                                      double ppm_ms2, 
                                      double rt_tol,
                                      bool useRT,
                                      bool useEntropy){
    
    // prepare spectrum/spectra data
    NumericVector thisclss = thisSpectrum(_,2);
    NumericVector uni_thisClss = unique(thisclss);
    //cout << "now this uni_thisClss is --> " << uni_thisClss << endl;
    vector<NumericMatrix> main_exp_spectrum(uni_thisClss.size());
    if(uni_thisClss.size() == 1){
      main_exp_spectrum[0] = thisSpectrum;
    } else {
      int n=0, min_bd, max_bd;
      for(double u : uni_thisClss){
        //cout << "n --> " << n << endl;
        IntegerVector idxVec = whichTrue(u == thisclss);
        min_bd = min(idxVec); // lower limit
        max_bd = max(idxVec); // upper limit
        NumericMatrix m2 = thisSpectrum(Range(min_bd,max_bd) ,_);
        main_exp_spectrum[n] = m2;
        n++;
      }
    }
    
    // we extract contamination reference spectrum based on chimeric precurosors (mz and/or rt)
    // columns of precrsr_mtx: 0 mzs, 1 rts, 2 ints, 3 orphs[0, not orphan isotope; 1, is an orphan isotope]
    //cout << "OK, now we are extracting chimeras ref spectrum from db -- precrsr_mtx : " << precrsr_mtx.nrow() << endl;
    NumericVector allMzs = precrsr_mtx(_,0);
    NumericVector allRts = precrsr_mtx(_,1);
    NumericVector allInts = precrsr_mtx(_,2);
    NumericVector allOrphan = precrsr_mtx(_,3);
    
    // I. Summarize (unique) all precursors (m/z ppm + intensity + rt + orphan iso)
    // m/z <-- intensity (highest); rt, mean; 
    // 1), for the ones not orphan isotopes
    IntegerVector idx_nOrph = whichTrue(allOrphan == 0);
    NumericVector allMzs_nOrphan = allMzs[idx_nOrph];
    NumericVector allInts_nOrphan = allInts[idx_nOrph];
    NumericVector allRts_nOrphan = allRts[idx_nOrph];
    //cout << "allMzs_nOrphan --> " << allMzs_nOrphan << endl;
    
    double current_max_int = max(allInts_nOrphan), current_mz, current_rt;
    vector <double> new_mzs, new_rts;
    int idxx;
    bool bl1;
    while(allInts_nOrphan.size()>0){
      //cout << "allInts_nOrphan -> " << allInts_nOrphan << endl;
      //cout << "new_mzs -> " << new_mzs.size() << " | --> " << new_rts.size() << endl;
      if(new_mzs.size() == 0){
        idxx = whichTrue1(allInts_nOrphan == current_max_int);
        current_mz = allMzs_nOrphan[idxx];
        current_rt = allRts_nOrphan[idxx];
        new_mzs.push_back(current_mz);
        new_rts.push_back(current_rt);
      } else {
        bl1 = false;
        idxx = whichTrue1(allInts_nOrphan == current_max_int);
        current_mz = allMzs_nOrphan[idxx];
        for(double lmz : new_mzs){
          if(abs(lmz - current_mz)/current_mz < ppm_ms1*1e-6){
            bl1 = true;
            break;
          }
        }
        if(!bl1){
          // not included for this centroid (within the ppm)
          new_mzs.push_back(current_mz);
          new_rts.push_back(current_rt);
        }
      }
      // erase
      allInts_nOrphan.erase(idxx);
      allRts_nOrphan.erase(idxx);
      allMzs_nOrphan.erase(idxx);
      // find a new maximum
      current_max_int = max(allInts_nOrphan);
    }
    
    // 2). for the ones are orphan isotopes
    IntegerVector idx_Orph = whichTrue(allOrphan == 1);
    vector <double> new_mzs_iso, new_rts_iso;
    if(idx_Orph.size()>0){
      // code will work only orphan isotopes exit
      NumericVector allMzs_Orphan = allMzs[idx_Orph];
      NumericVector allInts_Orphan = allInts[idx_Orph];
      NumericVector allRts_Orphan = allRts[idx_Orph];
      
      double current_max_int = max(allInts_Orphan), current_mz, current_rt;
      while(allInts_Orphan.size()>0){
        //cout << "allInts_Orphan -> " << allInts_Orphan << endl;
        //cout << "new_mzs -> " << new_mzs.size() << " | --> " << new_rts_iso.size() << endl;
        if(new_mzs_iso.size() == 0){
          idxx = whichTrue1(allInts_Orphan == current_max_int);
          current_mz = allMzs_Orphan[idxx];
          current_rt = allRts_Orphan[idxx];
          new_mzs_iso.push_back(current_mz);
          new_rts_iso.push_back(current_rt);
        } else {
          bl1 = false;
          idxx = whichTrue1(allInts_Orphan == current_max_int);
          current_mz = allMzs_Orphan[idxx];
          for(double lmz : new_mzs_iso){
            if(abs(lmz - current_mz)/current_mz < ppm_ms1*1e-6){
              bl1 = true;
              break;
            }
          }
          if(!bl1){
            // not included for this centroid (within the ppm)
            new_mzs_iso.push_back(current_mz);
            new_rts_iso.push_back(current_rt);
          }
        }
        // erase
        allInts_Orphan.erase(idxx);
        allRts_Orphan.erase(idxx);
        allMzs_Orphan.erase(idxx);
        // find a new maximum
        current_max_int = max(allInts_Orphan);
      }
    }
    
    // mz and rt summary done [new_mzs, new_rts] + [new_mzs_iso, new_rts_iso]
    // II. Extract reference spectrum from database
    // 1). Handle non-isotope precursor 
    //     (extract -> similarity matching -> return contamCandidateList)
    List contamCandidateList(new_mzs.size()+new_mzs_iso.size());
    int k =0;
    string ms2_ref_str;
    NumericMatrix ms2_ref_mtx;
    vector<int> allIDs;
    vector<string> allMS2refs;
    vector<NumericMatrix> contamCandidates;
    for(int i=0; i<new_mzs.size(); i++){
      double this_mz = new_mzs[i];
      double this_rt = new_rts[i];
      allIDs.clear();
      allMS2refs.clear();
      
      double ms_error = this_mz*ppm_ms1*1e-6;
      int x;
      if(useRT){
        x = SQLiteObj.extractIDMS2_with_mzrtRange(this_mz - ms_error, this_mz + ms_error, 
                                                  this_rt - rt_tol, this_rt + rt_tol);
      } else {
        x = SQLiteObj.extractIDMS2_with_mzRange_entireDB(this_mz - ms_error, this_mz + ms_error);
      }
      
      if(x == 0){
        warning("Database searching failed for m/z = " + std::to_string(this_mz) + " !\n");
      }
      
      allIDs = SQLiteObj.getIDsVec();
      allMS2refs = SQLiteObj.getMS2PeaksVec();
      
      if(allIDs.size()>0){
        //cout << this_mz << " <- this_mz | allIDs is over 0 now\n";
        for(int p=0; p<main_exp_spectrum.size();p++){
          NumericMatrix best_mtx;
          double curr_simi = 0.0, last_simi = -0.1;
          for(string ms2_ref_str : allMS2refs){
            //ms2_ref_str = allMS2refs[0];
            ms2_ref_mtx = ms2peak_parse(ms2_ref_str);
            // cout << "extracting " << precrsr_mz << ", and got num of ID is -> " << allIDs.size() << endl;
            // NumericVector t1 = ms2_ref_mtx(_,0);
            // NumericVector t2 = ms2_ref_mtx(_,1);
            // cout << "mzs:-->" << t1 << endl;
            // cout << "ints:-->" << t2 << endl;
            if(useEntropy){
              curr_simi = entropySimilarity(main_exp_spectrum[p], ms2_ref_mtx, ppm_ms2);
            } else {
              curr_simi = spectrumSimilarity(main_exp_spectrum[p], ms2_ref_mtx, ppm_ms2);
            }
            
            if(curr_simi > last_simi){
              last_simi = curr_simi;
              best_mtx = ms2_ref_mtx;
            }
            //cout << "ok, now the sim num is --> " << curr_simi << endl;
          }
          contamCandidates.push_back(best_mtx);
        }
      } else {
        // cout << this_mz << " <- this_mz | allIDs is 0 now\n";
      }
      contamCandidateList[k] = contamCandidates;
      k++;
    }
    
    // 2). Handle orphan isotope precursor 
    //     (calculate M+0 -> extract -> spectrum prediction -> similarity matching -> return contamCandidateList)
    for(int i=0; i < new_mzs_iso.size(); i++){
      double this_mz_iso = new_mzs_iso[i];
      double this_rt_iso = new_rts_iso[i];
      double this_mz0 = this_mz_iso-1.003; //1.003 is diff between c13 and c12
      allIDs.clear();
      allMS2refs.clear();
      
      double ms_error = this_mz0*ppm_ms1*1e-6;
      int x;
      if(useRT){
        x = SQLiteObj.extractIDMS2_with_mzrtRange(this_mz0 - ms_error, this_mz0 + ms_error, 
                                                  this_rt_iso - rt_tol, this_rt_iso + rt_tol);
      } else {
        x = SQLiteObj.extractIDMS2_with_mzRange_entireDB(this_mz0 - ms_error, this_mz0 + ms_error); //_entireDB
      }
      
      if(x == 0){
        warning("Database searching failed for m/z = " + std::to_string(this_mz0) + " !\n");
      }
      
      allIDs = SQLiteObj.getIDsVec();
      allMS2refs = SQLiteObj.getMS2PeaksVec();
      if(allIDs.size()>0){
        // found some ref spectra, predict orphan spectrum -> similarity matching -> attached to contamCandidateList
        // cout << "found isotopes from db --> " << this_mz0 << " | allMS2refs-> " << allMS2refs.size() << endl;
        
        for(int p=0; p<main_exp_spectrum.size();p++){
          NumericMatrix best_mtx, new_ms2_ref_mtx;
          double curr_simi = 0.0, last_simi = -0.1;
          for(string ms2_ref_str : allMS2refs){
            ms2_ref_mtx = ms2peak_parse(ms2_ref_str);
            // NumericVector t1 = ms2_ref_mtx(_,0);
            // NumericVector t2 = ms2_ref_mtx(_,1);
            // cout << "mzs:-->" << t1 << endl;
            // cout << "ints:-->" << t2 << endl;
            new_ms2_ref_mtx = orphanIsotopePredict(ms2_ref_mtx, this_mz0, ppm_ms1, ppm_ms2);
            //cout << new_ms2_ref_mtx.nrow() << " <--- new_ms2_ref_mtx | ms2_ref_mtx--> " << ms2_ref_mtx.nrow() << endl;
            if(useEntropy){
              curr_simi = entropySimilarity(main_exp_spectrum[p], new_ms2_ref_mtx, ppm_ms2);
            } else {
              curr_simi = spectrumSimilarity(main_exp_spectrum[p], new_ms2_ref_mtx, ppm_ms2);
            }
            if(curr_simi > last_simi){
              last_simi = curr_simi;
              best_mtx = ms2_ref_mtx;
            }
            //cout << "ok, now the sim num is --> " << curr_simi << endl;
          }
          contamCandidates.push_back(best_mtx);
        }
        contamCandidateList[k] = contamCandidates;
        k++;
      } else {
        //cout << "NOT found isotopes from db --> " << this_mz0 << endl;
      }
      
    }
    
    // cout << "new_mzs.size()     -> " << new_mzs.size() << endl;
    // cout << "new_mzs_iso.size() -> " << new_mzs_iso.size() << endl;
    // 
    
    // Clean empty List component
    List contamCandidateList_clean;
    for(int g=0; g< contamCandidateList.size(); g++){
      vector<NumericMatrix> vec_num = contamCandidateList[g];
      if(vec_num.size() !=0){
        contamCandidateList_clean.push_back(vec_num);
      }
    }
    
    return contamCandidateList_clean;
  }
  
  // int isomericDeco(){
  //   // m/z is 100% same, but RT may be different
  //   // If EIC of ms2 exits, use the similarity of EIC-MS2 peaks
  //   // Using clustering to distinguish belonging of MS2 fragments
  //   
  //   // This function maybe used later to consider EICs of MS2. 
  //   // Details: extract EICs of MS2 and select model peaks and deconvolute shared EICs by isomerics
  //   // TODO later.
  //   return 1;
  // }
  
  NumericMatrix static reconstructSpectrum(NumericVector exp_vec, NumericVector mzs_vec,
                                           NumericMatrix exp_mtx, NumericMatrix candidate_Spec, 
                                           NumericVector coef_Vec) {
    NumericMatrix reconsctedSpe;
    if(coef_Vec[0] == 0) {
      // substract contaminating ions and return remaining as the deconvolution results
      double thisCoef;
      for(int cf=1; cf<coef_Vec.size();cf++){
        NumericVector ref_ints = candidate_Spec(_,cf);
        if(coef_Vec[cf] > 0){
          thisCoef = coef_Vec[cf];
          NumericVector cofed_ints = thisCoef*ref_ints;
          for(int exp_idx=0; exp_idx<exp_vec.size(); exp_idx++){
            exp_vec[exp_idx] = exp_vec[exp_idx] - cofed_ints[exp_idx]; // simply remove the contamination contribution
          }
        }
      }
      NumericVector new_exp_mz, new_exp_int;
      bool blv1;
      for(int exp_idx=0; exp_idx<exp_vec.size(); exp_idx++){
        blv1 = false;
        for(int j=0; j<exp_mtx.nrow(); j++){
          if(abs(exp_mtx(j,0) - mzs_vec[exp_idx])/(mzs_vec[exp_idx]) < 2e-6){
            blv1 = true;
          }
        }
        if((exp_vec[exp_idx]>0) & (blv1)){
          new_exp_mz.push_back(mzs_vec[exp_idx]);
          new_exp_int.push_back(exp_vec[exp_idx]);
        }
      }
      reconsctedSpe = cbind(new_exp_mz, new_exp_int);
      
    } else {
      // extract main spectrum and return 
      // A deconvoluted spectrum includes
      NumericVector new_exp_mz, new_exp_int;
      double residue = 0.0, main_ints, thisCoef, candidates_ints, final_ints, max_ints;
      bool blv1;
      // calculate residue
      for(int exp_idx=0; exp_idx<exp_vec.size(); exp_idx++) {
        candidates_ints = 0.0;
        for(int cf=0;cf<coef_Vec.size();cf++){
          thisCoef = coef_Vec[cf];
          candidates_ints = candidates_ints + candidate_Spec(exp_idx, cf)*thisCoef;
        }
        residue = residue + (exp_vec[exp_idx] - candidates_ints);
      }
      residue = residue / ((double) exp_vec.size());
      
      for(int exp_idx=0; exp_idx<exp_vec.size(); exp_idx++) {
        main_ints = 0.0;
        final_ints = 0.0;
        blv1 = false;
        // calculate main part
        main_ints = candidate_Spec(exp_idx, 0);
        main_ints = main_ints*coef_Vec[0];
        final_ints = main_ints + residue;
        for(int j=0; j<exp_mtx.nrow(); j++){
          if(abs(exp_mtx(j,0) - mzs_vec[exp_idx])/(mzs_vec[exp_idx]) < 2e-6){
            blv1 = true;
          }
        }
        if((final_ints>0.0) & (blv1)){
          new_exp_mz.push_back(mzs_vec[exp_idx]);
          new_exp_int.push_back(final_ints);
        }
      }
      // intensity data to be re-normalized
      max_ints = max(new_exp_int);
      new_exp_int = new_exp_int / max_ints;
      reconsctedSpe = cbind(new_exp_mz, new_exp_int);
    }
    
    if(reconsctedSpe.nrow() > exp_mtx.nrow()){
      stop("Error: Unexpected condition occurred when doing spectra reconstruction -> code: error2026DDA\n");
    }
    
    return reconsctedSpe;
  }
  
  List static predictRef_Spectrum(SqliteDriver SQLiteObj, 
                                  NumericMatrix thisSpectrum, 
                                  NumericMatrix precrsr_mtx, 
                                  double ppm_ms1,
                                  double ppm_ms2, 
                                  double rt_tol,
                                  bool useRT,
                                  int rule_opt,
                                  bool useEntropy) {
    // This function will be called only if there is zero or one contam Candidate member
    // This function we only consider non-isotope precursors
    // cout << "Starting to predict spectrum...\n";
    
    // Prepare spectrum/spectra data
    NumericVector thisclss = thisSpectrum(_,2);
    NumericVector uni_thisClss = unique(thisclss);
    //cout << "now this uni_thisClss is --> " << uni_thisClss << endl;
    vector<NumericMatrix> main_exp_spectrum(uni_thisClss.size());
    if(uni_thisClss.size() == 1){
      main_exp_spectrum[0] = thisSpectrum;
    } else {
      int n=0, min_bd, max_bd;
      for(double u : uni_thisClss){
        //cout << "n --> " << n << endl;
        IntegerVector idxVec = whichTrue(u == thisclss);
        min_bd = min(idxVec); // lower limit
        max_bd = max(idxVec); // upper limit
        NumericMatrix m2 = thisSpectrum(Range(min_bd,max_bd) ,_);
        main_exp_spectrum[n] = m2;
        n++;
      }
    }
    
    // Process data for extraction
    NumericVector iso_status = precrsr_mtx(_,3);
    NumericVector all_mzs = precrsr_mtx(_,0);
    NumericVector all_rts = precrsr_mtx(_,1);
    NumericVector all_ints = precrsr_mtx(_,2);
    
    IntegerVector nonIso_idx = whichTrue(iso_status == 0);
    all_mzs = all_mzs[nonIso_idx];
    all_rts = all_rts[nonIso_idx];
    all_ints = all_ints[nonIso_idx];
    double min_rt = min(all_rts);
    double max_rt = max(all_rts);
    
    // merge all mzs within ppm range
    double current_max_int = max(all_ints), current_mz, curr_simi = -1.0, best_simi = -1.0;
    vector <double> new_mzs;
    int idxx;
    bool bl1;
    while(all_mzs.size()>0){
      if(new_mzs.size() == 0){
        idxx = whichTrue1(all_ints == current_max_int);
        current_mz = all_mzs[idxx];
        new_mzs.push_back(current_mz);
      } else {
        bl1 = false;
        idxx = whichTrue1(all_ints == current_max_int);
        current_mz = all_mzs[idxx];
        for(double lmz : new_mzs){
          if(abs(lmz - current_mz)/current_mz < ppm_ms1*1e-6){
            bl1 = true;
            break;
          }
        }
        if(!bl1){
          // not included for this centroid (within the ppm)
          new_mzs.push_back(current_mz);
        }
      }
      // erase
      all_mzs.erase(idxx);
      all_ints.erase(idxx);
      // find a new maximum
      current_max_int = max(all_ints);
    }
    
    // run the data extraction from db
    all_mzs = new_mzs;
    vector<int> all_dirs;
    vector<double> all_rule_ms;
    vector<int> * all_rule_formulas;
    if(rule_opt == 0) {
      all_dirs = getAll_dir_change();
      all_rule_ms = getAll_ms_change();
      all_rule_formulas = getAll_formulas();
    } else if(rule_opt == 1) {
      // Bio-transformation
      all_dirs = getBio_dir_change();
      all_rule_ms = getBio_ms_change();
      all_rule_formulas = getBio_formulas();
    } else if(rule_opt == 2) {
      // adducts
      all_dirs = getAdc_dir_change();
      all_rule_ms = getAdc_ms_change();
      all_rule_formulas = getAdc_formulas();
    } else if(rule_opt == 3) {
      // fragments
      all_dirs = getFrg_dir_change();
      all_rule_ms = getFrg_ms_change();
      all_rule_formulas = getFrg_formulas();
    } else {
      all_dirs = getFrg_dir_change();
      all_rule_ms = getFrg_ms_change();
      all_rule_formulas = getFrg_formulas();
    }
    
    vector<double> predicted_prec_mzs;
    vector<int> all_rules_idx, all_mzs_idx, allIDs;
    double thismz, thisrule_mz, thisprdc_mz, ms_error;
    int thisdir, x;
    
    NumericMatrix ms2_ref_mtx, ms2_ref_mtx_best;
    vector<string> allMS2refs;
    List allMS2refs_List;
    IntegerVector formular_res;
    
    for(int j=0; j<all_mzs.size(); j++) {
      for(int i=0; i<all_dirs.size();i++) {
        
        thisdir = all_dirs[i];
        thismz = all_mzs[j];
        thisrule_mz = all_rule_ms[i];
        allIDs.clear();
        allMS2refs.clear();
        if(thisdir == 0){// 1 and -1
          // 1
          thisprdc_mz = thismz + thisrule_mz;
          ms_error = thisprdc_mz*ppm_ms1*1e-6;
          if(useRT){
            x = SQLiteObj.extractIDMS2_with_mzrtRange(thisprdc_mz - ms_error, thisprdc_mz + ms_error,
                                                      min_rt - rt_tol, max_rt + rt_tol);
          } else {
            x = SQLiteObj.extractIDMS2_with_mzRange_expDB(thisprdc_mz - ms_error, thisprdc_mz + ms_error);
          }
          
          allIDs = SQLiteObj.getIDsVec();
          allMS2refs = SQLiteObj.getMS2PeaksVec();
          
          if(allIDs.size()>0){
            all_rules_idx.push_back(i);
            all_mzs_idx.push_back(j);
            predicted_prec_mzs.push_back(thisprdc_mz);
            allMS2refs_List.push_back(allMS2refs);
          }
          
          // -1
          thisprdc_mz = thismz - thisrule_mz;
          if(thisprdc_mz < 50.0) {
            continue;
          } else {
            ms_error = thisprdc_mz*ppm_ms1*1e-6;
            if(useRT){
              x = SQLiteObj.extractIDMS2_with_mzrtRange(thisprdc_mz - ms_error, thisprdc_mz + ms_error,
                                                        min_rt - rt_tol, max_rt + rt_tol);
            } else {
              x = SQLiteObj.extractIDMS2_with_mzRange_expDB(thisprdc_mz - ms_error, thisprdc_mz + ms_error);
            }
            allIDs = SQLiteObj.getIDsVec();
            allMS2refs = SQLiteObj.getMS2PeaksVec();
            
            if(allIDs.size()>0){
              all_rules_idx.push_back(i);
              all_mzs_idx.push_back(j);
              predicted_prec_mzs.push_back(thisprdc_mz);
              allMS2refs_List.push_back(allMS2refs);
            }
          }
          
        } else {
          // 1 or -1
          thisprdc_mz = thismz + thisrule_mz*thisdir;
          if(thisprdc_mz < 50.0){
            continue;
          }
          ms_error = thisprdc_mz*ppm_ms1*1e-6;
          if(useRT){
            x = SQLiteObj.extractIDMS2_with_mzrtRange(thisprdc_mz - ms_error, thisprdc_mz + ms_error,
                                                      min_rt - rt_tol, max_rt + rt_tol);
          } else {
            x = SQLiteObj.extractIDMS2_with_mzRange_expDB(thisprdc_mz - ms_error, thisprdc_mz + ms_error);
          }
          allIDs = SQLiteObj.getIDsVec();
          allMS2refs = SQLiteObj.getMS2PeaksVec();
          
          if(allIDs.size()>0){
            all_rules_idx.push_back(i);
            all_mzs_idx.push_back(j);
            predicted_prec_mzs.push_back(thisprdc_mz);
            allMS2refs_List.push_back(allMS2refs);
          }
        }
      }
    }
    
    // cout << "Now size -- > " << predicted_prec_mzs.size() << " | -- > "
    //      << all_rules_idx.size() << " | --> " << all_mzs_idx.size() <<
    // " | allMS2refs_List-> " << allMS2refs_List.size() << endl;
    List contamCandidateList_clean;
    allMS2refs.clear();
    int this_idx0, this_idx1;
    
    if(predicted_prec_mzs.size() != 0){
      for(int k=0;k<allMS2refs_List.size();k++) {
        vector<string> allMS2refs = allMS2refs_List[k];
        // cout << "vector<string> allMS2refs--> " << allMS2refs.size() << endl;
        for(int r=0; r<allMS2refs.size();r++){
          string ms2_ref_str = allMS2refs[r];
          ms2_ref_mtx = ms2peak_parse(ms2_ref_str);
          for(int p=0; p<main_exp_spectrum.size();p++){
            if(useEntropy){
              curr_simi = entropySimilarity(main_exp_spectrum[p], ms2_ref_mtx, ppm_ms2);
            } else {
              curr_simi = spectrumSimilarity(main_exp_spectrum[p], ms2_ref_mtx, ppm_ms2);
            }
            
            if(curr_simi > best_simi){
              best_simi = curr_simi;
              ms2_ref_mtx_best = ms2_ref_mtx;
              this_idx0 = k;
              this_idx1 = all_mzs_idx[k];
            }
          }
        }
      }
      
      // already found most similar spectrum --> starting generating a pseudo spectrum
      int rule_idx;
      vector<int> this_formula, this_sub_formula;
      List sub_formula_list;
      NumericVector all_ref_mzs, all_ref_ints, new_ref_mzs, new_ref_ints;
      NumericMatrix new_ref_mtx;
      if(best_simi>0.001){
        // cout << best_simi << " <--- best_simi | this_idx ---> " << this_idx0 << " | --> " << ms2_ref_mtx_best.nrow() << endl;
        // vector<string> allMS2refs = allMS2refs_List[this_idx0];
        rule_idx = all_rules_idx[this_idx0];
        //this_formula = all_rule_formulas[rule_idx];
        all_ref_mzs = ms2_ref_mtx_best(_,0);
        all_ref_ints = ms2_ref_mtx_best(_,1);
        //cout << "all_ref_mzs --> " << all_ref_mzs.size() << endl;
        //cout << "all_mzs[this_idx1] --> " << all_mzs[this_idx1] << endl;
        this_formula = assignFullFormula_accurate(all_mzs[this_idx1], ppm_ms1);
        //cout << "C"<< this_formula[0] << "H" << this_formula[1] << "O" << this_formula[2]<< "N" << this_formula[3] << "P" << this_formula[4] <<"S" << this_formula[5] << endl;
        
        double ms2mz; bool keep;
        IntegerVector keepIdx;
        for(int m=0; m<all_ref_mzs.size(); m++){
          keep=false;
          ms2mz = all_ref_mzs[m];
          sub_formula_list = assignSubFormula_all(ms2mz, ppm_ms2);
          //cout << "sub_formula_list size --> " << sub_formula_list.size() << endl;
          for(int n=0; n<sub_formula_list.size();n++){
            //cout << "line 2451 <-- \n";
            vector<int> this_sub_formula = sub_formula_list[n];
            if((this_sub_formula[0] <= this_formula[0]) &
               (this_sub_formula[1] <= this_formula[1]) &
               (this_sub_formula[2] <= this_formula[2]) &
               (this_sub_formula[3] <= this_formula[3]) &
               (this_sub_formula[4] <= this_formula[4]) &
               (this_sub_formula[5] <= this_formula[5])){
              //cout << "FOUND good fragment -> " << m << endl;
              keep=true;
              break;
            }
          }
          if(keep){
            keepIdx.push_back(m);
          }
          
          // if((this_sub_formula[0] > this_formula[0])){
          //   cout << "FOUND bad fragment C -> " << m << endl;
          // }
          // if((this_sub_formula[1] > this_formula[1])){
          //   cout << "FOUND bad fragment H -> " << m << endl;
          // }
          // if((this_sub_formula[2] > this_formula[2])){
          //   cout << "FOUND bad fragment O -> " << m << endl;
          // }
          // if((this_sub_formula[3] > this_formula[3])){
          //   cout << "FOUND bad fragment N -> " << m << endl;
          // }
          // if((this_sub_formula[4] > this_formula[4])){
          //   cout << "FOUND bad fragment P -> " << m << endl;
          // }
          // if((this_sub_formula[5] > this_formula[5])){
          //   cout << "FOUND bad fragment S -> " << m << endl;
          // }
        }
        //cout << "Useful frgs num -> " << keepIdx.size() << endl;
        
        if(keepIdx.size() > 0){
          new_ref_mzs = all_ref_mzs[keepIdx];
          new_ref_ints = all_ref_ints[keepIdx];
          new_ref_mtx = cbind(new_ref_mzs, new_ref_ints);
        }
        
        contamCandidateList_clean.push_back(new_ref_mtx);
      }
    }
    
    return contamCandidateList_clean;
  }
  
  NumericMatrix static extractSingleSpectrum(NumericMatrix thisSpectra, int idx){
    NumericMatrix resx;
    NumericVector thisclss = thisSpectra(_,2);
    NumericVector uni_thisClss = unique(thisclss);
    int n=0, min_bd, max_bd;
    for(double u : uni_thisClss){
      if(u == (double) idx){
        IntegerVector idxVec = whichTrue(u == thisclss);
        min_bd = min(idxVec); // lower limit
        max_bd = max(idxVec); // upper limit
        NumericMatrix res = thisSpectra(Range(min_bd,max_bd) ,_);
        resx = res;
        break;
      }
    }
    return resx;
  }
  
  //// no need to consider isobaric ions [already seperated by HRMS]
  // int isobaricDeco(){
  //   // this function is designed to deconvolve the isobarics based on the RT difference 
  //   // maybe a little bit eariler or later compared to the contaminated (main, chimeric) ion.
  //   // At the same time, the m/z is also different.
  //   
  //   return 1;
  // }
  
public:
  
  void setDDA_Arguments(NumericMatrix pm, 
                        NumericVector scant1,
                        NumericVector scant2,
                        List scanms1,
                        List scanms2,
                        NumericMatrix prec_mzs,
                        double win_size,
                        double ppm1,
                        double ppm2,
                        double sn,
                        double filt,
                        bool show_Output,
                        int thread_i,
                        string fileNM,
                        bool useEntpy,
                        double rt_s = 0,
                        double inclu_thres = 1e4,
                        int ionmode = 0,
                        string db_path = "") {
    
    /*
     * 
     Explanation of all arguments in this function
     1. peak_matrix : pm, targeted peak matrix of MS1 level, three columns 
     (m/z, rt_min, rt_max) or four columns (mz_min, mz_max, rt_min, rt_max);
     2. scanrt_ms1: scant1, RT of all MS1 scans (rt information, sec);
     3. scanrt_ms2: scant2, RT of all MS2 scans (rt information, sec);
     4. scan_ms1: scanms1, all MS1 scans data (two columns: mz + intensity);
     5. scan_ms2: scanms2, all MS2 scans data (two columns: mz + intensity);
     6. prec_mzs, precursors_mzs, all mz values of precursors;
     7. isolation_window_size: win_size, the windows size of MS instrument;
     8. ppm_ms1: ppm for ms1 processing
     9. ppm_ms2: ppm for ms2 processing
     10. sn_thres: signal to noise threshold;
     11. filter_thres: threshold for ms1 peak filtration;
     12. rt_size, RT windows size, seconds
     13. inclu_thres: threshold of intensity inclusion threshold
     14. ion_mode: ESI mode: 0, is negative; 1, is positive.
     15. useEntropy: use spectral entropy similarity method (true) or dot-product (false);
     */
    peak_matrix = pm;
    scanrt_ms1 = scant1;
    scanrt_ms2 = scant2;
    scan_ms1 = scanms1;
    scan_ms2 = scanms2;
    precursors_mzs = prec_mzs;
    ppm_ms1 = ppm1;
    ppm_ms2 = ppm2;
    sn_thres = sn;
    isolation_window_size = win_size;
    filter_thres = filt;
    file_name = fileNM;
    thread_int = thread_i;
    rt_size = rt_s;
    inclusion_inten_thre = inclu_thres;
    ion_mode = ionmode;
    database_path = db_path;
    showOutput = show_Output;
    useEntropy = useEntpy;
    if(showOutput){
      cout << "Arguments configuration done!" << "\n";
    }
  }
  
  void PerformDDAProcess_core(bool decoOn, bool showOutput) {
    if(showOutput){
      cout << "PerformDDAProcess starting..." << "\n";
    }
    int res_stp0 = formatPeakMatrix();
    int res_stp1 = precursorsGrouping();
    int res_stp2 = 0;
    if(decoOn){
      if(showOutput){
        //cout << "Deconvolution is going to be executed .. \n";
      }
      res_stp2 = chimericSpectraDetection();
    }
    
    // If need to run deconvolution
    int res_stp3, res_stp4, res_stp5;
    if(!decoOn) {
      MS2listSummarize_noDeco();
      //MS2listSummarize_noMerge();
    } else if(ContmGroupIdxVec.size() == 0){
      if(showOutput){
        //cout << "No convoluted MS2 peaks found, skipped deconvolution" << endl;
      }
      MS2listSummarize_noDeco();
    } else {
      findContaminationList();
      List decoSpecRes = performDeco();
      tmpmtx = decoSpecRes;
      //cout << "Deconvolution of this thread executed successfully\n";
      MS2listSummarize();
    }
    
    // cout << "============== Status ============== " << endl;
    // if(res_stp0 == 1){
    //   cout << "Format Peaks Matrix Succeeded!" << endl;
    // }
    // if(res_stp1 == 1){
    //   cout << "Precursors Grouping Succeeded!" << endl;
    // }
    // if(res_stp2 == 1){
    //   cout << "Chimeras Detection Succeeded!" << endl;
    // }

  }

  // Results Getters
  List getResultsList(){
    return results;
  }
  IntegerVector getCleanGroup(){
    return CleanGroupIdxVec;
  }
  IntegerVector getContmGroup(){
    return ContmGroupIdxVec;
  }
  List getconvolvedPrcList(){
    return convolved_prec_List;
  }
  
  void setUseRT(bool use_rt){
    useRT = use_rt;
  }
  IntegerVector getFeatureIdxVec(){
    return featureIdxVec;
  }
  
  List gettmpmtx(){
    return tmpmtx;
  }
  
  NumericMatrix get_peak_matrix(){
    return peak_matrix;
  }
  
};

List PerformDDA_main(NumericMatrix pm, 
                     NumericVector scant1,
                     NumericVector scant2,
                     List scanms1,
                     List scanms2,
                     NumericMatrix prec_mzs,
                     double win_size,
                     double ppm1,
                     double ppm2,
                     double sn,
                     double filt,
                     double intensity_thresh,
                     int ionmode,
                     string db_path, 
                     bool decoOn,
                     bool useEntropy,
                     bool showOutput,
                     int thread_num,
                     string file_nm);


#endif