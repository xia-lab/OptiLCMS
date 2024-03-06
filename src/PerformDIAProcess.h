#ifndef DIA_H
#define DIA_H

#include <RcppArmadillo.h>
#include "DecoSpectra.h"
#include "utilities.h"

using namespace Rcpp;
using namespace std;


class diaprocess {
  
private:
  List peak_list;
  NumericMatrix rt_info;
  NumericMatrix swath_info;
  NumericVector scantimes_ms1;
  NumericVector scantimes_ms2;
  List scan_ms1;
  List scan_ms2;
  double peakwidth_min;
  double ppm_ms2;
  double sn_thres;
  double smooth_span;
  double filter_thres;
  
  List results;
  
public:
  
  void setDIA_Arguments(List pm, 
                        NumericMatrix swath,
                        NumericVector scanrt1,
                        NumericVector scanrt2,
                        List scanms1,
                        List scanms2,
                        double pkw_min,
                        double ppm2,
                        double sn,
                        double sm_span,
                        double filt) {
    
    /*
     * 
     Explanation of all arguments in this function
     1. peak_matrix : the peak grouping matrix of MS1 level;
     2. rt_info: retention time information of all MS1 peaks;
     3. swath_info: SWATH design information;
     4. scantimes_ms1: RT of all MS1 scans;
     5. scantimes_ms2: RT of all MS2 scans;
     6. scan_ms1: all MS1 scans data (two columns: mz + intensity);
     7. scan_ms2: all MS2 scans data (two columns: mz + intensity);
     8. peakwidth_min: minimum peakwidth;
     9. ppm_ms2: ppm for ms2 processing
     10. sn_thres: signal to noise threshold;
     11. smooth_span: span value for smoothing;
     12. filter_thres: threshold for ms1 peak filtration
     */
    peak_list = pm;
    swath_info = swath;
    scantimes_ms1 = scanrt1;
    scantimes_ms2 = scanrt2;
    scan_ms1 = scanms1;
    scan_ms2 = scanms2;
    peakwidth_min = pkw_min;
    ppm_ms2 = ppm2;
    sn_thres = sn;
    smooth_span = sm_span;
    filter_thres = filt;
    cout << "Arguments configuration done!" << "\n";
  }
  
  void PerformDIAProcess_core() {
    cout << "PerformDIAProcess starting..." << "\n";
    cout << "now, peak width min is -->" << peakwidth_min << "\n";
    cout << "Total number of MS1 features is " << peak_list.size() << endl;
    cout << scan_ms1.size() << " MS1 scans and " << scan_ms2.size() << " MS2 scans are imported!\n";
    // 1) prepare all parameters (swath info and index of scans)
    IntegerVector swath_idx_ms2(scan_ms2.size()); //swath_idx_ms2 is corresponding to the rows of "swath_info"
    double tl, tu; // scan time bounds (lower and upper);
    int swath_idx, k_count = 0;
    if(scantimes_ms1[0] > scantimes_ms2[0]){
      // for the case starting with MS2 acquisition
      tl = 0;
      tu = scantimes_ms1[0];
      for(int s=k_count; s<scan_ms2.size(); s++){
        if((tl<scantimes_ms2[s]) & (tu>scantimes_ms2[s])){
          swath_idx_ms2[s] = swath_idx;
          swath_idx++;
          k_count++;
        } else {
          break;
        }
      }
    }
    for(int t=1; t<scan_ms1.size(); t++){
      tl = scantimes_ms1[t-1];
      tu = scantimes_ms1[t];
      swath_idx = 0;
      for(int s=k_count; s<scan_ms2.size(); s++){
        if((tl<scantimes_ms2[s]) & (tu>scantimes_ms2[s])){
          swath_idx_ms2[s] = swath_idx;
          swath_idx++;
          k_count++;
        } else {
          break;
        }
      }
    }
    
    if(k_count < scantimes_ms2.size()){
      // there are still some scans left not indexed
      swath_idx = 0;
      for(int s=k_count; s<scan_ms2.size(); s++){
        swath_idx_ms2[s] = swath_idx;
        swath_idx++;
      }
    }

    // 2) start deconvolution for all features
    int idx_swath;
    double this_rt, this_rt_min, this_rt_max, rt_expansion, this_rt_min_ext, this_rt_max_ext;
    double this_mz, this_mz_min, this_mz_max; // no need to expand mz, because this has already been a range.
    double baseline_val, maxo_val, sn_val;
    vector<NumericMatrix> deco_results(peak_list.size());
    cout << "===== Start deconvolution. Progress shown below ===== \n";
    for(int f=0; f<peak_list.size() ; f++){ // 477 | f<peak_list.size() 6198 // this for-loop can be parallel, later
      NumericVector this_ft_list = peak_list[f];
      this_mz = this_ft_list[0];
      this_mz_min = this_ft_list[1];
      this_mz_max = this_ft_list[2];
      this_rt = this_ft_list[3];
      this_rt_min = this_ft_list[4];
      this_rt_max = this_ft_list[5];
      baseline_val = this_ft_list[6];
      // maxo_val = this_ft_list[8];
      // sn_val = this_ft_list[9];
      
      rt_expansion = (this_rt_max - this_rt_min)*2.5;
      if(rt_expansion < 1){rt_expansion = 2.5;}
      this_rt_min_ext = this_rt - rt_expansion;
      if(this_rt_min_ext <0){this_rt_min_ext = 0;}
      this_rt_max_ext = this_rt + rt_expansion;
      
      for(int sw=0; sw<swath_info.nrow();sw++){
        if((this_mz >= swath_info(sw,0)) & (this_mz <= swath_info(sw,1))){
          idx_swath = sw;
          break;
        }
      }
      // baseline_val = maxo_val/sn_val;
      
      // 2.1) prepare MS1 data and EICs
      //cout << "this_rt_min_ext ---> " << this_rt_min_ext << endl;
      //cout << "this_rt_max_ext ---> " << this_rt_max_ext << endl;
      IntegerVector idx_ms1_ext = whichTrue((scantimes_ms1 >= this_rt_min_ext) & (scantimes_ms1 <= this_rt_max_ext));
      //cout << "idx_ms1_ext size ---> " << idx_ms1_ext.size() << endl;
      IntegerVector idx_ms1 = whichTrue((scantimes_ms1 >= this_rt_min) & (scantimes_ms1 <= this_rt_max));
      IntegerVector is_keep_idx(idx_ms1.size());
      for(int m=0; m<idx_ms1.size(); m++){
        is_keep_idx[m] = whichTrue1(idx_ms1[m] == idx_ms1_ext);
      }
      
      List scan_ms1_this = scan_ms1[idx_ms1_ext];
      NumericVector tmpVec = {this_mz_min, this_mz_max};
      tmpVec.attr("dim") = Dimension(1, 2);
      NumericMatrix ms1_mz_range = as<NumericMatrix>(tmpVec);
      NumericVector mzVec1 = {this_mz};
      
      List ms1_eic0 = extractEIC(scan_ms1_this, ms1_mz_range, mzVec1);
      NumericMatrix ms1_eic = ms1_eic0[0];
      
      NumericVector mzVecs = ms1_eic(_,0);
      NumericVector intVecs = ms1_eic(_,1);
      NumericVector rtVecs = scantimes_ms1[idx_ms1_ext];
      // cout << "rtVecs size  ==> " << rtVecs.size() << endl;
      // cout << "mzVecs size  ==> " << mzVecs.size() << endl;
      // cout << "intVecs size ==> " << intVecs.size() << endl;
      NumericMatrix peaks_ms1_ext = cbind((NumericVector)idx_ms1_ext, rtVecs, mzVecs, intVecs);
      NumericVector peak_ms1_ext_smooth = SmoothLoess(peaks_ms1_ext, 0.1);
      NumericVector pk_ms1_sm_int = peak_ms1_ext_smooth[is_keep_idx];
      
      if(is_true(all(pk_ms1_sm_int == 0))){
        // need to give a value to the final result set : todo later
        continue;
      }

      int idx_apex_ms1 = whichTrue1(peak_ms1_ext_smooth == max(pk_ms1_sm_int));
      double idx_rt_apex_ms1 = rtVecs[idx_apex_ms1];
      // cout << "pk_ms1_sm_int       ===> " << pk_ms1_sm_int << endl;
      // cout << "idx_apex_ms1        ===> " << idx_apex_ms1 << endl;
      // cout << "peak_ms1_ext_smooth ===> " << peak_ms1_ext_smooth << endl;
      // cout << "idx_rt_apex_ms1     ===> " << idx_rt_apex_ms1 << endl;
      // cout << "this_rt_min_ext     ===>" << this_rt_min_ext << endl;
      // cout << "this_rt_max_ext     ===>" << this_rt_max_ext << endl;
      // cout << "idx_swath           ===>" << idx_swath << endl;
      // cout << "swath_idx_ms2 size  ===>" << swath_idx_ms2.size() << endl;
      //cout << "swath_idx_ms2==>" << swath_idx_ms2 << endl;
      // 2.2) prepare MS2 data and EICs
      IntegerVector idx_ms2_ext_all = whichTrue((scantimes_ms2 >= this_rt_min_ext) & 
                                                (scantimes_ms2 <= this_rt_max_ext) & 
                                                (swath_idx_ms2 == idx_swath));
      IntegerVector idx_ms2_all = whichTrue((scantimes_ms2 >= this_rt_min) & 
                                            (scantimes_ms2 <= this_rt_max) & 
                                            (swath_idx_ms2 == idx_swath));
      //cout << "line 192 --> " << idx_ms2_ext_all << "\n";
      List spec_exp_ext_ms2 = scan_ms2[idx_ms2_ext_all];
      
      NumericMatrix spec_apex_ms2;
      
      if(idx_apex_ms1<spec_exp_ext_ms2.size()){
        NumericMatrix spec_apex_ms20 = spec_exp_ext_ms2[idx_apex_ms1];
        spec_apex_ms2 = spec_apex_ms20;
      } else if(idx_apex_ms1 !=0){
        // cout << "line 201 --> " << idx_apex_ms1 << "\n";
        // cout << "line 203 --> " << spec_exp_ext_ms2.size() << "\n";
        NumericMatrix spec_apex_ms20 = spec_exp_ext_ms2[idx_apex_ms1-1];
        spec_apex_ms2 = spec_apex_ms20;
      } else{
        NumericMatrix spec_apex_ms20 = spec_exp_ext_ms2[0];
        spec_apex_ms2 = spec_apex_ms20;
      }

      //cout << "spec_apex_ms2 nrow -> " << spec_apex_ms2.nrow() << endl;
      int m = 0;
      for(int n=0; n<spec_apex_ms2.nrow(); n++){
        if((spec_apex_ms2(n,1) > filter_thres) & 
           (spec_apex_ms2(n,0) <= this_mz)){
          m++;
        }
      }
      NumericVector mzs4ms2(m);
      m = 0;
      for(int n=0; n<spec_apex_ms2.nrow(); n++){
        if((spec_apex_ms2(n,1) > filter_thres) & 
           (spec_apex_ms2(n,0) <= this_mz)){
          mzs4ms2[m] = spec_apex_ms2(n,0);
          m++;
        }
      }
      
      NumericMatrix ms2_mz_range (mzs4ms2.size(), 2);
      for(int n=0; n<ms2_mz_range.nrow(); n++){
        ms2_mz_range(n,0) = mzs4ms2[n] - mzs4ms2[n]*ppm_ms2*1e-6;
        ms2_mz_range(n,1) = mzs4ms2[n] + mzs4ms2[n]*ppm_ms2*1e-6;
      }
      
      List ms2_eic = extractEIC(spec_exp_ext_ms2, ms2_mz_range, mzs4ms2);

      for(int n=0; n<ms2_eic.size(); n++){
        NumericMatrix num_lst_this_ms2 = ms2_eic[n];
        NumericVector mz_ms2_val = num_lst_this_ms2(_,0);
        NumericVector int_ms2_val = num_lst_this_ms2(_,1);
        NumericVector sc_ms2_times = scantimes_ms2[idx_ms2_ext_all];
        ms2_eic[n] = cbind((NumericVector)idx_ms2_ext_all, sc_ms2_times, mz_ms2_val, int_ms2_val);
      }
      
      NumericVector info_peak_ms1 = {(double)is_keep_idx[0], 
                                     (double)idx_apex_ms1,
                                     (double)idx_rt_apex_ms1,
                                     (double)is_keep_idx[is_keep_idx.size()-1],
                                     baseline_val,
                                     0.0};
      
      // 2.3) run the decoSpectra
      //cout << "DecoSpectra --> " << f << "\n";
      List spec_decon = DecoSpectra(f,
                                    ms2_eic,
                                    pk_ms1_sm_int,
                                    idx_ms2_ext_all.size(),
                                    idx_apex_ms1,
                                    info_peak_ms1,
                                    peakwidth_min,
                                    sn_thres,
                                    false);
      //results = spec_decon;
      NumericMatrix spec_mtx_res = spec_decon[0];
      deco_results[f] = spec_mtx_res;
      //cout << " -> done! " << endl;
      cout << ".";
    }
    cout << endl;
    // 3) summarize results for consensus
    NumericMatrix tmp_mtx4cout;
    int icc = 0;
    for(int ic=0; ic<peak_list.size(); ic++){
      tmp_mtx4cout = deco_results[ic];
      if(tmp_mtx4cout.nrow()>0){
        icc++;
      }
    }
    
    // Initialize results variables
    // Spectra List -> contains a lot of lists ( number is icc here) : ms/ms spectrum NumericMatrix
    // Indicator List -> contains a lost of lists ( number is icc here) : indicator to show they spectrum are "deconvoluted" (all 1 [deconvoluted] here)
    List specList(icc), indicList(icc);
    vector<int> ftVec(icc);

    vector<NumericMatrix> resMtx(1);
    vector<int> indVec(1);
    indVec[0] = 1;
    int jcc = 0;
    for(int i=0; i<peak_list.size(); i++){
      //resMtx.clear();
      NumericMatrix tmp_mtx4return = deco_results[i];
      if(tmp_mtx4return.nrow() > 0){
        resMtx[0] = tmp_mtx4return;
        specList[jcc] = resMtx;
        indicList[jcc] = indVec;
        ftVec[jcc] = i;
        jcc++;
      }
    }

    results = List::create(Named("Spectra") = specList,
                           Named("Indicator") = indicList,
                           Named("FeatureIdx") = ftVec);
  }
  
  
  List getResults(){
    return results;
  }
  
  static int findMaxIdx2(double d[], int len) {
    int iMin = 0;
    for (int i = 1; i < len; i++){
      if (d[iMin] < d[i]) {
        iMin = i;
      }
    }
    return iMin;
  }
  
  static int count_increaser(int & x){
    return x + 1;
  }
  
  static List extractEIC(List specExp, NumericMatrix mzRange, NumericVector mz) {
    int lenMZ = mz.length();
    int lenSpec = specExp.size();
    List result(lenMZ);
    //arma::vec Pos = {};
    int iInRange = 0;
    for (int iMZ=0; iMZ<lenMZ; iMZ++) {
      NumericVector mzSpec(lenSpec), intSpec(lenSpec);
      for (int iSpec = 0; iSpec < lenSpec; iSpec++) {
        
        NumericMatrix spec = specExp[iSpec];
        double mz1 = mz(iMZ);
        double mzStart = mzRange(iMZ, 0);
        double mzEnd = mzRange(iMZ, 1);
        IntegerVector idxInRange;
        int nr = spec.nrow();
        iInRange = 0;
        for (int iRec = 0; iRec<nr; iRec++) {
          bool bl1 = (spec(iRec, 0) >= mzStart);
          bool bl2 = (spec(iRec, 0) <= mzEnd);
          if (bl1 & bl2) {
            idxInRange.push_back(iRec); 
            iInRange++;
          }
        }
        
        if (iInRange > 1) {
          double intInRange[iInRange];
          for (int i = 0; i < iInRange; i++) {
            intInRange[i] = spec(idxInRange[i], 1);
          }
          int iMax = findMaxIdx2(intInRange, iInRange);
          mzSpec(iSpec) = spec(idxInRange[iMax], 0);
          intSpec(iSpec) = spec(idxInRange[iMax], 1);
        } else if (iInRange == 1){
          mzSpec(iSpec) = spec(idxInRange[0], 0);
          intSpec(iSpec) = spec(idxInRange[0], 1);
        } else {
          mzSpec(iSpec) = mz1;
          intSpec(iSpec) = 0;
        }
        
      }
      
      // result[iMZ] = List::create(Rcpp::Named("mz") = mzSpec,
      //                            Rcpp::Named("intensity") = intSpec);
      NumericMatrix mtx = cbind(mzSpec, intSpec);
      result[iMZ] = mtx;
    }
    
    return result;
  }
  
  
};


List PerformDIA_main(List pm, 
                     NumericMatrix swath,
                     NumericVector scanrt1,
                     NumericVector scanrt2,
                     List scanms1,
                     List scanms2,
                     double pkw_min,
                     double ppm2,
                     double sn,
                     double sm_span,
                     double filt);
#endif
