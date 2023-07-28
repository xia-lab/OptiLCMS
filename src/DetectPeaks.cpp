#include "DetectPeaks.h"

List DetectPeaks(NumericMatrix eic, double peakwidth_min, int num_scantime, int idx_apex_eic,
                 double snthr, bool is_smooth, int n_skip_max) {
  
  Rcpp::List DetectedPeaks;
  
  bool is_peak = false;
  NumericVector eic_ints, eic_int;
  NumericVector eic_raw_int = eic( _ , 3);
  if(is_smooth){
    eic_int = SmoothLoess(eic, 0.1);
    eic_ints = clone(eic_int);
  } else {
    eic_int = eic_raw_int;
  }
  //cout << "eic_raw_int --> " << eic_raw_int << endl;
  NumericVector is_roi = getContinuousPtsAboveThrIdx(eic_raw_int, 0, peakwidth_min, 0.0, 1);
  //cout << "ios_roi --> " << is_roi << endl;
  IntegerVector idx_fr_roi = GetRoi(is_roi, idx_apex_eic);
  //cout << "idx_fr_roi --> " << idx_fr_roi << endl;
  if(idx_fr_roi.size() == 0)
    return DetectedPeaks;
  
  NumericVector eic_idxs = eic( _ , 0);
  IntegerVector eic_idx = as<IntegerVector>(eic_idxs);
  //cout << "eic_idxs --> " << eic_idxs << endl;
  //cout << "eic_idx --> " << eic_idx << endl;
  //cout << "eic_raw_int 1 --> " << eic_raw_int << endl;
  NumericVector x_tmp = clone(eic_raw_int);
  //cout << "x_tmp --> " << x_tmp << endl;
  double noise = EstimateChromNoise(x_tmp, 0.05, 12);
  //cout << "Arriving here --> 1 \n";
  // cout << "noise --> " << noise << endl;
  NumericVector noise_local = GetLocalNoiseEstimate(eic_raw_int, 
                                                    idx_fr_roi, 
                                                    12, 42, 
                                                    num_scantime, 
                                                    noise, 4);
  
  // cout << "noise_local --> " << noise_local << endl;

  NumericVector tmpVec = {noise, noise_local[0]};
  NumericVector baselineVec = {1, min(tmpVec)};
  double baseline = max(baselineVec);
  NumericVector sdnoiseVec = {1, noise_local[1]};
  double sdnoise = max(sdnoiseVec);
  //cout << "sdnoise --> " << sdnoise << endl;
  NumericVector v =
    NumericVector::create(sdnoise);
  LogicalVector l1 = is_na(v);
  if(l1[0]) sdnoise = 0;
  double intthr = sdnoise*snthr;
  //cout << "baseline --> " << baseline << endl;
  NumericVector not_noise = getContinuousPtsAboveThrIdx(eic_raw_int, 0, 4, baseline, 1);
  IntegerVector idx_not_noise = GetRoi(not_noise, idx_apex_eic);
  //cout << "not_noise --> " << not_noise << endl;
  //cout << "idx_apex_eic --> " << idx_apex_eic << endl;

  // cout << "idx_not_noise --> " << idx_not_noise << endl;
  if(idx_not_noise.size() == 0) return DetectedPeaks;
  // cout << "Arriving here --> 2.5 \n";
  IntegerVector idx_fr_roi_skip  = clone(idx_fr_roi);
  
  //cout << "idx_not_noise --> " << idx_not_noise << endl;
  //cout << " idx_fr_roi_skip --> " << idx_fr_roi_skip << endl;
  IntegerVector idx_fr_roi_denoise;
  for(int i = 0; i < idx_fr_roi_skip.size(); i++){
    LogicalVector tmpLog = idx_fr_roi_skip[i] == idx_not_noise;
    bool dt = is_true(any(tmpLog));
    if(dt){
     idx_fr_roi_denoise.push_back(idx_fr_roi_skip[i]);
    }
    tmpLog.erase(tmpLog.begin(), tmpLog.end());
  }
  // cout << " idx_fr_roi_denoise --> " << idx_fr_roi_denoise << endl;
  if(idx_fr_roi_denoise.size() == 0){
    return DetectedPeaks;
  }
  // cout << "Arriving here --> 3 \n";
  // cout << "eic_raw_int 2 --> " << eic_raw_int << endl;
  // cout << "baseline 2 --> " << baseline << endl;
  // cout << "intthr 2 --> " << intthr << endl;
  
  bool nonabovethre = is_false(any(eic_raw_int - baseline >= intthr));
  // cout << "nonabovethre ---> " << nonabovethre << endl;
  
  int lb, rb;
  IntegerVector idx_apex_frroi;
  List info_pk;
  
  if(nonabovethre) {
    lb = 0;
    rb = eic_raw_int.size() - 1;
    idx_apex_frroi = idx_apex_eic;
    
    info_pk = List::create(_["lb"] = lb, 
                           _["apex"] = idx_apex_frroi,
                           _["rb"] = rb, 
                           _["bl"] = baseline,
                           _["is_peak"] = is_peak);
  } else {
    NumericVector eicss = eic_ints[idx_fr_roi_denoise];
    IntegerVector idx_apex_roi = FindLocalMax(eicss, 3, baseline);
    IntegerVector idx_bnd_roi = FindLocalMin(eicss, 3);

    if(idx_apex_roi.size() == 0){
      lb = 0;
      rb = eic_raw_int.size() - 1;
      idx_apex_frroi = idx_apex_eic;
      
      info_pk = List::create(_["lb"] = lb, 
                             _["apex"] = idx_apex_frroi,
                             _["rb"] = rb, 
                             _["bl"] = baseline,
                             _["is_peak"] = is_peak);
      
    } else {
      if(idx_bnd_roi.size() == 0){
        IntegerVector idx_bnd_frroi = {min(idx_fr_roi_denoise), max(idx_fr_roi_denoise)};;
      }

      idx_apex_frroi = idx_fr_roi_denoise[idx_apex_roi];
      IntegerVector idx_bnd_ffroi = idx_fr_roi_denoise[idx_bnd_roi];
      
      // cout << "idx_fr_roi_denoise --> " << idx_fr_roi_denoise << endl;
      // cout << "idx_apex_frroi --> " << idx_apex_frroi << endl;
      // cout << "idx_bnd_ffroi -->" << idx_bnd_ffroi << endl;
      
      idx_bnd_ffroi.push_back(min(idx_fr_roi_denoise));
      idx_bnd_ffroi.push_back(max(idx_fr_roi_denoise));
      IntegerVector idx_bnd_can = sort_unique(idx_bnd_ffroi);
      //cout << "idx_bnd_can --> " << idx_bnd_can << endl;
      //cout << "Arriving here --> 3.5 \n";
      IntegerVector idx_diff, s_idx_diff, idx_sign_change, lbVec, rbVec, sharpness, i_rm;
      for(int idx = 0; idx < idx_apex_frroi.size(); idx++){
        idx_diff = idx_bnd_can - idx_apex_frroi[idx];
        s_idx_diff = sign(idx_diff);
        //cout << "Arriving here ------------>  x "<< idx << " <----\n";
        //cout << "\nidx_diff --> " << idx_diff << endl;
        //cout << "\ns_idx_diff --> " << s_idx_diff << endl;
        IntegerVector idx_sign_change = whichTrue(diff(s_idx_diff) > 0);
        IntegerVector lb_idx0 = idx_bnd_can[idx_sign_change];
        //cout << "lb_idx0 ---> " << lb_idx0 << endl;
        lbVec.push_back(lb_idx0[0]);
        //cout << idx << " <--- lbVec ---> " << lbVec << endl;

        if(idx_sign_change[0] == (idx_bnd_can.size() -1)){
          IntegerVector rb_idx0 = idx_bnd_can[idx_sign_change];
          //cout << "rb_idx0 -if--> " << rb_idx0 << endl;
          rbVec.push_back(rb_idx0[0]);
        } else {
          IntegerVector rb_idx0 = idx_bnd_can[idx_sign_change + 1];
          //cout << "rb_idx0 -else--> " << rb_idx0 << endl;
          rbVec.push_back(rb_idx0[0]);
        }

        if(idx > 0){
          if(lbVec[idx] < idx_apex_frroi[idx- 1]){
            //cout << "\nidx_apex_frroi[idx - 1] --> " << idx_apex_frroi[idx - 1] << "\n";
            //cout << "idx_apex_frroi[idx] --> " << idx_apex_frroi[idx] << "\n";
            int loval0 = idx_apex_frroi[idx-1];
            int hival0 = idx_apex_frroi[idx];
            if(loval0 >  hival0){
              loval0 = idx_apex_frroi[idx];
              hival0 = idx_apex_frroi[idx-1];
            }
            //cout << "loval 0--> " << loval0 << "\n";
            //cout << "hival 0--> " << hival0 << "\n";
            IntegerVector idx_vec1 = seq(loval0, hival0);
            NumericVector eics_apex_roi = eic_int[idx_vec1];
            lbVec[idx] = idx_apex_frroi[idx - 1] + which_min(eics_apex_roi);
            // cout << idx << " <---- idx_apex_frroi[idx - 1] ---> " <<idx_apex_frroi[idx - 1] << endl;
            // cout << idx << " <---- which_min(eics_apex_roi) ---> " << which_min(eics_apex_roi) << endl;
            // cout << idx << " <---- lbVec[idx] ---> " << lbVec[idx] << endl;
          }
        }
        //cout << "idx_apex_frroi.size --> " << idx_apex_frroi.size() << endl;
        if(idx < idx_apex_frroi.size() - 1){
          if(rbVec[idx] > idx_apex_frroi[idx+1]){
            // cout << "\nidx_apex_frroi[idx] --> " << idx_apex_frroi[idx] << "\n";
            // cout << "idx_apex_frroi[idx+1] --> " << idx_apex_frroi[idx+1] << "\n\n";
            int loval = idx_apex_frroi[idx];
            int hival = idx_apex_frroi[idx+1];
            if(loval > hival){
              loval = idx_apex_frroi[idx+1];
              hival = idx_apex_frroi[idx];
            }
            // cout << "idx 1--> " << idx << "\n";
            // cout << "idx_apex_frroi --> " << idx_apex_frroi << "\n";
            // cout << "loval 1--> " << loval << "\n";
            // cout << "hival 1--> " << hival << "\n";
            IntegerVector idx_vec2 = seq(loval, hival);
            NumericVector eics_apex_roi2 = eic_int[idx_vec2];
            rbVec[idx] = idx_apex_frroi[idx] + which_min(eics_apex_roi2);
          }
        }
        
        if((lbVec[idx] == idx_apex_frroi[idx]) | (rbVec[idx] == idx_apex_frroi[idx])){
          i_rm.push_back(idx);
          continue;
        }
      }
      // cout << "Arriving here --> lbVec --> " << lbVec << endl;
      // cout << "Arriving here --> rbVec --> " << rbVec << endl;
      // cout << "i_rm -------> " << i_rm << endl;
      
      if(i_rm.size() > 0){

        int ii = 0;
        for(IntegerVector::iterator it = lbVec.begin(); it < lbVec.end(); it++){
          if(is_true(any(ii == i_rm))){
            it = lbVec.erase(it);
            it--;
          }
          ii++;
        };
        
        ii = 0;
        for(IntegerVector::iterator it = rbVec.begin(); it < rbVec.end(); it++){
          if(is_true(any(ii == i_rm))){
            it = rbVec.erase(it);
            it--;
          }
          ii++;
        };
        
        ii = 0;
        for(IntegerVector::iterator it = idx_apex_frroi.begin(); it < idx_apex_frroi.end(); it++){
          if(is_true(any(ii == i_rm))){
            it = idx_apex_frroi.erase(it);
            it--;
          }
          ii++;
        };

      }
      
      if(idx_apex_frroi.size() > 0){
        is_peak = true;
      }
      //cout << "Arriving here --> 5 \n";
      info_pk = List::create(_["lb"] = lbVec, 
                             _["apex"] = idx_apex_frroi,
                             _["rb"] = rbVec, 
                             _["bl"] = baseline,
                             _["is_peak"] = is_peak);
      
    }
  }
  


  DetectedPeaks = List::create(_["info_pk"] = info_pk,
                               _["info_eic"] = eic_ints);

  return DetectedPeaks;
}


