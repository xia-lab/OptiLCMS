#include "utilities.h"

  
vector<double> lowessCpp(IntegerVector x, NumericVector y, double spanVal){
  const vector<double> xs = as<vector<double>>(x);
  const vector<double> ys = as<vector<double>>(y);
  vector<double> res;
  lowess(xs, ys, spanVal, res);
  return res;
}

NumericVector SmoothLoess(NumericMatrix &eic, double span) {
  //cout << "Running into function --- > SmoothLoess < ----" << endl;
  //cout << "eic_size: --> " << eic.size() << endl;
  NumericVector ress;
  int n_spec = eic.nrow();
  
  NumericVector ints = eic( _ , 3); // get all intensity values
  
  if(n_spec > 5) {
    vector<double> res;
    double tmpVal = 4 / (double)n_spec;
    double spanVal = (tmpVal > span) ? tmpVal : span;
    IntegerVector idxs = seq(1, n_spec);
    res = lowessCpp(idxs, ints, spanVal);
    ress = res;
  } else {
    ress = ints;
  }
  //cout << "ress: -> " << ress << endl;
  return ress;
}

bool checkContinuousPtsAboveThr(NumericVector v, int iStart, double num, double thr, int nSkipMax) {
  int cnt = 0;
  bool res = false;
  int nSkip = 0;
  int nSkipPre = 0;
  for (int i = iStart; i < v.length(); i++) {
    if (v[i] > thr) {
      cnt++;
      nSkip = 0;
      nSkipPre = 0;
    } else {
      if (cnt > 0) {
        nSkipPre = nSkip;
        nSkip++;
      } else {
        nSkip = nSkipMax + 1;
      }
      if (nSkip > nSkipMax) {
        cnt = 0;
        nSkip = 0;
        nSkipPre = 0;
      } else {
        if (nSkipPre < nSkipMax) {
          cnt++;
        } else {
          cnt = cnt - nSkip + 1;
          nSkipPre = 0;
        }
      }
    }
    
    if (cnt >= num) {
      return(true);
    }
  }
  
  return(res);
}

NumericVector getContinuousPtsAboveThrIdx(NumericVector v, int iStart, int num, double thr, int nSkipMax) {
  int cnt = 0;
  int stidx = 0;
  int enidx = 0;
  int nv = v.length();
  int nSkip = 0;
  int nSkipPre = 0;
  NumericVector res(nv);
  for (int i = iStart; i < nv; i++) {
    res[i] = false;
    if (v[i] > thr) {
      cnt++;
      nSkip = 0;
      nSkipPre = 0;
      if (cnt == 1) {
        stidx = i;
      } else {
        enidx = i;
      }
    } else {
      if (cnt > 0) {
        nSkipPre = nSkip;
        nSkip++;
      } else {
        nSkip = nSkipMax + 1;
      }
      if (nSkip > nSkipMax) {
        cnt = 0;
        nSkip = 0;
        nSkipPre = 0;
      } else {
        if (nSkipPre < nSkipMax) {
          cnt++;
        } else {
          cnt = cnt - nSkip + 1;
          nSkipPre = 0;
        }
      }
    }
    
    if ((cnt == 0 || i == nv - 1) && (enidx - stidx + 1) >= num) {
      for (int j = stidx; j <= enidx; j++) {
        res[j] = true;
      }
      stidx = 0;
      enidx = 0;
    }
  }
  return(res);
}

IntegerVector whichTrue(LogicalVector vecValues){
  IntegerVector trueVec;
  for(int i = 0; i < vecValues.size(); i++) {
    if(vecValues[i]){
      trueVec.push_back(i);
    }
  }
  return trueVec;
}

int whichTrue1(LogicalVector vecValues){
  int trueVec = -1;
  for(int i = 0; i < vecValues.size(); i++) {
    if(vecValues[i]){
      trueVec = i;
      return trueVec;
    }
  }
  return trueVec;
}

IntegerVector GetRoi(NumericVector is_roi, int idx_apex_eic) {
   // cout << "Running into the function GetROI <<---- \n";
   // cout << "is_roi --> " << is_roi << endl;
   // cout << "idx_apex_eic --> " << idx_apex_eic << endl;
  
  IntegerVector idx_fr_roi;
  
  LogicalVector res0 = (is_roi == 0);
  //cout << "res0 --> " << res0 << endl;
  bool res = is_true(all(res0 == TRUE));
  if(res){
    return idx_fr_roi;
  }
  
  IntegerVector i_s = whichTrue(diff(is_roi) == 1);
  IntegerVector i_e = whichTrue(diff(is_roi) == -1);
  //for (double i : diff(is_roi)) cout << "diff(is_roi) i --> " << i << "" << endl;
  
   // cout << "i_s --> " << i_s << endl;
   // cout << "i_e --> " << i_e << endl;
  
  if(is_roi[0] == 1) i_s.push_front(-1);
  if(is_roi[is_roi.size() - 1] == 1) i_e.push_back(is_roi.size() - 1);
  
  IntegerVector i_s_new;
  IntegerVector i_e_new;
  int maxSize = 1;
  bool bool1 = false;
  bool bool2 = false;
  
  if(i_s.size() == 0) {
    i_s_new.push_back(0);
  } else {
    i_s_new = i_s + 1;
    maxSize = i_s_new.size();
    bool1 = true;
  }
  
  if(i_e.size() == 0) {
    i_e_new = is_roi.size() - 1;
  } else {
    i_e_new = i_e;
    maxSize = i_e_new.size();
    bool2 = true;
  }
  
  // cout << "i_s_new --> " << i_s_new << endl;
  // cout << "i_e_new --> " << i_e_new << endl;
  // cout << "maxSize --> " << maxSize << endl;
  // cout << "bool1 --> " << bool1 << endl;
  // cout << "bool2 --> " << bool2 << endl;
  
  if(bool1 & !bool2){
    for(int j = 0; j < maxSize; j++){
      if((i_s_new[j] <= idx_apex_eic) & (i_e_new[0] >= idx_apex_eic)){
        idx_fr_roi = seq(i_s_new[j], i_e_new[0]);
        break;
      }
    }
  } 
  if(bool2 & !bool1){
    for(int j = 0; j < maxSize; j++){
      if((i_s_new[0] <= idx_apex_eic) & (i_e_new[j] >= idx_apex_eic)){
        idx_fr_roi = seq(i_s_new[0], i_e_new[j]);
        break;
      }
    }
  }
  if(bool2 & bool1){
    for(int j = 0; j < maxSize; j++){
      if((i_s_new[j] <= idx_apex_eic) & (i_e_new[j] >= idx_apex_eic)){
        idx_fr_roi = seq(i_s_new[j], i_e_new[j]);
        break;
      }
    }
  }
  
  return idx_fr_roi;
}

double EstimateChromNoise(NumericVector &x, double trim = 0.05, int min_pts = 20) {
  int countVal = 0;
  double res;
  for(int j = 0; j < x.size(); j++){
    if(x[j] > 0) countVal++;
  }
  if(countVal < min_pts) {
    res = mean(x);
    return res;
  }
  NumericVector newX = x[x > 0.0];
  int rm_len = floor(newX.size()*trim);
  newX = newX.sort();
  newX.erase(newX.size() - (rm_len), newX.size());
  newX.erase(0, rm_len);
  res = mean(newX);
  return res;
}

NumericVector GetLocalNoiseEstimate(NumericVector d, IntegerVector idx_fr_roi, double noiserange_min, double noiserange_max, 
                                    int Nscantime, double threshold, int num){
  NumericVector res;
  IntegerVector drange = idx_fr_roi;
  
  if(drange.size() < Nscantime) {
    //cout << "RUnning into 1st option --> GetLocalNoiseEstimate <--- \n";
    res = CalculateBL(d, drange, threshold, num, 0, noiserange_min);
    //cout << "CalculateBL RES ---> " << res << "\n";
  } else {
    //cout << "RUnning into 2nd option --> GetLocalNoiseEstimate <--- \n";
    int mid = ceil((double) d.size()/2);
    //cout << "mid --> " << mid << endl;
    NumericVector tmpVec = {1, mid - noiserange_min};
    int dr1 = max(tmpVec);
    tmpVec = {mid + noiserange_min, (double)d.size()};
    int dr2 = min(tmpVec);
    drange = seq(dr1 - 1, dr2 - 1);
    res = CalculateBL(d, drange, threshold, num, 0, noiserange_min);
    // cout << "drange ---> " << drange << "\n";
    // cout << "res ---> " << res << "\n";
    tmpVec = {1, mid - noiserange_max};
    dr1 = max(tmpVec);
    tmpVec = {mid + noiserange_max, (double)d.size()};
    dr2 = min(tmpVec);
    drange = seq(dr1 -1, dr2-1);
    // cout << "drange.size() --> " << drange.size() << endl;
    // cout << "drange --> " << drange << endl;
    // cout << "Nscantime --> " << Nscantime << endl;
    if(drange.size() < Nscantime) {
      //cout << "RUnning into 3rd option --> GetLocalNoiseEstimate <--- \n";
      NumericVector res2 = CalculateBL(d, drange, threshold, num, 0, noiserange_min);
      NumericVector resFinal = {min(res[0], res2[0]), max(res[1], res2[1])};
      //cout << "resFinal --> " << resFinal << "\n";
      return resFinal;
    }
  }
  return res;
}

NumericVector CalculateBL(NumericVector d, IntegerVector drange, double threshold, int num, int n_skip_max, double noiserange_min){
  //cout << "\n --------------- \n Running into the function --> CalculateBL <--- \n" << endl;
  //cout << "d --> " << d << endl;
  //cout << "drange --> " << drange << endl;
  //cout << "min drange --> " << min(drange) << endl;
  //cout << "max drange --> " << max(drange) << endl;
  NumericVector d0 = clone(d);
  
  double min_drange = min(drange);
  double max_drange = max(drange);
  
  NumericVector BLres;
  
  d.erase(min_drange, max_drange + 1);
  //cout << "nd -----> " << d << endl;
  NumericVector n1_cp = getContinuousPtsAboveThrIdx(d, 0, num, threshold, n_skip_max);

  NumericVector::iterator it = d.begin();
  for(int i = 0; i < n1_cp.size(); i++){
    if(n1_cp[i] != 0){
      it = d.erase(it);
    } else {
      it++;
    }
  }
  //cout << "d.size() --> " << d.size() << endl;
  double baseline1 = 1.0;
  double sdnoise1 = 1.0;
  if(d.size() > 1){
    baseline1 = mean(d);
    sdnoise1 = sd(d);
  }
  
  int d1 = drange[0];
  int d2 = drange[drange.size() - 1];
  
  //cout << "d1 --> " << d1 << endl;
  //cout << "d2 --> " << d2 << endl;
  
  NumericVector nrange2;
  NumericVector nrange2_2;
  
  NumericVector tNumVec = {0, d1 - noiserange_min};
  nrange2 = seq(max(tNumVec), d1);
  tNumVec = {(double)d0.size() - 1, d2 + noiserange_min};
  nrange2_2 = seq(d2, min(tNumVec));

  double tmpVal;
  for(int k = 0; k < nrange2_2.size(); k++){
    tmpVal = nrange2_2[k];
    nrange2.push_back(tmpVal);
  }
  
  //cout << "nrange2 --> " << nrange2 << endl;
  NumericVector n2 = d0[nrange2];
  //cout << "n2 --> " << n2 << endl;
  NumericVector n2_cp = getContinuousPtsAboveThrIdx(n2, 0, num, threshold, n_skip_max);
  //cout << "n2_cp --> " << n2_cp << endl;
  
  NumericVector::iterator it2 = n2.begin();
  for(int i = 0; i < n2_cp.size(); i++){
    if(n2_cp[i] != 0){
      it2 = n2.erase(it2);
    } else {
      it2++;
    }
  }
  //cout << "n2 new --> " << n2 << endl;
  double baseline2 = 1.0;
  double sdnoise2 = 1.0;
  if(n2.size() > 1){
    baseline2 = mean(n2);
    sdnoise2 = sd(n2);
  }
  
  BLres = {min(baseline1, baseline2), min(sdnoise1, sdnoise2)};
  return BLres;
}

IntegerVector FindLocalMax(NumericVector x, int m = 3, double v = 20) {
  IntegerVector pksf;
  IntegerVector shape = diff(sign(diff(x)));
  //cout << "shape --> " << shape << endl;
  IntegerVector idx_true = whichTrue(shape < 0);
  //cout << "idx_true --> " << idx_true << endl;
  int i, z, w;
  IntegerVector tmpIntVec1, tmpIntVec2;
  
  for(int j = 0; j < idx_true.size(); j++){
    i = idx_true[j];
    z = i - m + 1;
    if(z <= 0) z = 0;
    w = i + m + 1;
    if(w >= (x.size() - 1)) w = x.size() - 1;
    //IntegerVector tmpIntVec1 = Rcpp::Range(z, i);
    if(z < i){
      tmpIntVec1 = seq(z, i);
    } else {
      tmpIntVec1 = seq(i, z);
    }
    
    if(i+2 < w){
      tmpIntVec2 = seq(i+2, w);
    } else {
      tmpIntVec2 = seq(w, i+2);
    }

    NumericVector t1 = x[tmpIntVec1];
    NumericVector t2 = x[tmpIntVec2];
    LogicalVector tmpLogVec1 = (t1 <= x[i + 1]);
    LogicalVector tmpLogVec2 = (t2 <= x[i + 1]);
    //cout << "x[i+1] " << x[i+1] << endl;
    //cout << "tmpLogVec1 --> " << tmpLogVec1 << endl;
    //cout << "tmpLogVec2 --> " << tmpLogVec2 << endl;
    if(is_true(all(tmpLogVec1)) & is_true(all(tmpLogVec2))){
      pksf.push_back(i + 1);
    }
  }
  //cout << "pksf --> " << pksf << endl;
  //cout << "x --> " << x << endl;
  NumericVector FnumV = x[pksf];
  LogicalVector FlogV =  FnumV >= v;
  //cout << "FnumV --> " << FnumV << endl;
  //cout << "FlogV --> " << FlogV << endl;
  return pksf[FlogV];
}

IntegerVector FindLocalMin(NumericVector x, int m = 3) {
  IntegerVector pksf;
  IntegerVector shape = diff(sign(diff(x)));
  IntegerVector idx_true = whichTrue(shape > 0);
  int i, z, w;
  IntegerVector tmpIntVec1, tmpIntVec2;
  
  for(int j = 0; j < idx_true.size(); j++){
    i = idx_true[j];
    z = i - m + 1;
    if(z <= 0) z = 0;
    w = i + m + 1;
    if(w >= (x.size() - 1)) w = x.size() - 1;
    if(z < i){
      tmpIntVec1 = seq(z, i);
    } else {
      tmpIntVec1 = seq(i, z);
    }
    
    if(i+2 < w){
      tmpIntVec2 = seq(i+2, w);
    } else {
      tmpIntVec2 = seq(w, i+2);
    }
    
    NumericVector t1 = x[tmpIntVec1];
    NumericVector t2 = x[tmpIntVec2];
    LogicalVector tmpLogVec1 = (t1 >= x[i + 1]);
    LogicalVector tmpLogVec2 = (t2 >= x[i + 1]);

    if(is_true(all(tmpLogVec1)) & is_true(all(tmpLogVec2))){
      pksf.push_back(i + 1);
    }
  }

  return pksf;
}

NumericMatrix mergeEIC(NumericMatrix x, NumericMatrix y) {
  int xr = x.nrow()-1;
  int yr = y.nrow()-1;
  int is, ie;
  if (x(0, 0) < y(0, 0)) {
    is = x(0, 0);
  } else {
    is = y(0, 0);
  }
  
  if (x(xr, 0) > y(yr, 0)) {
    ie = x(xr, 0);
  } else {
    ie = y(yr, 0);
  }
  NumericMatrix dm(ie-is+1, 3);
  for (int i=0; i<=xr; i++) {
    int idx = x(0, 0) - is + i;
    dm(idx ,0) = x(i, 0);
    dm(idx, 1) = x(i, 1);
  }
  for (int i=0; i<=yr; i++) {
    int idx = y(0, 0) - is + i;
    
    dm(idx ,0) = y(i, 0);
    dm(idx, 2) = y(i, 1);
  }
  return dm;
}

double Gauss(int x, int h, int mu, int sigma){
  double a = pow(x-mu, 2);
  double b = 2*pow(sigma, 2);
  return h*exp(-(a/b));
}

double cor_fast (NumericVector x, NumericVector y){
  
  double sumx = 0; 
  double sumx2 = 0;
  double sumy = 0; 
  double sumy2 = 0;
  
  for (int r = 0; r < x.size(); r++) {
    double d = x[r];
    sumx += d;
    sumx2 += pow(d,2);
  }
  for (int r = 0; r < y.size(); r++) {
    double d = y[r];
    sumy += d;
    sumy2 += pow(d,2);
  }
  
  double res_stdev_x = sqrt((x.size()) * sumx2 - pow(sumx, 2));
  double res_stdev_y = sqrt((y.size()) * sumy2 - pow(sumy, 2));
  
  double sXY = 0;
  for(int r = 0; r < x.size(); r++){
    sXY += x[r] * y[r];
  }
  
  return (x.size() * sXY - sumx * sumy) / (res_stdev_x * res_stdev_y);
  
  //return 0.0;
}

double GetDistantP(NumericMatrix peak1, 
                   NumericMatrix peak2){
  
  NumericMatrix dt = mergeEIC(peak1, peak2);
  // NumericVector dt_1 = dt( _ , 1);
  // cout << " dt_1" << which_max(dt( _ , 1)) << endl;
  // cout << " dt_2" << which_max(dt( _ , 2)) << endl;
  double c_g = Gauss(which_max(dt( _ , 1)) - which_max(dt( _ , 2)),
                     1, 0, 12);
  if(c_g < 0.7) return 1.0;
  double c_p = cor_fast(dt( _ , 1), dt( _ , 2)) * c_g;
  
  if(isnan(c_p)){
    c_p = 0;
  }
  // cout << "c_p is --> " << c_p << " | " << c_g << endl;
  return 1-c_p;
}

// double funOptimc(NumericVector x, NumericVector M, NumericMatrix S) {
//   int Sr = S.nrow();
//   int Sc = S.ncol();
//   double  res = 0;
//   NumericVector m(Sr);
//   for (int i=0; i<Sr; i++) {
//     m(i) = 0;
//   }
//   for (int i=0; i<Sr; i++) {
//     for (int j=0; j<Sc; j++) {
//       m(i) = m(i) + x(j) * S(i, j);
//     }
//   }
//   for (int i=0; i<Sr; i++) {
//     res = res + pow(M(i) - m(i), 2);
//   }
//   return res;
// }
// 
// NumericVector optim_real(NumericMatrix mpk_mtx, NumericVector eic){
//   NumericVector res;
//   Rcpp::Environment stats("package:stats");
//   Rcpp::Function optimR = stats["optim"];
//   
//   NumericVector init_val (mpk_mtx.ncol());
//   Rcpp::List opt_results = optimR(Rcpp::_["par"]  = init_val,
//                                   Rcpp::_["fn"]   = Rcpp::InternalFunction(&funOptimc),
//                                   Rcpp::_["lower"] = 0,
//                                   Rcpp::_["method"] = "L-BFGS-B",
//                                   Rcpp::_["M"] = eic,
//                                   Rcpp::_["S"] = mpk_mtx);
//   res = opt_results[0];
//   return res;
// }
// 
