#include "dda_utilities.h"
#include "entropy.h"
#include <regex>

double static lambdaEstimator(double x, NumericVector y){
  // This function is used to estimate the lamda value of 
  // y = λ*e^(-λx)
  // log(y) = log(λ) -λx*log(e)
  // log(y1) - log(y2) == λ*log(e)*(x2 - x1)
  // λ = (log(y1) - log(y2))/(log(e)*(x2 - x1))
  double y0 = y[0];
  double y2 = y[y.size()-1];
  double lambda = (log(y2) - log(y0))/(x);
  return lambda;
}

NumericMatrix row_erase (NumericMatrix& x, IntegerVector& rowID) {
  rowID = rowID.sort();
  
  NumericMatrix x2(Dimension(x.nrow()- rowID.size(), x.ncol()));
  int iter = 0; 
  int del = 1; // to count deleted elements
  for (int i = 0; i < x.nrow(); i++) {
    if (i != rowID[del - 1]) {
      x2.row(iter) = x.row(i);
      iter++;
    } else {
      del++;
    }
  }
  return x2;
}

// static int calculateTrue(LogicalVector logcVec){
//   int num = 0;
//   for(bool bl : logcVec){
//     if(bl) num++;
//   }
//   return num;
// }

List MSCentroidsGrouping(NumericVector allMZs){
  // This function is used to group all mz centroids into groups
  List res;
  if(allMZs.size() < 100000){
    warning("Too few MS centroids included (as least 10k is required)!");
    return res;
  }
  //cout << "RUnning into function MSCentroidsGrouping --> " << allMZs.size() << endl;
  
  NumericVector mzdiffVec(allMZs.size()-1);
  int k = 0;
  for(int i = 1; i < allMZs.size(); i++){
    mzdiffVec[k] = allMZs[i] - allMZs[i-1];
    k++;
  }
  
  NumericVector mzdiffVec_sorted = clone(mzdiffVec);
  mzdiffVec_sorted.sort();
  double bin_size;
  int bin_num;
  if(mzdiffVec_sorted.size() > 3000000){
    bin_size = (max(mzdiffVec) - min(mzdiffVec))/500;
    bin_num = 500;
  } else {
    bin_size = (max(mzdiffVec) - min(mzdiffVec))/200;
    bin_num = 200;
  }

  IntegerVector statsVec(bin_num);
  k = 0;
  for(double indiNum : mzdiffVec_sorted){
    //cout << "indiNum --> " << indiNum << endl;
    if(indiNum >= (k+1)*bin_size){
      //cout << k << " <<-- (k+1)*bin_size -> " << (k+1)*bin_size << endl;
      while(indiNum > (k+1)*bin_size){
        k++;
      }
    }
    statsVec[k] = statsVec[k] + 1;
  }

  // cout << "statsVec --> " << statsVec << endl;
  // res.push_back(allMZs);
  // return res;
  
  int size4lambd = bin_num;//*0.8;
  NumericVector lambdVec(size4lambd);
  for(k = 0; k < size4lambd-1; k++){
    int idx0, idx1, idx2;
    idx0 = bin_num - k;
    idx1 = bin_num - (k+1);
    idx2 = bin_num - (k+2);
    double y0, y1, y2;
    y0 = statsVec[idx0];
    y1 = statsVec[idx1];
    y2 = statsVec[idx2];
    NumericVector y_vals = NumericVector::create(y2,y1,y0);
    lambdVec[k] = lambdaEstimator(-2, y_vals);
  }
// 
  // cout << "lambdVec      --> " << lambdVec << endl;
  // cout << "lambdVec size --> " << lambdVec.size() << endl;
  k = 0;
  NumericVector lambdVec_real;
  IntegerVector binIdx_vec;
  for (double lmv : lambdVec){
    k++;
    if((lmv != R_PosInf) & (lmv != R_NegInf) & (lmv != R_NaN) & (lmv > 0)){
      binIdx_vec.push_back(k);
      lambdVec_real.push_back(lmv);
    }
  }
  // cout << "lambdVec_real is --> " << lambdVec_real << endl;
  // cout << "bindIdx_vec is   --> " << binIdx_vec << endl;

  // remove 10% head and 10% tail because of unstable lambd
  int rmNum = floor(lambdVec_real.size()*0.1);
  NumericVector clean_lambd_vec;
  for(int i = rmNum; i < lambdVec_real.size() - rmNum; i++){
    clean_lambd_vec.push_back(lambdVec_real[i]);
  }
  
  double cutoff = min(clean_lambd_vec);
  cutoff = cutoff/2.0;
  // cout << "cutoff is   --> " << cutoff << endl;
  // cout << "bin_size is   --> " << bin_size << endl;
  int cutoff_idx = max(binIdx_vec) - 1;
  for(int i = lambdVec_real.size()-1; i >= 0; i--){
    if(lambdVec_real[i] < cutoff){
      cutoff_idx = binIdx_vec[i];
      break;
    }
  }
  
  cutoff_idx = bin_num - cutoff_idx + 1;
  
  double mz_cutoff = bin_size*cutoff_idx;
  //cout << cutoff_idx << " <-- mz_cutoff -> " << mz_cutoff << endl;
  
  double mzdiff;
  vector<double> thisMZcluster;
  NumericVector allMZsU = Rcpp::unique(allMZs);
  allMZsU = Rcpp::signif(allMZsU, 8);
  allMZsU = Rcpp::unique(allMZsU);
  allMZsU.sort();
  // cout << "allMZsU size: " << allMZsU.size() << endl;
  for(int i = 1; i < allMZsU.size(); i++){
    //cout << "i is --> " << i << endl;
    mzdiff = allMZsU[i] - allMZsU[i-1];
    if(mzdiff < mz_cutoff){
      thisMZcluster.push_back(allMZsU[i-1]);
    } else {
      thisMZcluster.push_back(allMZsU[i-1]);
      if(thisMZcluster.size() > 5){
        res.push_back(thisMZcluster);
      }
      thisMZcluster.clear();
    }
  }

  //res.push_back(lambdVec_real);
  // res.push_back(allMZs);
  // res.push_back(allMZsU);
  return res;
}

float dot(NumericVector a, NumericVector b, bool norm){
  
  arma::vec ay = Rcpp::as<arma::vec>(a);
  arma::vec by = Rcpp::as<arma::vec>(b);
  if(norm){
    return arma::norm_dot(ay, by);
  }
  return arma::dot(ay, by);
}


double entropy(NumericMatrix peaks_a, NumericMatrix peaks_b){
  double res = 0;
  res = r_calculate_entropy_similarity(peaks_a, peaks_b,
                                       -1.0, 2, // ppm / da
                                       true, // clean or not
                                       10, 1500, // mz range
                                       0.01, // noise thre
                                       100); // max peak
  return res;
}

NumericMatrix cosineSimilarity(NumericMatrix Xr) {
  //Rcpp::NumericMatrix Xr(Xs);  // creates Rcpp matrix from SEXP
  int n = Xr.nrow(), k = Xr.ncol();
  arma::mat X(Xr.begin(), n, k, false); // reuses memory and avoids extra copy
  arma::mat Y = arma::trans(X) * X; // matrix product
  arma::mat res = (1 - Y / (arma::sqrt(arma::diagvec(Y)) * arma::trans(arma::sqrt(arma::diagvec(Y)))));
  return Rcpp::wrap(res);
}

NumericMatrix ms2peak_parse(string text){
  
  //string text = "70.098\t13000\n80.099\t19084\n450.910\t78921\n";
  string spliter = "\n";
  string spliter2 = "\t";
  vector<string> items{};
  NumericVector mz_vec;
  NumericVector int_vec;
  
  size_t pos = 0, pos2 = 0;
  string mzstr, intstr;
  double mz_num, int_num;
  while ((pos = text.find(spliter)) != string::npos) {
    items.push_back(text.substr(0, pos));
    text.erase(0, pos + 1);
  }
  for (string &str : items) {
    pos2 = str.find(spliter2);
    mzstr = str.substr(0, pos2);
    intstr = str.substr(pos2+1, str.size());
    mz_num = std::stod(mzstr);
    int_num = std::stod(intstr);
    //cout << "mz_num - " << mz_num << " int_num - " << int_num << endl;
    mz_vec.push_back(mz_num);
    int_vec.push_back(int_num);
  }
  
  NumericMatrix res = cbind(mz_vec, int_vec);
  return res;
}

double spectrumSimilarity(NumericMatrix mtx1, NumericMatrix mtx2, double ppm_ms2){
  // raw matrix would be ok, will be re-sorted as two vectors for 'dot' similarity
  NumericVector vec_r, vec_s;
  IntegerVector idx_r, idx_s;

  // NumericVector nmv1 = mtx1(_,0);
  // NumericVector nmv2 = mtx1(_,1);
  // cout << "mtx1(_,0)" <<  nmv1 << endl;
  // cout << "mtx1(_,1)" <<  nmv2 << endl;
  // 
  // NumericVector nmv21 = mtx2(_,0);
  // NumericVector nmv22 = mtx2(_,1);
  // cout << "mtx2(_,0)" <<  nmv21 << endl;
  // cout << "mtx2(_,1)" <<  nmv22 << endl;
  
  for(int r=0; r<mtx1.nrow(); r++) {
    for(int s=0; s<mtx2.nrow(); s++) {
      if(abs(mtx1(r,0)-mtx2(s,0))/(mtx1(r,0)) < (ppm_ms2*1e-6)) {
        idx_r.push_back(r);
        idx_s.push_back(s);
        vec_r.push_back(mtx1(r,1));
        vec_s.push_back(mtx2(s,1));
        break;
      }
    }
  }
  for(int r=0; r<mtx1.nrow(); r++){
    if(is_false(any(r == idx_r))){
      vec_r.push_back(mtx1(r,1));
      vec_s.push_back(0.0);
    }
  }
  for(int s=0; s<mtx2.nrow(); s++){
    if(is_false(any(s == idx_s))){
      vec_r.push_back(0.0);
      vec_s.push_back(mtx2(s,1));
    }
  }
  
  float resf = dot(vec_r, vec_s, true);
  double resd = (double) resf;
  return resd;
}

double entropySimilarity(NumericMatrix mtx1, NumericMatrix mtx2, double ppm_ms2){
  // raw matrix would be ok, calculate entropy similarity
  cout << "We are using entropySimilarity now " << endl;
  double res = 0;
  res = r_calculate_entropy_similarity(mtx1, mtx2,
                                       -1.0, ppm_ms2, // da, ppm
                                       true, // clean or not
                                       10, 1500, // mz range
                                       0.01, // noise thre
                                       100); // max peak
  return res;
}


double get_mass_sodium(){
  return 22.98976928;
}

double get_mass_potassium(){
  return 38.963707;
}

vector<int> parse_formula(std::string formula_txt){
  // C H O N P S
  // cout << "formula_txt -> " << formula_txt << endl;
  if((formula_txt.find("[") != -1) | (formula_txt.find("+") != -1) | (formula_txt.find("-") != -1) | (formula_txt.find("(") != -1)){
    //std::replace(formula_txt.begin(), formula_txt.end(), "[", "");
    formula_txt = std::regex_replace(formula_txt, std::regex("\\[|\\]|\\+|\\-|\\(|\\)"), "");
    //formula_txt.replace(formula_txt.find("["), 1, "");
  }
  //cout << "formula_txt -> " << formula_txt << endl;
  vector<int> formula = {0,0,0,0,0,0};
  
  int pos_c = formula_txt.find("C");
  if(pos_c != -1){
    string tmp_c = formula_txt.substr(pos_c+1, 1);
    bool bl_c = std::regex_match(tmp_c, std::regex("^[A-Za-z]+$"));
    if(!bl_c & (tmp_c != "")){
      tmp_c = formula_txt.substr(pos_c+1, 3);
      int num_c = std::stoi(tmp_c);
      formula[0] = num_c;
    } else {
      formula[0] = 1;
    }
  }
  
  int pos_h = formula_txt.find("H");
  if(pos_h != -1){
    string tmp_h = formula_txt.substr(pos_h+1, 1);
    bool bl_h = std::regex_match(tmp_h, std::regex("^[A-Za-z]+$"));
    if(!bl_h & (tmp_h != "")){
      tmp_h = formula_txt.substr(pos_h+1, 3);
      int num_h = std::stoi(tmp_h);
      formula[1] = num_h;
    } else {
      formula[1] = 1;
    }
  }
  
  int pos_o = formula_txt.find("O");
  if(pos_o != -1){
    string tmp_o = formula_txt.substr(pos_o+1, 1);
    bool bl_o = std::regex_match(tmp_o, std::regex("^[A-Za-z]+$"));
    if(!bl_o & (tmp_o != "")){
      tmp_o = formula_txt.substr(pos_o+1, 3);
      int num_o = std::stoi(tmp_o);
      formula[2] = num_o;
    } else {
      formula[2] = 1;
    }
  }
  int pos_n = formula_txt.find("N");
  if(pos_n != -1){
    string tmp_n = formula_txt.substr(pos_n+1, 1);
    bool bl_n = std::regex_match(tmp_n, std::regex("^[A-Za-z]+$"));
    if(!bl_n & (tmp_n != "")){
      tmp_n = formula_txt.substr(pos_n+1, 3);
      int num_n = std::stoi(tmp_n);
      formula[3] = num_n;
    } else {
      formula[3] = 1;
    }
  }
  int pos_p = formula_txt.find("P");
  if(pos_p != -1){
    string tmp_p = formula_txt.substr(pos_p+1, 1);
    bool bl_p = std::regex_match(tmp_p, std::regex("^[A-Za-z]+$"));
    if(!bl_p & (tmp_p != "")){
      tmp_p = formula_txt.substr(pos_p+1, 3);
      int num_p = std::stoi(tmp_p);
      formula[4] = num_p;
    } else {
      formula[4] = 1;
    }
  }
  int pos_s = formula_txt.find("S");
  if(pos_s != -1){
    string tmp_s = formula_txt.substr(pos_s+1, 1);
    bool bl_s = std::regex_match(tmp_s, std::regex("^[A-Za-z]+$"));
    if(!bl_s & (tmp_s != "")){
      tmp_s = formula_txt.substr(pos_s+1, 3);
      int num_s = std::stoi(tmp_s);
      formula[5] = num_s;
    } else {
      formula[5] = 1;
    }
  }
  
  // cout << "H -> " << pos_h << endl;
  // cout << "O -> " << pos_o << endl;
  // cout << "N -> " << pos_n << endl;
  // cout << "P -> " << pos_p << endl;
  // 
  
  return formula;
}


double neutral_loss_similarity(NumericVector exp_mzs, NumericVector exp_ints, 
                               NumericVector ref_mzs, NumericVector ref_ints,
                               double ppm2) {
  
  NumericVector all_mzs = clone(exp_mzs);
  for(int m=0; m<ref_mzs.size(); m++){
    if(is_true(all((abs(ref_mzs[m] - exp_mzs)/ref_mzs[m] > ppm2*1e-6)))){
      all_mzs.push_back(ref_mzs[m]);
    }
  }
  
  all_mzs.sort();
  NumericVector exp_ints_new(all_mzs.size()), ref_ints_new(all_mzs.size());
  for(int n=0; n<all_mzs.size(); n++){
    for(int exp_idx=0; exp_idx < exp_mzs.size(); exp_idx++){
      if(abs(exp_mzs[exp_idx]- all_mzs[n])/all_mzs[n] < ppm2*1e-6){
        exp_ints_new[n] = exp_ints[exp_idx];
      }
    }
    for(int ref_idx=0; ref_idx < ref_mzs.size(); ref_idx++){
      if(abs(ref_mzs[ref_idx] - all_mzs[n])/all_mzs[n] < ppm2*1e-6){
        ref_ints_new[n] = ref_ints[ref_idx];
      }
    }
  }
  exp_ints_new = exp_ints_new/max(exp_ints_new);
  ref_ints_new = ref_ints_new/max(ref_ints_new);
  double res = dot(exp_ints_new, ref_ints_new, true);

  // cout << "exp_mzs -> " << exp_mzs << endl;
  // cout << "exp_ints-> " << exp_ints << endl;
  // cout << "ref_mzs -> " << ref_mzs << endl;
  // cout << "ref_ints-> " << ref_ints << endl;
  // cout << "----- " << endl;
  // cout << "all_mzs      -> " << all_mzs << endl;
  // cout << "exp_ints_new -> " << exp_ints_new << endl;
  // cout << "ref_ints_new -> " << ref_ints_new << endl;
  // cout << "======================== " << endl;
  
  return res;
}

NumericVector unique_num(NumericVector x){
  NumericVector res(1);
  res[0] = x[0];
  bool notfound;
  for(int i=1; i<x.size(); i++){
    notfound = true;
    for(int j=0; j<res.size(); j++){
      if(res[j] == x[i]){
        notfound = false;
      }
    }
    if(notfound){
      res.push_back(x[i]);
    }
  }
  return res;
}

