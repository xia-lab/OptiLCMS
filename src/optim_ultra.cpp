#include "optim_ultra.h"


class OptimU : public Functor {
public:
  arma::mat mpkMatrix;
  arma::vec eicVec;
  
  double operator()(const arma::vec &x) {
    double res = 0;
    int Sr = mpkMatrix.n_rows;
    int Sc = mpkMatrix.n_cols;
    arma::vec m = arma::ones<arma::vec>(Sr) * 0;
    for(int i = 0; i < Sr; i++){
      for(int j =0; j < Sc; j++){
        m(i) = m(i) + x(j) * mpkMatrix(i, j);
      }
    }
    
    for(int i = 0; i < Sr; i++){
      res = res+ pow(eicVec(i) - m(i), 2);
    }
    return res;
  }
};


double optim_ultra(NumericMatrix mpkmtx, NumericVector vec_eic, int main_idx) {
  
  arma::vec eic = as<arma::vec>(wrap(vec_eic));
  arma::mat mpk_mtx = as<arma::mat>(wrap(mpkmtx));
  
  // arma::vec eic = arma::ones<arma::vec>(vec_eic.size()) * 0;
  // arma::mat mpk_mtx = arma::zeros<arma::mat>(mpkmtx.nrow(), mpkmtx.ncol());
  // //arma::mat mpk_mtx;
  // for(int i = 0; i < vec_eic.size(); i++){
  //   eic[i] = vec_eic[i];
  //   for(int j = 0; j < mpk_mtx.n_cols; j++){
  //     mpk_mtx(i,j) = mpkmtx(i,j);
  //   }
  // }
  
  OptimU f;
  f.eicVec = eic;
  f.mpkMatrix = mpk_mtx;
  arma::vec lower = arma::ones<arma::vec>(mpk_mtx.n_cols) * 0;
  Roptim<OptimU> opt("L-BFGS-B");
  opt.set_lower(lower);
  opt.control.trace = 0;
  arma::vec x = arma::ones<arma::vec>(mpk_mtx.n_cols) * 0;
  opt.minimize(f, x);
  
  arma::vec res = opt.par();
  return res[main_idx];
  //Rcpp::Rcout << "------res --> " << res << std::endl;
}