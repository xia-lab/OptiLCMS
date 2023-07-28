// This script is a wrapper to provide a port to access Linear Regression Model from pense package
// To provide an interface for function "LsEnRegressionDispatch"
// Here, we use the interface from pense package directly
// CItation: The methods are proposed in Cohen Freue, G. V., Kepplinger, D., Salibián-Barrera, M., and Smucler, E. (2019) <https://projecteuclid.org/euclid.aoas/1574910036>. The package implements the extensions and algorithms described in Kepplinger, D. (2020) <doi:10.14288/1.0392915>.
#include "linear_regression.h"

// COV function added
double static cov_run(NumericMatrix x0, NumericVector y0){
  // step 1: format matrix
  arma::mat x(x0.nrow(), x0.ncol());
  for(int i=0; i<x0.ncol(); i++){
    NumericVector valc = x0(_,i);
    for(int j=0; j<valc.size();j++){
      x(j,i) = x0(j,i);
    }
  }
  // step 2: rank y and convert numericvector as arma::vec
  arma::vec y = as<arma::vec>(y0);
  // step 3: calculate cov
  arma::vec C = arma::cov(x,y);
  // step 4: summarize results and return
  double res = arma::max(C);
  return res;
}

NumericVector static generateLambdaVec(NumericMatrix X, NumericVector Y){
  // This function is originated from PENSE package to prepare a series of lambda value
  // we are generating 10 different lambd values for further usage
  // The meaning of all input parameters are same as the function PerformLinearRegress
  NumericVector lambdVecRes(10);
  
  double max_lambda;
  double max_cov = cov_run(X, Y);
  double nr = (double) X.nrow();
  max_lambda = fabs(((nr-1.0)/nr)*max_cov/0.0001);// 1% the minimum of alpha

  //cout << max_cov << " <- max_cov || max_lambda --> " << max_lambda << endl;
  double min_lambda = max_cov*1e-4; //1% of the min value of alpha
  if(log(min_lambda) < 0){min_lambda = 1.001;}
  if(log(max_lambda) < 0){max_lambda = 1e6;}
  
  double diff = (log(max_lambda) - log(min_lambda))/10.0;
  min_lambda = log(min_lambda);
  
  for(int i=0;i<10;i++){
    lambdVecRes[i] = min_lambda + diff*i;
  }
  return lambdVecRes;
}


// [[Rcpp::plugins("cpp14")]]
// /*/ // [[Rcpp::export]]
NumericVector PerformLinearRegress(NumericMatrix X, NumericVector Y, NumericVector penalty_loadings){
  // This function is used to perform linear regression analysis
  // The input includes two parameters:
  // 1. X, a NumericMatrix, including different components as different columns;
  // 2. Y, a NumericVector. This is a vector including response variables.
  
  // This function will return a NumericVector, which includes the coefficients of the linear regression model
  // In this function, we are going to include 11 alpha values (0~1);
  // Lambda values also have 10 members, which is evaluated by function generateLambdaVec
  
  // 1) prepare alpha and lambda value vectors
  NumericVector alphaVec = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  NumericVector lambdaVec = generateLambdaVec(X,Y);
  
  // 2) prepare parameter List
  List L(110), List_param, Lx, en_options;
  int num = 0;
  for(int i=0; i<11; i++){
    for(int j=0; j<10; j++){
      L[num] = List::create(_["lambda"] = lambdaVec[j], 
                            _["alpha"] = alphaVec[i]);
      num++;
    }
  }
  
  // 3) prepare include_intercept parameter
  LogicalVector v_intercept = {1}; // include include_intercept
  
  // 4) prepare other parameters
  en_options = List::create(_["algorithm"] = 5, 
                            _["eps"] = 1e-6,
                            _["sparse"] = false);
  if(penalty_loadings.size()==0){
    List_param = List::create(_["en_options"] = en_options);
  } else {
    List_param = List::create(_["en_options"] = en_options, 
                              _["pen_loadings"] = penalty_loadings);
  }
  
  // 5) Run linear regression with all parameters
  Lx = LsEnRegression(X, Y, L, v_intercept, List_param);
  
  // 6) Summarize all results to select the coefficient with minimal residue
  List estimatedResList = Lx["estimates"];
  
  List tmpL;
  NumericVector theseBetas, residueVec(110);
  double thisbeta, thisResidue;
  NumericVector currentCol, mCol;
  for(int n=0; n<110; n++){
    thisResidue = 0.0;
    std::vector<double> sumCol(Y.size());
    tmpL = estimatedResList[n];
    theseBetas = tmpL["beta"];
    // calculate residue
    for(int m=0; m<theseBetas.size();m++){
      thisbeta = theseBetas[m];
      currentCol = X(_,m);
      mCol = currentCol*thisbeta;
      for(int s=0; s<mCol.size(); s++){
        sumCol[s] = sumCol[s] + mCol[s];
      }
    }
    
    for(int s=0; s<sumCol.size();s++){
      thisResidue = thisResidue + fabs(sumCol[s] - Y[s]);
    }
    residueVec[n] = thisResidue;
    
  }
  // cout << "lambdaVec -> " << lambdaVec << endl;
  int idx_min = which_min(residueVec);
  // cout << "Now the min residue is --> " << min(residueVec) << " | idx-> " << idx_min << endl;
  tmpL = estimatedResList[idx_min];
  theseBetas = tmpL["beta"];
  
  // NumericVector tlmb = tmpL["lambda"];
  // cout <<  " lambdaVec[idx_min]-> " << tlmb << endl;
  return theseBetas;
}
