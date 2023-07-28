#include <RcppArmadillo.h>

#include "lowess.h"
#include "DetectPeaks.h"

using namespace Rcpp;
using namespace std;


class funcSets {

public:
  void testFUN() {
    cout << "RUnning testFUN in Qiangclass" << "\n";
  }
  
  std::vector<double> lowessCpp(NumericVector x, NumericVector y, double spanVal){
    const std::vector<double> xs = as<std::vector<double>>(x);
    const std::vector<double> ys = as<std::vector<double>>(y);
    std::vector<double> res;
    lowess(xs, ys, spanVal, res);
    return res;
  }
  

};

// [[Rcpp::export]]
List SpectraDeconvCore(int idx, 
                       List spectra_eics, 
                       NumericMatrix ms1Peak, int ScanNum, NumericVector ms1PeakInfo, int idx_apex_ms1,
                       double min_peakwidth, double max_peakwidth, double snthr) {
  
  /*
   * This function is used to speed up --> DecoSpectra
   * Arguments Explainations 
   */
  
  // idx : is the index of the ms1 peak
  // spectra_eics : contains all ms2 spectra info within the MS1 RT range
  // ms1Peak: is the EIC info of the MS1 peak
  // ScanNum: is the number of the MS2 scan of all MS2 spectra within the RT range
  // ms1PeakInfo: is the ms1 peak summary vector
  // idx_apex_ms1: is the index of the apex mz point
  // min_peakwidth: is the minimum peak width for ms2 peak deconvolution
  // max_peakwidth: is the maximum peak width for ms2 peak deconvolution
  // snthr: is the signal to noise threshold
  
  
  
  Rcpp::List spec_decon(idx);
  Rcout << "Running into the FUNCTION ---> SpectraDeconvCore" << endl;
  
  
  return spec_decon;
}





std::vector<double> LowessFun(NumericVector x, NumericVector y, double spanVal) {
  funcSets testObj; // TODO: to delete when everything is done
  return testObj.lowessCpp(x,y,spanVal);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
res <- SpectraDeconvCore(3, ms2.eic.ext, peak.ms1.smooth, length(idx.ms2.ext),
                  as.numeric(info.pk.ms1), 11, 5, 20, 10)
*/
