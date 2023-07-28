#include <RcppArmadillo.h>
#include "PerformDDAProcess.h"
#include "PerformDIAProcess.h"

using namespace Rcpp;

/* 
 * This script is designed to export DDA + DIA performing cpp function into R to accept all input parameters
 */

// [[Rcpp::export]]
List PerformDDADeco(NumericMatrix pm, 
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
                    std::string db_path,
                    bool decoOn,
                    bool useEntropy,
                    bool show_output,
                    int thread_num,
                    std::string file_nm){
  
  List res0 = PerformDDA_main(pm, 
                              scant1, scant2, 
                              scanms1, scanms2, 
                              prec_mzs, win_size, 
                              ppm1, ppm2, 
                              sn, filt, 
                              intensity_thresh,
                              ionmode, 
                              db_path,
                              decoOn,
                              useEntropy,
                              show_output,
                              thread_num,
                              file_nm);
  List res = res0;//[1];
  
  return res;
}

// [[Rcpp::export]]
List PerformDIADeco(List pm, 
                    NumericMatrix swath,
                    NumericVector scant1,
                    NumericVector scant2,
                    List scanms1,
                    List scanms2,
                    double pkw_min,
                    double ppm2,
                    double sn,
                    double span,
                    double filt){
  
  List res0 = PerformDIA_main(pm, 
                              swath,
                              scant1,
                              scant2,
                              scanms1,
                              scanms2,
                              pkw_min,
                              ppm2,
                              sn,
                              span,
                              filt);
  List res = res0;
  
  return res;
}



/*** R
# rm(list = ls()[ls() != "PerformDDADeco"])
# if(file.exists("~/Github/OptiLCMS2ID/data/dda_input_data.rda")){
#   load("~/Github/OptiLCMS2ID/data/dda_input_data_hilic.rda")
# } else if(file.exists("~/../Data/Github/OptiLCMS2ID/data/dda_input_data.rda")){
#   load("~/../Data/Github/OptiLCMS2ID/data/dda_input_data_hilic.rda")
# }
# # "dda_input_data_hilic"" is a hilic urine example, positive mode
# # system.time(res <- PerformDDA_main(matrix(), scanrts_ms1, scanrts_ms2, scan_ms1, scan_ms2,
# #                                    prec_mzs, 0.4, 10, 20, 12, 1000))
# # peak_mtx[,2] <- peak_mtx[,2]-1
# # peak_mtx[,3] <- peak_mtx[,3]+1
# system.time(res <- PerformDDADeco(peak_mtx, scanrts_ms1, scanrts_ms2, scan_ms1, scan_ms2,
#                                            prec_mzs, 1, 10, 20, 12, 1000))
*/