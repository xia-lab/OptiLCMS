#include "DecoSpectra.h"
#include "PerformDIAProcess.h"

// [[Rcpp::export]]
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
                     double filt){
  cout << "PerformDIAProcess Preparation starting..." << "\n";
  diaprocess runObj;
  runObj.setDIA_Arguments(pm, swath, scanrt1, scanrt2, 
                       scanms1, scanms2, 
                       pkw_min, ppm2, sn, sm_span, filt);
  runObj.PerformDIAProcess_core();
  List res = runObj.getResults();
  cout << "====== DIA data deconvolution done! ====== \n";
  return res;
}

// [[Rcpp::export]]
List dia_feature_preparation(NumericMatrix groupPkMtx, NumericMatrix chromPeaks, List peakidx){
  List ft_list(peakidx.size());
  for(int i=0; i<peakidx.size(); i++){ //peakidx.size()
    NumericVector pkidxs = peakidx[i];
    vector<double> this_ft_info(7);
    this_ft_info[0] = groupPkMtx(i,0); //mz
    this_ft_info[1] = groupPkMtx(i,1); //mz min
    this_ft_info[2] = groupPkMtx(i,2); //mz max
    
    this_ft_info[3] = groupPkMtx(i,3); //rt
    this_ft_info[4] = groupPkMtx(i,4); //rt min
    this_ft_info[5] = groupPkMtx(i,5); //rt max
    
    double baseline = 0, maxo, sn, tmpv;
    int kc = 0;
    for(int j=0; j<pkidxs.size(); j++){
      double thisidx = pkidxs[j];
      maxo = chromPeaks(thisidx-1,8);
      sn = chromPeaks(thisidx-1,9);
      tmpv = maxo/sn;
      if(isnan(tmpv)){
        continue;
      }
      baseline = baseline + maxo/sn;
      kc++;
    }
    this_ft_info[6] = baseline/kc;
    ft_list[i] = this_ft_info;
  }
  return ft_list;
}


/*** R
load("~/Github/OptiLCMS2ID/data/all_diaProcess_params.rda");
system.time(
asss <- PerformDIA_main(pm = pg.raw,
                swath = as.matrix(info.swath[[2]]),
                scanrt1 = scantime.ms1,
                scanrt2 = scantime.ms2,
                scanms1 = scan.ms1,
                scanms2 = scan.ms2,
                pkw_min = peakwidth[1],
                ppm = ppm.ms2.mtp, 
                sn = snthr,
                sm_span = span,
                filt = int.filter)
)

*/
