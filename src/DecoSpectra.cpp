#include "DecoSpectra.h"
#include "hclust_ultrafast.h"
#include "optim_ultra.h"

// [[Rcpp::export]]
List DecoSpectra(int idx_pg, 
                 List spectra_eics,
                 NumericVector peak_ms1,
                 int num_scantime,
                 int idx_apex_eic,
                 NumericVector info_pk_ms1,
                 double peakwidth_min,
                 double snthr,
                 bool is_dec_smoothed) {

  // Initialize the results set
  NumericMatrix spec_null;
  Rcpp::List spec_decon = List::create(_["spec"] = spec_null);
  
  List peaks_ms2_info(spectra_eics.size());
  List detectedEIC;
  IntegerVector peaks_ms2_info_names;
  List sharpnessList, sd_smoothList, lbList, apexList, apextList, rbList;
  LogicalVector is_peakVec;
  NumericVector blVec;
  
  int j = 0;
  NumericVector eic_int, idx_a, idx_s, idx_e, eic_pk_alli, eic_pk_rt;
  IntegerVector idx_eic_pk, idx_eic_pk_keep, idx_l, idx_r;
  
  // Detect EIC peak and calculate the characteristics (like sharpness..)
  for(int i = 0; i < spectra_eics.size(); i++){
    NumericMatrix eic = spectra_eics[i];
    NumericVector sharpnessVector, sd_smoothVector;
    detectedEIC = DetectPeaks(eic, peakwidth_min, num_scantime, idx_apex_eic, snthr, true, 1);
   
    if(detectedEIC.size() != 0){
      peaks_ms2_info[j] = detectedEIC;
      j++;
      peaks_ms2_info_names.push_back(i);

      List info_pk = detectedEIC[0]; // Information of the EIC
      eic_int = detectedEIC[1]; // this is the smoothed intensity of the EIC

      idx_a = info_pk[1]; // apex of the EIC
      idx_s = info_pk[0]; // lbvec
      idx_e = info_pk[2]; // rb_vec
      bool is_peak = info_pk[4]; // is a peak or not
      double bl = info_pk[3]; // baseline
      eic_pk_alli = eic( _ , 3); // intentisy ori EIC
      eic_pk_rt = eic( _ , 1); // rt of ori EIC
      
      // this for-loop is used to calculate all sub-eics based on apex vector
      for(int idx = 0; idx < idx_a.size(); idx++){
        idx_eic_pk = seq(idx_s[idx], idx_e[idx]);
        NumericMatrix eic_pk = eic( Range(idx_s[idx], idx_e[idx]) , _ );
        NumericVector eic_pk_i = eic_pk( _ , 3);
        idx_eic_pk_keep = idx_eic_pk[whichTrue(eic_pk_i > 0)];

        if(idx_a[idx] - 1 < idx_s[idx]){
          idx_l = idx_s[idx];
        } else {
          idx_l = seq(idx_s[idx], (idx_a[idx] - 1));
        }
        
        if(idx_e[idx] < (idx_a[idx] + 1)){
          idx_r = seq(idx_e[idx], idx_a[idx] + 1);
        } else {
          idx_r = seq(idx_a[idx] + 1, idx_e[idx]);
        }
        
        int iddxx = idx_a[idx];
        NumericVector sharpl(idx_l.size());
        NumericVector sharp2(idx_r.size());
        
        int ji = 0;
        for(int idxl : idx_l){
          double tmpval = floor((eic_int[iddxx] - eic_int[idxl])/(idx_a[idx] - idxl));
          sharpl[ji] = tmpval;
          ji++;
        }
        
        // Calculate sharpness of the peak (with the specific apex)
        double sharpness_l = median(sharpl);
        ji = 0;
        for(int idx2 : idx_r) {
          double tmpval2 = floor((eic_int[iddxx] - eic_int[idx2])/(idx2 - idx_a[idx]));
          sharp2[ji] = tmpval2;
          ji++;
        }
        double sharpness_r = median(sharp2);
        NumericVector sharpVec = {sharpness_l, sharpness_r};
        double sharpness = mean(sharpVec);
        
        // calculate smoothness
        NumericVector sd_smoth_val1 =  eic_pk_alli[idx_eic_pk_keep];
        
        NumericVector sd_smoth_val2 =  eic_int[idx_eic_pk_keep];
        double sd_smoth_val3 = max(sd_smoth_val2) - min(sd_smoth_val2);
        double sd_smooth = sd(sd_smoth_val1 - sd_smoth_val2)/sd_smoth_val3;

        sharpnessVector.push_back(sharpness);
        sd_smoothVector.push_back(sd_smooth);
      }
      
      // Summarize all results
      sharpnessList.push_back(sharpnessVector);
      sd_smoothList.push_back(sd_smoothVector);
      lbList.push_back(idx_s);
      apexList.push_back(idx_a);
      apextList.push_back(eic_pk_rt[idx_a]); // rt of apex
      rbList.push_back(idx_e);
      blVec.push_back(bl);
      is_peakVec.push_back(is_peak);
    }
  }
  
  // a commented results summary - for debugging
  // List info_pk_ms2_raw = List::create(_["sharpness"] = sharpnessList,
  //                                 _["sd_smooth"] = sd_smoothList,
  //                                 _["lb"] = lbList,
  //                                 _["apex"] = apexList,
  //                                 _["apex_t"] = apextList,
  //                                 _["rb"] = rbList,
  //                                 _["bl"] = blVec,
  //                                 _["is_peak"] = is_peakVec);
  // return info_pk_ms2_raw;

  //// Filtering out sd_smooth > 0.35 or not a real peak (is_peak == false)
  List nr_rmList;
  LogicalVector logVec, is_simpVec;
  
  for(int i = 0; i < is_peakVec.size(); i++){
    NumericVector sd_smoothVec = sd_smoothList[i];
    LogicalVector nr_rm_vec(sd_smoothVec.size());
    int k = 0; // a counter of peak number 
    for(int j = 0; j < sd_smoothVec.size(); j++){
      if(sd_smoothVec[j] > 0.35){
        nr_rm_vec[j] = true; // will remove
      } else {
        k++;
        nr_rm_vec[j] = false;
      }
    }
    nr_rmList.push_back(nr_rm_vec);
    if(is_true(all(nr_rm_vec)) | !is_peakVec[i]){
      logVec.push_back(true);
    } else {
      logVec.push_back(false);
    }
    if(k>1){
      is_simpVec.push_back(false);
    } else {
      is_simpVec.push_back(true);
    }
  }
  
  if(is_true(all(logVec))) {
    // if all peaks are not real or all peaks should be removed and return nothing
    return spec_decon;
  }
  
  // results summary : characteristics of all sub-eics in each ms2 EIC peaks
  List info_pk_ms2 = List::create(_["sharpness"] = sharpnessList,
                                  _["sd_smooth"] = sd_smoothList,
                                  _["lb"] = lbList,
                                  _["apex"] = apexList,
                                  _["apex_t"] = apextList,
                                  _["rb"] = rbList,
                                  _["is_simp"] = is_simpVec,
                                  _["nr_rm"] = nr_rmList,
                                  _["bl"] = blVec,
                                  _["is_peak"] = is_peakVec);
  
  // initialize variables for cleaning
  NumericVector rt_vec, tmprtvec, tmprtvec2;
  List info_pk_ms2_clean;
  LogicalVector tmpLogVec;
  CharacterVector names_vec;
  IntegerVector names_vec_val1;
  List eic_peak, tmpList, tmpList2, tmpList3;
  // to generate "info_pk_ms2_clean" which includes all clean ms2 peaks (sub_eics) of all EICs
  for(int i = 0; i < is_peakVec.size(); i++){
    tmpLogVec = nr_rmList[i]; // a logical vector containing information on removing or not
    if(!is_peakVec[i]) continue;
    for(int j = 0; j < tmpLogVec.size(); j++){
      if(!tmpLogVec[j]){
        IntegerVector idx_eic, tmpidvec1, tmpidvec2;
        // Process peak names generation
        tmprtvec = apextList[i];
        rt_vec.push_back(tmprtvec[j]);
        names_vec.push_back(to_string(peaks_ms2_info_names[i]) + "_" + to_string(j));
        names_vec_val1.push_back(peaks_ms2_info_names[i]);
        // Process all eic peaks info
        tmpidvec1 = lbList[i];
        tmpidvec2 = rbList[i];
        idx_eic = seq(tmpidvec1[j], tmpidvec2[j]);
        
        tmpList2 = peaks_ms2_info[i];
        tmprtvec2 = tmpList2[1];
        
        tmpList = List::create(_["rt"] = idx_eic, // idxs of sub-eic
                               _["i_s"] = tmprtvec2[idx_eic]);
        eic_peak.push_back(tmpList);
        
        // Process info_pk_ms2_clean List [including sharpness + is_simp info]
        NumericVector shrpItems = sharpnessList[i];
        NumericVector rbItems = rbList[i];
        NumericVector lbItems = lbList[i];
        tmpList3 = List::create(_["sharpness"] = shrpItems[j],
                                _["is_simp"] = is_simpVec[i],
                                _["bl"] = blVec[i],
                                _["rb"] = rbItems[j],
                                _["lb"] = lbItems[j]);
        info_pk_ms2_clean.push_back(tmpList3);
      }
    }
  }
  // return info_pk_ms2_clean; // used for debugging
  
  // handle ms 1 peaks
  IntegerVector idx_eic1 = seq(info_pk_ms1[0], info_pk_ms1[3]); // two boundaries of ms1 pk -> as index
  tmprtvec2 = peak_ms1;//(_,4);
  tmpList = List::create(_["rt"] = idx_eic1,
                         _["i_s"] = tmprtvec2);
  eic_peak.push_back(tmpList);

  // need to add rt of MS1 peak into rt_vec and names_vec
  rt_vec.push_back(info_pk_ms1[2]); // retention time
  names_vec.push_back("0"); // ms1 peak marked as "0"
  IntegerVector names_vec_idx = seq(0, names_vec.size()-1); // idxes of all names

  // h-clustering of retention time (and unique the clusters' name)
  NumericVector r_cl_rt = auto_hclust(rt_vec);
  NumericVector r_cl_rt_uni = unique(r_cl_rt);

  // check the member of the cluster to identify the ones with only 1 peak (to remove)
  List cl_rt_nm, cl_rt_nm_idx;
  int idx_cl_rt_main = -1;
  CharacterVector names_vec2rm;
  IntegerVector names_vec_idx2rm;
  for(int num = 0; num < r_cl_rt_uni.size(); num++){
    CharacterVector thisCluster;
    IntegerVector thisCluster_idx;
    for(int j = 0; j < r_cl_rt.size(); j++){
      if(r_cl_rt[j] == r_cl_rt_uni[num]){
        thisCluster.push_back(names_vec[j]);
        thisCluster_idx.push_back(names_vec_idx[j]);
        if(names_vec[j] == "0") idx_cl_rt_main = num;
      }
    }
    if(thisCluster.size()>1){
      cl_rt_nm.push_back(thisCluster); // List, "names" of rt-based clusters | one member in the list is a cluster (charactervector)
      cl_rt_nm_idx.push_back(thisCluster_idx); // List, idx of the "names"
    } else {
      names_vec2rm.push_back(thisCluster[0]);
      names_vec_idx2rm.push_back(thisCluster_idx[0]);
    }
  }
  
  // SHOULD NOT REMOVE THEM DIRECTLY, BECAUSE IT WILL AFFECT THE INDEX
  // Use this if condition to trim the data and vectors for following processing
  IntegerVector names_vec_idx0 = clone(names_vec_idx);
  if(names_vec2rm.size() > 0){
    for(String nm2rm : names_vec2rm){
      if(nm2rm == "0") {continue;}
      for(int nmidx : names_vec_idx){
        if(nm2rm == names_vec[nmidx]){
          names_vec.erase(nmidx);
          names_vec_idx0.erase(nmidx);
          break;
        }
      }
    }
    
    if(names_vec.size() == 0) {
      // nothing left, stop processing
      return spec_decon;
    }
    
    names_vec_idx = seq(0, names_vec.size()-1); // create a new names_vec_idx

    // Correct cl_rt_nm_idx
    for(int ci = 0; ci < cl_rt_nm_idx.size(); ci++){
      IntegerVector tmpCi = cl_rt_nm_idx[ci];
      for(int tci = 0; tci < tmpCi.size(); tci++){
        int tciVal = tmpCi[tci];
        tmpCi[tci] = whichTrue1(tciVal == names_vec_idx0);
      }
      cl_rt_nm_idx[ci] = tmpCi;
    }
    
    //correct eic_peak + info_pk_ms2_clean + names_vec_val1
    for(int ww : names_vec_idx2rm.sort(true)) {
      eic_peak.erase(ww);
      if(ww > (info_pk_ms2_clean.size()-1)){
        continue;
      }
      info_pk_ms2_clean.erase(ww);
      names_vec_val1.erase(ww);
    }
  }
  
  // no any rt-based clusters found, neither "idx_cl_rt_main" (ms1 peak) remains
  if((cl_rt_nm.size() == 0) | (idx_cl_rt_main == -1)){
    return spec_decon;
  }
  
  // starting check clusters of components  
  List cl_components_idx, cl_components;
  
  // for-loop with retention time-based name cluster idx
  for(int i = 0; i < cl_rt_nm_idx.length(); i++){
    IntegerVector thisCluster_idx = cl_rt_nm_idx[i];
    if(thisCluster_idx.size() < 2) continue;
    NumericMatrix distant_mtx( thisCluster_idx.size() , thisCluster_idx.size() ); // initialize distance matrix

    // calculate peak corrrelation similarity within this rt-cluster
    int idxp1, idxp2;
    for(idxp1 = 0; idxp1 < thisCluster_idx.size(); idxp1++){
      for(idxp2 = idxp1+1; idxp2 < thisCluster_idx.size(); idxp2++){
        // extract the idxes of sub-eic peaks from this specific cluster
        int idx1 = thisCluster_idx[idxp1];
        int idx2 = thisCluster_idx[idxp2];

        List eic_p1 = eic_peak[idx1];
        List eic_p2 = eic_peak[idx2];
        
        NumericVector nv1, nv2, nv3, nv4;
        nv1 = eic_p1[0]; // rt
        nv2 = eic_p1[1]; // intensity
        nv3 = eic_p2[0]; // rt
        nv4 = eic_p2[1]; // intensity
        
        // construct new eic peak
        NumericMatrix eic_peak1 = cbind(nv1, nv2);
        NumericMatrix eic_peak2 = cbind(nv3, nv4);
        
        // calculate similarity distance of two sub-eics of correlations
        double distant = GetDistantP(eic_peak1, eic_peak2);
        if(distant > 1){
          distant_mtx(idxp1, idxp2) = 1.0;
        } else {
          distant_mtx(idxp1, idxp2) = distant;
        }

      }
    }
    
    // further h-clustering the sub-eics based on peak similarity (correlations from the steps above)
    NumericVector r_cl_pk = matrix_hclust(distant_mtx);

    // Names_vec AND names_vec_idx should be split to make sure that 
    // the idx and names are correct for multiple "cl_rt_nm_idx" Lists
    CharacterVector new_names_vec = names_vec[thisCluster_idx];
    IntegerVector new_names_vec_idx = names_vec_idx[thisCluster_idx];

    // Unique and sort clusters
    NumericVector r_cl_pk_uni = unique(r_cl_pk);
    r_cl_pk_uni = r_cl_pk_uni.sort();
    
    // generate multiple componenets within the same rt-cluster (pushback names and idxs)
    List cl_pk_nm, cl_pk_nm_idx;
    for(int num=0; num < r_cl_pk_uni.size(); num++){
      CharacterVector thispkCluster;
      IntegerVector thispkCluster_idx;
      for(int j=0; j < r_cl_pk.size(); j++){
        if(r_cl_pk[j] == r_cl_pk_uni[num]){
          thispkCluster.push_back(new_names_vec[j]);
          thispkCluster_idx.push_back(new_names_vec_idx[j]);
        }
      }
      if(thispkCluster.size()>1){
        cl_pk_nm.push_back(thispkCluster);
        //cl_pk_nm_idx.push_back(thispkCluster_idx);
        cl_components.push_back(thispkCluster); // names of this component
        cl_components_idx.push_back(thispkCluster_idx); //idx of the names in this component
      }
    }

    if((cl_pk_nm.size() == 0)){
      continue;
    }
    // cl_components.push_back(cl_pk_nm_idx);
  }

  // Initiate variables to select model/main peaks from components
  IntegerVector nm_pk_main_idx, nm_mpk_idx;
  NumericVector dist_main;
  CharacterVector nm_mpk, dist_main_names;
  int nm_mpk_main_idx = -1, ms1_pkidx = -1; // nm_mpk_main_idx, idx of main model peak
  
  // for-loop to circulate all components
  for(int i = 0; i< cl_components.size(); i++){
    IntegerVector nm_pk_idx = cl_components_idx[i];
    CharacterVector nm_pk = cl_components[i];

    bool bl1 = false; // ms1 peak included (true) or not (false)
    bool bl2 = false; // there is only one peak in this component (true) or not (false)
    
    // to see if ms1 peak included in this component
    for(int j =0; j < nm_pk.size(); j++){
      if(nm_pk[j] == "0") {
        bl1 = true;
        ms1_pkidx = nm_pk_idx[j];
        break;
      }
    }
    
    if(bl1 & (nm_pk.size() == 1)){bl2 = true;}
    if(bl1){
      if(!bl2){
        // size > 1, there is some other member(s) in nm_pk
        for(int jj =0; jj < nm_pk.size(); jj++){
          if(nm_pk[jj] != "0"){
            nm_pk_main_idx.push_back(nm_pk_idx[jj]);
            // extract ms1 peak and other ms2 peak to see which is the most similar one to ms1.
            List eic_p0 = eic_peak[ms1_pkidx];
            List eic_px = eic_peak[nm_pk_idx[jj]];
            
            NumericVector nv1, nv2, nv3, nv4;
            nv1 = eic_p0[0];
            nv2 = eic_p0[1];
            nv3 = eic_px[0];
            nv4 = eic_px[1];
            NumericMatrix eic_peak0 = cbind(nv1, nv2);
            NumericMatrix eic_peakx = cbind(nv3, nv4);
            
            double distant = GetDistantP(eic_peak0, eic_peakx);
            dist_main.push_back(distant);
            dist_main_names.push_back(nm_pk[jj]);
          }
        }
        int pk_min_idx = which_min(dist_main);
        nm_pk_idx = nm_pk_main_idx;
        nm_mpk_main_idx = nm_pk_idx[pk_min_idx];
      } else {
        // there is only one member ("0") in nm_pk
        nm_mpk_main_idx = nm_pk_idx[0];
        nm_pk_main_idx.push_back(nm_pk_idx[0]);
      }
      
      nm_mpk_idx.push_back(nm_mpk_main_idx);
      
    } else {
      // cout << "Other option <---\n";
      // no ms1 peak included, need simp and sharp info here 
      NumericVector shrpVec0, shrpVec;
      IntegerVector idxVec0, idxVec;
      
      for(int d : nm_pk_idx){
        List thisList = info_pk_ms2_clean[d];
        int isSimp = thisList[1];
        shrpVec0.push_back(thisList[0]);
        idxVec0.push_back(d);
        if(isSimp == 1){
          shrpVec.push_back(thisList[0]);
          idxVec.push_back(d);
        }
      }
      // select the one with the highest sharpness as the model peak of this component
      if(shrpVec.size() > 0){
        int val1 = idxVec[which_max(shrpVec)];
        nm_mpk_idx.push_back(val1);
      } else {
        int val1 = idxVec0[which_max(shrpVec0)];
        nm_mpk_idx.push_back(val1);
      }
    }
  }
  
  // organize name | mpk <- model peak
  for(int nm : nm_mpk_idx){
    nm_mpk.push_back(names_vec[nm]);
  }

  if(nm_pk_main_idx.size() == 0){
    // no main peak found, use the MS1 peak directly
    for(int kw : names_vec_idx){
      if(names_vec[kw] == "0"){
        nm_mpk_main_idx = kw;
        nm_pk_main_idx = kw;
        ms1_pkidx = kw;
        cl_components.push_back("0");
        cl_components_idx.push_back(kw);
        nm_mpk.push_back("0");
        nm_mpk_idx.push_back(kw);
        
        break;
      }
    }
  }
  
  // summarize simple ms2 peaks here | no need to do deconvolution
  IntegerVector nm_pk_simp_idx; // idx of all simple ms2 peaks
  for(int i : nm_pk_main_idx){
    if(i == info_pk_ms2_clean.size()){
      continue;
    }
    List l1 = info_pk_ms2_clean[i];
    int simp = l1[1];
    if(simp == 1){nm_pk_simp_idx.push_back(i);}
  }
  
  NumericVector spec_smp_mzmzvec;
  NumericVector spec_smp_intsvec;
  NumericVector spec_smp_distvec;
  NumericMatrix spec_smp;
  if(nm_pk_simp_idx.size() > 0){
    for(int nm : nm_pk_simp_idx){
      int nm_eic0 = names_vec_val1[nm];
      int nm_eic = whichTrue1(peaks_ms2_info_names == nm_eic0);

      double mz, intensity, dis;
      NumericMatrix thiseic = spectra_eics[nm_eic0];
      // mz
      mz = thiseic( idx_apex_eic , 2);
      // intensity
      List thisL1 = peaks_ms2_info[nm_eic];
      List thisL10 = thisL1[0];
      double bl0 = thisL10[3];
      intensity = thiseic( idx_apex_eic , 3) - bl0;
      // distant
      String nm_eic_name = names_vec[nm];
      int res_dis_idx = -1;
      for(int dmnm_idx = 0; dmnm_idx < dist_main_names.size(); dmnm_idx++){
        if(dist_main_names[dmnm_idx] == nm_eic_name) res_dis_idx = dmnm_idx;
      }

      if(res_dis_idx != -1)
        dis = dist_main[res_dis_idx];
      else
        dis = 0.0;
      
      // SUMMARY --->
      spec_smp_mzmzvec.push_back(mz);
      spec_smp_intsvec.push_back(intensity);
      spec_smp_distvec.push_back(dis);
    }
    spec_smp = cbind(spec_smp_mzmzvec,
                     spec_smp_intsvec,
                     spec_smp_distvec);
  }
  
  // prepare model ms2 peak eic
  // List eic_model_pk;
  vector<NumericVector> eic_model_pk;
  int k = -1;
  List l2 = peaks_ms2_info[0];
  NumericVector nmvec1 = l2[1];
  int lengthVal = nmvec1.size(); // nrow of model peak matrix (number of centroids of all model peaks)
  
  // circulate all model peaks
  for(int np : nm_mpk_idx){
    NumericVector eic;
    k++;
    if(nm_mpk[k] == "0"){
      eic = rep(0, lengthVal);
      int intval1 = (int)info_pk_ms1[0];
      int intval2 = (int)info_pk_ms1[3];
      double blval = (double)info_pk_ms1[4];
      if(intval2 >= eic.size()){
        NumericVector ms1_eic_vec_x = peak_ms1 - blval;
        eic[seq(intval1, intval2)] = ms1_eic_vec_x[seq(0,eic.size()-1)];
      } else {
        eic[seq(intval1, intval2)] = peak_ms1 - blval;
      }
    } else {
      if(names_vec_val1.size() <= np){continue;}
      int intvl0 = names_vec_val1[np];
      int intvl1 = whichTrue1(intvl0 == peaks_ms2_info_names);

      List thisL1 = peaks_ms2_info[intvl1];
      NumericVector eic0 = thisL1[1];
      eic = rep(0, eic0.size());

      List thisL2 = info_pk_ms2_clean[np];
      int intval1 = thisL2[4];
      int intval2 = thisL2[3];
      double blval = thisL2[2];
      eic[seq(intval1, intval2)] = eic0[seq(intval1, intval2)] - blval;
    }
    // remove negative values
    for(int i = 0; i < eic.size(); i++){
      if(eic[i] < 0) eic[i] = 0;
    }

    // zero eic
    if(max(eic) == 0) {continue;}; //return spec_decon; // Not sure should skip this eic or simply give up
    NumericVector eic_norm = eic/max(eic);
    eic_model_pk.push_back(eic_norm);
  }
  
  // remove all idx (centroids) with the sum of all idx (rows) is 0
  IntegerVector idx_sep;
  double rowsum;
  for(int j = 0; j < lengthVal; j++){
    rowsum = 0;
    for(int i = 0; i < eic_model_pk.size(); i++){
      NumericVector subL = eic_model_pk[i];
      double subL_val = subL[j];
      rowsum = rowsum + subL_val;
      if(rowsum > 0){
        break;
      }
    }
    
    if(rowsum == 0){
      idx_sep.push_back(j);
    }
  }
  
  int i_start, i_end;
  if(idx_sep.size() >0 ){
    idx_sep.push_back(-1);
    idx_sep.push_back(lengthVal);
    idx_sep = unique(idx_sep);
    idx_sep.sort();

    IntegerVector numVec1 = idx_sep[idx_sep < idx_apex_eic];
    IntegerVector numVec2 = idx_sep[idx_sep > idx_apex_eic];

    i_start = max(numVec1) + 1;
    i_end = min(numVec2) - 1;
  } else {
    i_start = 0;
    i_end = lengthVal-1;
  }
  
  // model peak matrix construction
  NumericMatrix mpk_mtx( (i_end + 1 - i_start), 0);
  IntegerVector mpk_mtx_idx;
  for(int i = 0; i < eic_model_pk.size(); i++){
    NumericVector subL = eic_model_pk[i];
    subL = subL[seq(i_start, i_end)];
    bool bl3 = is_true(any(subL > 0));
    if(bl3){
      mpk_mtx = cbind(mpk_mtx, subL);
      mpk_mtx_idx.push_back(nm_mpk_idx[i]);
    }
  }
  // no model peak to deconvolve EIC, return retults
  if(mpk_mtx.ncol() < 2){
    spec_decon = List::create(_["spec"] = spec_smp,
                              _["spec_smp"] = spec_smp,
                              _["spec_cpl"] = R_NilValue);
    return spec_decon;
  }
  
  // used for spec matrix later [simple components]
  NumericVector spec_nv0;
  NumericVector spec_nv1;
  NumericVector spec_nv2;
  if(spec_smp.ncol() > 0){
    spec_nv0 = spec_smp( _, 0);
    spec_nv1 = spec_smp( _, 1);
    spec_nv2 = spec_smp( _, 2);
  }
// 
//   cout << "mpk_mtx_idx     -> "<< mpk_mtx_idx << endl;
//   cout << "nm_mpk_main_idx -> "<< nm_mpk_main_idx << endl;
  IntegerVector idx_mpk_main = whichTrue(mpk_mtx_idx == nm_mpk_main_idx);
  
  int new_ms1pk_apex = info_pk_ms1[1] - i_start;// + 1;

  IntegerVector nm_eic_dec;
  
  for(int i = 0; i < info_pk_ms2_clean.size(); i++){
    List tmpL4 = info_pk_ms2_clean[i];
    int tmpi = tmpL4[1];
    if(tmpi == 0){
      nm_eic_dec.push_back(names_vec_val1[i]);
    }
  }
  nm_eic_dec = unique(nm_eic_dec);
  nm_eic_dec.sort();

  NumericMatrix spec;
  NumericMatrix spec_cpl(3, 0);
  
  if(nm_eic_dec.length() > 0){
    List eic_pre_dec;
    List info_decomposite;
    for(int idx_eic : nm_eic_dec) {
      int int_idx = whichTrue1(idx_eic == peaks_ms2_info_names);
      NumericVector lbVec = lbList[int_idx];
      NumericVector rbVec = rbList[int_idx];
      IntegerVector i_scan = seq(min(lbVec), max(rbVec));
      
      NumericVector eicVec (lengthVal), eicVec0;
      if(is_dec_smoothed){
        List tmpL5 = peaks_ms2_info[int_idx];
        eicVec0 = tmpL5[1];
      } else {
        NumericMatrix tmpL6 = spectra_eics[idx_eic];
        eicVec0 = tmpL6( _ , 3);
      }
      
      double blVal0 = blVec[int_idx];
      eicVec[i_scan] = eicVec0[i_scan];
      eicVec = eicVec - blVal0;
      
      for(int u = 0; u< eicVec.size(); u++){
        if(!is_true(any(u == i_scan))){
          eicVec[u] = 0;
        }
        if(eicVec[u] < 0){
          eicVec[u] = 0;
        }
      }
      
      eicVec = eicVec[seq(i_start, i_end)];
      eic_pre_dec.push_back(eicVec);

      // cout << "line 697 --> idx_mpk_main -> " << idx_mpk_main << endl;
      // cout << "eicVec -> " << eicVec << endl;
      // cout << "mpk_mtx --> " << mpk_mtx.ncol() << endl;
      // cout << "mpk_mtx --> " << mpk_mtx.nrow() << endl;
      
      if(idx_mpk_main.size() == 0){
        return spec_decon;
      }
      double optiRes = optim_ultra(mpk_mtx, eicVec, idx_mpk_main[0]);
      
      NumericVector tmpMPKVec = mpk_mtx( _ , idx_mpk_main[0] ) * optiRes;
      NumericMatrix dtm = spectra_eics[idx_eic];
      NumericVector infoCompVec = 
        NumericVector::create(_["mz"]=dtm(new_ms1pk_apex, 2), 
                              _["intensity"]= tmpMPKVec[new_ms1pk_apex],//tmpMPK(new_ms1pk_apex, idx_mpk_main[0]),
                              _["dist"]=0);
      spec_cpl = cbind(spec_cpl, infoCompVec);
    }
    
    spec_cpl = transpose(spec_cpl);
    LogicalVector is_keep_col =  spec_cpl( _ ,1) > 0;
    CharacterVector ch = colnames(spec_cpl);
    
    if(is_true(any(is_keep_col))){
      // IntegerVector invecs = whichTrue(is_keep_col);
      NumericVector mzVecx;
      NumericVector intensityVecx;
      for(int i = 0; i < is_keep_col.size(); i++){
        if(is_keep_col[i]){
          mzVecx.push_back(spec_cpl(i,0));
          intensityVecx.push_back(spec_cpl(i,1));
          // used for spec matrix
          spec_nv0.push_back(spec_cpl(i,0));
          spec_nv1.push_back(spec_cpl(i,1));
          spec_nv2.push_back(0);
        }
      }
      
      NumericVector distVecx (intensityVecx.size());
      spec_cpl = cbind(mzVecx, intensityVecx, distVecx);
      colnames(spec_cpl) = ch;
      
      // used for spec_cpl
      // spec = cbind(spec_nv0, spec_nv1, spec_nv2);
      // sort spec
      NumericVector spec_nv0_sorted;
      NumericVector spec_nv1_sorted;
      NumericVector spec_nv2_sorted;
      
      spec_nv0_sorted = clone(spec_nv0); //spec_nv0.sort();
      spec_nv0_sorted.sort();
      int idxVal;
      for(int mz_idx = 0; mz_idx < spec_nv0_sorted.size(); mz_idx++){
        idxVal = whichTrue1(spec_nv0_sorted[mz_idx] == spec_nv0);
        //cout << "idxVal --> " << idxVal << endl;
        spec_nv1_sorted.push_back(spec_nv1[idxVal]);
        spec_nv2_sorted.push_back(spec_nv2[idxVal]);
      }
      spec = cbind(spec_nv0_sorted, spec_nv1_sorted, spec_nv2_sorted);
    } else {
      spec_cpl = NumericMatrix (0,0);
    }
  } else {
    spec = spec_smp;
    spec_cpl = NumericMatrix (0,0);
  }
  
  if((spec_smp.ncol() == 0) & (spec_cpl.ncol() == 0)){
    return spec_decon;
  } else if((spec_smp.ncol() == 0) & (spec_cpl.ncol() != 0)){
    spec = spec_cpl;
  } else if((spec_smp.ncol() != 0) & (spec_cpl.ncol() == 0)){
    spec = spec_smp;
  }
  
  spec_decon = List::create(_["spec"] = spec,
                            _["spec_smp"] = spec_smp,
                            _["spec_cpl"] = spec_cpl);
  return spec_decon;
}




/*** R
# info.pk.ms1_new <- info.pk.ms1;
# info.pk.ms1_new[1,2] <- info.pk.ms1[1,2] - 1;
# info.pk.ms1_new[1,1] <- info.pk.ms1_new[1,1]-1
# info.pk.ms1_new[1,4] <- info.pk.ms1_new[1,4]-1
# system.time(
# spec.decon <- DecoSpectra(idx.pg,
#                           ms2.eic.ext,
#                           peak.ms1.smooth[,5],
#                           length(idx.ms2.ext),
#                           idx_apex_eic = idx.apex.ms1 - 1,
#                           info_pk_ms1 = as.numeric(info.pk.ms1_new),
#                           peakwidth = peakwidth[1],
#                           snthr = snthr,
#                           is_dec_smoothed = is.dec.smoothed)
#)
*/
