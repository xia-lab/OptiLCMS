#' PerformInfusionImport
#'
#' @param filesPath 
#'
#' @return mSet Object
#' @export
#'
#' @examples todo
PerformInfusionImport <- function(filesPath = NA){
   
   require(mzR)
   
   msFiles_list <- lapply(filesPath, function(x){mzR::openMSfile(x, backend = "pwiz")})
   pks_list <- lapply(msFiles_list, mzR::spectra)
   hdr_list <- lapply(msFiles_list, mzR::header)
   
   # Process ion mode (1 is positive, 0 is negative)
   all_polarities <- sapply(1:length(filesPath), function(x){
     unique(hdr_list[[x]][["polarity"]])
   })
   ion_mode <- unique(unlist(all_polarities))
   if(length(ion_mode)>1){
     stop("Scans with mutiple polarities are not allowed! Please keep one polarity at a time!")
   }
   
   # Split scanrts_ms1, scanrts_ms2, scan_ms1, scan_ms2 based on files' information
   scanrts_ms1_List = 
     scan_ms1_List = list();
   
   # Process the data of MS1 level
   scanrts_ms1_List <- lapply(1:length(filesPath), function(x){
     hdr_list[[x]][["retentionTime"]][hdr_list[[x]][["msLevel"]] == 1]
   })
   scan_ms1_idx <- lapply(1:length(filesPath), function(x){
     hdr_list[[x]][["msLevel"]] == 1
   })
   scan_ms1_List <- lapply(1:length(filesPath), function(x){
     pks_list[[x]][scan_ms1_idx[[x]]]
   })
   
   MSn_data <- list(scanrts_ms1 = scanrts_ms1_List,
                    scan_ms1 = scan_ms1_List,
                    fileNames = basename(filesPath))
   
   #if(is.null(mSet)){
   mSet <- new("mSet")
   mSet@MSnData <- MSn_data
   return(mSet)
 }
 

#' PerformInusionFeatureDection
#'
#' @param mSet 
#' @param method 
#' @param noise 
#' @param snthresh 
#' @param ppm 
#' @param scales 
#' @param minFraction 
#' @param minSamples 
#' @param ncores 
#'
#' @return mSet Object
#' @export
#'
#' @examples todo
PerformInusionFeatureDection <- function(mSet, method = "cwt", noise = NA, 
                                         snthresh = 3, ppm = 10, 
                                         scales = c(1,4,9), 
                                         groups = NULL,
                                         ncores = 1L){
  
  mSet@MSnData$scan_ms1 -> scan_ms1
  mSet@MSnData[["fileNames"]] -> spec_files;
  
  if(ncores == 1L){
    require("MassSpecWavelet")
    peaks_list <- rep(list(NA), length(scan_ms1))
    for(i in 1:length(scan_ms1)){
      message("Processing spectra file == ", spec_files[i], " == ", appendLF = F)
      this_file <- scan_ms1[[i]];
      res <- lapply(1:length(this_file), function(x){
        #cat(x, "\n")
        y<-this_file[[x]]
        ints <- y[,2]
        mzs <- y[,1]
        if(method =="cwt"){
          peaks <- infusionDectectCWT(ints, mzs, scales = scales, snthresh = snthresh)
        } else {
          mzdiff <- 500*ppm*1e-6
          peaks <- infusionDectectKDE(ints, mzs, noise = noise, snthresh = snthresh, ppm = ppm, mzdiff = mzdiff)
        }
        
        return(cbind(peaks, sample = x))
      })
      
      peaks_table <- do.call(rbind, res)
      if(length(which(is.na(peaks_table[,1])))!=0){
        peaks_table <- peaks_table[-which(is.na(peaks_table[,1])), ]
      }
      
      leng_log <- nrow(peaks_table)
      if(length(leng_log) == 0){
        peaks_table <- as.data.frame(t(as.matrix(peaks_table)))
      } else {
        peaks_table <- as.data.frame(peaks_table)
      }
      
      message("...", appendLF = F)
      # need to do a merging of detected peaks within same file - mzCluster
      res_merged <- do_groupPeaks_mzClustx(peaks_table,
                                           sampleGroups = rep(1, length(this_file)),
                                           ppm = ppm,
                                           absMz = 0,
                                           minFraction = 0,
                                           minSamples = 0)
      
      res_merged[["featureDefinitions"]] -> merged_peak0
      res_merged[["peakIndex"]] -> merged_peak_idx
      col_idx <- which(colnames(peaks_table) == "sn")
      
      merged_peak_list <- lapply(merged_peak_idx, function(x){
        peaks_table_sub<- peaks_table[x,]
        row_idx <- which.max(peaks_table_sub[,col_idx])[1]
        return(peaks_table_sub[row_idx,])
      })
      merged_peak1 <- do.call(rbind, merged_peak_list)
      merged_peak1$sample <- i
      rownames(merged_peak1) <- NULL
      peaks_list[[i]] <- merged_peak1
      message(" ... Done!", appendLF = T)
    }
  } else {
    # use multiple cores
    require("MassSpecWavelet")
    peaks_list <- rep(list(NA), length(scan_ms1))
    for(i in 1:length(scan_ms1)){
      message("Processing spectra file == ", spec_files[i], " == ", appendLF = F)
      this_file <- scan_ms1[[i]];
      # res <- bplapply(1:length(this_file), function(x){
      #   y<-this_file[[x]]
      #   ints <- y[,2]
      #   mzs <- y[,1]
      #   peaks <- infusionDectectCWT(ints, mzs)
      #   return(cbind(peaks, sample = x))
      # })
      require(parallel);
      cl <- makeCluster(getOption("cl.cores", ncores))
      clusterExport(cl, c("infusionDectectCWT", "this_file", "infusionDectectKDE",
                          "noise", "snthresh", "ppm", "scales",
                          "sumIntos", "maxIntos","do_groupPeaks_mzClustx", ".fix_mz_clust_peaks",
                          "mzClustGeneric", "mzClust_hclust_cpp", "R_mzClust_hclust_rcpp"), envir = environment())
      res <- list()
      
      if(method =="cwt"){
        res <- parLapply(cl, 
                         1:length(this_file), 
                         function(x){
                           y<-this_file[[x]]
                           ints <- y[,2]
                           mzs <- y[,1]
                           peaks <- infusionDectectCWT(ints, mzs, scales = scales, snthresh = snthresh)
                           return(cbind(peaks, sample = x))
                         })
      } else {
        res <- parLapply(cl, 
                         1:length(this_file), 
                         function(x){
                           y<-this_file[[x]]
                           ints <- y[,2]
                           mzs <- y[,1]
                           mzdiff <- 500*ppm*1e-6
                           peaks <- infusionDectectKDE(ints, mzs, noise = noise, snthresh = snthresh, ppm = ppm, mzdiff = mzdiff)
                           return(cbind(peaks, sample = x))
                         })
      }
      stopCluster(cl)
      
      peaks_table <- do.call(rbind, res)
      if(length(which(is.na(peaks_table[,1])))!=0){
        peaks_table <- peaks_table[-which(is.na(peaks_table[,1])), ]
      }
      if(nrow(peaks_table)==1){
        peaks_table <- as.data.frame(t(as.matrix(peaks_table)))
      } else {
        peaks_table <- as.data.frame(peaks_table)
      }
      
      message("...", appendLF = F)
      # need to do a merging of detected peaks within same file - mzCluster
      res_merged <- do_groupPeaks_mzClustx(peaks_table,
                                           sampleGroups = rep(1, length(this_file)),
                                           ppm = ppm,
                                           absMz = 0,
                                           minFraction = 0,
                                           minSamples = 0)
      
      res_merged[["featureDefinitions"]] -> merged_peak0
      res_merged[["peakIndex"]] -> merged_peak_idx
      col_idx <- which(colnames(peaks_table) == "sn")
      
      merged_peak_list <- lapply(merged_peak_idx, function(x){
        peaks_table_sub<- peaks_table[x,]
        row_idx <- which.max(peaks_table_sub[,col_idx])[1]
        return(peaks_table_sub[row_idx,])
      })
      merged_peak1 <- do.call(rbind, merged_peak_list)
      merged_peak1$sample <- i
      rownames(merged_peak1) <- NULL
      peaks_list[[i]] <- merged_peak1
      message(" ... Done!", appendLF = T)
    }
  }
  
  peaks_list_all <- do.call(rbind, peaks_list)
  if(length(which(is.na(peaks_list_all$mz)))>0){
    peaks_list_all <- peaks_list_all[!is.na(peaks_list_all$mz), ]
  }
  peaks_list_all <- peaks_list_all[peaks_list_all$sn>= snthresh,]

  # peak grouping across different samples and groups
  message("MZ alignment across samples .. ", appendLF = F)
  if(is.null(groups)){
    sampleGroups_idx <- unique(peaks_list_all$sample)
  } else {
    
    spec_files_nms <- basename(spec_files)
    sampleGroups_idx <- sapply(peaks_list_all$sample, function(x){
      sfn <- spec_files_nms[x]
      idx <- sfn == sampleGroups[,1]
      sampleGroups[idx, 2]
    })
    if(is.list(sampleGroups_idx)){
      sampleGroups_idx <- unlist(sampleGroups_idx)
    }
  }
  res_merged_files <- do_groupPeaks_mzClustx(peaks_list_all,
                                       sampleGroups = sampleGroups_idx,
                                       ppm = ppm,
                                       absMz = 0,
                                       minFraction = 0,
                                       minSamples = 0)
  message(" Done! ", appendLF = T)
  
  peak_mzs <- res_merged_files[["featureDefinitions"]][,c(1:3)]
  
  peak_ints_list <- lapply(1:nrow(peak_mzs), function(x){
    row_idx <- res_merged_files$peakIndex[[x]]
    peak_sub <- peaks_list_all[row_idx, ]
    df <- data.frame(matrix(NA, ncol = length(spec_files)))
    if(any(table(peak_sub$sample)>1)){
      unique_col_idx <- unique(peak_sub$sample)
      unique_into_val <- vapply(unique_col_idx, FUN = function(x){
        vals <- peak_sub$into[x == peak_sub$sample]
        vals[which.max(vals)]
        }, double(1L))
      df[1, unique_col_idx] <- unique_into_val
    } else {
      df[1, peak_sub$sample] <- peak_sub$into
    }
    df
  })
  
  peak_ints_dt <- do.call(rbind, peak_ints_list)
  colnames(peak_ints_dt) <- spec_files
  
  peak_mzs_dt <- cbind(peak_mzs, peak_ints_dt)
  mSet@MSnData[["peak_matrix"]] <- peak_mzs_dt;
  return(mSet)
}



#' PerformIndividualInusionFeatureDection
#'
#' @param mSet 
#' @param method 
#' @param noise 
#' @param snthresh 
#' @param ppm 
#' @param scales 
#' @param ncores 
#'
#' @return NA
#' @export
#'
#' @examples todo
PerformIndividualInusionFeatureDection <- function(mSet, method = "cwt", noise = NA, 
                                         snthresh = 3, ppm = 10, 
                                         scales = c(1,4,9), 
                                         ncores = 1L){
  
  mSet@MSnData$scan_ms1 -> scan_ms1
  mSet@MSnData[["fileNames"]] -> spec_files;
  if(length(scan_ms1)>1){
    stop("This function cannot be used to perform multiple spectra files.")
  }
  if(ncores == 1L){
    require("MassSpecWavelet")
    peaks_list <- rep(list(NA), length(scan_ms1))
    #for(i in 1:length(scan_ms1)){
    i = 1;
      message("Processing spectra file == ", spec_files[i], " == ", appendLF = F)
      this_file <- scan_ms1[[i]];
      res <- lapply(1:length(this_file), function(x){
        #cat(x, "\n")
        y<-this_file[[x]]
        ints <- y[,2]
        mzs <- y[,1]
        if(method =="cwt"){
          peaks <- infusionDectectCWT(ints, mzs, scales = scales, snthresh = snthresh)
        } else {
          mzdiff <- 500*ppm*1e-6
          peaks <- infusionDectectKDE(ints, mzs, noise = noise, snthresh = snthresh, ppm = ppm, mzdiff = mzdiff)
        }
        
        return(cbind(peaks, sample = x))
      })
      
      peaks_table <- do.call(rbind, res)
      if(length(which(is.na(peaks_table[,1])))!=0){
        peaks_table <- peaks_table[-which(is.na(peaks_table[,1])), ]
      }
      
      if(nrow(peaks_table)==1){
        peaks_table <- as.data.frame(t(as.matrix(peaks_table)))
      } else {
        peaks_table <- as.data.frame(peaks_table)
      }
      
      message("...", appendLF = F)
      # need to do a merging of detected peaks within same file - mzCluster
      res_merged <- do_groupPeaks_mzClustx(peaks_table,
                                           sampleGroups = rep(1, length(this_file)),
                                           ppm = ppm,
                                           absMz = 0,
                                           minFraction = 0,
                                           minSamples = 0)
      
      res_merged[["featureDefinitions"]] -> merged_peak0
      res_merged[["peakIndex"]] -> merged_peak_idx
      col_idx <- which(colnames(peaks_table) == "sn")
      
      merged_peak_list <- lapply(merged_peak_idx, function(x){
        peaks_table_sub<- peaks_table[x,]
        row_idx <- which.max(peaks_table_sub[,col_idx])[1]
        return(peaks_table_sub[row_idx,])
      })
      merged_peak1 <- do.call(rbind, merged_peak_list)
      merged_peak1$sample <- i
      rownames(merged_peak1) <- NULL
      peaks_list[[i]] <- merged_peak1
      message(" ... Done!", appendLF = T)
    #}
  } else {
    # use multiple cores
    require("MassSpecWavelet")
    peaks_list <- rep(list(NA), length(scan_ms1))
    #for(i in 1:length(scan_ms1)){
    i = 1;
      message("Processing spectra file == ", spec_files[i], " == ", appendLF = F)
      this_file <- scan_ms1[[i]];
      # res <- bplapply(1:length(this_file), function(x){
      #   y<-this_file[[x]]
      #   ints <- y[,2]
      #   mzs <- y[,1]
      #   peaks <- infusionDectectCWT(ints, mzs)
      #   return(cbind(peaks, sample = x))
      # })
      require(parallel);
      cl <- makeCluster(getOption("cl.cores", ncores))
      clusterExport(cl, c("infusionDectectCWT", "this_file", "infusionDectectKDE",
                          "noise", "snthresh", "ppm", "scales",
                          "sumIntos", "maxIntos","do_groupPeaks_mzClustx", ".fix_mz_clust_peaks",
                          "mzClustGeneric", "mzClust_hclust_cpp", "R_mzClust_hclust_rcpp"), envir = environment())
      res <- list()
      
      if(method =="cwt"){
        res <- parLapply(cl, 
                         1:length(this_file), 
                         function(x){
                           y<-this_file[[x]]
                           ints <- y[,2]
                           mzs <- y[,1]
                           peaks <- infusionDectectCWT(ints, mzs, scales = scales, snthresh = snthresh)
                           return(cbind(peaks, sample = x))
                         })
      } else {
        res <- parLapply(cl, 
                         1:length(this_file), 
                         function(x){
                           y<-this_file[[x]]
                           ints <- y[,2]
                           mzs <- y[,1]
                           mzdiff <- 500*ppm*1e-6
                           peaks <- infusionDectectKDE(ints, mzs, noise = noise, snthresh = snthresh, ppm = ppm, mzdiff = mzdiff)
                           return(cbind(peaks, sample = x))
                         })
      }
      stopCluster(cl)
      
      peaks_table <- do.call(rbind, res)
      if(length(which(is.na(peaks_table[,1])))!=0){
        peaks_table <- peaks_table[-which(is.na(peaks_table[,1])), ]
      }
      if(nrow(peaks_table)==1){
        peaks_table <- as.data.frame(t(as.matrix(peaks_table)))
      } else {
        peaks_table <- as.data.frame(peaks_table)
      }
      
      message("...", appendLF = F)
      # need to do a merging of detected peaks within same file - mzCluster
      res_merged <- do_groupPeaks_mzClustx(peaks_table,
                                           sampleGroups = rep(1, length(this_file)),
                                           ppm = ppm,
                                           absMz = 0,
                                           minFraction = 0,
                                           minSamples = 0)
      
      res_merged[["featureDefinitions"]] -> merged_peak0
      res_merged[["peakIndex"]] -> merged_peak_idx
      col_idx <- which(colnames(peaks_table) == "sn")
      
      merged_peak_list <- lapply(merged_peak_idx, function(x){
        peaks_table_sub<- peaks_table[x,]
        row_idx <- which.max(peaks_table_sub[,col_idx])[1]
        return(peaks_table_sub[row_idx,])
      })
      merged_peak1 <- do.call(rbind, merged_peak_list)
      merged_peak1$sample <- i
      rownames(merged_peak1) <- NULL
      peaks_list[[i]] <- merged_peak1
      message(" Done! ", appendLF = T)
      
  }
  qs::qsave(peaks_list, file = paste0(sub(".mzML|.mzXML|.cdf|.CDF|.mzData", "", basename(spec_files[i])), "_peakTable.qs"))

}



#' FormatInfusionFeatureTable
#'
#' @param mSet 
#' @param minFraction 
#' @param minSamples 
#'
#' @return mSet Object
#' @export
#'
#' @examples todo
FormatInfusionFeatureTable <- function(mSet, minFraction = 0.1,  minSamples = 0){
  
  
}

infusionDectectCWT <- function(int, mz, snthresh = 3, scales = c(1, seq(2, 30, 2))){
  require("MassSpecWavelet")
  ## MassSpecWavelet Calls
  peakInfo <- peakDetectionCWT(int, SNR.Th = snthresh, scales = scales)
  if(length(peakInfo[["majorPeakInfo"]][["peakIndex"]]) == 0){
    df <- matrix(NA, ncol = 11)
    colnames(df) <- c( "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "maxo", "sn", "intf", "maxf")
    return(df)
  }
  majorPeakInfo <- peakInfo$majorPeakInfo
  

  betterPeakInfo <- tuneInPeakInfo(int,
                                   majorPeakInfo)
  peakIndex <- betterPeakInfo$peakIndex
  nPeaks <- length(peakIndex)
  
  ## sum and max of raw values, sum and max of filter-response
  rints <- numeric(nPeaks)
  fints <- NA
  maxRints <- numeric(nPeaks)
  maxFints <- NA
  
  for (a in 1:nPeaks) {
    rints[a] <- sumIntos(int, peakIndex[a],
                         betterPeakInfo$peakScale[a])
    maxRints[a] <- maxIntos(int, peakIndex[a],
                            betterPeakInfo$peakScale[a])
  }
  ## filter-response is not summed here, the maxF-value is the one
  ## which was "xcmsRaw$into" earlier
  
  ## Assemble result
  basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax",
                 "into","maxo","sn","intf","maxf")
  
  peaklist <- matrix(-1, nrow = nPeaks, ncol = length(basenames))
  colnames(peaklist) <- c(basenames)
  
  peaklist[,"mz"] <- mz[peakIndex]
  peaklist[,"mzmin"] <- mz[(peakIndex - betterPeakInfo$peakScale)]
  peaklist[,"mzmax"] <- mz[(peakIndex + betterPeakInfo$peakScale)]
  
  ## peaklist[,"rt"]    <- rep(-1, length(peakIndex))
  ## peaklist[,"rtmin"] <- rep(-1, length(peakIndex))
  ## peaklist[,"rtmax"] <- rep(-1, length(peakIndex))
  
  peaklist[,"into"] <- rints ## sum of raw-intensities
  peaklist[,"maxo"] <- maxRints
  peaklist[,"intf"] <- rep(NA, nPeaks)
  peaklist[,"maxf"] <- betterPeakInfo$peakValue
  
  peaklist[,"sn"]   <- betterPeakInfo$peakSNR
  
  return(peaklist)
}

infusionDectectKDE <- function(ints, mzs, snthresh = 3, noise = 100, ppm = 10, mzdiff = 0.001){
  
  # there are two case 
  # 1). The ions are already well controided, minimum mz_diff > ppm
  # 2). The ions are not well-controided, we need to identify mzs from same ion from others
  
  ion_idx <- ints>noise

  ints_filtered <- ints[ion_idx];
  mzs_filtered <- mzs[ion_idx]
  
  mz_max_num <- length(mzs)
  sn_values <- vapply(1:length(mzs_filtered), function(x){
    
    idx_mz <- which(mzs_filtered[x] == mzs)
    snr = 1
    if(idx_mz == 1){
      if(ints[2]>noise){
        snr <- ints[1]/ints[2]
      } else {
        snr <- ints[1]/noise
      }
    } else if(idx_mz == mz_max_num) {
      if(ints[mz_max_num-1]>noise){
        snr <- ints[mz_max_num]/ints[mz_max_num-1]
      } else {
        snr <- ints[mz_max_num]/noise
      }
    } else {
      if(mean(ints[idx_mz-1], ints[idx_mz+1])>noise){
        snr <- ints[idx_mz]/mean(ints[idx_mz-1], ints[idx_mz+1])
      } else {
        snr <- ints[idx_mz]/noise
      }
    }
    return(snr)
  }, double(length = 1L))
  
  # ion_idx2 <- sn_values>snthresh
  # 
  # ints_filtered <- ints_filtered[ion_idx2];
  # mzs_filtered <- mzs_filtered[ion_idx2]
  # sn_values <- sn_values[ion_idx2]
  
  
  mz_diff <- diff(mzs_filtered)
  if(length(mz_diff)==0){
    min_mz_diff <- 10
  } else {
    min_mz_diff <- min(mz_diff, na.rm = T)
  }
  
  
  # 1er case
  {
    if(min_mz_diff>mzdiff){
      
      res_table <- data.frame(mz= mzs_filtered, 
                              mzmin = mzs_filtered,
                              mzmax = mzs_filtered, 
                              rt = -1,
                              rtmin = -1,
                              rtmax = -1, 
                              into = ints_filtered,
                              maxo = ints_filtered,
                              sn = sn_values,
                              intf = NA, 
                              maxf = ints_filtered)
      
    }
  }

  # 2me case
  {
    if(min_mz_diff<=mzdiff){
      
      mz_diff_sorted <- sort(mz_diff)
      # Kernel density estimate
      dens <- density(mz_diff_sorted, bw = 500*ppm*1e-6)
      #plot(dens)
      dense_log_vals <- log(dens[["y"]])
      idx_cutoff0 <- which(dense_log_vals<0)[1]
      
      if(length(idx_cutoff0)!=0){
        idx_cutoff <- idx_cutoff0-1
        dens$x <- dens$x[1:idx_cutoff]
        dens$y <- log(dens$y[1:idx_cutoff])

      } else{
        # to do later
        stop("todo 1")
      }
      
      leng_num <- length(dens$x)
      if(leng_num > 15){
        pos_cutoff <- ceiling(leng_num*33/100) # last 1/3 
        
        # fit a linear model
        lineardata <- data.frame(
          x = dens$x[(leng_num-pos_cutoff):leng_num],
          y = dens$y[(leng_num-pos_cutoff):leng_num]
        )
        
        model <- lm(lineardata$y ~ lineardata$x)
        
        # View model summary
        coef_res <- coef(model)
        coef_val <- unname(coef_res[2])
        
        y_smooth <- predict(loess(dens$y ~ dens$x, span=0.5))
        dy <- diff(y_smooth) / diff(dens$x)
        
        coef_val_thresh <- coef_val*1.5
        
        dy_res <- dy < coef_val_thresh
        groups_mtx <- suppressWarnings(matrix(dy_res, ncol = 3, byrow = TRUE))
        idx_bool <- apply(groups_mtx, 1, function(x){
          length(which(x))>1
        })
        idx<- (which(!idx_bool)[1]-1)*3
        mz_diff_val <- mz_diff_sorted[idx]
        
        idx_new_peaks <- c(1, which(mz_diff>=mz_diff_val)+1)
        
        res_df_list <- lapply(1:length(idx_new_peaks), function(x){
          if(x == length(idx_new_peaks)) {
            # the last ion
            mzs_this_peak <- mzs_filtered[idx_new_peaks[x]]
            ints_this_peak <- ints_filtered[idx_new_peaks[x]]
            
            sn_ratio <- sn_values[idx_new_peaks[x]]
            res_df1 <- data.frame(mz = mzs_this_peak, mzmin = mzs_this_peak, mzmax = mzs_this_peak, 
                                  into = ints_this_peak, maxo = ints_this_peak, sn = sn_ratio)
          } else if((idx_new_peaks[x+1] - idx_new_peaks[x])>1){
            # multiple ion
            all_mzs_this_peak <- mzs_filtered[c(idx_new_peaks[x]:(idx_new_peaks[x+1]-1))]
            all_ints_this_peak <- ints_filtered[c(idx_new_peaks[x]:(idx_new_peaks[x+1]-1))]
            
            mz <- mean(all_mzs_this_peak)
            mzmin <- min(all_mzs_this_peak)
            mzmax <- max(all_mzs_this_peak)
            into <- sum(all_ints_this_peak)
            maxo <- max(all_ints_this_peak)
            sn_ratio <- mean(sn_values[c(idx_new_peaks[x]:(idx_new_peaks[x+1]-1))])
            res_df1 <- data.frame(mz = mz, mzmin = mzmin, mzmax = mzmax, into = into, maxo = maxo, sn = sn_ratio)
          } else {
            # single ion
            mzs_this_peak <- mzs_filtered[idx_new_peaks[x]]
            ints_this_peak <- ints_filtered[idx_new_peaks[x]]
            
            sn_ratio <- sn_values[idx_new_peaks[x]]
            res_df1 <- data.frame(mz = mzs_this_peak, mzmin = mzs_this_peak, mzmax = mzs_this_peak, 
                                  into = ints_this_peak, maxo = ints_this_peak, sn = sn_ratio)
          }
          return(res_df1)
        })
        
        res_df_table <- do.call(rbind, res_df_list)
        
        
        res_table <- data.frame(mz= res_df_table$mz, 
                                mzmin = res_df_table$mzmin,
                                mzmax = res_df_table$mzmax, 
                                rt = -1,
                                rtmin = -1,
                                rtmax = -1, 
                                into = res_df_table$into,
                                maxo = res_df_table$maxo,
                                sn = res_df_table$sn,
                                intf = NA, 
                                maxf = res_df_table$into)
      } else {
        res_table <- data.frame(mz= mzs_filtered, 
                                mzmin = mzs_filtered,
                                mzmax = mzs_filtered, 
                                rt = -1,
                                rtmin = -1,
                                rtmax = -1, 
                                into = ints_filtered,
                                maxo = ints_filtered,
                                sn = sn_values,
                                intf = NA, 
                                maxf = ints_filtered)
      }

      
      # plot(dens$x, y_smooth)
      # coef_val_all <- vapply(2:(leng_num-pos_cutoff), function(x){
      #   lmdata <- data.frame(
      #     x = dens$x[1:x],
      #     y = dens$y[1:x]
      #   )
      #   
      #   modelx <- lm(lmdata$y ~ lmdata$x)
      #   coef_valz <- unname(abs(coef(modelx)[2]))
      #   coef_valz
      # }, FUN.VALUE = double(1L))
      # 
      # plot(dens)
      # abline(model, col="red", lwd=2)
      
      
    }
    
  }

  return(res_table)
  
}

sumIntos <- function(into, inpos, scale){
  scale = floor(scale)
  sum(into[(inpos-scale):(inpos+scale)])
}
maxIntos <- function(into, inpos, scale){
  scale = floor(scale)
  max(into[(inpos-scale):(inpos+scale)])
}

do_groupPeaks_mzClustx <- function (peaks, sampleGroups, ppm = 20, absMz = 0, minFraction = 0.5, 
          minSamples = 1) 
{
  if (missing(sampleGroups)) 
    stop("Parameter 'sampleGroups' is missing! This should be a vector of ", 
         "length equal to the number of samples specifying the group ", 
         "assignment of the samples.")
  if (missing(peaks)) 
    stop("Parameter 'peaks' is missing!")
  if (!(is.matrix(peaks) | is.data.frame(peaks))) 
    stop("Peaks has to be a 'matrix' or a 'data.frame'!")
  .reqCols <- c("mz", "sample")
  if (!all(.reqCols %in% colnames(peaks))) 
    stop("Required columns ", paste0("'", .reqCols[!.reqCols %in% 
                                                     colnames(peaks)], "'", collapse = ", "), " not found in 'peaks' parameter")
  if (!is.factor(sampleGroups)) 
    sampleGroups <- factor(sampleGroups, levels = unique(sampleGroups))
  sampleGroupNames <- levels(sampleGroups)
  sampleGroupTable <- table(sampleGroups)
  nSampleGroups <- length(sampleGroupTable)
  if (max(peaks[, "sample"]) > length(sampleGroups)) 
    stop("Sample indices in 'peaks' are larger than there are sample", 
         " groups specified with 'sampleGroups'!")
  peaks <- .fix_mz_clust_peaks(peaks)
  grps <- mzClustGeneric(peaks[, .reqCols, drop = FALSE], 
                         sampclass = sampleGroups, mzppm = ppm, mzabs = absMz, 
                         minsamp = minSamples, minfrac = minFraction)
  grpmat <- grps$mat
  if (is.null(nrow(grpmat))) {
    matColNames <- names(grpmat)
    grpmat <- matrix(grpmat, ncol = length(grpmat), byrow = FALSE)
    colnames(grpmat) <- matColNames
  }
  rts <- rep(-1, nrow(grpmat))
  cns <- colnames(grpmat)
  grpmat <- cbind(grpmat[, 1:3, drop = FALSE], rts, rts, rts, 
                  grpmat[, 4:ncol(grpmat), drop = FALSE])
  colnames(grpmat) <- c(cns[1:3], c("rtmed", "rtmin", "rtmax"), 
                        cns[4:length(cns)])
  return(list(featureDefinitions = grpmat, peakIndex = grps$idx))
}

.fix_mz_clust_peaks <- function (x) 
{
  nas <- is.na(x[, "mz"])
  if (any(nas)) {
    if (all(c("mzmin", "mzmax") %in% colnames(x)) && !any(is.na(x[nas, 
                                                                  c("mzmin", "mzmax")]))) {
      warning("Got ", sum(nas), " peaks with missing values in column ", 
              "'mz'. Replaced them with the mean of values in columns ", 
              "'mzmin' and 'mzmax' values.")
      x[nas, "mz"] <- rowMeans(x[nas, c("mzmin", "mzmax")])
    }
    else {
      stop("Got ", sum(nas), " peaks with missing values in column 'mz'.")
    }
  }
  x
}

mzClustGeneric <- function (p, sampclass = NULL, mzppm = 20, mzabs = 0, minsamp = 1, 
          minfrac = 0.5) 
{
  makeBin <- function(pos) {
    if (pos > numpeaks) 
      return(list(pos = pos, bin = c(-1)))
    bin <- pord[pos]
    pos <- pos + 1
    basepeak <- p[bin[1], 1]
    error_range <- c(basepeak, basepeak * error_window + 
                       basepeak + 2 * mzabs)
    while (pos < numpeaks && p[pord[pos], 1] <= error_range[2]) {
      bin <- c(bin, pord[pos])
      pos <- pos + 1
    }
    if (pos%%(numpeaks%/%100 + 1) == 0) {
      # cat(format(((pos - 1)/numpeaks * 100), digits = 1, 
      #            nsmall = 2), " ")
      flush.console()
    }
    lst <- list(pos = pos, bin = bin)
    lst
  }
  meanDeviationOverLimit <- function(bin) {
    bin_mz <- p[bin, 1]
    m <- mean(bin_mz)
    error_range <- c(m - ppm_error * m - mzabs, ppm_error * 
                       m + m + mzabs)
    if (length(bin_mz[(bin_mz > error_range[2]) | (bin_mz < 
                                                   error_range[1])]) > 0) {
      return(TRUE)
    }
    else {
      FALSE
    }
  }
  bin2output <- function(bin) {
    gcount <- integer(length(classnum))
    if (length(gcount) != 0) {
      for (i in seq(along = bin)) {
        class_idx <- sampclass[p[bin[i], 2]]
        gcount[class_idx] <- gcount[class_idx] + 1
      }
    }
    if (length(bin) < minsamp || (!any(gcount >= classnum * 
                                       minfrac) && length(gcount) > 0)) 
      return(list())
    groupvec <- c(rep(NA, 4 + length(gcount)))
    groupvec[1] <- mean(p[bin, 1])
    groupvec[2:3] <- range(p[bin, 1])
    groupvec[4] <- length(bin)
    sorted <- order(p[bin, 1])
    grp_members <- bin[sorted]
    groupvec[4 + seq(along = gcount)] <- gcount
    lst <- list(stat = groupvec, members = grp_members)
    lst
  }
  ppm_error <- mzppm/1e+06
  error_window <- 2 * ppm_error
  if (is.null(sampclass)) {
    classnum <- integer(0)
    classnames <- seq(along = classnum)
  }
  else {
    classnames <- levels(sampclass)
    sampclass <- as.vector(unclass(sampclass))
    classnum <- integer(max(sampclass))
  }
  for (i in seq(along = classnum)) classnum[i] <- sum(sampclass == 
                                                        i)
  numpeaks <- nrow(p)
  groupmat <- matrix(nrow = 512, ncol = 4 + length(classnum))
  groupindex <- vector("list", 512)
  pord <- order(p[, 1])
  pos <- c(1)
  binNumber <- 1
  newbin <- makeBin(pos)
  binA <- newbin$bin
  pos <- newbin$pos
  while (1) {
    if (binNumber + 4 > nrow(groupmat)) {
      groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), 
                                         ncol = ncol(groupmat)))
      groupindex <- c(groupindex, vector("list", length(groupindex)))
    }
    newbin <- makeBin(pos)
    binB <- newbin$bin
    pos <- newbin$pos
    if (binB[1] < 0) {
      out <- bin2output(binA)
      if (length(out) != 0) {
        groupmat[binNumber, ] <- out$stat
        groupindex[[binNumber]] <- out$members
        binNumber <- binNumber + 1
      }
      break
    }
    max_binA <- max(p[binA, 1])
    min_binB <- min(p[binB, 1])
    binclust <- 0
    if (max_binA + max_binA * error_window + 2 * mzabs >= 
        min_binB && min_binB - min_binB * error_window - 
        2 * mzabs <= max_binA) {
      binC <- c(binA, binB)
      binclust <- 1
    }
    else {
      if (meanDeviationOverLimit(binA)) {
        binC <- binA
        binclust <- 2
      }
    }
    if (binclust != 0) {
      groups <- mzClust_hclust_cpp(p[binC, 1], ppm_error, 
                               mzabs)
      last_group <- groups[which.max(p[binC, 1])]
      binA <- binC[which(groups == last_group)]
      if (max(groups) > 1) {
        for (c in 1:max(groups)) {
          if (c == last_group) {
            next
          }
          tmp_grp <- which(groups == c)
          tmp_c <- binC[tmp_grp]
          out <- bin2output(tmp_c)
          if (length(out) != 0) {
            groupmat[binNumber, ] <- out$stat
            groupindex[[binNumber]] <- out$members
            binNumber <- binNumber + 1
          }
        }
      }
    }
    if (binclust != 1) {
      out <- bin2output(binA)
      if (length(out) != 0) {
        groupmat[binNumber, ] <- out$stat
        groupindex[[binNumber]] <- out$members
        binNumber <- binNumber + 1
      }
      binA <- binB
    }
  }
  colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "npeaks", 
                          classnames)
  binNumber <- binNumber - 1
  groupmat <- groupmat[seq(length = binNumber), ]
  groupindex <- groupindex[seq(length = binNumber)]
  # cat("\n")
  flush.console()
  return(list(mat = groupmat, idx = groupindex))
}

mzClust_hclust_cpp <- function(x, eppm, eabs)
{
  N <- length(x)
  d <- dist(x)
  g <- R_mzClust_hclust_rcpp(
    x = as.double(x),
    num = N,
    d = as.double(d),
    eppm = as.double(eppm),
    eabs = as.double(eabs))
  return(g)
}

# Rcpp::sourceCpp("src/mzClust.cpp")