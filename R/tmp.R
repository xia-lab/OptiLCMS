
#' @export
controller.modifier <- function(new_command_set, last_command_set, plan){
  
  ###-------------Operators definition ------------//
  
  #' operators_1 : PerformDataTrimming ~ PerformParamsOptimization;
  #' operators_2 : PerformParamsOptimization ~ peak_profiling;  #' 
  #' operators_3 : PerformPeakProfiling ~ PerformPeakAnnotation;  #' 
  #' operators_4 : ;
  
  
  #' others_1 c1 : Control Optimization or not;
  #' others_1 c2 : ;
  #' others_1 c3 : ;
  #' others_1 c4 : .
  
  ##----------------------------------------------//
  
  if(class(new_command_set) == "OptiCommandSet"){
    
    # 1. ROIExtraction: -----
    # 1.1 Note on controller: c1, read; c2, trim; c3, write; c4, plot; C5, rmConts.
    new_ROIExtraction <- new_command_set@ROIExtraction[[3]];
    last_ROIExtraction <- last_command_set@ROIExtraction[[3]];
    ChangedArugsArray <- NULL;
    
    for (i in 2:length(new_ROIExtraction)){
      for (j in 2:length(last_ROIExtraction)){
        if (names(new_ROIExtraction[i]) == names(last_ROIExtraction[j]) & 
            new_ROIExtraction[[i]] != last_ROIExtraction[[j]]) {
          ChangedArugsArray <- c(ChangedArugsArray, names(new_ROIExtraction[i]))
        } 
      }
    }

    newArugsNMs <- names(new_ROIExtraction)[-1];
    lastArugsNMs <- names(last_ROIExtraction)[-1];
    NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
    RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
    
    ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray)
    
    if(is.null(ChangedArugsArray) | length(ChangedArugsArray) == 0){
      # Not changed !
      plan@running.controller@ROI_extract[c(1:5)] <- rep(FALSE, 5);
    }
    
    if("datapath" %in% ChangedArugsArray){
      plan@running.controller@ROI_extract[c(1:5)] <- rep(TRUE,5);
      plan@running.controller@operators[1] <- TRUE;
    } else {
      plan@running.controller@ROI_extract[1] <- FALSE;
    }
    
    if("rmConts" %in% ChangedArugsArray){
      plan@running.controller@ROI_extract[5] <- TRUE;
      plan@running.controller@operators[1] <- TRUE;
    }else {
      plan@running.controller@ROI_extract[5] <- FALSE;
    }
    
    if(any(c("mode", "mz", "mzdiff", "rt", "rtdiff", "rt.idx") %in% ChangedArugsArray)){
      # if any of these params changed, need to run the following steps
      plan@running.controller@ROI_extract[c(2:4)] <- rep(TRUE,3);
      plan@running.controller@operators[1] <- TRUE;
    } else if (any(plan@running.controller@ROI_extract[c(1,5)])){
      # if datapath changed, still need to run the following steps
      plan@running.controller@ROI_extract[c(2:4)] <- rep(TRUE,3);
      plan@running.controller@operators[1] <- TRUE;
    } else {
      # if all these params and datapath did not changed, skipped trimming step
      plan@running.controller@ROI_extract[2] <- FALSE;
    }
    
    if("write" %in% ChangedArugsArray){
      plan@running.controller@ROI_extract[3] <- TRUE;
    }
    
    if("plot" %in% ChangedArugsArray){
      plan@running.controller@ROI_extract[4] <- TRUE;
    }
    
    # 2. ParamsOptimization: -----
    # 2.1 Note on controller: others c1, run or not.
    new_ParamsOptimization <- new_command_set@ParamsOptimization[[3]];
    last_ParamsOptimization <- last_command_set@ParamsOptimization[[3]];
    ChangedArugsArray <- NULL;
    
    for (i in 2:length(new_ParamsOptimization)){
      for (j in 2:length(last_ParamsOptimization)){
        if (names(new_ParamsOptimization[i]) == names(last_ParamsOptimization[j]) & 
            new_ParamsOptimization[[i]] != last_ParamsOptimization[[j]]) {
          ChangedArugsArray <- c(ChangedArugsArray, names(new_ParamsOptimization[i]))
        } 
      }
    }
    
    newArugsNMs <- names(new_ParamsOptimization)[-1];
    lastArugsNMs <- names(last_ParamsOptimization)[-1];
    NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
    RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
    
    ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray)
    
    if(is.null(ChangedArugsArray) | length(ChangedArugsArray) == 0){
      # Not changed !
      plan@running.controller@others_1[1] <- FALSE;
    }
    
    if("param" %in% ChangedArugsArray){
      
      if(length(RmArugsArray == "param") != 0){
        if(RmArugsArray == "param"){
          # refer to: the omitted Argus is param
          if(last_ParamsOptimization$param == "SetPeakParam()"){
            # handle: the omitted param is SetPeakParam()
            plan@running.controller@others_1[1] <- FALSE;
          } else {
            # handle: the omitted param is not SetPeakParam()
            plan@running.controller@others_1[1] <- TRUE;
          }
        }
      }

      if (length(NewArugsArray == "param") != 0) {
        if (NewArugsArray == "param") {
          # refer to: the added Argus is param
          if (new_ParamsOptimization$param == "SetPeakParam()") {
            # handle: the added param is SetPeakParam()
            plan@running.controller@others_1[1] <- FALSE;
          } else {
            # handle: the added param is not SetPeakParam()
            plan@running.controller@others_1[1] <- TRUE;
          }
        }
      }
      
      if (length(NewArugsArray == "param") == 0 & 
          length(RmArugsArray == "param") == 0) {
        plan@running.controller@others_1[1] <- TRUE;
      }
      
    }
    
    if(plan@running.controller@operators[1]){
      # Data ROI changed ! Optimization has to be re-run !
      plan@running.controller@others_1[1] <- TRUE;
    }
    
    if(plan@running.controller@others_1[1]){
      # the peak profiling also need to be re-run
      plan@running.controller@operators[2] <- TRUE;
    }
  }
  
  # 3. ImportRawMSData: -----
  # 3.1 Note on controller: others c2, data import. c3, plot.
  
  new_ImportRawMSData <- new_command_set@ImportRawMSData[[3]];
  last_ImportRawMSData <- last_command_set@ImportRawMSData[[3]];
  ChangedArugsArray <- NULL;
  
  for (i in 2:length(new_ImportRawMSData)){
    for (j in 2:length(last_ImportRawMSData)){
      if (names(new_ImportRawMSData[i]) == names(last_ImportRawMSData[j]) & 
          new_ImportRawMSData[[i]] != last_ImportRawMSData[[j]]) {
        ChangedArugsArray <- c(ChangedArugsArray, names(new_ImportRawMSData[i]))
      } 
    }
  }
  
  newArugsNMs <- names(new_ImportRawMSData)[-1];
  lastArugsNMs <- names(last_ImportRawMSData)[-1];
  NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
  RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
  
  ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray);
  
  if(is.null(ChangedArugsArray) | length(ChangedArugsArray) == 0){
    # Not changed !
    plan@running.controller@data_import[c(1:2)] <- rep(TRUE, 2);
  }
  
  if("foldername" %in% ChangedArugsArray){
    plan@running.controller@data_import[1] <- TRUE;
  }
  
  if("plotSettings" %in% ChangedArugsArray){
    plan@running.controller@data_import[2] <- TRUE;
  }
  
  # 4. PeakProfiling: -----
  # 4.1 Note on controller peak_profiling: c1, peak picking; c2, peak alignment; c3, peak filing; c4, plotting.
  
  new_PeakProfiling <- new_command_set@PeakProfiling[[3]];
  last_PeakProfiling <- last_command_set@PeakProfiling[[3]];
  ChangedArugsArray <- NULL;
  
  for (i in 2:length(new_PeakProfiling)){
    for (j in 2:length(last_PeakProfiling)){
      if (names(new_PeakProfiling[i]) == names(last_PeakProfiling[j])) {
        if(new_PeakProfiling[[i]] != last_PeakProfiling[[j]]){
          ChangedArugsArray <- c(ChangedArugsArray, names(new_PeakProfiling[i]))
        }
      } 
    }
  }
  
  newArugsNMs <- names(new_PeakProfiling)[-1];
  lastArugsNMs <- names(last_PeakProfiling)[-1];
  NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
  RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
  
  ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray);
  
  if(is.null(ChangedArugsArray) | length(ChangedArugsArray) == 0){
    # Not changed !
    plan@running.controller@peak_profiling <- rep(FALSE, 4);
    names(plan@running.controller@peak_profiling) <- c("c1","c2","c3","c4");
  }
  
  if(plan@running.controller@operators[2]){
    # OPtimized parameters changed ! re-run everything in peak profiling.
    plan@running.controller@peak_profiling <- rep(TRUE, 4);
    names(plan@running.controller@peak_profiling) <- c("c1","c2","c3","c4");
    
    plan@running.controller@operators[3] <- TRUE;
  }
  
  if("plotSettings" %in% ChangedArugsArray){
    plan@running.controller@peak_profiling[4] <- TRUE;
  }

  if(.on.public.web){
    # load params.rda and envir.rds and do a comparison
    envir.path <- paste0(getwd(),"/temp/envir");
    envir_tmp <- readRDS(paste0(envir.path,"/envir.rds"));
    last_param <- envir_tmp[["param"]];
    load("params.rda");
    new_param <- peakParams;
    # TODO: need to configure with the web; handle the annotate params at the same time
    
  } else {
    
    if("Params" %in% ChangedArugsArray){
      
      newParamArgus <- new_PeakProfiling$Params;
      lastParamArgus <- last_PeakProfiling$Params;
      
      ChangedParamArgus <- ParamsChanged(lastParamArgus, newParamArgus);

      if(is.null(ChangedParamArgus)){
        plan@running.controller@peak_profiling[c(1:3)] <- rep(FALSE,3);
      } else if (any(ChangedParamArgus %in% c("min_peakwidth","max_peakwidth","mzdiff","ppm","noise","prefilter","value_of_prefilter",
                                              "Peak_method","snthresh","fwhm","sigma","steps"))){
        plan@running.controller@peak_profiling[c(1:3)] <- rep(TRUE,3);
        plan@running.controller@operators[3] <- TRUE;
      } else if (any(ChangedParamArgus %in% c("bw","RT_method","minFraction","minSamples","maxFeatures","family","smooth",
                                              "span","integrate","mzCenterFun","verbose.columns","fitgauss"))){
        plan@running.controller@peak_profiling[c(1:3)] <- c(FALSE, TRUE, TRUE);
        plan@running.controller@operators[3] <- TRUE;
      } else if (ChangedParamArgus == "all"){
        plan@running.controller@peak_profiling[c(1:3)] <- rep(TRUE,3);
        plan@running.controller@operators[3] <- TRUE;
      }
    }
  }
  
  # 5. PeakAnnotation: -----
  # 5.1 Note on controller PeakAnnotation: c1, peak annotation; others not used for now.

  new_PeakAnnotation <- new_command_set@PeakAnnotation[[3]];
  last_PeakAnnotation <- last_command_set@PeakAnnotation[[3]];
  ChangedArugsArray <- NULL;
  
  for (i in 2:length(new_PeakAnnotation)){
    for (j in 2:length(last_PeakAnnotation)){
      if (names(new_PeakAnnotation[i]) == names(last_PeakAnnotation[j])) {
        if(new_PeakAnnotation[[i]] != last_PeakAnnotation[[j]]){
          ChangedArugsArray <- c(ChangedArugsArray, names(new_PeakAnnotation[i]))
        }
      } 
    }
  }
  
  newArugsNMs <- names(new_PeakAnnotation)[-1];
  lastArugsNMs <- names(last_PeakAnnotation)[-1];
  NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
  RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
  
  ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray);
  
  if(.on.public.web){
    # load params.rda and envir.rds and do a comparison
    # TODO: need to configure with the web and profiling
    
  } else {
    
    if("Params" %in% ChangedArugsArray){
      
      newAnnoParamArgus <- new_PeakAnnotation$annotaParam;
      lastAnnoParamArgus <- last_PeakAnnotation$annotaParam;
      
      ChangedParamArgus <- ParamsChanged(lastAnnoParamArgus, newAnnoParamArgus);
      
      if(is.null(ChangedParamArgus)){
        plan@running.controller@peak_annotation[1] <- FALSE;
      } else {
        plan@running.controller@peak_annotation[1] <- TRUE;
      }
    }
    
  }
  
  if(plan@running.controller@operators[3]){
    plan@running.controller@peak_annotation[1] <- TRUE;
  }
  
  # 6. FormatPeakList: -----
  # 6.1 Note on controller FormatPeakList.
  # No need to resume for this function
  
  

  return(plan)
  
  # 
}
